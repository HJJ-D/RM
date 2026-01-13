"""
Simple, robust single-cell tracking for microscopy time-lapse data.

Based on:
- Frame-wise global assignment using Hungarian matching
- IoU overlap matching between instance masks
- Centroid distance gating
- Optional gap closing (allow missing one frame with strict criteria)

No deep learning. No external tracking libraries.
Uses only: numpy, scikit-image, scipy
"""

import numpy as np
from skimage.measure import regionprops
from scipy.optimize import linear_sum_assignment
from typing import List, Dict, Tuple


def extract_instances(label_img: np.ndarray) -> List[Dict]:
    """
    Extract instance properties from a label image.
    
    Parameters
    ----------
    label_img : np.ndarray
        2D label image where 0 = background, 1..N = cell instances
        
    Returns
    -------
    List[Dict]
        List of dictionaries, one per instance, containing:
        - label: instance ID in this frame
        - bbox: (minr, minc, maxr, maxc)
        - centroid: (r, c)
        - area: number of pixels
        - major_axis_length: major axis length
        - minor_axis_length: minor axis length
        - eccentricity: eccentricity
        - crop_mask: boolean array of the mask cropped to bbox
    """
    instances = []
    
    # Extract region properties
    regions = regionprops(label_img)
    
    for region in regions:
        # Get bounding box: (min_row, min_col, max_row, max_col)
        bbox = region.bbox  # returns (minr, minc, maxr, maxc)
        
        # Extract crop mask: boolean array of the region within its bbox
        crop_mask = region.image  # binary mask of the region
        
        # Build instance dictionary
        instance = {
            'label': region.label,
            'bbox': bbox,
            'centroid': region.centroid,  # (row, col)
            'area': region.area,
            'major_axis_length': region.major_axis_length,
            'minor_axis_length': region.minor_axis_length,
            'eccentricity': region.eccentricity,
            'crop_mask': crop_mask,
        }
        
        instances.append(instance)
    
    return instances


def iou_instance(a: Dict, b: Dict) -> float:
    """
    Compute IoU (Intersection over Union) between two instances.
    
    Uses bounding box intersection and crop masks for efficiency.
    Does NOT compute IoU over the whole image.
    
    Parameters
    ----------
    a, b : Dict
        Instance dictionaries with 'bbox' and 'crop_mask' keys
        
    Returns
    -------
    float
        IoU value in [0, 1]
    """
    # Unpack bounding boxes
    minr_a, minc_a, maxr_a, maxc_a = a['bbox']
    minr_b, minc_b, maxr_b, maxc_b = b['bbox']
    
    # Compute intersection of bounding boxes
    inter_minr = max(minr_a, minr_b)
    inter_minc = max(minc_a, minc_b)
    inter_maxr = min(maxr_a, maxr_b)
    inter_maxc = min(maxc_a, maxc_b)
    
    # Check if bboxes overlap
    if inter_minr >= inter_maxr or inter_minc >= inter_maxc:
        return 0.0
    
    # Compute the offset of the intersection bbox relative to each instance's bbox
    # For instance a
    offset_r_a = inter_minr - minr_a
    offset_c_a = inter_minc - minc_a
    # For instance b
    offset_r_b = inter_minr - minr_b
    offset_c_b = inter_minc - minc_b
    
    # Intersection dimensions
    inter_h = inter_maxr - inter_minr
    inter_w = inter_maxc - inter_minc
    
    # Extract the intersection region from each crop mask
    crop_a_inter = a['crop_mask'][
        offset_r_a:offset_r_a + inter_h,
        offset_c_a:offset_c_a + inter_w
    ]
    crop_b_inter = b['crop_mask'][
        offset_r_b:offset_r_b + inter_h,
        offset_c_b:offset_c_b + inter_w
    ]
    
    # Compute intersection and union
    intersection = np.sum(crop_a_inter & crop_b_inter)
    union = a['area'] + b['area'] - intersection
    
    # Avoid division by zero
    if union == 0:
        return 0.0
    
    return intersection / union


def _euclidean_distance(centroid_a: Tuple[float, float], 
                        centroid_b: Tuple[float, float]) -> float:
    """Compute Euclidean distance between two centroids."""
    return np.sqrt((centroid_a[0] - centroid_b[0])**2 + 
                   (centroid_a[1] - centroid_b[1])**2)


def match_frame_pair(
    cells_t: List[Dict],
    cells_t1: List[Dict],
    max_dist: float = 25.0,
    min_iou: float = 0.2,
    len_ratio_range: Tuple[float, float] = (0.75, 1.25),
    area_ratio_range: Tuple[float, float] = (0.7, 1.3),
    w_iou: float = 0.6,
    w_dist: float = 0.3,
    w_shape: float = 0.1
) -> Tuple[Dict[int, int], List[float], List[float]]:
    """
    Match instances between two consecutive frames using Hungarian algorithm.
    
    Builds a cost matrix based on:
    - IoU overlap (lower cost for higher IoU)
    - Centroid distance (normalized by max_dist)
    - Shape similarity (major axis length difference)
    
    Then applies Hungarian matching and filters results based on thresholds.
    
    Parameters
    ----------
    cells_t : List[Dict]
        Instances at frame t
    cells_t1 : List[Dict]
        Instances at frame t+1
    max_dist : float
        Maximum centroid distance for valid matches
    min_iou : float
        Minimum IoU threshold for accepting a match
    len_ratio_range : Tuple[float, float]
        Valid range for major_axis_length ratio (min, max)
    area_ratio_range : Tuple[float, float]
        Valid range for area ratio (min, max)
    w_iou : float
        Weight for IoU cost component
    w_dist : float
        Weight for distance cost component
    w_shape : float
        Weight for shape cost component
        
    Returns
    -------
    mapping : Dict[int, int]
        Dictionary mapping index in cells_t to index in cells_t1
    matched_ious : List[float]
        IoU values for accepted matches
    matched_dists : List[float]
        Distance values for accepted matches
    """
    Nt = len(cells_t)
    Nt1 = len(cells_t1)
    
    # Handle empty frames
    if Nt == 0 or Nt1 == 0:
        return {}, [], []
    
    # Build cost matrix
    cost_matrix = np.zeros((Nt, Nt1))
    
    for i, cell_a in enumerate(cells_t):
        for j, cell_b in enumerate(cells_t1):
            # Compute centroid distance
            dist = _euclidean_distance(cell_a['centroid'], cell_b['centroid'])
            
            # If distance too large, assign high cost
            if dist > max_dist:
                cost_matrix[i, j] = 1e6
                continue
            
            # Compute IoU
            iou = iou_instance(cell_a, cell_b)
            
            # Normalize distance
            dist_norm = dist / max_dist
            
            # Compute shape cost (major axis length difference)
            major_a = cell_a.get('major_axis_length', 1e-6)
            major_b = cell_b.get('major_axis_length', 1e-6)
            shape_cost = min(1.0, abs(major_a - major_b) / max(major_a, 1e-6))
            
            # Combine costs
            cost = w_iou * (1.0 - iou) + w_dist * dist_norm + w_shape * shape_cost
            cost_matrix[i, j] = cost
    
    # Apply Hungarian algorithm
    row_ind, col_ind = linear_sum_assignment(cost_matrix)
    
    # Filter matches based on thresholds
    mapping = {}
    matched_ious = []
    matched_dists = []
    
    for i, j in zip(row_ind, col_ind):
        # Skip high-cost assignments
        if cost_matrix[i, j] >= 1e6:
            continue
        
        cell_a = cells_t[i]
        cell_b = cells_t1[j]
        
        # Compute actual metrics
        iou = iou_instance(cell_a, cell_b)
        dist = _euclidean_distance(cell_a['centroid'], cell_b['centroid'])
        
        # Check IoU threshold
        if iou < min_iou:
            continue
        
        # Check distance threshold
        if dist > max_dist:
            continue
        
        # Check shape consistency (major axis length OR area)
        major_a = cell_a.get('major_axis_length', 1e-6)
        major_b = cell_b.get('major_axis_length', 1e-6)
        len_ratio = major_b / max(major_a, 1e-6)
        
        area_a = cell_a.get('area', 1)
        area_b = cell_b.get('area', 1)
        area_ratio = area_b / max(area_a, 1)
        
        # Accept if either length ratio OR area ratio is within range
        len_ok = len_ratio_range[0] <= len_ratio <= len_ratio_range[1]
        area_ok = area_ratio_range[0] <= area_ratio <= area_ratio_range[1]
        
        if not (len_ok or area_ok):
            continue
        
        # Accept this match
        mapping[i] = j
        matched_ious.append(iou)
        matched_dists.append(dist)
    
    return mapping, matched_ious, matched_dists


def track_cells(
    cells_by_frame: List[List[Dict]],
    T_early: int = 60,
    max_dist: float = 25.0,
    min_iou: float = 0.2,
    allow_gap: bool = False
) -> Tuple[List[Dict], Dict]:
    """
    Track cells across frames using frame-wise Hungarian matching.
    
    Algorithm:
    - Only use frames 0..T_early-1
    - For each consecutive frame pair (t, t+1):
        - Use Hungarian algorithm to find global optimal assignment
        - Filter matches based on IoU, distance, and shape criteria
    - Build tracks by chaining frame-to-frame matches
    - Optionally allow gap closing with strict criteria
    
    Parameters
    ----------
    cells_by_frame : List[List[Dict]]
        List where cells_by_frame[t] contains instances at frame t
    T_early : int
        Number of early frames to use for tracking
    max_dist : float
        Maximum centroid distance for candidate matching
    min_iou : float
        Minimum IoU threshold for accepting a match (default: 0.2)
    allow_gap : bool
        If True, allow one missing frame with strict criteria
        (iou >= 0.3, dist <= max_dist/2, tight shape bounds)
        
    Returns
    -------
    tracks : List[Dict]
        List of tracks, where each track is:
        {"frames": [t0, t1, ...], "cells": [cell_dict0, cell_dict1, ...]}
    diagnostics : Dict
        Dictionary containing:
        - "all_matched_ious": list of all matched IoU values
        - "all_matched_dists": list of all matched distances
        - "matches_per_frame": list of number of matches per frame pair
    """
    # Ensure we don't go beyond available frames
    T_early = min(T_early, len(cells_by_frame))
    
    # Store all diagnostics
    all_matched_ious = []
    all_matched_dists = []
    matches_per_frame = []
    
    # Store forward mappings for each frame pair
    # forward_map[t] = dict mapping index in frame t -> index in frame t+1
    forward_map = {}
    
    # Compute frame-to-frame matches
    for t in range(T_early - 1):
        mapping, ious, dists = match_frame_pair(
            cells_by_frame[t],
            cells_by_frame[t + 1],
            max_dist=max_dist,
            min_iou=min_iou
        )
        forward_map[t] = mapping
        all_matched_ious.extend(ious)
        all_matched_dists.extend(dists)
        matches_per_frame.append(len(mapping))
    
    # Build reverse mapping for tracking which instances are matched
    # reverse_map[t] = set of indices in frame t that are matched FROM previous frame
    reverse_map = {t: set() for t in range(T_early)}
    for t, mapping in forward_map.items():
        for j in mapping.values():
            reverse_map[t + 1].add(j)
    
    # Build tracks by following forward chains
    tracks = []
    used_starts = set()  # Track which (frame, index) have been used as track starts
    
    # Start tracks from frame 0
    for i in range(len(cells_by_frame[0])):
        if (0, i) in used_starts:
            continue
        
        track_frames = [0]
        track_cells = [cells_by_frame[0][i]]
        used_starts.add((0, i))
        
        # Extend forward
        current_t = 0
        current_i = i
        
        while current_t < T_early - 1:
            next_t = current_t + 1
            
            # Check if there's a match to next frame
            if current_t in forward_map and current_i in forward_map[current_t]:
                next_i = forward_map[current_t][current_i]
                track_frames.append(next_t)
                track_cells.append(cells_by_frame[next_t][next_i])
                current_t = next_t
                current_i = next_i
                continue
            
            # Try gap closing if enabled
            if allow_gap and (current_t + 2) < T_early:
                gap_t = current_t + 2
                # Use stricter criteria for gap closing
                gap_mapping, gap_ious, gap_dists = match_frame_pair(
                    [cells_by_frame[current_t][current_i]],
                    cells_by_frame[gap_t],
                    max_dist=max_dist / 2,  # Stricter distance
                    min_iou=0.3,  # Stricter IoU
                    len_ratio_range=(0.85, 1.15),  # Tighter shape bounds
                    area_ratio_range=(0.8, 1.2)
                )
                
                # Check if we found a gap match
                if 0 in gap_mapping:
                    gap_i = gap_mapping[0]
                    # Make sure this instance isn't already used
                    if (gap_t, gap_i) not in used_starts and gap_i not in reverse_map.get(gap_t, set()):
                        track_frames.append(gap_t)
                        track_cells.append(cells_by_frame[gap_t][gap_i])
                        used_starts.add((gap_t, gap_i))
                        current_t = gap_t
                        current_i = gap_i
                        continue
            
            # No match found, end this track
            break
        
        tracks.append({
            'frames': track_frames,
            'cells': track_cells
        })
    
    # Start new tracks from unmatched instances in later frames
    for t in range(1, T_early):
        for i in range(len(cells_by_frame[t])):
            # Skip if already matched from previous frame
            if i in reverse_map[t]:
                continue
            
            # Skip if already used as a track start
            if (t, i) in used_starts:
                continue
            
            # Start new track
            track_frames = [t]
            track_cells = [cells_by_frame[t][i]]
            used_starts.add((t, i))
            
            # Extend forward
            current_t = t
            current_i = i
            
            while current_t < T_early - 1:
                next_t = current_t + 1
                
                # Check if there's a match to next frame
                if current_t in forward_map and current_i in forward_map[current_t]:
                    next_i = forward_map[current_t][current_i]
                    track_frames.append(next_t)
                    track_cells.append(cells_by_frame[next_t][next_i])
                    current_t = next_t
                    current_i = next_i
                    continue
                
                # Try gap closing if enabled
                if allow_gap and (current_t + 2) < T_early:
                    gap_t = current_t + 2
                    gap_mapping, gap_ious, gap_dists = match_frame_pair(
                        [cells_by_frame[current_t][current_i]],
                        cells_by_frame[gap_t],
                        max_dist=max_dist / 2,
                        min_iou=0.3,
                        len_ratio_range=(0.85, 1.15),
                        area_ratio_range=(0.8, 1.2)
                    )
                    
                    if 0 in gap_mapping:
                        gap_i = gap_mapping[0]
                        if (gap_t, gap_i) not in used_starts and gap_i not in reverse_map.get(gap_t, set()):
                            track_frames.append(gap_t)
                            track_cells.append(cells_by_frame[gap_t][gap_i])
                            used_starts.add((gap_t, gap_i))
                            current_t = gap_t
                            current_i = gap_i
                            continue
                
                # No match found, end this track
                break
            
            tracks.append({
                'frames': track_frames,
                'cells': track_cells
            })
    
    # Prepare diagnostics
    diagnostics = {
        'all_matched_ious': all_matched_ious,
        'all_matched_dists': all_matched_dists,
        'matches_per_frame': matches_per_frame
    }
    
    return tracks, diagnostics


def filter_tracks(
    tracks: List[Dict],
    min_len: int = 20,
    key: str = "major_axis_length",
    max_jump_ratio: float = 0.30
) -> List[Dict]:
    """
    Filter tracks based on quality control criteria.
    
    Parameters
    ----------
    tracks : List[Dict]
        List of tracks to filter
    min_len : int
        Minimum track length to keep
    key : str
        Morphology key to check for consistency (default: "major_axis_length")
    max_jump_ratio : float
        Maximum relative change allowed between consecutive frames
        
    Returns
    -------
    List[Dict]
        Filtered list of tracks
    """
    filtered = []
    
    for track in tracks:
        # Check minimum length
        if len(track['frames']) < min_len:
            continue
        
        # Check morphology consistency
        if not _is_morphology_consistent(track['cells'], key, max_jump_ratio):
            continue
        
        filtered.append(track)
    
    return filtered


def _is_morphology_consistent(
    cells: List[Dict],
    key: str,
    max_jump_ratio: float
) -> bool:
    """
    Check if morphology is consistent across a track.
    
    Rejects track if relative change between consecutive frames exceeds threshold.
    Falls back to 'area' if key is missing or zero.
    """
    for i in range(len(cells) - 1):
        # Get values
        val_curr = cells[i].get(key, 0)
        val_next = cells[i + 1].get(key, 0)
        
        # Fallback to area if key is zero or missing
        if val_curr == 0 or val_next == 0:
            val_curr = cells[i].get('area', 0)
            val_next = cells[i + 1].get('area', 0)
        
        # Skip if still zero (shouldn't happen for area)
        if val_curr == 0 or val_next == 0:
            continue
        
        # Compute relative change
        rel_change = abs(val_next - val_curr) / val_curr
        
        # Check threshold
        if rel_change > max_jump_ratio:
            return False
    
    return True


# ============================================================================
# DEMO USAGE
# ============================================================================

def demo_tracking(masks: np.ndarray, T_early: int = 60):
    """
    Minimal demo showing how to use the tracking pipeline.
    
    Parameters
    ----------
    masks : np.ndarray
        3D array of shape (T, H, W) containing label images
    T_early : int
        Number of early frames to track
        
    Example
    -------
    >>> # Assuming you have loaded your masks as a 3D numpy array
    >>> masks = ...  # shape (120, H, W)
    >>> demo_tracking(masks, T_early=60)
    """
    print("=" * 60)
    print("Single-Cell Tracking Demo (Hungarian Matching)")
    print("=" * 60)
    
    # Step 1: Extract instances per frame
    print(f"\nExtracting instances from frames 0..{T_early-1}...")
    cells_by_frame = [extract_instances(masks[t]) for t in range(T_early)]
    
    # Print summary
    num_instances_per_frame = [len(cells) for cells in cells_by_frame]
    print(f"  Total instances across all frames: {sum(num_instances_per_frame)}")
    print(f"  Average instances per frame: {np.mean(num_instances_per_frame):.1f}")
    
    # Step 2: Track cells
    print("\nTracking cells with Hungarian matching...")
    tracks, diagnostics = track_cells(
        cells_by_frame,
        T_early=T_early,
        max_dist=25.0,
        min_iou=0.2,
        allow_gap=False
    )
    print(f"  Total tracks found: {len(tracks)}")
    
    # Print matching diagnostics
    if diagnostics['all_matched_ious']:
        print(f"\nMatching quality diagnostics:")
        print(f"  Total matches: {len(diagnostics['all_matched_ious'])}")
        print(f"  IoU - median: {np.median(diagnostics['all_matched_ious']):.3f}, "
              f"min: {np.min(diagnostics['all_matched_ious']):.3f}, "
              f"max: {np.max(diagnostics['all_matched_ious']):.3f}")
        print(f"  Distance - median: {np.median(diagnostics['all_matched_dists']):.2f}, "
              f"min: {np.min(diagnostics['all_matched_dists']):.2f}, "
              f"max: {np.max(diagnostics['all_matched_dists']):.2f}")
        print(f"  Avg matches per frame: {np.mean(diagnostics['matches_per_frame']):.1f}")
    
    # Step 3: Filter tracks
    print("\nFiltering tracks (min_len=20, max_jump_ratio=0.30)...")
    filtered_tracks = filter_tracks(
        tracks,
        min_len=20,
        key="major_axis_length",
        max_jump_ratio=0.30
    )
    print(f"  Tracks after filtering: {len(filtered_tracks)}")
    
    # Step 4: Summary statistics
    if len(filtered_tracks) > 0:
        track_lengths = [len(track['frames']) for track in filtered_tracks]
        print(f"\nTrack statistics:")
        print(f"  Average track length: {np.mean(track_lengths):.1f} frames")
        print(f"  Min track length: {np.min(track_lengths)} frames")
        print(f"  Max track length: {np.max(track_lengths)} frames")
        print(f"  Median track length: {np.median(track_lengths):.1f} frames")
    else:
        print("\nNo tracks passed filtering criteria.")
    
    print("\n" + "=" * 60)
    
    return filtered_tracks


if __name__ == "__main__":
    # Example with synthetic data for testing
    print("Testing tracking with synthetic data...")
    
    # Create a simple synthetic example: single moving cell
    H, W = 100, 100
    T = 30
    masks = np.zeros((T, H, W), dtype=np.int32)
    
    # Cell moves diagonally
    for t in range(T):
        y = 20 + t
        x = 20 + t
        # Draw a small circle
        rr, cc = np.ogrid[:H, :W]
        circle = ((rr - y)**2 + (cc - x)**2) <= 25
        masks[t][circle] = 1
    
    # Run demo
    tracks = demo_tracking(masks, T_early=T)
    
    print(f"\nSynthetic test result:")
    print(f"  Expected: 1 track of length ~{T}")
    print(f"  Got: {len(tracks)} tracks")
    if len(tracks) > 0:
        print(f"  First track length: {len(tracks[0]['frames'])}")
