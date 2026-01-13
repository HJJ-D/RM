"""
Run single-cell tracking analysis.
Track REF (control, no drug) and RIF (treatment, with drug) data.
"""

import numpy as np
from imageio import imread
from pathlib import Path
from tracking import extract_instances, track_cells, filter_tracks
import pickle


def load_masks_for_position(mask_dir: Path, num_frames: int = 120) -> np.ndarray:
    """
    Load all mask files for a position.
    
    Parameters
    ----------
    mask_dir : Path
        Directory containing MASK_img_*.tif files
    num_frames : int
        Number of frames to load (default 120)
        
    Returns
    -------
    np.ndarray
        3D array of shape (T, H, W)
    """
    # Read first frame to get image dimensions
    first_mask_path = mask_dir / f"MASK_img_{0:09d}.tif"
    if not first_mask_path.exists():
        raise FileNotFoundError(f"Mask file not found: {first_mask_path}")
    
    first_mask = imread(first_mask_path)
    H, W = first_mask.shape
    
    # Preallocate array
    masks = np.zeros((num_frames, H, W), dtype=first_mask.dtype)
    masks[0] = first_mask
    
    # Load remaining frames
    print(f"  Loading {num_frames} mask images...")
    for t in range(1, num_frames):
        mask_path = mask_dir / f"MASK_img_{t:09d}.tif"
        if mask_path.exists():
            masks[t] = imread(mask_path)
        else:
            print(f"  Warning: frame {t} file does not exist, using blank mask")
    
    return masks


def run_tracking_for_position(
    position_name: str,
    mask_dir: Path,
    T_early: int = 60,
    max_dist: float = 25.0,
    min_iou: float = 0.10,
    allow_gap: bool = True,
    min_track_len: int = 20,
    max_jump_ratio: float = 0.30
):
    """
    Run tracking for a single position.
    
    Parameters
    ----------
    position_name : str
        Position name (e.g., "Pos101")
    mask_dir : Path
        Directory containing PreprocessedPhaseMasks
    T_early : int
        Number of early frames to use
    max_dist : float
        Maximum centroid distance
    min_iou : float
        Minimum IoU threshold
    allow_gap : bool
        Whether to allow gap closing
    min_track_len : int
        Minimum track length
    max_jump_ratio : float
        Maximum morphology change ratio
        
    Returns
    -------
    dict
        Dictionary containing tracking results
    """
    print(f"\n{'='*60}")
    print(f"Processing {position_name}")
    print(f"{'='*60}")
    
    # Load masks
    masks = load_masks_for_position(mask_dir, num_frames=120)
    print(f"  Mask shape: {masks.shape}")
    print(f"  Data type: {masks.dtype}")
    
    # Check instance count
    unique_labels = np.unique(masks[0])
    num_instances_frame0 = len(unique_labels) - 1  # Subtract background
    print(f"  Frame 0 instances: {num_instances_frame0}")
    
    # Extract instances
    print(f"\nExtracting instances from first {T_early} frames...")
    cells_by_frame = []
    for t in range(T_early):
        cells = extract_instances(masks[t])
        cells_by_frame.append(cells)
        if t == 0 or t == T_early - 1:
            print(f"    Frame {t}: {len(cells)} instances")
    
    total_instances = sum(len(cells) for cells in cells_by_frame)
    avg_instances = total_instances / T_early
    print(f"  Total instances: {total_instances}")
    print(f"  Average instances per frame: {avg_instances:.1f}")
    
    # Track
    print(f"\nTracking cells (max_dist={max_dist}, min_iou={min_iou}, allow_gap={allow_gap})...")
    tracks, diagnostics = track_cells(
        cells_by_frame,
        T_early=T_early,
        max_dist=max_dist,
        min_iou=min_iou,
        allow_gap=allow_gap
    )
    print(f"  Raw tracks: {len(tracks)}")
    
    # Print matching quality diagnostics
    if diagnostics['all_matched_ious']:
        print(f"\nMatching quality diagnostics:")
        print(f"  Total matches: {len(diagnostics['all_matched_ious'])}")
        print(f"  IoU - median: {np.median(diagnostics['all_matched_ious']):.3f}, "
              f"min: {np.min(diagnostics['all_matched_ious']):.3f}, "
              f"max: {np.max(diagnostics['all_matched_ious']):.3f}")
        print(f"  Distance - median: {np.median(diagnostics['all_matched_dists']):.2f}, "
              f"min: {np.min(diagnostics['all_matched_dists']):.2f}, "
              f"max: {np.max(diagnostics['all_matched_dists']):.2f}")
        print(f"  Average matches per frame: {np.mean(diagnostics['matches_per_frame']):.1f}")
    
    # Filter
    print(f"\nFiltering tracks (min_len={min_track_len}, max_jump_ratio={max_jump_ratio})...")
    filtered_tracks = filter_tracks(
        tracks,
        min_len=min_track_len,
        key="major_axis_length",
        max_jump_ratio=max_jump_ratio
    )
    print(f"  Filtered tracks: {len(filtered_tracks)}")
    
    # Statistics
    if len(filtered_tracks) > 0:
        track_lengths = [len(track['frames']) for track in filtered_tracks]
        print(f"\nTrack statistics:")
        print(f"  Average length: {np.mean(track_lengths):.1f} frames")
        print(f"  Shortest: {np.min(track_lengths)} frames")
        print(f"  Longest: {np.max(track_lengths)} frames")
        print(f"  Median: {np.median(track_lengths):.1f} frames")
    
    return {
        'position': position_name,
        'tracks': filtered_tracks,
        'num_tracks': len(filtered_tracks),
        'params': {
            'T_early': T_early,
            'max_dist': max_dist,
            'min_iou': min_iou,
            'allow_gap': allow_gap,
            'min_track_len': min_track_len,
            'max_jump_ratio': max_jump_ratio
        }
    }


def run_tracking_for_group(
    group_name: str,
    base_dir: Path,
    positions: list,
    output_file: str = None,
    **tracking_params
):
    """
    Run tracking for a group of positions.
    
    Parameters
    ----------
    group_name : str
        Group name ("REF" or "RIF")
    base_dir : Path
        Base data directory
    positions : list
        List of positions (e.g., ["Pos101", "Pos102", ...])
    output_file : str
        Output pickle filename (optional)
    **tracking_params
        Parameters passed to run_tracking_for_position
        
    Returns
    -------
    dict
        Dictionary containing results for all positions
    """
    print(f"\n{'#'*60}")
    print(f"# Processing group: {group_name}")
    print(f"# Positions: {len(positions)}")
    print(f"{'#'*60}")
    
    results = {}
    
    for pos in positions:
        mask_dir = base_dir / pos / "PreprocessedPhaseMasks"
        
        if not mask_dir.exists():
            print(f"\nWarning: mask directory for {pos} does not exist, skipping")
            continue
        
        try:
            result = run_tracking_for_position(pos, mask_dir, **tracking_params)
            results[pos] = result
        except Exception as e:
            print(f"\nError: failed to process {pos}: {e}")
            continue
    
    # Summary statistics
    print(f"\n{'='*60}")
    print(f"Group {group_name} Summary")
    print(f"{'='*60}")
    
    total_tracks = sum(r['num_tracks'] for r in results.values())
    print(f"Total positions: {len(results)}")
    print(f"Total tracks: {total_tracks}")
    
    if len(results) > 0:
        avg_tracks_per_pos = total_tracks / len(results)
        print(f"Average tracks per position: {avg_tracks_per_pos:.1f}")
        
        # Track length distribution for all tracks
        all_track_lengths = []
        for result in results.values():
            for track in result['tracks']:
                all_track_lengths.append(len(track['frames']))
        
        if len(all_track_lengths) > 0:
            print(f"\nTrack length statistics for all tracks:")
            print(f"  Average: {np.mean(all_track_lengths):.1f} frames")
            print(f"  Median: {np.median(all_track_lengths):.1f} frames")
            print(f"  Shortest: {np.min(all_track_lengths)} frames")
            print(f"  Longest: {np.max(all_track_lengths)} frames")
    
    # Save results
    if output_file:
        output_path = Path(output_file)
        with open(output_path, 'wb') as f:
            pickle.dump(results, f)
        print(f"\nResults saved to: {output_path}")
    
    return results


if __name__ == "__main__":
    # Data directory (relative to script location)
    script_dir = Path(__file__).parent
    data_root = script_dir / "data"
    
    # REF group (control, no drug)
    ref_base = data_root / "REF" / "1" / "HR_REF_masks"
    ref_positions = [f"Pos{i}" for i in range(101, 111)]  # Pos101-Pos110
    
    # RIF group (treatment, with drug)
    rif_base = data_root / "RIF" / "1" / "HR_RIF10_masks"
    rif_positions = [f"Pos{i}" for i in range(201, 211)]  # Pos201-Pos210
    
    # Tracking parameters
    tracking_params = {
        'T_early': 60,           # Use only first 60 frames (early frames are cleaner)
        'max_dist': 25.0,        # Maximum centroid distance (pixels)
        'min_iou': 0.2,          # Minimum IoU threshold (raised to 0.2 for Hungarian method)
        'allow_gap': False,      # Disable gap closing by default (Hungarian method is more stable)
        'min_track_len': 20,     # Minimum track length
        'max_jump_ratio': 0.30   # Maximum morphology change ratio
    }
    
    # Run REF group
    print("\n" + "="*60)
    print("Processing REF group (control, no drug)")
    print("="*60)
    ref_results = run_tracking_for_group(
        group_name="REF",
        base_dir=ref_base,
        positions=ref_positions,
        output_file=script_dir / "tracking_results_REF.pickle",
        **tracking_params
    )
    
    # Run RIF group
    print("\n" + "="*60)
    print("Processing RIF group (treatment, with drug)")
    print("="*60)
    rif_results = run_tracking_for_group(
        group_name="RIF",
        base_dir=rif_base,
        positions=rif_positions,
        output_file=script_dir / "tracking_results_RIF.pickle",
        **tracking_params
    )
    
    # Final comparison
    print("\n" + "#"*60)
    print("# Final Comparison: REF vs RIF")
    print("#"*60)
    
    ref_total = sum(r['num_tracks'] for r in ref_results.values())
    rif_total = sum(r['num_tracks'] for r in rif_results.values())
    
    print(f"\nREF group (control, no drug):")
    print(f"  Positions: {len(ref_results)}")
    print(f"  Total tracks: {ref_total}")
    if len(ref_results) > 0:
        print(f"  Average tracks/position: {ref_total/len(ref_results):.1f}")
    
    print(f"\nRIF group (treatment, with drug):")
    print(f"  Positions: {len(rif_results)}")
    print(f"  Total tracks: {rif_total}")
    if len(rif_results) > 0:
        print(f"  Average tracks/position: {rif_total/len(rif_results):.1f}")
    
    print(f"\nComplete!")
