"""
Visualize tracking results: overlay tracked cells on mask images.
Supports exporting to video or viewing frame by frame.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib import colors as mcolors
from imageio import imread
from pathlib import Path
from tracking import extract_instances, track_cells
import pickle


def generate_colors(n_tracks):
    """
    Generate n distinct colors for visualization.
    """
    # Use tab20 and tab20b/c colormaps
    if n_tracks <= 20:
        cmap = plt.cm.tab20
        return [cmap(i / 20) for i in range(n_tracks)]
    else:
        # Mix multiple colormaps for more colors
        colors = []
        cmap1 = plt.cm.tab20
        cmap2 = plt.cm.tab20b
        cmap3 = plt.cm.tab20c
        
        for i in range(n_tracks):
            if i < 20:
                colors.append(cmap1(i / 20))
            elif i < 40:
                colors.append(cmap2((i - 20) / 20))
            else:
                colors.append(cmap3((i - 40) / 20))
        
        return colors


def visualize_tracks_on_frames(
    masks,
    tracks,
    start_frame=0,
    end_frame=30,
    output_dir=None,
    show_id=True,
    show_bbox=True,
    show_centroid=True,
    alpha_mask=0.3,
    save_individual_frames=True
):
    """
    Overlay tracking results on mask images.
    
    Parameters
    ----------
    masks : np.ndarray
        3D array (T, H, W) of label images
    tracks : list
        Tracking results
    start_frame : int
        Start frame
    end_frame : int
        End frame
    output_dir : Path or str
        Output directory (if saving)
    show_id : bool
        Whether to show track ID
    show_bbox : bool
        Whether to show bounding box
    show_centroid : bool
        Whether to show centroid
    alpha_mask : float
        Mask transparency
    save_individual_frames : bool
        Whether to save individual frame images
    """
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(exist_ok=True, parents=True)
    
    # Filter tracks visible in display range
    visible_tracks = []
    for track in tracks:
        track_frames = track['frames']
        # Include track if it appears in the display range
        if any(start_frame <= f <= end_frame for f in track_frames):
            visible_tracks.append(track)
    
    print(f"Tracks in display range: {len(visible_tracks)}")
    
    # Assign colors to each track
    track_colors = generate_colors(len(visible_tracks))
    
    # Build frame -> track mapping
    # frame_to_tracks[t] = [(track_idx, cell_dict), ...]
    frame_to_tracks = {t: [] for t in range(start_frame, end_frame + 1)}
    
    for track_idx, track in enumerate(visible_tracks):
        for frame, cell in zip(track['frames'], track['cells']):
            if start_frame <= frame <= end_frame:
                frame_to_tracks[frame].append((track_idx, cell))
    
    # Create figure
    fig, ax = plt.subplots(1, 1, figsize=(12, 10))
    
    # Process each frame
    for t in range(start_frame, end_frame + 1):
        ax.clear()
        
        # Display mask image as background
        mask_img = masks[t]
        
        # Create RGB image for display
        H, W = mask_img.shape
        rgb_img = np.zeros((H, W, 3), dtype=np.float32)
        
        # Gray background for cells
        rgb_img[mask_img > 0] = [0.2, 0.2, 0.2]
        
        # Overlay tracked cells with track colors
        for track_idx, cell in frame_to_tracks[t]:
            label = cell['label']
            color = track_colors[track_idx][:3]  # RGB part
            
            # Color this cell
            rgb_img[mask_img == label] = color
        
        ax.imshow(rgb_img)
        
        # Draw annotations
        for track_idx, cell in frame_to_tracks[t]:
            color = track_colors[track_idx]
            centroid = cell['centroid']
            bbox = cell['bbox']  # (minr, minc, maxr, maxc)
            
            # Draw bounding box
            if show_bbox:
                minr, minc, maxr, maxc = bbox
                width = maxc - minc
                height = maxr - minr
                rect = Rectangle(
                    (minc, minr), width, height,
                    linewidth=1.5,
                    edgecolor=color,
                    facecolor='none',
                    linestyle='--'
                )
                ax.add_patch(rect)
            
            # Draw centroid
            if show_centroid:
                ax.plot(centroid[1], centroid[0], 'o', 
                       color=color, markersize=4, 
                       markeredgecolor='white', markeredgewidth=0.5)
            
            # Show track ID
            if show_id:
                ax.text(
                    centroid[1], centroid[0] - 8,
                    f'T{track_idx}',
                    color='white',
                    fontsize=8,
                    fontweight='bold',
                    ha='center',
                    va='bottom',
                    bbox=dict(boxstyle='round,pad=0.3', 
                             facecolor=color, alpha=0.7, edgecolor='none')
                )
        
        # Set title and info
        n_cells_in_frame = len(frame_to_tracks[t])
        ax.set_title(
            f'Frame {t} | Tracked Cells: {n_cells_in_frame} / {len(visible_tracks)} tracks',
            fontsize=14, fontweight='bold'
        )
        ax.axis('off')
        
        # Save or display
        if save_individual_frames and output_dir:
            output_path = output_dir / f"frame_{t:04d}.png"
            plt.savefig(output_path, dpi=100, bbox_inches='tight')
            if t == start_frame:
                print(f"Saving frame images to: {output_dir}")
        
        if t == start_frame:
            plt.pause(0.1)
    
    plt.close()
    
    return visible_tracks


def create_track_summary_figure(tracks, start_frame=0, end_frame=30, output_path=None):
    """
    Create track summary figure: show temporal distribution of all tracks.
    """
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Filter tracks
    visible_tracks = []
    for track in tracks:
        if any(start_frame <= f <= end_frame for f in track['frames']):
            visible_tracks.append(track)
    
    # Sort by start frame
    visible_tracks.sort(key=lambda t: t['frames'][0])
    
    # Generate colors
    track_colors = generate_colors(len(visible_tracks))
    
    # Draw timeline for each track
    for i, track in enumerate(visible_tracks):
        frames = track['frames']
        # Only show frames in range
        frames_in_range = [f for f in frames if start_frame <= f <= end_frame]
        
        if not frames_in_range:
            continue
        
        color = track_colors[i]
        
        # Draw continuous timeline
        for j in range(len(frames_in_range) - 1):
            f1, f2 = frames_in_range[j], frames_in_range[j + 1]
            if f2 - f1 == 1:
                # Continuous
                ax.plot([f1, f2], [i, i], '-', color=color, linewidth=3, alpha=0.8)
            else:
                # Has gap
                ax.plot([f1, f2], [i, i], ':', color=color, linewidth=2, alpha=0.5)
        
        # Mark start and end points
        ax.plot(frames_in_range[0], i, 'o', color=color, markersize=6, 
               markeredgecolor='black', markeredgewidth=0.5)
        ax.plot(frames_in_range[-1], i, 's', color=color, markersize=6,
               markeredgecolor='black', markeredgewidth=0.5)
        
        # Show track ID
        ax.text(start_frame - 1, i, f'T{i}', 
               fontsize=7, va='center', ha='right', color='gray')
    
    ax.set_xlabel('Frame', fontsize=12, fontweight='bold')
    ax.set_ylabel('Track ID', fontsize=12, fontweight='bold')
    ax.set_title(f'Track Timeline (Frames {start_frame}-{end_frame})', 
                fontsize=14, fontweight='bold')
    ax.set_xlim(start_frame - 2, end_frame + 1)
    ax.set_ylim(-1, len(visible_tracks))
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"Track summary saved to: {output_path}")
    
    plt.tight_layout()
    plt.show()
    
    return fig


def visualize_position_tracking(
    position_name,
    mask_dir,
    start_frame=0,
    end_frame=30,
    output_base_dir=None,
    tracking_params=None
):
    """
    Complete visualization pipeline: run tracking on a single position and visualize.
    
    Parameters
    ----------
    position_name : str
        Position name
    mask_dir : Path
        Mask directory
    start_frame : int
        Start frame
    end_frame : int
        End frame
    output_base_dir : str
        Output base directory (default: ./track_visualization)
    tracking_params : dict
        Tracking parameters
    """
    # Use script directory as default
    script_dir = Path(__file__).parent
    if output_base_dir is None:
        output_base_dir = script_dir / "track_visualization"
    
    print("="*60)
    print(f"Visualizing Tracking Results: {position_name}")
    print("="*60)
    
    # Default parameters
    if tracking_params is None:
        tracking_params = {
            'T_early': max(end_frame + 1, 60),
            'max_dist': 25.0,
            'min_iou': 0.2,
            'allow_gap': False
        }
    
    # Load masks
    print(f"\nLoading masks (frames 0-{tracking_params['T_early']-1})...")
    masks = []
    for t in range(tracking_params['T_early']):
        mask_path = mask_dir / f"MASK_img_{t:09d}.tif"
        if mask_path.exists():
            masks.append(imread(mask_path))
        else:
            print(f"Warning: frame {t} does not exist")
            break
    
    masks = np.array(masks)
    print(f"  Loaded {len(masks)} frames, shape: {masks.shape}")
    
    # Extract instances
    print("\nExtracting instances...")
    cells_by_frame = []
    for t in range(len(masks)):
        cells = extract_instances(masks[t])
        cells_by_frame.append(cells)
    
    total_instances = sum(len(cells) for cells in cells_by_frame)
    print(f"  Total instances: {total_instances}")
    
    # Run tracking
    print("\nRunning tracking...")
    tracks, diagnostics = track_cells(cells_by_frame, **tracking_params)
    print(f"  Raw tracks: {len(tracks)}")
    
    # Show diagnostic info
    if diagnostics['all_matched_ious']:
        print(f"\nMatching quality:")
        print(f"  IoU median: {np.median(diagnostics['all_matched_ious']):.3f}")
        print(f"  Distance median: {np.median(diagnostics['all_matched_dists']):.2f} pixels")
    
    # Create output directory
    output_dir = Path(output_base_dir) / position_name
    
    # Visualize
    print(f"\nGenerating visualization images...")
    visible_tracks = visualize_tracks_on_frames(
        masks,
        tracks,
        start_frame=start_frame,
        end_frame=end_frame,
        output_dir=output_dir,
        show_id=True,
        show_bbox=True,
        show_centroid=True,
        save_individual_frames=True
    )
    
    # Create summary figure
    print("\nGenerating track timeline summary...")
    summary_path = output_dir / f"track_summary_{start_frame}_{end_frame}.png"
    create_track_summary_figure(
        visible_tracks,
        start_frame=start_frame,
        end_frame=end_frame,
        output_path=summary_path
    )
    
    # Statistics
    print("\n" + "="*60)
    print("Statistics:")
    print("="*60)
    print(f"Display range: frames {start_frame}-{end_frame}")
    print(f"Tracks in range: {len(visible_tracks)}")
    
    track_lengths = [len(t['frames']) for t in visible_tracks]
    if track_lengths:
        print(f"Track lengths: mean={np.mean(track_lengths):.1f}, "
              f"median={np.median(track_lengths):.1f}, "
              f"range=[{np.min(track_lengths)}, {np.max(track_lengths)}]")
    
    print(f"\nOutput directory: {output_dir}")
    print("="*60)
    
    return tracks, visible_tracks


if __name__ == "__main__":
    # Get script directory for relative paths
    script_dir = Path(__file__).parent
    
    # Example 1: REF group Pos102
    print("\n" + "#"*60)
    print("# Example: Visualize tracking results for REF Pos102")
    print("#"*60)
    
    position = "Pos102"
    mask_dir = script_dir / "data" / "REF" / "1" / "HR_REF_masks" / position / "PreprocessedPhaseMasks"
    
    if mask_dir.exists():
        tracks, visible_tracks = visualize_position_tracking(
            position_name=position,
            mask_dir=mask_dir,
            start_frame=0,
            end_frame=30,
            output_base_dir=script_dir / "track_visualization" / "REF",
            tracking_params={
                'T_early': 60,
                'max_dist': 25.0,
                'min_iou': 0.2,
                'allow_gap': False
            }
        )
    else:
        print(f"Error: directory does not exist {mask_dir}")
    
    # Example 2: RIF group Pos201
    print("\n" + "#"*60)
    print("# Example: Visualize tracking results for RIF Pos201")
    print("#"*60)
    
    position = "Pos201"
    mask_dir = script_dir / "data" / "RIF" / "1" / "HR_RIF10_masks" / position / "PreprocessedPhaseMasks"
    
    if mask_dir.exists():
        tracks_rif, visible_tracks_rif = visualize_position_tracking(
            position_name=position,
            mask_dir=mask_dir,
            start_frame=0,
            end_frame=30,
            output_base_dir=script_dir / "track_visualization" / "RIF",
            tracking_params={
                'T_early': 60,
                'max_dist': 25.0,
                'min_iou': 0.2,
                'allow_gap': False
            }
        )
    else:
        print(f"Error: directory does not exist {mask_dir}")
   