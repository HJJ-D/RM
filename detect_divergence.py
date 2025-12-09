"""
Bacterial Growth Divergence Detection
======================================
Detects when RIF-treated bacteria growth curves diverge from untreated controls.

Based on Lecture 3 requirements:
- Growth rate from sliding window (30 min / 15 frames at 2 min intervals)
- Detect earliest time where treated differs from control (p < 0.05)
- Use Omnipose segmentation masks provided
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # 使用非交互式后端，避免Qt冲突
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats
from imageio import imread
from scipy.optimize import curve_fit
import pickle
import warnings
warnings.filterwarnings('ignore')


class GrowthAnalyzer:
    """Analyze bacterial growth from segmentation masks."""
    
    def __init__(self, time_interval=2, window_minutes=30, pixel_size=0.0733):
        """
        Parameters:
        -----------
        time_interval : float
            Minutes between frames (default: 2)
        window_minutes : float
            Sliding window size in minutes (default: 30)
        pixel_size : float
            Pixel size in micrometers (default: 0.0733 μm at 150x)
        """
        self.time_interval = time_interval
        self.window_size = int(window_minutes / time_interval)  # Convert to frames
        self.pixel_size = pixel_size
        
    def load_areas_from_pickle_or_masks(self, folder_path):
        """Load area data from pickle file (preferred) or calculate from masks."""
        folder = Path(folder_path)
        pickle_file = folder / "growth_areas.pickle"
        
        # Try to load from pickle first
        if pickle_file.exists():
            with open(pickle_file, 'rb') as f:
                areas = pickle.load(f)
            return np.array(areas)
        
        # Fallback: calculate from mask files
        return self.load_masks_from_folder(folder)
    
    def load_masks_from_folder(self, folder_path):
        """Calculate total area from mask files (fallback method)."""
        folder = Path(folder_path)
        mask_files = sorted(folder.glob("MASK_*.tif"))
        
        areas = []
        for mask_file in mask_files:
            mask = imread(mask_file)
            # Convert to binary (0/1) in case mask has other values
            mask = (mask > 0).astype(np.uint8)
            # Total area = count of bacterial pixels (pixels with value 1)
            total_area = np.sum(mask)
            areas.append(total_area)
        
        return np.array(areas)
    
    def load_all_positions(self, base_folder, positions):
        """Load data from multiple positions (replicates)."""
        all_areas = []
        for pos in positions:
            pos_folder = Path(base_folder) / pos / "PreprocessedPhaseMasks"  
            if pos_folder.exists():
                # Use pickle if available, otherwise calculate from masks
                areas = self.load_areas_from_pickle_or_masks(pos_folder)
                all_areas.append(areas)
                print(f"Loaded {pos}: {len(areas)} timepoints")
        
        return all_areas
    
    def smooth_sliding_window(self, data, window_size=None):
        """Apply sliding window smoothing."""
        if window_size is None:
            window_size = self.window_size
        
        smoothed = np.convolve(data, np.ones(window_size)/window_size, mode='valid')
        return smoothed
    
    def calculate_growth_rate_derivative(self, areas):
        """Calculate growth rate as derivative of smoothed area."""
        # Smooth the area curve
        smoothed = self.smooth_sliding_window(areas)
        
        # Calculate derivative (change per frame)
        growth_rate = np.gradient(smoothed)
        
        return growth_rate
    
    def calculate_growth_rate_exponential(self, areas, window_size=None):
        """Calculate growth rate by fitting exponential in sliding windows."""
        if window_size is None:
            window_size = self.window_size
        
        growth_rates = []
        
        for i in range(len(areas) - window_size + 1):
            window_data = areas[i:i+window_size]
            time_points = np.arange(len(window_data))
            
            # Fit exponential: y = a * e^(b*t), where b is the growth rate
            # Use log transform: log(y) = log(a) + b*t
            try:
                log_data = np.log(window_data + 1)  # +1 to avoid log(0)
                # Linear fit in log space
                coeffs = np.polyfit(time_points, log_data, 1)
                growth_rate = coeffs[0]  # Slope = growth rate
                growth_rates.append(growth_rate)
            except:
                growth_rates.append(0)
        
        return np.array(growth_rates)
    
    def normalize_to_control(self, treated_data, control_data):
        """Normalize treated data by control mean."""
        control_mean = np.mean(control_data, axis=0)
        normalized = treated_data / (control_mean + 1e-10)
        return normalized, control_mean
    
    def detect_divergence_ttest(self, treated_growth_rates, control_growth_rates, 
                                 alpha=0.05, min_consecutive=3):
        """
        Detect first timepoint where treated significantly differs from control.
        
        Parameters:
        -----------
        treated_growth_rates : list of arrays
            Growth rates for each treated replicate
        control_growth_rates : list of arrays
            Growth rates for each control replicate
        alpha : float
            Significance level (default: 0.05)
        min_consecutive : int
            Minimum consecutive significant timepoints required
        
        Returns:
        --------
        divergence_time : int or None
            First timepoint of divergence (frame number)
        p_values : array
            P-values at each timepoint
        """
        # Find minimum length across all replicates
        min_length = min([len(gr) for gr in treated_growth_rates + control_growth_rates])
        
        p_values = []
        
        for t in range(min_length):
            # Get growth rates at this timepoint from all replicates
            treated_at_t = [gr[t] for gr in treated_growth_rates]
            control_at_t = [gr[t] for gr in control_growth_rates]
            
            # Independent t-test
            t_stat, p_val = stats.ttest_ind(treated_at_t, control_at_t)
            p_values.append(p_val)
        
        p_values = np.array(p_values)
        
        # Find first significant point with consecutive significance
        divergence_time = None
        for i in range(len(p_values) - min_consecutive + 1):
            if np.all(p_values[i:i+min_consecutive] < alpha):
                divergence_time = i
                break
        
        return divergence_time, p_values
    
    def get_time_in_minutes(self, frame_idx):
        """Convert frame index to time in minutes."""
        return frame_idx * self.time_interval


def main():
    """Main analysis pipeline."""
    
    print("=" * 70)
    print("Bacterial Growth Divergence Detection")
    print("=" * 70)
    print()
    
    # Configuration
    REF_BASE = Path("F:/RM/REF_masks101_110")
    RIF_BASE = Path("F:/RM/RIF10_masks201_210")
    
    # Positions to analyze (replicates)
    ref_positions = [f"Pos{i}" for i in range(101, 111)]  # Pos101-110
    rif_positions = [f"Pos{i}" for i in range(201, 211)]  # Pos201-210
    
    # Initialize analyzer
    analyzer = GrowthAnalyzer(time_interval=2, window_minutes=30)
    
    print("Step 1: Loading REF (untreated control) data...")
    print("-" * 70)
    ref_areas = analyzer.load_all_positions(REF_BASE, ref_positions)
    print()
    
    print("Step 2: Loading RIF10 (rifampicin treated) data...")
    print("-" * 70)
    rif_areas = analyzer.load_all_positions(RIF_BASE, rif_positions)
    print()
    
    print("Step 3: Calculating growth rates with sliding window...")
    print("-" * 70)
    
    # Calculate growth rates for each replicate
    ref_growth_rates = []
    for i, areas in enumerate(ref_areas):
        gr = analyzer.calculate_growth_rate_exponential(areas)
        ref_growth_rates.append(gr)
        print(f"REF replicate {i+1}: {len(gr)} timepoints")
    
    rif_growth_rates = []
    for i, areas in enumerate(rif_areas):
        gr = analyzer.calculate_growth_rate_exponential(areas)
        rif_growth_rates.append(gr)
        print(f"RIF replicate {i+1}: {len(gr)} timepoints")
    
    print()
    
    print("Step 4: Detecting divergence (t-test, p < 0.05)...")
    print("-" * 70)
    
    # Detect divergence
    divergence_frame, p_values = analyzer.detect_divergence_ttest(
        rif_growth_rates, 
        ref_growth_rates,
        alpha=0.05,
        min_consecutive=3
    )
    
    if divergence_frame is not None:
        # Account for window offset
        actual_frame = divergence_frame + analyzer.window_size - 1
        divergence_minutes = analyzer.get_time_in_minutes(actual_frame)
        print(f"✓ DIVERGENCE DETECTED at frame {actual_frame}")
        print(f"  Time: {divergence_minutes:.1f} minutes ({divergence_minutes/60:.2f} hours)")
        print(f"  P-value: {p_values[divergence_frame]:.4f}")
    else:
        print("✗ No significant divergence detected")
        divergence_minutes = None
        actual_frame = None
    
    print()
    
    print("Step 5: Creating visualizations...")
    print("-" * 70)
    
    # Create comprehensive figure
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    # Plot 1: Total Area over time
    ax1 = axes[0]
    min_length = min([len(a) for a in ref_areas + rif_areas])
    times = np.arange(min_length) * analyzer.time_interval / 60  # Convert to hours
    
    # Convert to micrometers squared
    pixel_to_um2 = analyzer.pixel_size ** 2
    
    for i, areas in enumerate(ref_areas):
        ax1.plot(times, areas[:min_length] * pixel_to_um2, 'b-', alpha=0.3, linewidth=1)
    for i, areas in enumerate(rif_areas):
        ax1.plot(times, areas[:min_length] * pixel_to_um2, 'r-', alpha=0.3, linewidth=1)
    
    # Plot means
    ref_mean = np.mean([a[:min_length] for a in ref_areas], axis=0) * pixel_to_um2
    rif_mean = np.mean([a[:min_length] for a in rif_areas], axis=0) * pixel_to_um2
    ax1.plot(times, ref_mean, 'b-', linewidth=2, label='REF (Untreated)')
    ax1.plot(times, rif_mean, 'r-', linewidth=2, label='RIF10 (Treated)')
    
    if actual_frame is not None and actual_frame < min_length:
        ax1.axvline(x=divergence_minutes/60, color='green', linestyle='--', 
                   linewidth=2, label=f'Divergence ({divergence_minutes:.1f} min)')
    
    ax1.set_xlabel('Time (hours)', fontsize=12)
    ax1.set_ylabel(r'Total Bacterial Area ($\mu m^2$)', fontsize=12)
    ax1.set_title('Bacterial Growth: Total Area Over Time', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Growth Rate over time
    ax2 = axes[1]
    min_gr_length = min([len(gr) for gr in ref_growth_rates + rif_growth_rates])
    gr_times = (np.arange(min_gr_length) + analyzer.window_size - 1) * analyzer.time_interval / 60
    
    for i, gr in enumerate(ref_growth_rates):
        ax2.plot(gr_times, gr[:min_gr_length], 'b-', alpha=0.3, linewidth=1)
    for i, gr in enumerate(rif_growth_rates):
        ax2.plot(gr_times, gr[:min_gr_length], 'r-', alpha=0.3, linewidth=1)
    
    # Plot means
    ref_gr_mean = np.mean([gr[:min_gr_length] for gr in ref_growth_rates], axis=0)
    rif_gr_mean = np.mean([gr[:min_gr_length] for gr in rif_growth_rates], axis=0)
    ax2.plot(gr_times, ref_gr_mean, 'b-', linewidth=2, label='REF (Untreated)')
    ax2.plot(gr_times, rif_gr_mean, 'r-', linewidth=2, label='RIF10 (Treated)')
    
    if actual_frame is not None and actual_frame < min_length:
        ax2.axvline(x=divergence_minutes/60, color='green', linestyle='--', 
                   linewidth=2, label=f'Divergence ({divergence_minutes:.1f} min)')
    
    ax2.set_xlabel('Time (hours)', fontsize=12)
    ax2.set_ylabel('Growth Rate (exponential fit)', fontsize=12)
    ax2.set_title('Growth Rate Comparison (30-min sliding window)', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: P-values over time
    ax3 = axes[2]
    ax3.plot(gr_times[:len(p_values)], p_values, 'k-', linewidth=2, label='P-value')
    ax3.axhline(y=0.05, color='orange', linestyle='--', linewidth=2, label='α = 0.05')
    
    if divergence_frame is not None:
        ax3.axvline(x=divergence_minutes/60, color='green', linestyle='--', 
                   linewidth=2, label=f'Divergence ({divergence_minutes:.1f} min)')
    
    ax3.set_xlabel('Time (hours)', fontsize=12)
    ax3.set_ylabel('P-value (t-test)', fontsize=12)
    ax3.set_title('Statistical Significance: Treated vs Control', fontsize=14, fontweight='bold')
    ax3.set_yscale('log')
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('F:/RM/divergence_analysis.png', dpi=300, bbox_inches='tight')
    print("✓ Saved: divergence_analysis.png")
    plt.close()
    
    # Save results
    results = {
        'divergence_frame': actual_frame,
        'divergence_minutes': divergence_minutes,
        'divergence_hours': divergence_minutes / 60 if divergence_minutes else None,
        'p_values': p_values,
        'ref_areas': ref_areas,
        'rif_areas': rif_areas,
        'ref_growth_rates': ref_growth_rates,
        'rif_growth_rates': rif_growth_rates,
        'parameters': {
            'time_interval': analyzer.time_interval,
            'window_size': analyzer.window_size,
            'window_minutes': analyzer.window_size * analyzer.time_interval,
            'alpha': 0.05,
            'min_consecutive': 3
        }
    }
    
    with open('F:/RM/divergence_results.pkl', 'wb') as f:
        pickle.dump(results, f)
    print("✓ Saved: divergence_results.pkl")
    
    # Create summary report
    with open('F:/RM/divergence_report.txt', 'w', encoding='utf-8') as f:
        f.write("=" * 70 + "\n")
        f.write("BACTERIAL GROWTH DIVERGENCE ANALYSIS REPORT\n")
        f.write("=" * 70 + "\n\n")
        
        f.write("DATASET:\n")
        f.write(f"  - Control (REF): {len(ref_areas)} replicates (Pos101-110)\n")
        f.write(f"  - Treated (RIF10): {len(rif_areas)} replicates (Pos201-210)\n")
        f.write(f"  - Imaging interval: {analyzer.time_interval} minutes\n\n")
        
        f.write("METHODOLOGY:\n")
        f.write(f"  - Sliding window: {analyzer.window_size * analyzer.time_interval} minutes\n")
        f.write(f"    ({analyzer.window_size} frames)\n")
        f.write(f"  - Growth rate: Exponential fit (log-linear)\n")
        f.write(f"  - Statistical test: Independent t-test\n")
        f.write(f"  - Significance level: alpha = 0.05\n")
        f.write(f"  - Minimum consecutive significant points: 3\n\n")
        
        f.write("RESULTS:\n")
        if divergence_frame is not None:
            f.write(f"  [+] DIVERGENCE DETECTED\n")
            f.write(f"  - Frame: {actual_frame}\n")
            f.write(f"  - Time: {divergence_minutes:.1f} minutes\n")
            f.write(f"  - Time: {divergence_minutes/60:.2f} hours\n")
            f.write(f"  - P-value at divergence: {p_values[divergence_frame]:.6f}\n")
        else:
            f.write(f"  [X] NO SIGNIFICANT DIVERGENCE DETECTED\n")
        
        f.write("\n" + "=" * 70 + "\n")
    
    print("✓ Saved: divergence_report.txt")
    print()
    
    print("=" * 70)
    print("ANALYSIS COMPLETE!")
    print("=" * 70)
    print()
    
    if divergence_frame is not None:
        print(f"[ANSWER] Treated and control diverge at {divergence_minutes:.1f} minutes")
        print(f"         ({divergence_minutes/60:.2f} hours after start)")
    else:
        print("[WARNING] No significant divergence detected in this dataset")
    
    print()
    print("Output files:")
    print("  - divergence_analysis.png  (visualization)")
    print("  - divergence_results.pkl   (full results data)")
    print("  - divergence_report.txt    (summary report)")
    print()


if __name__ == "__main__":
    main()
