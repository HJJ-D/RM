"""
Complete analysis pipeline: from tracking results to growth rate analysis.

Using real REF (control) and RIF (treatment) tracking data.
"""

import pickle
from pathlib import Path
from growth_analysis import analyze_drug_response
import matplotlib.pyplot as plt
import numpy as np


def load_tracking_results():
    """
    Load tracking results.
    
    Returns
    -------
    tracks_ref : list
        All tracks from REF group (control)
    tracks_rif : list
        All tracks from RIF group (treatment)
    """
    # Check if files exist (relative to script location)
    script_dir = Path(__file__).parent
    ref_file = script_dir / "tracking_results_REF.pickle"
    rif_file = script_dir / "tracking_results_RIF.pickle"
    
    if not ref_file.exists():
        print(f"Error: {ref_file} does not exist")
        print("Please run first: python run_tracking.py")
        return None, None
    
    if not rif_file.exists():
        print(f"Error: {rif_file} does not exist")
        print("Please run first: python run_tracking.py")
        return None, None
    
    # Load data
    print("Loading tracking results...")
    with open(ref_file, 'rb') as f:
        ref_results = pickle.load(f)
    print(f"  REF group: {len(ref_results)} positions")
    
    with open(rif_file, 'rb') as f:
        rif_results = pickle.load(f)
    print(f"  RIF group: {len(rif_results)} positions")
    
    # Merge all position tracks
    tracks_ref = []
    for pos_name, pos_data in ref_results.items():
        tracks_ref.extend(pos_data['tracks'])
        print(f"    {pos_name}: {len(pos_data['tracks'])} tracks")
    
    print(f"  REF total tracks: {len(tracks_ref)}")
    
    tracks_rif = []
    for pos_name, pos_data in rif_results.items():
        tracks_rif.extend(pos_data['tracks'])
        print(f"    {pos_name}: {len(pos_data['tracks'])} tracks")
    
    print(f"  RIF total tracks: {len(tracks_rif)}")
    
    return tracks_ref, tracks_rif


def visualize_growth_analysis(results, save_path=None, output_dir=None):
    """
    Visualize growth rate analysis results.
    
    Parameters
    ----------
    results : dict
        Results returned by analyze_drug_response
    save_path : str
        Path to save combined figure (default: ./growth_analysis_results.png)
    output_dir : str
        Directory to save individual figures (default: current directory)
    """
    # Use script directory as default
    script_dir = Path(__file__).parent
    if save_path is None:
        save_path = script_dir / "growth_analysis_results.png"
    if output_dir is None:
        output_dir = script_dir
    output_dir = Path(output_dir)
    
    # Extract data
    times_ctrl = results['control']['times']
    g_ctrl = results['control']['growth_rates']
    times_treat = results['treatment']['times']
    g_treat = results['treatment']['growth_rates']
    times_norm = results['normalized']['times']
    g_norm = results['normalized']['growth_rates']
    t_star = results['normalized']['response_time']
    stats_ctrl = results['control']['statistics']
    stats_treat = results['treatment']['statistics']
    
    # ========== Fig 1: Growth Rates Time Course ==========
    fig1, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_ctrl, g_ctrl, 'b-', linewidth=2.5, label='Control (REF)', alpha=0.7)
    ax.plot(times_treat, g_treat, 'r-', linewidth=2.5, label='Treatment (RIF)', alpha=0.7)
    if t_star is not None:
        ax.axvline(t_star, color='green', linestyle='--', linewidth=2.5, 
                   label=f'Response time: frame {t_star}')
    ax.set_xlabel('Frame', fontsize=14)
    ax.set_ylabel('Growth Rate (pixels/frame)', fontsize=14)
    ax.set_title('Growth Rates: Control vs Treatment', fontsize=16, fontweight='bold')
    ax.legend(fontsize=12)
    ax.tick_params(axis='both', labelsize=12)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    fig1.savefig(f"{output_dir}/fig1_growth_timecourse.png", dpi=150, bbox_inches='tight')
    print(f"  Saved: fig1_growth_timecourse.png")
    plt.close(fig1)
    
    # ========== Fig 2: Normalized Growth Rate ==========
    fig2, ax = plt.subplots(figsize=(8, 4))
    ax.plot(times_norm, g_norm, 'purple', linewidth=2.5, label='Normalized Growth')
    ax.axhline(1.0, color='gray', linestyle='--', linewidth=1.5, label='No effect (=1.0)')
    ax.axhline(0.8, color='orange', linestyle='--', linewidth=1.5, label='20% reduction (=0.8)')
    if t_star is not None:
        ax.axvline(t_star, color='green', linestyle='--', linewidth=2.5,
                   label=f'Response at frame {t_star}')
    ax.set_xlabel('Frame', fontsize=14)
    ax.set_ylabel('Normalized Growth (Treatment/Control)', fontsize=14)
    ax.set_title('Normalized Growth Rate', fontsize=16, fontweight='bold')
    ax.legend(fontsize=12)
    ax.tick_params(axis='both', labelsize=12)
    ax.grid(alpha=0.3)
    ax.set_ylim([0, max(1.5, g_norm.max() * 1.1)])
    plt.tight_layout()
    fig2.savefig(f"{output_dir}/fig2_normalized_growth.png", dpi=150, bbox_inches='tight')
    print(f"  Saved: fig2_normalized_growth.png")
    plt.close(fig2)
    
    # ========== Fig 3: Growth Rate Distribution ==========
    fig3, ax = plt.subplots(figsize=(8, 4))
    data = [g_ctrl, g_treat]
    bp = ax.boxplot(data, labels=['Control', 'Treatment'], patch_artist=True,
                    showmeans=True, meanline=True)
    bp['boxes'][0].set_facecolor('lightblue')
    bp['boxes'][1].set_facecolor('lightcoral')
    ax.set_ylabel('Growth Rate (pixels/frame)', fontsize=14)
    ax.set_title('Growth Rate Distribution', fontsize=16, fontweight='bold')
    ax.tick_params(axis='both', labelsize=12)
    ax.grid(alpha=0.3, axis='y')
    # Add statistics info
    textstr = f"Control:\n  Median: {stats_ctrl['median']:.4f}\n  IQR: [{stats_ctrl['q25']:.4f}, {stats_ctrl['q75']:.4f}]\n\n"
    textstr += f"Treatment:\n  Median: {stats_treat['median']:.4f}\n  IQR: [{stats_treat['q25']:.4f}, {stats_treat['q75']:.4f}]"
    ax.text(0.98, 0.98, textstr, fontsize=11, transform=ax.transAxes,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    plt.tight_layout()
    fig3.savefig(f"{output_dir}/fig3_distribution.png", dpi=150, bbox_inches='tight')
    print(f"  Saved: fig3_distribution.png")
    plt.close(fig3)
    
    # ========== Fig 4: Drug-Induced Inhibition ==========
    fig4, ax = plt.subplots(figsize=(8, 4))
    inhibition = (1 - g_norm) * 100  # Convert to percentage
    ax.plot(times_norm, inhibition, 'darkgreen', linewidth=2.5)
    ax.axhline(0, color='gray', linestyle='-', linewidth=1.5)
    ax.axhline(20, color='orange', linestyle='--', linewidth=1.5, label='20% inhibition')
    if t_star is not None:
        ax.axvline(t_star, color='green', linestyle='--', linewidth=2.5,
                   label=f'Response at frame {t_star}')
    ax.fill_between(times_norm, 0, inhibition, where=(inhibition > 0), 
                     color='green', alpha=0.3, label='Growth inhibition')
    ax.set_xlabel('Frame', fontsize=14)
    ax.set_ylabel('Growth Inhibition (%)', fontsize=14)
    ax.set_title('Drug-Induced Growth Inhibition', fontsize=16, fontweight='bold')
    ax.legend(fontsize=12)
    ax.tick_params(axis='both', labelsize=12)
    ax.grid(alpha=0.3)
    plt.tight_layout()
    fig4.savefig(f"{output_dir}/fig4_inhibition.png", dpi=150, bbox_inches='tight')
    print(f"  Saved: fig4_inhibition.png")
    plt.close(fig4)
    
    # ========== Combined figure (keep original functionality) ==========
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Subplot 1: Raw growth rates comparison
    ax = axes[0, 0]
    ax.plot(times_ctrl, g_ctrl, 'b-', linewidth=2, label='Control (REF)', alpha=0.7)
    ax.plot(times_treat, g_treat, 'r-', linewidth=2, label='Treatment (RIF)', alpha=0.7)
    if t_star is not None:
        ax.axvline(t_star, color='green', linestyle='--', linewidth=2, 
                   label=f'Response time: frame {t_star}')
    ax.set_xlabel('Frame', fontsize=12)
    ax.set_ylabel('Growth Rate (pixels/frame)', fontsize=12)
    ax.set_title('Growth Rates: Control vs Treatment', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    
    # Subplot 2: Normalized growth rate
    ax = axes[0, 1]
    ax.plot(times_norm, g_norm, 'purple', linewidth=2, label='Normalized Growth')
    ax.axhline(1.0, color='gray', linestyle='--', linewidth=1, label='No effect (=1.0)')
    ax.axhline(0.8, color='orange', linestyle='--', linewidth=1, label='20% reduction (=0.8)')
    if t_star is not None:
        ax.axvline(t_star, color='green', linestyle='--', linewidth=2,
                   label=f'Response at frame {t_star}')
    ax.set_xlabel('Frame', fontsize=12)
    ax.set_ylabel('Normalized Growth (Treatment/Control)', fontsize=12)
    ax.set_title('Normalized Growth Rate', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    ax.set_ylim([0, max(1.5, g_norm.max() * 1.1)])
    
    # Subplot 3: Growth rate distribution (boxplot)
    ax = axes[1, 0]
    data = [g_ctrl, g_treat]
    bp = ax.boxplot(data, labels=['Control', 'Treatment'], patch_artist=True,
                    showmeans=True, meanline=True)
    bp['boxes'][0].set_facecolor('lightblue')
    bp['boxes'][1].set_facecolor('lightcoral')
    ax.set_ylabel('Growth Rate (pixels/frame)', fontsize=12)
    ax.set_title('Growth Rate Distribution', fontsize=13, fontweight='bold')
    ax.grid(alpha=0.3, axis='y')
    textstr = f"Control:\n  Median: {stats_ctrl['median']:.4f}\n  IQR: [{stats_ctrl['q25']:.4f}, {stats_ctrl['q75']:.4f}]\n\n"
    textstr += f"Treatment:\n  Median: {stats_treat['median']:.4f}\n  IQR: [{stats_treat['q25']:.4f}, {stats_treat['q75']:.4f}]"
    ax.text(1.5, ax.get_ylim()[1] * 0.95, textstr, fontsize=9,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # Subplot 4: Inhibition over time
    ax = axes[1, 1]
    inhibition = (1 - g_norm) * 100
    ax.plot(times_norm, inhibition, 'darkgreen', linewidth=2)
    ax.axhline(0, color='gray', linestyle='-', linewidth=1)
    ax.axhline(20, color='orange', linestyle='--', linewidth=1, label='20% inhibition')
    if t_star is not None:
        ax.axvline(t_star, color='green', linestyle='--', linewidth=2,
                   label=f'Response at frame {t_star}')
    ax.fill_between(times_norm, 0, inhibition, where=(inhibition > 0), 
                     color='green', alpha=0.3, label='Growth inhibition')
    ax.set_xlabel('Frame', fontsize=12)
    ax.set_ylabel('Growth Inhibition (%)', fontsize=12)
    ax.set_title('Drug-Induced Growth Inhibition', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"  Saved combined figure: {save_path}")
    plt.show()


def main():
    """
    Main pipeline: Load data -> Analyze -> Visualize
    """
    print("="*60)
    print("Complete Analysis Pipeline: Real Data Analysis")
    print("="*60)
    
    # 1. Load tracking results
    tracks_ref, tracks_rif = load_tracking_results()
    
    if tracks_ref is None or tracks_rif is None:
        print("\nPlease run tracking first:")
        print("  python run_tracking.py")
        return
    
    # 2. Run growth rate analysis
    print("\n" + "="*60)
    print("Using improved normalization method (v2 - smooth mode)")
    print("="*60)
    results = analyze_drug_response(
        tracks_treat=tracks_rif,  # RIF is treatment group
        tracks_ctrl=tracks_ref,   # REF is control group
        use_v2=True,              # Use improved method
        norm_mode="smooth",       # Smooth control data
        smooth_window=7           # Smoothing window size
    )
    
    # 3. Save results
    script_dir = Path(__file__).parent
    result_file = script_dir / "growth_analysis_results.pickle"
    with open(result_file, 'wb') as f:
        pickle.dump(results, f)
    print(f"\nAnalysis results saved to: {result_file}")
    
    # 4. Visualize
    print("\nGenerating visualization figures...")
    visualize_growth_analysis(results)
    
    # 5. Summary report
    print("\n" + "="*60)
    print("Analysis Summary")
    print("="*60)
    
    stats_ctrl = results['control']['statistics']
    stats_treat = results['treatment']['statistics']
    t_star = results['normalized']['response_time']
    
    print(f"\nControl group (REF - no drug):")
    print(f"  Median growth rate: {stats_ctrl['median']:.4f} pixels/frame")
    print(f"  IQR: [{stats_ctrl['q25']:.4f}, {stats_ctrl['q75']:.4f}]")
    print(f"  Sample size: {stats_ctrl['n_points']} time points")
    
    print(f"\nTreatment group (RIF - with drug):")
    print(f"  Median growth rate: {stats_treat['median']:.4f} pixels/frame")
    print(f"  IQR: [{stats_treat['q25']:.4f}, {stats_treat['q75']:.4f}]")
    print(f"  Sample size: {stats_treat['n_points']} time points")
    
    if stats_ctrl['median'] > 0:
        reduction = (1 - stats_treat['median'] / stats_ctrl['median']) * 100
        print(f"\nDrug effect:")
        print(f"  Growth inhibition: {reduction:.1f}%")
    
    if t_star is not None:
        print(f"\nDrug response time:")
        print(f"  First sustained inhibition detected: frame {t_star}")
        print(f"  (Normalized growth rate sustained below 0.8, i.e., 20% inhibition)")
    else:
        print(f"\nDrug response time:")
        print(f"  No sustained growth inhibition detected")
        print(f"  (May need to adjust threshold or extend observation time)")
    
    print("\n" + "="*60)
    print("Analysis complete!")
    print("="*60)


if __name__ == "__main__":
    main()
