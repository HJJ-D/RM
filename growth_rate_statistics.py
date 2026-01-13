"""
Generate summary statistics of growth rates under control and treatment conditions.

Outputs:
- Mean, median, standard deviation, quartiles for each group
- Comparison table (suitable for publication)
- Optional statistical tests (t-test, Mann-Whitney U)
"""

import numpy as np
import pickle
from pathlib import Path
from scipy import stats
from growth_analysis import aggregate_growth_rates, compute_growth_statistics


def load_tracking_data():
    """Load tracking data."""
    script_dir = Path(__file__).parent
    ref_file = script_dir / "tracking_results_REF.pickle"
    rif_file = script_dir / "tracking_results_RIF.pickle"
    
    if not ref_file.exists() or not rif_file.exists():
        print("Error: Please run run_tracking.py first to generate tracking results")
        return None, None
    
    with open(ref_file, 'rb') as f:
        ref_results = pickle.load(f)
    with open(rif_file, 'rb') as f:
        rif_results = pickle.load(f)
    
    # Merge all positions
    tracks_ref = []
    for pos_data in ref_results.values():
        tracks_ref.extend(pos_data['tracks'])
    
    tracks_rif = []
    for pos_data in rif_results.values():
        tracks_rif.extend(pos_data['tracks'])
    
    return tracks_ref, tracks_rif


def compute_detailed_statistics(times, growth_rates, group_name):
    """Compute detailed statistics."""
    valid_rates = growth_rates[~np.isnan(growth_rates)]
    
    if len(valid_rates) == 0:
        return None
    
    stats_dict = {
        'group': group_name,
        'n_timepoints': len(valid_rates),
        'mean': np.mean(valid_rates),
        'std': np.std(valid_rates),
        'sem': np.std(valid_rates) / np.sqrt(len(valid_rates)),  # Standard Error of Mean
        'median': np.median(valid_rates),
        'q25': np.percentile(valid_rates, 25),
        'q75': np.percentile(valid_rates, 75),
        'iqr': np.percentile(valid_rates, 75) - np.percentile(valid_rates, 25),
        'min': np.min(valid_rates),
        'max': np.max(valid_rates),
        'range': np.max(valid_rates) - np.min(valid_rates),
        # Temporal analysis
        'early_mean': np.mean(valid_rates[:len(valid_rates)//3]) if len(valid_rates) >= 3 else np.nan,
        'mid_mean': np.mean(valid_rates[len(valid_rates)//3:2*len(valid_rates)//3]) if len(valid_rates) >= 3 else np.nan,
        'late_mean': np.mean(valid_rates[2*len(valid_rates)//3:]) if len(valid_rates) >= 3 else np.nan,
    }
    
    return stats_dict


def perform_statistical_tests(growth_ctrl, growth_treat):
    """Perform statistical tests."""
    # Remove NaN values
    ctrl_valid = growth_ctrl[~np.isnan(growth_ctrl)]
    treat_valid = growth_treat[~np.isnan(growth_treat)]
    
    results = {}
    
    # 1. Independent t-test (assumes normal distribution)
    t_stat, t_pvalue = stats.ttest_ind(ctrl_valid, treat_valid)
    results['t_test'] = {
        'statistic': t_stat,
        'p_value': t_pvalue,
        'significant_0.05': t_pvalue < 0.05,
        'significant_0.01': t_pvalue < 0.01
    }
    
    # 2. Mann-Whitney U test (non-parametric, no normality assumption)
    u_stat, u_pvalue = stats.mannwhitneyu(ctrl_valid, treat_valid, alternative='two-sided')
    results['mann_whitney'] = {
        'statistic': u_stat,
        'p_value': u_pvalue,
        'significant_0.05': u_pvalue < 0.05,
        'significant_0.01': u_pvalue < 0.01
    }
    
    # 3. Effect size (Cohen's d)
    pooled_std = np.sqrt((np.var(ctrl_valid) + np.var(treat_valid)) / 2)
    cohens_d = (np.mean(ctrl_valid) - np.mean(treat_valid)) / pooled_std if pooled_std > 0 else 0
    results['effect_size'] = {
        'cohens_d': cohens_d,
        'interpretation': interpret_cohens_d(cohens_d)
    }
    
    # 4. Percent change
    mean_ctrl = np.mean(ctrl_valid)
    mean_treat = np.mean(treat_valid)
    pct_change = ((mean_treat - mean_ctrl) / abs(mean_ctrl)) * 100 if mean_ctrl != 0 else 0
    results['percent_change'] = pct_change
    
    return results


def interpret_cohens_d(d):
    """Interpret Cohen's d effect size."""
    d = abs(d)
    if d < 0.2:
        return "negligible"
    elif d < 0.5:
        return "small"
    elif d < 0.8:
        return "medium"
    else:
        return "large"


def print_summary_table(stats_ctrl, stats_treat, test_results):
    """Print summary table (suitable for copying to paper)."""
    print("\n" + "="*80)
    print("SUMMARY STATISTICS OF GROWTH RATES")
    print("="*80)
    
    # Table header
    print(f"\n{'Metric':<25} {'Control (REF)':<20} {'Treatment (RIF)':<20}")
    print("-"*65)
    
    # Sample info
    print(f"{'N (time points)':<25} {stats_ctrl['n_timepoints']:<20} {stats_treat['n_timepoints']:<20}")
    
    # Central tendency
    print(f"\n--- Central Tendency ---")
    print(f"{'Mean ± SEM':<25} {stats_ctrl['mean']:.4f} ± {stats_ctrl['sem']:.4f}       {stats_treat['mean']:.4f} ± {stats_treat['sem']:.4f}")
    print(f"{'Median':<25} {stats_ctrl['median']:<20.4f} {stats_treat['median']:<20.4f}")
    
    # Dispersion
    print(f"\n--- Dispersion ---")
    print(f"{'Std Dev':<25} {stats_ctrl['std']:<20.4f} {stats_treat['std']:<20.4f}")
    print(f"{'IQR (Q75-Q25)':<25} {stats_ctrl['iqr']:<20.4f} {stats_treat['iqr']:<20.4f}")
    print(f"{'Range (Max-Min)':<25} {stats_ctrl['range']:<20.4f} {stats_treat['range']:<20.4f}")
    
    # Quartiles
    print(f"\n--- Quartiles ---")
    print(f"{'Q25 (25th percentile)':<25} {stats_ctrl['q25']:<20.4f} {stats_treat['q25']:<20.4f}")
    print(f"{'Q75 (75th percentile)':<25} {stats_ctrl['q75']:<20.4f} {stats_treat['q75']:<20.4f}")
    print(f"{'Min':<25} {stats_ctrl['min']:<20.4f} {stats_treat['min']:<20.4f}")
    print(f"{'Max':<25} {stats_ctrl['max']:<20.4f} {stats_treat['max']:<20.4f}")
    
    # Temporal analysis
    print(f"\n--- Temporal Analysis ---")
    print(f"{'Early phase mean':<25} {stats_ctrl['early_mean']:<20.4f} {stats_treat['early_mean']:<20.4f}")
    print(f"{'Mid phase mean':<25} {stats_ctrl['mid_mean']:<20.4f} {stats_treat['mid_mean']:<20.4f}")
    print(f"{'Late phase mean':<25} {stats_ctrl['late_mean']:<20.4f} {stats_treat['late_mean']:<20.4f}")
    
    # Statistical tests
    print(f"\n" + "="*80)
    print("STATISTICAL COMPARISONS")
    print("="*80)
    
    print(f"\n--- Percent Change ---")
    pct = test_results['percent_change']
    direction = "decrease" if pct < 0 else "increase"
    print(f"Treatment vs Control: {pct:+.2f}% ({direction})")
    
    print(f"\n--- Independent t-test ---")
    t = test_results['t_test']
    print(f"t-statistic: {t['statistic']:.4f}")
    print(f"p-value: {t['p_value']:.2e}")
    print(f"Significant at α=0.05: {'Yes' if t['significant_0.05'] else 'No'}")
    print(f"Significant at α=0.01: {'Yes' if t['significant_0.01'] else 'No'}")
    
    print(f"\n--- Mann-Whitney U test (non-parametric) ---")
    mw = test_results['mann_whitney']
    print(f"U-statistic: {mw['statistic']:.4f}")
    print(f"p-value: {mw['p_value']:.2e}")
    print(f"Significant at α=0.05: {'Yes' if mw['significant_0.05'] else 'No'}")
    print(f"Significant at α=0.01: {'Yes' if mw['significant_0.01'] else 'No'}")
    
    print(f"\n--- Effect Size ---")
    es = test_results['effect_size']
    print(f"Cohen's d: {es['cohens_d']:.4f} ({es['interpretation']})")
    
    print("\n" + "="*80)


def generate_latex_table(stats_ctrl, stats_treat, test_results):
    """Generate LaTeX table code."""
    print("\n" + "="*80)
    print("LaTeX TABLE CODE (for paper)")
    print("="*80)
    
    latex = r"""
\begin{table}[htbp]
\centering
\caption{Summary statistics of growth rates under control and treatment conditions}
\label{tab:growth_stats}
\begin{tabular}{lcc}
\toprule
\textbf{Metric} & \textbf{Control} & \textbf{Treatment} \\
\midrule
N (time points) & %d & %d \\
Mean $\pm$ SEM & %.4f $\pm$ %.4f & %.4f $\pm$ %.4f \\
Median & %.4f & %.4f \\
Std Dev & %.4f & %.4f \\
IQR & %.4f & %.4f \\
Range & %.4f & %.4f \\
\midrule
\multicolumn{3}{l}{\textit{Statistical comparison}} \\
Percent change & \multicolumn{2}{c}{%+.2f\%%} \\
t-test p-value & \multicolumn{2}{c}{%.2e} \\
Mann-Whitney p-value & \multicolumn{2}{c}{%.2e} \\
Cohen's d & \multicolumn{2}{c}{%.4f (%s)} \\
\bottomrule
\end{tabular}
\end{table}
""" % (
        stats_ctrl['n_timepoints'], stats_treat['n_timepoints'],
        stats_ctrl['mean'], stats_ctrl['sem'], stats_treat['mean'], stats_treat['sem'],
        stats_ctrl['median'], stats_treat['median'],
        stats_ctrl['std'], stats_treat['std'],
        stats_ctrl['iqr'], stats_treat['iqr'],
        stats_ctrl['range'], stats_treat['range'],
        test_results['percent_change'],
        test_results['t_test']['p_value'],
        test_results['mann_whitney']['p_value'],
        test_results['effect_size']['cohens_d'],
        test_results['effect_size']['interpretation']
    )
    
    print(latex)


def main():
    print("="*80)
    print("GROWTH RATE SUMMARY STATISTICS")
    print("Control (REF) vs Treatment (RIF)")
    print("="*80)
    
    # Load data
    print("\nLoading tracking data...")
    tracks_ref, tracks_rif = load_tracking_data()
    
    if tracks_ref is None:
        return
    
    print(f"  Control tracks: {len(tracks_ref)}")
    print(f"  Treatment tracks: {len(tracks_rif)}")
    
    # Compute population-level growth rates
    print("\nComputing population-level growth rates...")
    
    times_ctrl, growth_ctrl = aggregate_growth_rates(
        tracks_ref, method="median", window=9
    )
    print(f"  Control: {len(times_ctrl)} time points")
    
    times_treat, growth_treat = aggregate_growth_rates(
        tracks_rif, method="median", window=9
    )
    print(f"  Treatment: {len(times_treat)} time points")
    
    # Compute detailed statistics
    stats_ctrl = compute_detailed_statistics(times_ctrl, growth_ctrl, "Control")
    stats_treat = compute_detailed_statistics(times_treat, growth_treat, "Treatment")
    
    # Perform statistical tests
    test_results = perform_statistical_tests(growth_ctrl, growth_treat)
    
    # Output results
    print_summary_table(stats_ctrl, stats_treat, test_results)
    
    # Generate LaTeX table
    generate_latex_table(stats_ctrl, stats_treat, test_results)
    
    return {
        'control': stats_ctrl,
        'treatment': stats_treat,
        'tests': test_results,
        'raw_data': {
            'times_ctrl': times_ctrl,
            'growth_ctrl': growth_ctrl,
            'times_treat': times_treat,
            'growth_treat': growth_treat
        }
    }


if __name__ == "__main__":
    results = main()
