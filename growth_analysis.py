"""
Post-tracking growth rate analysis pipeline for single-cell microscopy data.

Computes:
- Single-cell growth rates from tracked cell elongation
- Aggregated population-level growth using robust statistics (median/quantile)
- Normalized treatment vs control responses
- Automatic detection of earliest drug response time

Designed for scientific analysis with robustness to outliers and resistant cells.
"""

import numpy as np
from scipy import stats
from typing import List, Dict, Tuple, Optional
import warnings


def extract_length_signal(track: Dict, key: str = "major_axis_length") -> Tuple[np.ndarray, np.ndarray]:
    """
    Extract temporal length signal from a single track.
    
    Parameters
    ----------
    track : Dict
        Track dictionary containing 'frames' and 'cells' keys
    key : str
        Morphology key to extract (default: "major_axis_length")
        Falls back to "area" if key is missing or zero
        
    Returns
    -------
    times : np.ndarray
        Frame indices (time points)
    values : np.ndarray
        Corresponding length/area values
        
    Notes
    -----
    - Filters out zero/missing values
    - Returns empty arrays if no valid data
    """
    frames = track['frames']
    cells = track['cells']
    
    times = []
    values = []
    
    for t, cell in zip(frames, cells):
        # Get primary key value
        val = cell.get(key, 0)
        
        # Fallback to area if key is missing or zero
        if val == 0 or val is None or np.isnan(val):
            val = cell.get('area', 0)
        
        # Skip invalid values
        if val > 0 and not np.isnan(val):
            times.append(t)
            values.append(val)
    
    return np.array(times), np.array(values)


def compute_local_growth_rate(
    times: np.ndarray,
    values: np.ndarray,
    window: int = 9
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute local growth rate using sliding window linear regression.
    
    For each window centered at time t:
    - Fit linear model: value = a + b * time
    - Growth rate = slope b
    
    Parameters
    ----------
    times : np.ndarray
        Time points (frame indices)
    values : np.ndarray
        Measured values (e.g., cell length)
    window : int
        Sliding window size (odd number preferred for symmetric windows)
        
    Returns
    -------
    center_times : np.ndarray
        Time points at window centers
    growth_rates : np.ndarray
        Local growth rates (slopes)
        
    Notes
    -----
    - Uses simple linear regression within each window
    - Skips windows with insufficient points (< 3)
    - Handles NaNs gracefully
    - Growth rate units: [value_units / frame]
    """
    n = len(times)
    
    # Ensure odd window for symmetric centering
    if window % 2 == 0:
        window += 1
    
    half_window = window // 2
    
    center_times = []
    growth_rates = []
    
    for i in range(n):
        # Define window bounds
        start = max(0, i - half_window)
        end = min(n, i + half_window + 1)
        
        # Extract window data
        t_window = times[start:end]
        v_window = values[start:end]
        
        # Filter out NaNs
        valid_mask = ~np.isnan(v_window)
        t_valid = t_window[valid_mask]
        v_valid = v_window[valid_mask]
        
        # Need at least 3 points for meaningful regression
        if len(t_valid) < 3:
            continue
        
        # Linear regression: value = a + b * time
        # Use scipy.stats.linregress for robustness
        try:
            slope, intercept, r_value, p_value, std_err = stats.linregress(t_valid, v_valid)
            
            center_times.append(times[i])
            growth_rates.append(slope)
            
        except Exception:
            # Skip if regression fails
            continue
    
    return np.array(center_times), np.array(growth_rates)


def aggregate_growth_rates(
    tracks: List[Dict],
    method: str = "median",
    q: float = 0.25,
    key: str = "major_axis_length",
    window: int = 9,
    min_cells: int = 5
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Aggregate growth rates across all tracks using robust statistics.
    
    Procedure:
    1. For each track: extract length signal and compute local growth rates
    2. For each time point: collect all growth rates from tracks covering that time
    3. Aggregate using median or quantile
    
    Parameters
    ----------
    tracks : List[Dict]
        List of track dictionaries
    method : str
        Aggregation method: "median" or "quantile"
    q : float
        Quantile value if method="quantile" (e.g., 0.25 for 25th percentile)
    key : str
        Morphology key for length extraction
    window : int
        Sliding window size for growth rate computation
    min_cells : int
        Minimum number of cells required at a time point for aggregation
        
    Returns
    -------
    times : np.ndarray
        Time points (sorted)
    aggregated_growth : np.ndarray
        Aggregated growth rates at each time point
        
    Notes
    -----
    - Uses median (default) or low quantile to suppress outliers/resistant cells
    - Ignores time points with too few samples (< min_cells)
    - Returns NaN for time points with insufficient data
    """
    # Step 1: Compute growth rates for all tracks
    all_growth_data = []
    
    for track in tracks:
        # Extract length signal
        times_track, values_track = extract_length_signal(track, key=key)
        
        if len(times_track) < window:
            # Track too short for sliding window
            continue
        
        # Compute local growth rates
        center_times, growth_rates = compute_local_growth_rate(
            times_track, values_track, window=window
        )
        
        if len(center_times) == 0:
            continue
        
        # Store as list of (time, growth_rate) tuples
        for t, gr in zip(center_times, growth_rates):
            all_growth_data.append((t, gr))
    
    if len(all_growth_data) == 0:
        return np.array([]), np.array([])
    
    # Step 2: Organize data by time point
    # Convert to structured array for easier grouping
    all_growth_data = np.array(all_growth_data, dtype=[('time', float), ('growth', float)])
    
    # Get unique time points
    unique_times = np.unique(all_growth_data['time'])
    unique_times = np.sort(unique_times)
    
    # Step 3: Aggregate at each time point
    aggregated_values = []
    valid_times = []
    
    for t in unique_times:
        # Get all growth rates at this time
        mask = all_growth_data['time'] == t
        growth_at_t = all_growth_data['growth'][mask]
        
        # Filter NaNs
        growth_at_t = growth_at_t[~np.isnan(growth_at_t)]
        
        # Check minimum sample size
        if len(growth_at_t) < min_cells:
            continue
        
        # Aggregate
        if method == "median":
            agg_value = np.median(growth_at_t)
        elif method == "quantile":
            agg_value = np.quantile(growth_at_t, q)
        else:
            raise ValueError(f"Unknown method: {method}. Use 'median' or 'quantile'.")
        
        valid_times.append(t)
        aggregated_values.append(agg_value)
    
    return np.array(valid_times), np.array(aggregated_values)


def smooth_series(
    times: np.ndarray,
    values: np.ndarray,
    method: str = "rolling_median",
    window: int = 7
) -> np.ndarray:
    """
    Temporally smooth a time series to reduce noise.
    
    Parameters
    ----------
    times : np.ndarray
        Time points (frame indices), must be sorted
    values : np.ndarray
        Values to smooth
    method : str
        Smoothing method: "rolling_median", "moving_average", or "savgol"
    window : int
        Window size (odd number preferred for symmetric windows)
        
    Returns
    -------
    smoothed_values : np.ndarray
        Smoothed values aligned with input
        
    Notes
    -----
    - rolling_median is most robust to outliers (recommended)
    - moving_average is simple but sensitive to outliers
    - savgol uses Savitzky-Golay filter (requires scipy)
    - Endpoints are handled by using smaller windows
    - NaNs are ignored in smoothing
    """
    if len(values) == 0:
        return values.copy()
    
    # Ensure odd window for symmetry
    if window % 2 == 0:
        window += 1
    
    half_window = window // 2
    smoothed = np.zeros_like(values, dtype=float)
    
    if method == "rolling_median":
        # Rolling median - most robust
        for i in range(len(values)):
            # Define window bounds
            start = max(0, i - half_window)
            end = min(len(values), i + half_window + 1)
            
            # Extract window and filter NaNs
            window_vals = values[start:end]
            valid_vals = window_vals[~np.isnan(window_vals)]
            
            if len(valid_vals) > 0:
                smoothed[i] = np.median(valid_vals)
            else:
                smoothed[i] = np.nan
    
    elif method == "moving_average":
        # Simple moving average
        for i in range(len(values)):
            start = max(0, i - half_window)
            end = min(len(values), i + half_window + 1)
            
            window_vals = values[start:end]
            valid_vals = window_vals[~np.isnan(window_vals)]
            
            if len(valid_vals) > 0:
                smoothed[i] = np.mean(valid_vals)
            else:
                smoothed[i] = np.nan
    
    elif method == "savgol":
        # Savitzky-Golay filter (requires scipy)
        try:
            from scipy.signal import savgol_filter
            
            # Handle NaNs by interpolation
            valid_mask = ~np.isnan(values)
            if np.sum(valid_mask) < window:
                warnings.warn("Too few valid points for savgol filter, using median instead")
                return smooth_series(times, values, method="rolling_median", window=window)
            
            # Only smooth valid regions
            if np.all(valid_mask):
                # No NaNs, apply directly
                polyorder = min(3, window - 2)
                smoothed = savgol_filter(values, window, polyorder)
            else:
                # Has NaNs, interpolate first
                valid_indices = np.where(valid_mask)[0]
                smoothed = values.copy()
                
                if len(valid_indices) >= window:
                    polyorder = min(3, window - 2)
                    smoothed[valid_mask] = savgol_filter(
                        values[valid_mask], window, polyorder
                    )
        
        except ImportError:
            warnings.warn("scipy not available for savgol, using rolling_median instead")
            return smooth_series(times, values, method="rolling_median", window=window)
    
    else:
        raise ValueError(f"Unknown smoothing method: {method}")
    
    return smoothed


def fit_control_baseline(
    times_ctrl: np.ndarray,
    g_ctrl: np.ndarray,
    method: str = "median"
) -> float:
    """
    Fit a single global baseline value for control growth rate.
    
    This provides a constant, stable baseline for normalization when
    control data is too noisy for temporal normalization.
    
    Parameters
    ----------
    times_ctrl : np.ndarray
        Time points for control group
    g_ctrl : np.ndarray
        Growth rates for control group
    method : str
        Method for computing baseline: "median", "mean", or "trimmed_mean"
        
    Returns
    -------
    baseline : float
        Scalar baseline value
        
    Notes
    -----
    - "median" is most robust (recommended)
    - "trimmed_mean" excludes top/bottom 10%
    - "mean" is sensitive to outliers
    """
    # Filter NaNs
    valid_vals = g_ctrl[~np.isnan(g_ctrl)]
    
    if len(valid_vals) == 0:
        return np.nan
    
    if method == "median":
        baseline = np.median(valid_vals)
    
    elif method == "mean":
        baseline = np.mean(valid_vals)
    
    elif method == "trimmed_mean":
        # Trim top and bottom 10%
        from scipy.stats import trim_mean
        baseline = trim_mean(valid_vals, proportiontocut=0.1)
    
    else:
        raise ValueError(f"Unknown baseline method: {method}")
    
    return baseline


def normalize_growth(
    treat_times: np.ndarray,
    g_treat: np.ndarray,
    ctrl_times: np.ndarray,
    g_ctrl: np.ndarray,
    eps: float = 1e-6
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Normalize treatment growth rates by control.
    
    Computes: g_norm(t) = g_treat(t) / (g_ctrl(t) + eps)
    
    Parameters
    ----------
    treat_times : np.ndarray
        Time points for treatment group
    g_treat : np.ndarray
        Growth rates for treatment group
    ctrl_times : np.ndarray
        Time points for control group
    g_ctrl : np.ndarray
        Growth rates for control group
    eps : float
        Small constant to avoid division by zero
        
    Returns
    -------
    times : np.ndarray
        Common time points (intersection of treat and ctrl)
    normalized_growth : np.ndarray
        Normalized growth rates
        
    Notes
    -----
    - Only computes normalized values for time points present in both groups
    - Normalized growth ~ 1.0 means no drug effect
    - Normalized growth < 1.0 means drug reduces growth
    - Normalized growth > 1.0 means drug increases growth (rare)
    """
    # Find common time points (intersection)
    common_times = np.intersect1d(treat_times, ctrl_times)
    
    if len(common_times) == 0:
        return np.array([]), np.array([])
    
    # Extract values at common times
    g_treat_common = []
    g_ctrl_common = []
    
    for t in common_times:
        # Find indices
        idx_treat = np.where(treat_times == t)[0]
        idx_ctrl = np.where(ctrl_times == t)[0]
        
        if len(idx_treat) > 0 and len(idx_ctrl) > 0:
            g_treat_common.append(g_treat[idx_treat[0]])
            g_ctrl_common.append(g_ctrl[idx_ctrl[0]])
    
    g_treat_common = np.array(g_treat_common)
    g_ctrl_common = np.array(g_ctrl_common)
    
    # Normalize: ratio of treatment to control
    normalized_growth = g_treat_common / (g_ctrl_common + eps)
    
    return common_times[:len(normalized_growth)], normalized_growth


def normalize_growth_v2(
    times_treat: np.ndarray,
    g_treat: np.ndarray,
    times_ctrl: np.ndarray,
    g_ctrl: np.ndarray,
    mode: str = "smooth",
    smooth_window: int = 7,
    baseline_method: str = "median",
    eps: float = 1e-6
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Improved normalization with stable control baseline.
    
    Provides two modes:
    - "smooth": temporally smooth control before normalization
    - "baseline": use a constant global baseline for control
    
    Parameters
    ----------
    times_treat : np.ndarray
        Time points for treatment group
    g_treat : np.ndarray
        Growth rates for treatment group
    times_ctrl : np.ndarray
        Time points for control group
    g_ctrl : np.ndarray
        Growth rates for control group
    mode : str
        Normalization mode: "smooth" or "baseline"
    smooth_window : int
        Window size for smoothing (if mode="smooth")
    baseline_method : str
        Method for baseline fitting (if mode="baseline"): "median", "mean", "trimmed_mean"
    eps : float
        Small constant to avoid division by zero
        
    Returns
    -------
    times : np.ndarray
        Common time points
    normalized_growth : np.ndarray
        Normalized growth rates
        
    Notes
    -----
    Mode "smooth":
    - Smooths noisy control data before normalization
    - Preserves temporal variation in control
    - Recommended when control shows trends
    
    Mode "baseline":
    - Uses constant baseline (median of control)
    - Most stable for very noisy control data
    - Assumes control growth is approximately constant
    """
    # Find common time points
    common_times = np.intersect1d(times_treat, times_ctrl)
    
    if len(common_times) == 0:
        return np.array([]), np.array([])
    
    # Extract values at common times
    g_treat_common = []
    g_ctrl_common = []
    
    for t in common_times:
        idx_treat = np.where(times_treat == t)[0]
        idx_ctrl = np.where(times_ctrl == t)[0]
        
        if len(idx_treat) > 0 and len(idx_ctrl) > 0:
            g_treat_common.append(g_treat[idx_treat[0]])
            g_ctrl_common.append(g_ctrl[idx_ctrl[0]])
    
    g_treat_common = np.array(g_treat_common)
    g_ctrl_common = np.array(g_ctrl_common)
    
    # Apply normalization based on mode
    if mode == "smooth":
        # Smooth control data before normalization
        g_ctrl_smooth = smooth_series(
            common_times,
            g_ctrl_common,
            method="rolling_median",
            window=smooth_window
        )
        normalized_growth = g_treat_common / (g_ctrl_smooth + eps)
    
    elif mode == "baseline":
        # Use constant baseline
        baseline = fit_control_baseline(times_ctrl, g_ctrl, method=baseline_method)
        normalized_growth = g_treat_common / (baseline + eps)
    
    else:
        raise ValueError(f"Unknown mode: {mode}. Use 'smooth' or 'baseline'.")
    
    return common_times, normalized_growth


def detect_response_time(
    times: np.ndarray,
    normalized_growth: np.ndarray,
    threshold: float = 0.8,
    consecutive: int = 3
) -> Optional[int]:
    """
    Detect earliest drug response time.
    
    Finds the first time t* where normalized growth drops below threshold
    and remains below for at least 'consecutive' frames.
    
    Parameters
    ----------
    times : np.ndarray
        Time points
    normalized_growth : np.ndarray
        Normalized growth rates (treatment / control)
    threshold : float
        Threshold for response detection (e.g., 0.8 = 20% reduction)
    consecutive : int
        Number of consecutive frames that must satisfy condition
        
    Returns
    -------
    t_star : int or None
        Earliest response time (frame index), or None if not detected
        
    Notes
    -----
    - Requires sustained reduction, not just a single dip
    - Conservative approach reduces false positives
    - Typical threshold: 0.7-0.9 (10-30% growth reduction)
    """
    if len(times) < consecutive:
        return None
    
    # Check each potential starting point
    for i in range(len(times) - consecutive + 1):
        # Check if next 'consecutive' points are all below threshold
        window_values = normalized_growth[i:i + consecutive]
        
        if np.all(window_values < threshold):
            # Found sustained response
            return int(times[i])
    
    # No sustained response detected
    return None


def detect_response_time_v2(
    times: np.ndarray,
    normalized_growth: np.ndarray,
    threshold: float = 0.8,
    consecutive: int = 3,
    min_after: int = 0,
    smooth_norm: bool = True,
    norm_smooth_window: int = 5
) -> Optional[int]:
    """
    Improved drug response time detection with optional smoothing.
    
    More stable than detect_response_time by optionally smoothing
    normalized growth before detection to avoid triggering on single dips.
    
    Parameters
    ----------
    times : np.ndarray
        Time points
    normalized_growth : np.ndarray
        Normalized growth rates (treatment / control)
    threshold : float
        Threshold for response detection (e.g., 0.8 = 20% reduction)
    consecutive : int
        Number of consecutive frames that must satisfy condition
    min_after : int
        Ignore frames before this time (e.g., to skip early artifacts)
    smooth_norm : bool
        If True, smooth normalized growth before detection
    norm_smooth_window : int
        Window size for smoothing normalized growth (if smooth_norm=True)
        
    Returns
    -------
    t_star : int or None
        Earliest response time (frame index), or None if not detected
        
    Notes
    -----
    - Smoothing normalized growth (default) makes detection more robust
    - min_after allows skipping initial frames (e.g., adjustment period)
    - Requires sustained reduction to avoid false positives
    """
    if len(times) < consecutive:
        return None
    
    # Optionally smooth normalized growth
    if smooth_norm:
        g_norm_smooth = smooth_series(
            times,
            normalized_growth,
            method="rolling_median",
            window=norm_smooth_window
        )
    else:
        g_norm_smooth = normalized_growth
    
    # Apply min_after filter
    valid_mask = times >= min_after
    times_valid = times[valid_mask]
    g_norm_valid = g_norm_smooth[valid_mask]
    
    if len(times_valid) < consecutive:
        return None
    
    # Check each potential starting point
    for i in range(len(times_valid) - consecutive + 1):
        # Check if next 'consecutive' points are all below threshold
        window_values = g_norm_valid[i:i + consecutive]
        
        if np.all(window_values < threshold):
            # Found sustained response
            return int(times_valid[i])
    
    # No sustained response detected
    return None


def compute_growth_statistics(
    times: np.ndarray,
    growth_rates: np.ndarray
) -> Dict:
    """
    Compute summary statistics for growth rates.
    
    Parameters
    ----------
    times : np.ndarray
        Time points
    growth_rates : np.ndarray
        Growth rate values
        
    Returns
    -------
    Dict
        Statistics including mean, median, std, quartiles, etc.
    """
    valid_rates = growth_rates[~np.isnan(growth_rates)]
    
    if len(valid_rates) == 0:
        return {
            'mean': np.nan,
            'median': np.nan,
            'std': np.nan,
            'q25': np.nan,
            'q75': np.nan,
            'min': np.nan,
            'max': np.nan,
            'n_points': 0
        }
    
    return {
        'mean': np.mean(valid_rates),
        'median': np.median(valid_rates),
        'std': np.std(valid_rates),
        'q25': np.percentile(valid_rates, 25),
        'q75': np.percentile(valid_rates, 75),
        'min': np.min(valid_rates),
        'max': np.max(valid_rates),
        'n_points': len(valid_rates)
    }


# ============================================================================
# EXAMPLE USAGE
# ============================================================================

def analyze_drug_response(
    tracks_treat: List[Dict],
    tracks_ctrl: List[Dict],
    use_v2: bool = True,
    norm_mode: str = "smooth",
    smooth_window: int = 7
):
    """
    Complete drug response analysis pipeline.
    
    Parameters
    ----------
    tracks_treat : List[Dict]
        Tracks from drug-treated group
    tracks_ctrl : List[Dict]
        Tracks from control group
    use_v2 : bool
        If True, use improved normalization and detection methods
    norm_mode : str
        Normalization mode: "smooth" or "baseline" (only if use_v2=True)
    smooth_window : int
        Window size for smoothing control (only if use_v2=True and norm_mode="smooth")
        
    Returns
    -------
    Dict
        Analysis results including growth rates, normalized response, and response time
    """
    print("=" * 60)
    print("Drug Response Analysis Pipeline")
    print("=" * 60)
    
    # Step 1: Aggregate growth rates for each group
    print("\nStep 1: Computing population-level growth rates...")
    print("  Treatment group...")
    times_treat, g_treat = aggregate_growth_rates(
        tracks_treat,
        method="median",  # Robust to outliers
        window=9
    )
    print(f"    {len(times_treat)} time points, median growth = {np.median(g_treat):.4f}")
    
    print("  Control group...")
    times_ctrl, g_ctrl = aggregate_growth_rates(
        tracks_ctrl,
        method="median",
        window=9
    )
    print(f"    {len(times_ctrl)} time points, median growth = {np.median(g_ctrl):.4f}")
    
    # Step 2: Normalize treatment by control
    print("\nStep 2: Normalizing treatment vs control...")
    if use_v2:
        print(f"  Using method: {norm_mode}")
        if norm_mode == "smooth":
            print(f"  Smoothing window: {smooth_window}")
        times_norm, g_norm = normalize_growth_v2(
            times_treat, g_treat,
            times_ctrl, g_ctrl,
            mode=norm_mode,
            smooth_window=smooth_window
        )
    else:
        times_norm, g_norm = normalize_growth(times_treat, g_treat, times_ctrl, g_ctrl)
    print(f"  {len(times_norm)} common time points")
    print(f"  Mean normalized growth: {np.mean(g_norm):.3f}")
    print(f"  Median normalized growth: {np.median(g_norm):.3f}")
    
    # Step 3: Detect drug response time
    print("\nStep 3: Detecting drug response time...")
    if use_v2:
        print("  Using improved detection with smoothing")
        t_star = detect_response_time_v2(
            times_norm,
            g_norm,
            threshold=0.8,         # 20% reduction
            consecutive=3,         # Sustained for 3 frames
            smooth_norm=True,      # Smooth before detection
            norm_smooth_window=5   # Small smoothing window
        )
    else:
        t_star = detect_response_time(
            times_norm,
            g_norm,
            threshold=0.8,
            consecutive=3
        )
    
    if t_star is not None:
        print(f"  ✓ Drug response detected at frame {t_star}")
        print(f"    (20% sustained growth reduction)")
    else:
        print(f"  ✗ No sustained drug response detected")
        print(f"    (No 20% reduction sustained for 3+ frames)")
    
    # Step 4: Compute statistics
    print("\nStep 4: Growth rate statistics...")
    
    stats_treat = compute_growth_statistics(times_treat, g_treat)
    print(f"  Treatment: median = {stats_treat['median']:.4f}, "
          f"IQR = [{stats_treat['q25']:.4f}, {stats_treat['q75']:.4f}]")
    
    stats_ctrl = compute_growth_statistics(times_ctrl, g_ctrl)
    print(f"  Control:   median = {stats_ctrl['median']:.4f}, "
          f"IQR = [{stats_ctrl['q25']:.4f}, {stats_ctrl['q75']:.4f}]")
    
    print("\n" + "=" * 60)
    
    # Return results
    return {
        'treatment': {
            'times': times_treat,
            'growth_rates': g_treat,
            'statistics': stats_treat
        },
        'control': {
            'times': times_ctrl,
            'growth_rates': g_ctrl,
            'statistics': stats_ctrl
        },
        'normalized': {
            'times': times_norm,
            'growth_rates': g_norm,
            'response_time': t_star
        }
    }


if __name__ == "__main__":
    """
    Minimal example with synthetic data
    """
    print("Testing growth analysis with synthetic tracks...")
    
    # Create synthetic tracks
    # Control: steady growth
    tracks_ctrl = []
    for i in range(20):
        frames = list(range(60))
        cells = []
        for t in frames:
            # Linear growth + small noise
            length = 10 + 0.5 * t + np.random.normal(0, 0.1)
            cells.append({
                'major_axis_length': length,
                'area': length * 3,
                'centroid': (50 + i, 50)
            })
        tracks_ctrl.append({'frames': frames, 'cells': cells})
    
    # Treatment: growth slows down after frame 20
    tracks_treat = []
    for i in range(20):
        frames = list(range(60))
        cells = []
        for t in frames:
            # Growth slows after t=20
            if t < 20:
                growth_rate = 0.5
            else:
                growth_rate = 0.2  # 60% reduction
            
            length = 10 + growth_rate * t + np.random.normal(0, 0.1)
            cells.append({
                'major_axis_length': length,
                'area': length * 3,
                'centroid': (50 + i, 50)
            })
        tracks_treat.append({'frames': frames, 'cells': cells})
    
    # Run analysis with improved methods (v2)
    print("\n" + "="*60)
    print("Testing with improved normalization (v2)")
    print("="*60)
    results = analyze_drug_response(
        tracks_treat, tracks_ctrl,
        use_v2=True,
        norm_mode="smooth",
        smooth_window=7
    )
    
    print("\nSynthetic test completed successfully!")
    print(f"Expected response time: ~20")
    print(f"Detected response time: {results['normalized']['response_time']}")
    
    # Also test baseline mode
    print("\n" + "="*60)
    print("Testing with baseline normalization")
    print("="*60)
    results_baseline = analyze_drug_response(
        tracks_treat, tracks_ctrl,
        use_v2=True,
        norm_mode="baseline"
    )
    print(f"Response time (baseline mode): {results_baseline['normalized']['response_time']}")
    
    # Demonstrate direct usage of new functions
    print("\n" + "="*60)
    print("Direct usage example")
    print("="*60)
    
    # Get aggregated growth rates
    times_treat, g_treat = aggregate_growth_rates(tracks_treat, window=9)
    times_ctrl, g_ctrl = aggregate_growth_rates(tracks_ctrl, window=9)
    
    # Method 1: Smooth normalization
    print("\nMethod 1: Smooth normalization")
    times, g_norm = normalize_growth_v2(
        times_treat, g_treat,
        times_ctrl, g_ctrl,
        mode="smooth",
        smooth_window=7
    )
    t_star = detect_response_time_v2(times, g_norm, threshold=0.8, consecutive=3)
    print(f"  Response time: {t_star}")
    
    # Method 2: Baseline normalization
    print("\nMethod 2: Baseline normalization")
    times2, g_norm2 = normalize_growth_v2(
        times_treat, g_treat,
        times_ctrl, g_ctrl,
        mode="baseline"
    )
    t_star2 = detect_response_time_v2(times2, g_norm2, threshold=0.8, consecutive=3)
    print(f"  Response time: {t_star2}")
    print(f"  Control baseline: {fit_control_baseline(times_ctrl, g_ctrl):.4f}")
