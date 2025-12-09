# Bacterial Growth Divergence Detection

**Course:** Lecture 3 - Antibiotic Effect Detection  
**Organism:** *M. smegmatis* NCTC 8159  
**Task:** Detect when RIF-treated bacteria diverge from untreated controls (p < 0.05)

---

## 📋 Overview

This code detects the **earliest time point** where rifampicin (RIF) treatment causes bacterial growth to significantly diverge from untreated controls, based on:

- Total bacterial area from Omnipose segmentation masks
- Growth rate calculated with 30-minute sliding window
- Statistical significance testing (independent t-test, α = 0.05)

## 🗂️ Data Structure

The code expects the following folder structure:

```
F:/RM/
├── REF_masks101_110/          # Untreated control
│   ├── Pos101/
│   │   └── PreprocessedPhaseMasks/
│   │       ├── growth_areas.pickle      # Pre-calculated areas (preferred)
│   │       ├── MASK_img_000000000.tif   # Segmentation masks (fallback)
│   │       ├── MASK_img_000000001.tif
│   │       └── ...
│   ├── Pos102/ ... Pos110/
│
├── RIF10_masks201_210/        # RIF-treated
│   ├── Pos201/
│   │   └── PreprocessedPhaseMasks/
│   │       ├── growth_areas.pickle      # Pre-calculated areas (preferred)
│   │       ├── MASK_img_000000000.tif   # Segmentation masks (fallback)
│   │       └── ...
│   ├── Pos202/ ... Pos210/
│
└── detect_divergence.py       # Analysis script
```

**Note:** The code automatically uses `growth_areas.pickle` if available (faster), otherwise calculates from `.tif` masks.


## 🚀 Quick Start

### Requirements

```bash
pip install numpy scipy matplotlib imageio
```

Or use conda:
```bash
conda install numpy scipy matplotlib imageio
```

### Running the Analysis

1. **Ensure data is in the correct location** (`F:/RM/` or modify paths in script)

2. **Run the script:**
   ```bash
   python detect_divergence.py
   ```

3. **Check outputs:**
   - `divergence_analysis.png` - Visualization of results
   - `divergence_results.pkl` - Full data for further analysis
   - `divergence_report.txt` - Summary report

**Note:** The script uses non-interactive matplotlib backend (Agg) to avoid Qt conflicts. Figures are saved directly without display.

## 📊 Output Example

```
======================================================================
BACTERIAL GROWTH DIVERGENCE ANALYSIS REPORT
======================================================================

DATASET:
  - Control (REF): 10 replicates (Pos101-110)
  - Treated (RIF10): 10 replicates (Pos201-210)
  - Imaging interval: 2 minutes

METHODOLOGY:
  - Sliding window: 30 minutes (15 frames)
  - Growth rate: Exponential fit (log-linear)
  - Statistical test: Independent t-test
  - Significance level: α = 0.05
  - Minimum consecutive significant points: 3

RESULTS:
  ✓ DIVERGENCE DETECTED
  - Frame: 45
  - Time: 90.0 minutes
  - Time: 1.50 hours
  - P-value at divergence: 0.032145

======================================================================
```

## 🔬 Methodology

### 1. Data Loading
- **Primary method**: Loads pre-calculated areas from `growth_areas.pickle` (faster, consistent with instructor's pipeline)
- **Fallback method**: Calculates from segmentation masks if pickle not found
- Supports multiple replicates (positions) for statistical robustness

### 2. Growth Rate Calculation
Uses **exponential growth model** with sliding window:
- Window size: 30 minutes (15 frames at 2-min intervals)
- Fits exponential: `A(t) = A₀ · e^(b·t)` 
- Log-linear transformation: `log(A) = log(A₀) + b·t`
- Growth rate `b` is the slope in log-space (standard microbiology method)
- Each window yields one growth rate value

### 3. Divergence Detection
- Performs **independent t-test** at each timepoint
- Compares treated vs control growth rates across all replicates
- Finds **first timepoint** where p < 0.05 for **3 consecutive frames**
- Conservative approach ensures robust detection (avoids false positives from noise)

### 4. Visualization
Creates comprehensive 3-panel figure:
1. **Total Area**: Raw bacterial area over time (all replicates + mean)
2. **Growth Rate**: Exponential growth rates for both conditions
3. **P-values**: Statistical significance over time (log scale)
   - Divergence point marked with vertical line
   - α = 0.05 threshold shown

## 🛠️ Customization

Edit these parameters in `detect_divergence.py`:

```python
# Timing parameters
analyzer = GrowthAnalyzer(
    time_interval=2,      # Minutes between frames
    window_minutes=30     # Sliding window size
)

# Detection parameters
divergence_frame, p_values = analyzer.detect_divergence_ttest(
    rif_growth_rates, 
    ref_growth_rates,
    alpha=0.05,           # Significance level
    min_consecutive=3     # Consecutive significant points
)
```

