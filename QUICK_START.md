# Quick Run Instructions

## For Evaluators - Fast Setup

### Step 1: Clone Repository
```bash
git clone https://github.com/YOUR_USERNAME/bacterial-growth-divergence.git
cd bacterial-growth-divergence
```

### Step 2: Install Dependencies
```bash
pip install -r requirements.txt
```

Or manually:
```bash
pip install numpy scipy matplotlib imageio
```

### Step 3: Get Data
The analysis requires segmentation masks from Lecture 3:
- **Control**: `REF_masks101_110/` folder (Positions 101-110)
- **Treated**: `RIF10_masks201_210/` folder (Positions 201-210)

**Data Location**: [Contact repository owner or check instructor's shared folder]

### Step 4: Update Paths
Edit `detect_divergence.py` (lines 187-188) to point to your data:
```python
REF_BASE = Path("YOUR_PATH/REF_masks101_110")
RIF_BASE = Path("YOUR_PATH/RIF10_masks201_210")
```

### Step 5: Run
```bash
python detect_divergence.py
```

### Step 6: Check Results
Look for:
- **divergence_report.txt** - The answer to "when do curves diverge?"
- **divergence_analysis.png** - Visual proof

## Expected Output

The script will print:
```
🎯 ANSWER: Treated and control diverge at XXX.X minutes
```

This is the answer to the lecture question: **"earliest time where growth rate differs from control with p < 0.05"**

## Troubleshooting

**Error: "No such file or directory"**
- Edit line 171-172 in `detect_divergence.py` to match your data location

**Error: "No module named 'scipy'"**
- Run: `pip install scipy numpy matplotlib`

**No divergence detected**
- Check that mask files exist and are readable
- Verify that both REF and RIF folders have data
- Try reducing `alpha` or `min_consecutive` parameters

## What the Code Does (Simple Explanation)

1. **Loads masks** → Counts bacterial pixels in each frame
2. **Calculates growth** → Fits exponential curve in 30-min windows
3. **Compares** → t-test between treated vs control at each timepoint
4. **Finds divergence** → First time p < 0.05

That's it! Simple and exactly what the lecture asked for.
