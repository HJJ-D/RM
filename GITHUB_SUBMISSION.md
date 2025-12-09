# GitHub Submission Guide

## 📦 What to Include in Repository

### Essential Files (Upload these)
```
✓ detect_divergence.py       # Main analysis script
✓ README.md                   # Full documentation
✓ QUICK_START.md             # Quick instructions for classmates
✓ .gitignore                 # Git ignore rules
✓ requirements.txt           # Python dependencies
✓ example_output/            # Sample results (optional)
    ├── divergence_analysis.png
    ├── divergence_report.txt
```

### Data Files (Do NOT upload - too large)
```
✗ REF_masks101_110/          # ~500MB+ (exclude with .gitignore)
✗ RIF10_masks201_210/        # ~500MB+ (exclude with .gitignore)
✗ *.tif files                # Individual mask images
```

### Pickle Files (Optional - small and useful)
```
? growth_areas.pickle         # Small files (~1KB each)
? growth_rate.pickle          # Can include for faster loading
```

---

## 🚀 Setup Instructions for GitHub

### Step 1: Initialize Repository
```bash
cd F:\RM
git init
git add detect_divergence.py README.md QUICK_START.md .gitignore requirements.txt
git commit -m "Initial commit: Bacterial growth divergence detection"
```

### Step 2: Create GitHub Repository
1. Go to https://github.com/new
2. Repository name: `bacterial-growth-divergence`
3. Description: "Detect when RIF-treated bacteria diverge from untreated controls (Lecture 3 Project)"
4. Keep it **Public** (for classmate access)
5. Do NOT initialize with README (you already have one)

### Step 3: Push to GitHub
```bash
git remote add origin https://github.com/YOUR_USERNAME/bacterial-growth-divergence.git
git branch -M main
git push -u origin main
```

---

## 📊 Data Sharing Strategy

Since data files are too large for GitHub, use one of these approaches:

### Option 1: Provide Data Location (Recommended)
Add to your README:
```markdown
## Data Access

The analysis requires preprocessed segmentation masks from Lecture 3:
- **Control (REF)**: Positions 101-110 from `REF_masks101_110` folder
- **Treated (RIF10)**: Positions 201-210 from `RIF10_masks201_210` folder

Data location: [Instructor's shared folder] or contact me for access.
```

### Option 2: Small Sample Dataset
Include just ONE position as example:
```bash
# Create sample data folder
mkdir -p sample_data/REF/Pos101/PreprocessedPhaseMasks
mkdir -p sample_data/RIF10/Pos201/PreprocessedPhaseMasks

# Copy only pickle files (small)
cp REF_masks101_110/Pos101/PreprocessedPhaseMasks/growth_areas.pickle sample_data/REF/Pos101/PreprocessedPhaseMasks/
cp RIF10_masks201_210/Pos201/PreprocessedPhaseMasks/growth_areas.pickle sample_data/RIF10/Pos201/PreprocessedPhaseMasks/

# Add to git
git add sample_data/
git commit -m "Add sample data for testing"
```

### Option 3: Google Drive / Dropbox
1. Upload full data to cloud storage
2. Create shareable link
3. Add link to README

---

## 📝 README Additions for Classmates

Add this section to your README.md:

```markdown
## 🎓 For Classmates (Peer Evaluation)

### My Results
- **Dataset Used**: 
  - REF (Control): Positions 101-110 (10 replicates)
  - RIF10 (Treated): Positions 201-210 (10 replicates)
- **Detected Divergence**: XX.X minutes (X.XX hours)
- **P-value**: 0.0XXX
- **Method**: Exponential growth rate + independent t-test

### How to Reproduce My Results

1. **Get the data** (see Data Access section)
2. **Install dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
3. **Update paths** in `detect_divergence.py` (lines 187-188):
   ```python
   REF_BASE = Path("YOUR_PATH/REF_masks101_110")
   RIF_BASE = Path("YOUR_PATH/RIF10_masks201_210")
   ```
4. **Run**:
   ```bash
   python detect_divergence.py
   ```
5. **Check outputs**: 
   - `divergence_analysis.png` 
   - `divergence_report.txt`

### Expected Runtime
- ~10-30 seconds (with pickle files)
- ~2-5 minutes (without pickle, calculating from masks)

### Troubleshooting
- **Qt error**: Code uses Agg backend, should work automatically
- **File not found**: Check data paths in lines 187-188
- **Module error**: Run `pip install -r requirements.txt`
```

---

## ✅ Final Checklist Before Submission

- [ ] Code runs successfully on your machine
- [ ] README.md includes clear instructions
- [ ] QUICK_START.md for fast setup
- [ ] requirements.txt lists all dependencies
- [ ] .gitignore excludes large data files
- [ ] Example output images included (optional)
- [ ] Repository is public
- [ ] Test: Clone to new folder and verify it works
- [ ] Data access method specified (shared folder / contact info)

---

## 🔗 Share with Classmates

Once uploaded, share:
```
Repository: https://github.com/YOUR_USERNAME/bacterial-growth-divergence
Quick Start: https://github.com/YOUR_USERNAME/bacterial-growth-divergence/blob/main/QUICK_START.md
```

---

## 📧 Contact Template for Classmates

```
Hi [Classmate Name],

For peer evaluation, here's my Lecture 3 project:

Repository: https://github.com/YOUR_USERNAME/bacterial-growth-divergence

Quick instructions:
1. Clone the repo
2. Get data from [instructor's folder / my Google Drive]
3. Update paths in detect_divergence.py (lines 187-188)
4. Run: python detect_divergence.py

My results: Divergence detected at XX.X minutes (p=0.0XX)

Let me know if you need help!

[Your Name]
```
