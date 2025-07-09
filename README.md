# Z-Score Analysis Tool for Docking Results
A Python-based tool for statistical post-processing of molecular docking results in virtual screening (VS) studies. This script computes Z-scores for binding energies and filters ligands based on user-defined significance thresholds. It also produces a summary PDF table and high-resolution visualizations suitable for publication or internal reporting.

# Description
The **Docking Score Z-Score Analysis Tool** processes docking results (in CSV format), calculates the Z-score of each compound’s binding affinity, and identifies statistically significant hits. The workflow includes:

- Calculation of Z-scores for each ligand, based on the mean and standard deviation of the complete docking score distribution.
- Selection of statistically significant compounds, defined by a user-adjustable Z-score threshold (default: Z ≤ -1.960).
- Automated generation of result files, including:
  - A CSV file containing filtered compounds and their corresponding Z-scores;
  - A structured PDF summary table for reporting or publication;
  - A high-resolution plot of the standard normal distribution with the significance threshold marked;
  - A docking score distribution plot with annotation of the final compound’s binding score.

# Suggested Z-score Thresholds
The script supports statistical filtering based on two-tailed Z-score significance levels. Depending on the desired confidence interval (CI), the corresponding absolute Z-score threshold (Z) can be selected as follows:

- For 50% confidence, use a threshold of Z ≤ -0.674
- For 75% confidence, use Z ≤ -1.150
- For 90% confidence, use Z ≤ -1.645
- For 95% confidence, use Z ≤ -1.960
- For 97% confidence, use Z ≤ -2.170
- For 99% confidence, use Z ≤ -2.576
- For 99.9% confidence, use Z ≤ -3.290

**These values correspond to the standard normal distribution under two-tailed statistical testing and can be adjusted in the script to meet specific requirements.

# Requirements & Installation

Python version:
- `Python ≥ 3.6`

Dependencies:
- `numpy`  
- `matplotlib`  
- `scipy`  
- `fpdf`

Install them via pip:

> pip install numpy matplotlib scipy fpdf

# Usage
Place the input CSV file (default name: FILE_NAME.csv) in the same directory as the script and run the following command in your terminal:

> python docking-zscore-analysis.py

**All generated artefacts - including the filtered CSV file, summary PDF report, and high-resolution PNG figures - will be automatically saved to the current working directory.

Note: The script uses a default Z-score threshold of -1.960 to select statistically significant compounds. You may edit the script to adjust this threshold as needed.

# Citation
If you use this tool in your research or publication, please cite it as follows:

Isaoğlu, M., & Durdağı, S. (2025). Z-Score Analysis Tool for Docking Results (Version 1.0) [Computer software]. https://github.com/DurdagiLab/docking-zscore-analysis
