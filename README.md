
# ASMOT: Auditing Synthetic Metagenomic Data with Multi-Scale Optimal Transport

**ASMOT** is an R package for auditing the fidelity of synthetic metagenomic data using **Multi-Scale Optimal Transport**. It provides a principled, geometry-aware framework to compare real and synthetic microbiome datasets across marginal, joint, and structural levels of dependency.

Unlike standard ecological metrics (Bray-Curtis, Alpha Diversity) that treat taxa as independent features, **ASMOT** quantifies realism as a **transport cost**:

> *How much effort is required to morph the synthetic distribution into the real one, given the phylogenetic and ecological structure?*

---

## Key Features

### Multi-Scale Auditing

ASMOT decomposes "realism" into three mathematically distinct layers:

| Level | Question | Method |
| --- | --- | --- |
| **Univariate** | Do individual taxon abundances match? | 1D Wasserstein () |
| **Joint** | Is the sample geometry (composition) preserved? | Unbalanced Optimal Transport (UOT) |
| **Structural** | Is the ecological network topology preserved? | Gromov-Wasserstein (GW) |

###  Adaptive Weighting 

Realism scores are not arbitrary. ASMOT learns component weights from the **null distribution variance**, automatically prioritizing audit levels that are most stable and discriminative for your specific dataset.

###  High-Order Interaction Scanning

A Monte-Carlo probe of -dimensional taxon subsets using **phylogenetically grounded OT**. This detects "mode collapse" or "spurious correlations" in complex motifs (e.g.,  triplets) that simple pairwise correlations miss.

###  Zero-Inflation Robustness

Native support for **Unbalanced Optimal Transport**, allowing the metric to handle the extreme sparsity of microbiome data without artificial imputation.

###  "Turing Test" via Classifiers

Includes an adversarial audit module (`asmot_classifier`) that trains a Random Forest to distinguish real vs. synthetic samples based on their OT geometry.

---

## Installation

Install directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("FloAI/ASMOT")

```

**Dependencies:** `transport`, `T4transport`, `randomForest`, `phangorn`, `ggplot2`.

---

## Quick Start

### 1. Preprocess & Initialize

ASMOT never applies hidden transformations. You explicitly choose your preprocessing (e.g., TSS for biology, CLR for correlations).

```r
library(ASMOT)

# Load your data (Rows = Samples, Cols = Taxa)
# Preprocess explicitly (optional but recommended)
real_data <- asmot_preprocess("real_counts.csv", method = "tss")
synth_data <- asmot_preprocess("synth_counts.csv", method = "tss")

# Create the Audit Object
# (Optionally pass a phylogenetic tree for ground-truth costs)
obj <- create_asmot(real = real_data, synth = synth_data)

```

### 2. Run the Statistical Audit

Use **Adaptive Weighting** (`weights="auto"`) to let the data determine the most important metric.

```r
# Run audit with 100 bootstraps
results <- asmot_audit(obj, B = 100, weights = "auto")

# Print S3 Summary
print(results)
# Output:
# Realism Score:   82.4%
# P-value:         0.042 (Significant deviation detected)
# Component Scores: ...

```

### 3. Visualize the "Realism Gap"

Generate a dashboard showing where the synthetic data falls relative to the natural variation of the real data.

```r
plots <- asmot_plot_dashboard(results)
print(plots)

```

### 4. Scan High-Order Interactions

Check if the generator fails at complex motifs (e.g.,  or ).

```r
# Scan dimensions k=1 to k=6
scan <- asmot_plot_scan(obj, max_k = 6)
print(scan)

```

---

## Methodology

ASMOT defines synthetic fidelity as:

> **The minimum transport cost required to morph the synthetic distribution into the real one, evaluated across increasing orders of dependency.**

### Why Optimal Transport?

Standard metrics (like Euclidean distance or KL-Divergence) fail in high-dimensional, sparse biological data because they do not account for the **geometry** of the space.

* If a model produces *E. coli* instead of *Salmonella*, standard metrics see a "mismatch."
* **ASMOT** uses the phylogenetic tree to see a "small shift," penalizing it less than producing a distinct phylum.

### Interpretation Guide

| Metric | Result | Interpretation |
| --- | --- | --- |
| **Realism Score** | **> 90%** | Indistinguishable from real data (Publication ready). |
| **P-value** | **< 0.05** | The synthetic data has statistically significant artifacts. |
| **Structural Error** | **High** | The model captures abundances but breaks the co-occurrence network. |
| **-Scan** | **Rising Error** | The model fails to capture complex () ecological rules. |

---

## Citation

If you use **ASMOT** in your research, please cite:

> [Your Name], et al. (2024). "Auditing Synthetic Metagenomic Data with Multi-Scale Optimal Transport." *[Journal Name/Preprint]*.

## License

This project is licensed under the MIT License.
