# ASMOT: Auditing Synthetic Metagenomic Data with Multi-Scale Optimal Transport

**ASMOT** is an R package for auditing the fidelity of synthetic metagenomic data using **Multi-Scale Optimal Transport**. It provides a principled, geometry-aware framework to compare real and synthetic microbiome datasets across marginal, joint, and structural levels of dependency.

Unlike standard ecological metrics (Bray-Curtis, Alpha Diversity) that treat taxa as independent features, **ASMOT** quantifies realism as a **transport cost**:

> *How much effort is required to morph the synthetic distribution into the real one, given the phylogenetic and ecological structure?*

---

## Key Features

### 🔹 Multi-Scale Auditing

ASMOT decomposes "realism" into three mathematically distinct layers:

| Level | Question | Method |
| --- | --- | --- |
| **Univariate** | Do individual taxon abundances match? | 1D Wasserstein () |
| **Joint** | Is the sample geometry (composition) preserved? | Unbalanced Optimal Transport (UOT) |
| **Structural** | Is the ecological network topology preserved? | Gromov-Wasserstein (GW) |

### 🔹 Adaptive Weighting

Realism scores are not arbitrary. ASMOT learns component weights from the **null distribution variance**, automatically prioritizing audit levels that are most stable and discriminative for your specific dataset.

### 🔹 High-Order Interaction Scanning

A Monte-Carlo probe of -dimensional taxon subsets. This detects "mode collapse" or "spurious correlations" in complex motifs (e.g.,  triplets) that simple pairwise correlations miss.

### 🔹 Benchmarking Ready

Includes built-in datasets (`synth_zinb`, `synth_htln`) to immediately benchmark different generative models against ground truth data (`amgut1.filt`).

---

## Installation

Install directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("FloAI/ASMOT")

```

**Dependencies:** `transport`, `T4transport`, `randomForest`, `phangorn`, `ggplot2`.

---

## Quick Start: Benchmarking Two Models

This example demonstrates how to compare a **ZINB** (Zero-Inflated Negative Binomial) baseline against a newer **HTLN** model using the package's built-in data.

### 1. Load Data & Create Objects

```r
library(ASMOT)

# Load built-in datasets
data("amgut1.filt")   # Real Ground Truth
data("synth_zinb")    # Baseline Model
data("synth_htln")    # New Model

# Create Audit Objects for both models
# We compare both against the same real data
obj_zinb <- create_asmot(real = amgut1.filt, synth = synth_zinb)
obj_htln <- create_asmot(real = amgut1.filt, synth = synth_htln)

```

### 2. Run the Statistical Audit

Use `weights="auto"` to let the data determine the most important metric, or fixed weights for a strict comparison.

```r
# Run Audits (B=50 bootstraps)
# Fixed weights allow for a fair, apples-to-apples comparison
weights <- c(0.2, 0.4, 0.4) 

res_zinb <- asmot_audit(obj_zinb, B = 50, weights = weights)
res_htln <- asmot_audit(obj_htln, B = 50, weights = weights)

# Print Scores
print(paste("ZINB Realism Score:", res_zinb$realism_score))
print(paste("HTLN Realism Score:", res_htln$realism_score))

```

### 3. Visual Benchmarking

Visualize the "Realism Gap" between the two models.

```r
library(ggplot2)

df <- data.frame(
  Model = c("ZINB", "HTLN"),
  Realism = c(res_zinb$realism_score, res_htln$realism_score)
)

ggplot(df, aes(x = Model, y = Realism, fill = Model)) +
  geom_bar(stat = "identity", width = 0.5) +
  ylim(0, 100) +
  labs(title = "Benchmark: ZINB vs HTLN", y = "ASMOT Realism Score (%)") +
  theme_minimal()

```

### 4. Scan High-Order Interactions

Check if the models fail at complex motifs (e.g.,  triplets).

```r
# Scan dimensions k=1 to k=5
scan <- asmot_plot_scan(obj_htln, max_k = 5)
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

##References
>

## License

This project is licensed under the MIT License.
