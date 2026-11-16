# Status Signals in the Manosphere — Reproducibility Materials

This repository contains the code, transcribed published summary statistics, and figures accompanying a methodological critique of Date Psychology's 2024 essay, [“Signaling Status in the Manosphere”](https://datepsychology.com/signaling-status-in-the-manosphere/).

**Scope:** These are *reconstructions from rounded, published statistics*, not an independent replication using participant-level data. The original survey microdata and analysis code were not available for this audit. The inferential assumptions, notably for the meta-analysis, should be interpreted accordingly.

## Contents

| File | Purpose |
| --- | --- |
| `reproduce_status_signals_figures_v4.py` | Python script that computes the statistics and regenerates the four charts. |
| `meta_reconstruction_inputs.csv` | Transcribed correlations and 95% confidence limits from the ten studies shown in the essay's forest plot, plus reconstructed random-effects weights. |
| `dunn_pairwise_recalculation.csv` | Alexander's published (rounded) z statistics and p-values alongside one-tail and conventional two-sided Bonferroni recalculations. |
| `reconstruction_results_v4.txt` | Reference output for the numerical calculations. |
| `figure_01_meta_analysis_v4.png` | Forest plot with reconstructed weights, pooled estimate, and prediction interval. |
| `figure_02_missing_displays_v4.png` | Hypothetical examples showing how the evidentiary value of *no public display* depends on base rates. These are **not observed data**. |
| `figure_03_relationship_display_v4.png` | Reconstructed differences and approximate confidence intervals from Krueger & Forest (2020), Study 2. |
| `figure_04_belief_correlations_v4.png` | Three reported/reconstructed correlation-scale estimates from the August 2024 survey. Dashed ±.18 lines indicate approximate 80%-power detectability, **not practical importance**. |
| `requirements.txt` | Python dependencies. |

## Reproduce the calculations and figures

Use Python 3.10 or newer:

```bash
python -m pip install -r requirements.txt
python reproduce_status_signals_figures_v4.py
```

The script prints the reconstructed numerical results and writes the PNGs and CSVs in the same folder. Files with the listed names will be overwritten, so run it in a working copy if you want to preserve the checked-in figures.

## Methods and interpretation

- The forest plot uses published rounded correlations and confidence intervals, Fisher-*z* transformations, and a DerSimonian–Laird random-effects variance estimate. Its prediction interval uses a *t* reference with `k - 2` degrees of freedom. These are descriptive reconstructions, not a systematic meta-analysis or an exact reanalysis of source datasets.
- The Dunn comparison checks use published rounded *z* statistics; they calculate one-tail and conventional two-sided p-values, with Bonferroni adjustment for three comparisons. The original statistical software and directional-testing intentions cannot be verified from the companion post alone.
- The experimental difference confidence intervals use reported group mean differences, test statistics, and degrees of freedom from Krueger & Forest (2020). These are reconstructed approximations, not author-provided confidence intervals.
- The missing-display graphic uses explicitly hypothetical rates to explain a likelihood-ratio concept. It is **not** an estimate of any influencer population.

## Principal sources

- Date Psychology, [“Signaling Status in the Manosphere” (2024)](https://datepsychology.com/signaling-status-in-the-manosphere/).
- Date Psychology, [“How Attractive Are Red Pill Influencers? Ratings From a Large Sample” (2024)](https://datepsychology.com/how-attractive-are-red-pill-influencers-ratings-from-a-large-sample/).
- Krueger, K. L., & Forest, A. L. (2020). [“Communicating Commitment: A Relationship-Protection Account of Dyadic Displays on Social Media.”](https://doi.org/10.1177/0146167219893998)

The companion critique provides the full reference list and detailed interpretation. No original copyrighted study PDFs, third-party webpage copies, or participant-level records are redistributed here.
