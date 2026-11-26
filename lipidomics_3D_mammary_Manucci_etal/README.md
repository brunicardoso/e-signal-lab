# Differential lipidomic profiles of 3D mammary cultures in response to the basement membrane

Code and input data for the analyses and figure panels of:

> Manucci AC, Yoshinaga MY, Inague A, Miyamoto S, Bruni-Cardoso A. *Differential lipidomic profiles of 3D mammary cultures in response to the basement membrane.* (submitted)

EpH4 mammary epithelial cells grown as 3D clusters were left non-treated (nT) or treated with a reconstituted basement membrane (rBM) for 48 h (lipidomics) or 0–72 h (western blots).

## Run

```bash
python3.12 -m pip install -r requirements.txt
bash run_all.sh            # writes everything to outputs/ (about 30 s)
```

## Figure → script → input

| Panel | Script | Input (`data/`) | Analysis |
|---|---|---|---|
| Fig. 1D | `fig1D_edu_ttest.py` | `edu_quantification.csv` | EdU-positive nuclei; two-tailed unpaired Student's t-test |
| Fig. 2B, Suppl. Table S1 | `fig2B_wb_stats.py`, `fig2B_wb_barplots.py`, `fig2B_wb_stats_table.py` | `western_blot/*_CellBioStats_version.xlsx` | z-scores; repeated-measures ANOVA, paired t-tests, Benjamini–Hochberg FDR |
| Fig. 3B | `fig3B_lipid_diversity.py` | `lipidomics_cells_DNA_normalized.csv` | species per class |
| Fig. 3C–D | `fig3CD_oplsda.py` | `lipidomics_cells_DNA_normalized.csv`, `metaboanalyst_oplsda_loadings.csv` | OPLS-DA score plot and S-plot; R²Y, Q², permutation test |
| Fig. 4A | `fig4A_class_log2FC.py` (`replicas` mode) | `lipidomics_cells_DNA_normalized.csv` | class sums; t-test on log10, BH FDR |
| Fig. 4B | `fig4B_fa_remodeling.py` | `lipidomics_cells_DNA_normalized.csv`, `fatty_acid_species_map.xlsx` | acyl-chain composition and abundance; t-tests, BH FDR |
| Fig. 5A | `fig5A_volcano.py` | `lipidomics_cells_DNA_normalized.csv` | species t-tests on log10, BH FDR; writes the statistics CSV |
| Fig. 5B | `fig5B_heatmap.py` | `lipidomics_cells_DNA_normalized.csv` + Fig. 5A CSV | z-scored heatmap, hierarchical clustering |

Fig. 6 (lipid ontology enrichment) was produced with the [LION/web](http://www.lipidontology.com) tool, not with code in this folder.

## Data

- `lipidomics_cells_DNA_normalized.csv`: 211 lipid species (columns) × 12 samples. Concentrations normalized to internal standards and DNA content.
- `fatty_acid_species_map.xlsx`: the sheet `Mapa espécies` maps each species to its acyl chains. The other sheets hold the original spreadsheet calculations of Supplementary Table S4.
- `metaboanalyst_oplsda_loadings.csv`: the original MetaboAnalyst S-plot export, used to check the Python OPLS-DA.
- `western_blot/`: densitometry normalized to Ponceau S (`Measument`), per condition and replicate. `Z_score` = (x − mean) / SD within each replicate (ddof = 1).
- `edu_quantification.csv`: % EdU-positive nuclei per field.

## Notes

- **Sample labels:** `lrECM` in the lipidomics data is the original label for rBM. `Pool` samples are pooled QCs and are excluded from all statistics.
- **OPLS-DA:** settings were recovered from the MetaboAnalyst export (log2, auto-scaling, 1 predictive + 1 orthogonal component). They reproduce R²X = 28.8 % / 24.5 % and the S-plot (206 of 211 lipids within 1e-4). MetaboAnalyst's random jitter on p[1] is not reproduced.
- **Western-blot output:** statistics files are timestamped. The plots and table use the newest file in `outputs/fig2B/`.

## Software

Python 3.12.16, NumPy 2.5.5, SciPy 1.18.1, pandas 3.1.1, Matplotlib 3.11.2, seaborn 0.13.2, openpyxl 3.1.5, statsmodels 0.15.0.

## License

Code: BSD 3-Clause (repository `LICENSE`). Data: CC BY 4.0.
