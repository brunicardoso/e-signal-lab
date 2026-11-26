#!/usr/bin/env bash
# Regenerate every analysis and figure panel into outputs/.
# Requires Python 3.12 with the packages in requirements.txt.
set -euo pipefail
cd "$(dirname "$0")"
D=data
O=outputs
mkdir -p "$O"/{fig1D,fig2B,fig3B,fig3CD,fig4A,fig4B,fig5A,fig5B}
PY=${PYTHON:-python}

echo "== Fig. 1D  EdU t-test";             $PY scripts/fig1D_edu_ttest.py "$D/edu_quantification.csv" "$O/fig1D/edu_ttest.txt"
echo "== Fig. 2B  western-blot statistics"; $PY scripts/fig2B_wb_stats.py --data-dir "$D/western_blot" --out-dir "$O/fig2B"
echo "== Fig. 2B  bar plots + table"
(cd scripts && WB_DATA_DIR="../$D/western_blot" WB_OUT_DIR="../$O/fig2B" $PY fig2B_wb_barplots.py \
            && WB_DATA_DIR="../$D/western_blot" WB_OUT_DIR="../$O/fig2B" $PY fig2B_wb_stats_table.py > /dev/null)
echo "== Fig. 3B  lipid diversity";        $PY scripts/fig3B_lipid_diversity.py "$D/lipidomics_cells_DNA_normalized.csv" "$O/fig3B/lipid_diversity"
echo "== Fig. 3C-D OPLS-DA";               $PY scripts/fig3CD_oplsda.py "$D/lipidomics_cells_DNA_normalized.csv" "$O/fig3CD/oplsda" "$D/metaboanalyst_oplsda_loadings.csv"
echo "== Fig. 4A  class sums";             $PY scripts/fig4A_class_log2FC.py "$D/lipidomics_cells_DNA_normalized.csv" "$O/fig4A/classes_log2FC" replicas
echo "== Fig. 4B  acyl-chain remodeling";  $PY scripts/fig4B_fa_remodeling.py "$D/lipidomics_cells_DNA_normalized.csv" "$D/fatty_acid_species_map.xlsx" "$O/fig4B/FA_remodeling.pdf"
echo "== Fig. 5A  volcano";                $PY scripts/fig5A_volcano.py "$D/lipidomics_cells_DNA_normalized.csv" "$O/fig5A/volcano_FDR010_logratio"
echo "== Fig. 5B  heatmap";                $PY scripts/fig5B_heatmap.py "$D/lipidomics_cells_DNA_normalized.csv" "$O/fig5A/volcano_FDR010_logratio.csv" "$O/fig5B/heatmap_significant.png"
echo "Done: outputs/"
