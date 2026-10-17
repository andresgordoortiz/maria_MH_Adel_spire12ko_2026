# Sources and attribution

Reference analysis: https://github.com/andresgordoortiz/oocyte_splicing_ssa_pladb_tub_2025/tree/e231983fd82736971aef7fc938cf34bed5ce3121

Paper: Al Shami et al., *The post-transcriptional RNA splicing landscape driving oocyte division*, bioRxiv v1, 20 August 2026, https://doi.org/10.64898/2026.08.19.745726 . Full JATS text was retrieved and read on 16 October 2026. New plots of the Spire and FMN2 datasets and their comparison are generated here. The upstream code is MIT licensed; its notice is preserved in LICENSE.reference.

Adaptations retain the reference thresholds, row scaling, miniheatmap clustering, feature tests, GO function and GC binning. Colours now follow the shared accessible palette in palettes.R; see ../docs/COLOURS.md. Drug-specific inputs and global R objects have been replaced with validated completed simulation tables and explicit inputs. The two-study comparison adds a common tested universe and distinguishes shared targets from concordant changes. The source tree and paper XML are cached locally under ignored reference/; they are not runtime dependencies.

Corrections are explicit: empirical FPR is labelled correctly; minus-strand flanking coordinates use the adjacent exon edges; undefined correlations of constant PSI rows are excluded; zero motif p-values are retained with a plotting floor; the GC gradient's reversed intron 5′/3′ slices are corrected. The main large-heatmap code uses Euclidean/complete clustering by default, whereas the miniheatmap explicitly uses Spearman/Ward.D2; both choices are preserved.

| Here | Pinned source scripts |
|---|---|
| splicing.R | figure2C_volcano_deltapsi.R; figure2D_number_significant_events.R; figure2F_miniheatmaps.R; figure_s1b_bigheatmaps.R; percetange_inclusionsvsskipping.R |
| annotations.R | exons_gc_length.R; introns_gc_length.R; Figure3A_SuppTable5_orf_cuantificaitn.R; figure3B__Supptable6_FigureS5_go.R |
| external.R | tableformatt.R; rmaps_tablecreation.R; Preprocessing/table_matt_rbp.R |
| matt_figures.R | Figure2G_matt_feature_summary.R; FiguresS2_S3_features_matt_plots.R; FigureS4G_H_length_gcc_review.R |
| gc_gradient.R | FigureS4_GC_gradient_metaplot_mm10.R (binning functions reused) |
| Retired individual RBP summary in legacy/rmaps3/development/ | Figure2H_S4A_B.R; FigureS4C_E.R; introns_rbps_rmaps.R |
| rbp_helpers.R; rbp_paper_figures.R; render_rmaps3.R | Figure2H_S4A_B.R; FigureS4C_E.R; panelB_boxplot_rmaps_exons.R; manuscript Figure 2H and S6A–E captions |

Matt is run from the exact container digest cited in the manuscript:
`andresgordoortiz/matt-container@sha256:6b729b82e312f06153f8faeb643402df4e1729485fba90359a6cf3928c24c20e` (Matt 1.3.2, checked locally).
Official Matt source/documentation: https://gitlab.com/aghr/matt .

Reference annotations are downloaded from VASTDB and GENCODE release M10, with expected SHA-256 hashes in prepare_references.py. URLs and hashes are exported to each result manifest. The available PROT_IMPACT resource is v3; the original script's v2.3 URL is unavailable (404). This prevents a claim of exact annotation-version reproduction.

rMAPS2: https://rmaps.cecsresearch.org/ ; method: https://pmc.ncbi.nlm.nih.gov/articles/PMC7319468/ . No data are uploaded automatically.

Local motif runs use the supplied rMAPS3 3.1.0 checkout at commit
`b9a9dce9ad8596fdef21c3a705fbc5f37d9fa45d`, in its existing Pixi environment.
The user authorized correcting a demonstrated motif-counting defect in an
isolated copy. See [RMAPS3.md](RMAPS3.md) for the regression tests, exact window
definition, Mann–Whitney setting, archived Fisher baselines and limitations.
The correction is now maintained in `../rMAPS3`, local branch
`fix/motif-window-counting`, commit `1076240`; GitHub publication is left to the
user. Its README distinguishes reproducible counting failures from separate
boundary and intermediate-output choices. Run manifests hash the actual modified source rather than identifying it as an
unmodified upstream release.

The bundled catalogue contains 114 motifs / 91 RBP labels. It does not contain
the paper's manually added ENCORI SF3B1 motif; that additional motif is not
reconstructed here. These local runs are an adaptation, not an exact replication
of the historical rMAPS2 server and its additional motif.
