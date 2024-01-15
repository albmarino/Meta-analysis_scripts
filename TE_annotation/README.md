## TE annotation and quantification

Scripts to annotate genome assemblies with EarlGrey and quantify TE content with dnaPipeTE.

dnaPipeTE quantifies TEs from unassembled short reads and is run in two rounds with the pipeline available at https://github.com/sigau/pipeline_dnapipe.
The full process from reads download to second clean output is done by `sradownload_dnapipe.py`.
`dnaPT_summary_overall_recent.sh` extracts the overall and recent TE content. It leverages `dnapt_recentTEs.R` which adapts part of `dnaPT_landscapes.sh` from https://github.com/clemgoub/dnaPT_utils.


## References

Baril, T., Imrie, R. M., & Hayward, A. (2022). Earl Grey: a fully automated user-friendly transposable element annotation and analysis pipeline [Preprint]. In Review. doi: 10.21203/rs.3.rs-1812599/v1

Goubert, C., Modolo, L., Vieira, C., ValienteMoro, C., Mavingui, P., & Boulesteix, M. (2015). De Novo Assembly and Annotation of the Asian Tiger Mosquito (Aedes albopictus) Repeatome with dnaPipeTE from Raw Genomic Reads and Comparative Analysis with the Yellow Fever Mosquito (Aedes aegypti). Genome Biology and Evolution, 7(4), 1192–1205. doi: 10.1093/gbe/evv050
