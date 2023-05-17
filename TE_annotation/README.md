## TE annotation of genome assemblies
Earl Grey 1.3 (Baril et al., 2022) was used to annotate 29 dipteran assemblies with the command below:
```
earlGrey -g GCA_001542645.1.fa -s Anopheles_gambiae -r metazoa -o ./Anopheles_gambiae_earlgrey -t 8
```

## TE annotation from unassembled reads
The pipeline available at https://github.com/sigau/pipeline_dnapipe is based on two rounds of dnaPipeTE (Goubert et al., 2015) and was used to quantify TEs from short reads. The full process from reads download to second clean output can be automatized with
```
python3 sradownload_dnapipe.py coevol_table.tsv
```
where `coevol_table.tsv` is Supplementary Table 2 (will work only on species with available reads experiment ID).
`dnaPT_summary_overall_recent.sh` was used to extract the overall and recent TE content. It leverages `dnapt_recentTEs.R` which adapts part of `dnaPT_landscapes.sh` from https://github.com/clemgoub/dnaPT_utils. The overall and recent TE content below 5% divergence were obtained with:
```
./dnaPT_summary_overall_recent.sh 5 recentTEs_5perc.tsv overallTEs.tsv
```
