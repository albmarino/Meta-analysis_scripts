## Traits correlation under phylogenetic non-independence with Coevol

Scripts to select GC3-poor and GC3-rich genes and to run Coevol.

#### Full geneset stats and filtering
`genespersp.sh` calculates the number of markers for each species.<br>
`gcstat.sh` reports GC3 level statistics per gene overall and at clade level.<br>
`gc3plots.R` plots GC3 levels per clade and outputs a list of species to remove.<br>

#### Extract genes post-filtering
`gcstat_filtgenespersp.sh` reports GC3 level statistics on the filtered dataset after species removal.<br>
`genespersp_filtgenespersp.sh` calculates the number of markers for each species on the filtered dataset.<br>
`gc3plots_gcrich.R` outputs lists of the top GC3-rich genes for every clade.<br>
`gc3plots.R` outputs lists of the top 50 GC3-poor genes for every clade.<br>

#### Run Coevol and read output
`concatenate_4coevol.sh` concatenates listed genes.<br>
`make_coevol_table.R` produces a matrix in input format for Coevol.<br>
`run_coevol.sh` and `read_coevol.sh` run Coevol and check output.<br>
`dnds_coevol.R` extracts Coevol dS and dN/dS.<br>



```
./genespersp.sh # outputs genespersp.tsv
./gcstat.sh # outputs gc3_report_allsp.tsv
Rscript gc3plots.R 50 50 95 gc3_report_allsp.tsv genespersp.tsv 50 # outputs sp2rm_50percgenespersp
./gcstat_filtgenespersp.sh sp2rm_50percgenespersp # outputs gc3_report_50percgenespersp.tsv
```

Only genes represented in at least 95% of the species were retained, the 50 GC3-poorest and GC3-richest genes extracted and concatenated for each clade:
```
./genespersp_filtgenespersp.sh sp2rm_50percgenespersp # outputs genespersp_50percgenespersp.tsv, a table with number of markers per species on the filtered geneset
Rscript ~/bin/gc3plots_gcrich.R 50 50 90 gc3_report_50percgenespersp.tsv genespersp_50percgenespersp.tsv # outputs lists of the top 50 GC3-rich genes for every clade
Rscript ~/bin/gc3plots.R 50 50 90 gc3_report_50percgenespersp.tsv genespersp_50percgenespersp.tsv # outputs lists of the top 50 GC3-poor genes for every clade
./concatenate_4coevol.sh gcpoorset_actinopteri_50percgenespersp_95_20_20.tsv sp2rm_50percgenespersp
```

Coevol was run for each clade with the GC3-poor and GC3-rich datasets:
```
Rscript make_coevol_table.R Actinopteri # outputs the traits matrix
~/bin/coevol -d concatenate_gcpoorset_actinopteri.fasta -t actinopteri_brlen_rooted.treefile -fixtimes -c Actinopteri_coevol.txt -dsom actinopteri_gcpoor_dsom
~/bin/readcoevol -x 400 +med +ci actinopteri_gcpoor_dsom
```

dS, dN/dS (on GC-rich and GC-poor genesets) and upstream branch lenghts of every species were extracted from `*postmeanbranchnonsynrate.tab` and `*postmeanbranchsynrate.tab` outputs:
```
Rscript dnds_coevol.R dnds_coevol.tsv
```
