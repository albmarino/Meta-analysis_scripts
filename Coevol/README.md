## Traits correlation under phylogenetic non-independence
Phylogenetic generalized Least Squares was performed with ape and nlme packages:
```
library(nlme)
library(ape)
df <- read.table("table_coevol.tsv", header=T, sep="\t")
fulltree <- read.tree("FULLTREE_brlengths_rooted.treefile")
rownames(df)<- df$Species
sp2rm <- fulltree$tip.label[! fulltree$tip.label %in% rownames(df)]
prunedtree <- drop.tip(fulltree, sp2rm)
df <-  df[prunedtree$tip.label, ]

pgls_gs_overallte <- gls(log(Unified_Cvalues_methods)~log(Overall_repeat_bp), correlation = corBrownian(phy = prunedtree, form =~Species), data = df, method = "ML", na.action=na.omit)

pgls_mass_dnds <- gls(Mass_g_log~dNdS, correlation = corBrownian(phy = prunedtree, form =~Species), data = df, method = "ML", na.action=na.omit)

pgls_longevity_dnds <- gls(MaxLongevity_y_log~dNdS, correlation = corBrownian(phy = prunedtree, form =~Species), data = df, method = "ML", na.action=na.omit)

pgls_gs_dnds <- gls(log(Unified_Cvalues_methods)~dNdS, correlation = corBrownian(phy = prunedtree, form =~Species), data = df, method = "ML", na.action=na.omit)

pgls_overallte_dnds <- gls(log(Overall_repeat_bp)~dNdS, correlation = corBrownian(phy = prunedtree, form =~Species), data = df, method = "ML", na.action=na.omit)

pgls_recentte_dnds <- gls(log(Overall_recent_bp)~dNdS, correlation = corBrownian(phy = prunedtree, form =~Species), data = df, method = "ML", na.action=na.omit)

pgls_gs_mass <- gls(log(Unified_Cvalues_methods)~Mass_g_log, correlation = corBrownian(phy = prunedtree, form =~Species), data = df, method = "ML", na.action=na.omit)

pgls_gs_longevity <- gls(log(Unified_Cvalues_methods)~MaxLongevity_y_log, correlation = corBrownian(phy = prunedtree, form =~Species), data = df, method = "ML", na.action=na.omit)

```

Bayesian inference of traits and substitution rates covariation was conducted with Coevol 1.6 (Lartillot & Poujol, 2011) on each clade separately, using a set of 50 clade-specific genes, both with a low and high GC content at the third codon position (GC3 content).
The average GC3 content for each gene and the number of markers available for each species were calculated, and species with less than 50% of single-copy markers were excluded:
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

Coevol was run for each clade with the GC3-poor and GC3-rich datasets as follows:
```
Rscript make_coevol_table.R Actinopteri # outputs the traits matrix
~/bin/coevol -d concatenate_gcpoorset_actinopteri.fasta -t actinopteri_brlen_rooted.treefile -fixtimes -c Actinopteri_coevol.txt -dsom actinopteri_gcpoor_dsom
~/bin/readcoevol -x 400 +med +ci actinopteri_gcpoor_dsom
```
