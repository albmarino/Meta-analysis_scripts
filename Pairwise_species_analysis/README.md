# Pairwise comparison of closely related species
The hypothesis that Ne differences can affect recent changes of GS and TE content is tested on independent species couples extracted from a subset of species with contrasted GS.
```
Rscript deltaGS_pairwise_test_plot.R FULLTREE_brlengths_rooted.treefile table_coevol_final.tsv 10 0.04 100 gcpoor_dsom.postmean1.tab
```

The script defines for each species a ΔdNdS, a ΔGS%, a ΔTE%, a ΔrecentTE%, and pairwise branch lengths from all the other species (brlen_pw):

- ΔdNdS is the simple pairwise difference between its dN/dS and the dN/dS of another species - the dN/dS estimated by Coevol for the GC3-poor geneset are used (extracted from `table_coevol_final.tsv`)
- ΔGS% corresponds to the pairwise difference between its GS and the GS of another species as a percentage of the GS inferred for their most recent common ancestor - the ancestral GS estimated by Coevol at internal nodes are used (extracted from Coevol outputs with suffix *gcpoor_dsom.postmean1.tab)
- ΔTE% corresponds to the pairwise difference between its overall TE content and the overall TE content of another species as a percentage of the GS inferred for their most recent common ancestor
- ΔrecentTE% corresponds to the pairwise difference between its recent TE content (sequence divergence from the consensus sequence lower than 5%) and the recent TE content of another species as a percentage of the GS inferred for their most recent common ancestor

The species with ΔGS% ≥ 10 (3rd argument) and brlen_pw ≤ 0.04 (4th argument) are retained as closely related taxa with sufficiently contrasted recent GS divergences.
From this set, independent species pairs are randomly chosen to test the association between ΔdNdS and ΔGS%, ΔTE% and ΔrecentTE% (chi-squared test). If ΔdNdS is positive, positive ΔGS%, ΔTE%, ΔrecentTE% should be expected: namely, an Ne lower than the one of a close species is expected to associate with a positive change in GS, in TE content or recent TE content compared to the species of comparison. The test is conducted on 100 (5th argument) randomized extractions.

......
