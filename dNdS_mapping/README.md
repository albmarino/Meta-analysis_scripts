## dN/dS estimation
Sequences with more than 10% of missing information are removed from the alignments. With preal being the unaligned fasta of a gene:
```
for preal in *prealign_species.fasta; do
java -jar ~/bin/macse_v2.06.jar -prog trimNonHomologousFragments -seq $preal -min_internal_homology_to_keep_seq 0.9 -out_trim_info "$clade"_"$gene"_stats.csv
done
ls *stats.csv > statlist
ls $1*aln > alnlist

for statfile in $(cat statlist); do
statgene=${statfile%_stats.csv}
if grep -q $statgene alnlist; then
grep true $statfile | cut -d';' -f1 > sp2keep.txt
seqkit grep -f sp2keep.txt "$statgene"_final_mask_align_NT.aln -o falserm_"$statgene"_final_mask_align_NT.aln
fi
done
```

bppml and mapnh (Dutheil et al., 2006; Guéguen et al., 2013; Romiguier et al., 2012; Guéguen & Duret, 2018) were run using the filtered alignments and the clade phylogenies:
```
for aln in falserm*NT.aln; do
# prune tree by alignment IDs
excludeSpeciesTree.R ${aln} ${clade}_out_107.fst.treefile

gene=${aln%_final_mask_align_NT.aln}
sed "s/aln_filename/$aln/g" bppML_V3.bpp > bppML_currentgene.bpp

# Fit a homogenous model over current gene for all branches
~/bin/MapNH_bppV3/bppML param=bppML_currentgene.bpp
sed "s/aln_filename/$aln/g" mapNH_V3.bpp > mapNH_currentgene.bpp

# Map substitutions for all categories
sed -i 's/curr_maptype/DnDs/g' mapNH_currentgene.bpp
~/bin/MapNH_bppV3/mapnh param=mapNH_currentgene.bpp
rename "s/^/${gene}_/" CM_splitcounts_*

done
```

The gene trees output by bppml are used as input to find deviant genes, such as: 
```
cat *ml_h.dnd_1 > trees
ls *ml_h.dnd_1 > treenames
```

For each clade outlier genes are identified with PhylteR:
```
library(phylter)
library(ape)
trees <- read.tree("trees")
names <- read.table("treenames", header=F)$V1
results_k2 <- phylter(trees, gene.names = names, k=2)
```
Outlier genes are later excluded from the dN/dS calculation of every species of the given clade, as well as genes with terminal branch lenghts shorter than 0.001, ${clade}_genes2keep being the list of non-outlier genes 
```
Rscript dndscalc_CM_noshortbranches_subset.R ${clade}_out_107.fst.treefile ${clade}_genes2keep ${clade}_dNdSCM_phylter_noshortbranches.tsv
```

## References

Dutheil, J., Gaillard, S., Bazin, E., Glémin, S., Ranwez, V., Galtier, N., & Belkhir, K. (2006). Bio++: a set of C++ libraries for sequence analysis, phylogenetics, molecular evolution and population genetics. BMC Bioinformatics, 7(1), 188. doi: 10.1186/1471-2105-7-188

Guéguen, L., & Duret, L. (2018). Unbiased Estimate of Synonymous and Nonsynonymous Substitution Rates with Nonstationary Base Composition. Molecular Biology and Evolution, 35(3), 734–742. doi: 10.1093/molbev/msx308

Guéguen, L., Gaillard, S., Boussau, B., Gouy, M., Groussin, M., Rochette, N. C., Bigot, T., Fournier, D., Pouyet, F., Cahais, V., Bernard, A., Scornavacca, C., Nabholz, B., Haudry, A., Dachary, L., Galtier, N., Belkhir, K., & Dutheil, J. Y. (2013). Bio++: Efficient Extensible Libraries and Tools for Computational Molecular Evolution. Molecular Biology and Evolution, 30(8), 1745–1750. doi: 10.1093/molbev/mst097

Lartillot, N., & Poujol, R. (2011). A Phylogenetic Model for Investigating Correlated Evolution of Substitution Rates and Continuous Phenotypic Characters. Molecular Biology and Evolution, 28(1), 729–744. doi: 10.1093/molbev/msq244

Ranwez, V., Douzery, E. J. P., Cambon, C., Chantret, N., & Delsuc, F. (2018). MACSE v2: Toolkit for the Alignment of Coding Sequences Accounting for Frameshifts and Stop Codons. Molecular Biology and Evolution, 35(10), 2582–2584. doi: 10.1093/molbev/msy159

Romiguier, J., Figuet, E., Galtier, N., Douzery, E. J. P., Boussau, B., Dutheil, J. Y., & Ranwez, V. (2012). Fast and Robust Characterization of Time-Heterogeneous Sequence Evolutionary Processes Using Substitution Mapping. PLOS ONE, 7(3), e33852. doi: 10.1371/journal.pone.0033852

Scornavacca, C., Belkhir, K., Lopez, J., Dernat, R., Delsuc, F., Douzery, E. J. P., & Ranwez, V. (2019). OrthoMaM v10: Scaling-Up Orthologous Coding Sequence and Exon Alignments with More than One Hundred Mammalian Genomes. Molecular Biology and Evolution, 36(4), 861–862. doi: 10.1093/molbev/msz015

Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE, 11(10), e0163962. doi: 10.1371/journal.pone.0163962
