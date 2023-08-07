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

