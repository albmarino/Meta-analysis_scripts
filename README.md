# Meta-analysis scripts
Commands and scripts for code reproducibility

## Assessing assembly quality
Annotate metazoan-conserved genes from BUSCO 5.2.2 database with Metaeuk gene predictor (Manni et al., 2021):
```
busco -i assemblies/ -o buscoout -l metazoa_odb10 --update-data -m genome --cpu 30
```
Calculate assembly metrics with Quast 5.0.2 (Gurevich et al., 2013):
```
python3 ~/softwares/quast-5.0.2/quast.py assembly.fa -e -s -o assembly.fa_quast_report
```

## Extracting single copy BUSCO genes
The nucleotide sequences of single copy genes were extracted from BUSCO outputs with `faatofa_singlecopy_busco_editgenecol.py`.

## Estimating expected C-values
table3 being Supplementary Table 3 and table1 being Supplementary Table 1, Weighted Least Square was calculated from C-values and assembly sizes in Supplementary Table 3 as follows:
```
library(dplyr)
m1 <- lm(Cvalue~Assembly_size, data=table3)
wt_m1 <- 1 / lm(abs(m1$residuals) ~ m1$fitted.values)$fitted.values^2
wls_m1 <- lm(Cvalue ~ Assembly_size, data = table3, weights=wt_m1)

#Predict values from assembly sizes with Quast_ContigN50 >= 50kb in table1
table1 <- filter(table1, Quast_ContigN50 >= 50kb)
as <- data.frame(Assembly_size = table1$Assembly_size)
corrected_as <- data.frame(Expected_Cvalue = predict(wls_m1,newdata = as))
```

## Testing reads effect
Method 1: z-test of the difference between the regression coefficients of WLS based on only SR and only LR
```
df_SR <- filter(table3, Reads_assembly %in% "SR")
df_LR <- filter(table3, Reads_assembly %in% c("LR", "LR-SR"))
m1_SR <- lm(Cvalue~Assembly_size, data=df_SR)
m1_LR <- lm(Cvalue~Assembly_size, data=df_LR)
wt_m1_SR <- 1 / lm(abs(m1_SR$residuals) ~ m1_SR$fitted.values)$fitted.values^2
wt_m1_LR <- 1 / lm(abs(m1_LR$residuals) ~ m1_LR$fitted.values)$fitted.values^2

wls_m1_SR <- lm(Cvalue ~ Assembly_size, data = df_SR, weights=wt_m1_SR)
wls_m1_LR <- lm(Cvalue ~ Assembly_size, data = df_LR, weights=wt_m1_LR)

compare.coeff <- function(coef_SR,se_SR,coef_LR,se_LR){
 return((coef_SR-coef_LR)/sqrt(se_SR^2+se_LR^2))
 }
  
coef_SR <- summary(wls_m1_SR)$coefficients[2,1]
se_SR <- summary(wls_m1_SR)$coefficients[2,2]
coef_LR <- summary(wls_m1_LR)$coefficients[2,1]
se_LR <- summary(wls_m1_LR)$coefficients[2,2]
  
pval <- 2*pnorm(-abs(compare.coeff(coef_SR,se_SR,coef_LR,se_LR)))
pval # 0.8127543
```

Method 2: test the effect of reads category (LR and SR) on the WLS slope
```
newheader <- "Reads_category"
newdf <- data.frame(matrix(nrow=0, ncol=length(newheader)))
for (i in table3$Reads_assembly) { 
 if (i=="LR" || i == "LR-SR") {newdf <- rbind(newdf, "LR")}
 else if (i=="SR") {newdf <- rbind(newdf, "SR")}
 }
colnames(newdf) <- newheader
table3 <- cbind(table3, newdf)
m1_readscombined <- lm(Cvalue ~ Assembly_size + Reads_category + Assembly_size:Reads_category, data=table3))
wt_m1_readscombined <- 1 / lm(abs(m1_readscombined$residuals) ~ m1_readscombined$fitted.values)$fitted.values^2
wls_m1_readscombined <- lm(Cvalue ~ Assembly_size + Reads_category + Assembly_size:Reads_category, data = table3, weights=wt_m1_readscombined)
summary(wls_m1_readscombined)

######################
Coefficients:
    Estimate   Std. Error t value  Pr(>|t|)    
(Intercept)                    -1.655e+07  1.281e+07  -1.292   0.1971    
Assembly_size                   1.214e+00  2.159e-02  56.222   <2e-16 ***
Reads_categorySR                5.468e+07  2.064e+07   2.649   0.0084 ** 
Assembly_size:Reads_categorySR -5.914e-03  3.246e-02  -0.182   0.8555
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Residual standard error: 1.441 on 384 degrees of freedom
Multiple R-squared:  0.9368, Adjusted R-squared:  0.9363 
F-statistic:  1898 on 3 and 384 DF,  p-value: < 2.2e-16
####################
```

## Gene alignments and phylogeny
The 954 BUSCO genes were grouped according to the taxa Aves, Actinopteri, Insecta, Mammalia, Mollusca. Alignments and phylogenies were first computed clade-wise with the OMM_MACSE pipeline (version 11.05) employing MACSE 2.06 (Ranwez et al., 2018; Scornavacca et al., 2019) and IQ-TREE 1.6.12 (Nguyen et al., 2015), and then progressively merged into a whole phylogeny.

## Aligning BUSCO genes
954 gene alignments were computed for each clade as below:
```
singularity run ~/bin/omm_macse_v11.05b.sif --out_dir ./actinopteri_${gene}_aligned --out_file_prefix actinopteri_${gene} --in_seq_file ${gene}_prealign_species.fasta --genetic_code_number 1 --java_mem 20000m
```

## Computing clade phylogenies
The output amino-acid alignments ${clade}_${gene}_final_align_AA.aln were used to compute the phylogeny for each clade. The same set of 107 BUSCO amino-acid sequences was concatenated and used to compute the phylogeny of all clades, except Mammalia for which another set of 100 genes was used. Phylogenies were calculated as below:
```
iqtree -s ${clade}_out_107.fst -st AA -m JTT+F+R10 -nt AUTO -ntmax 20 -bb 1000
```

## Merging clade trees into a full phylogeny
Here the example relative to Mammalia and Aves is reported, where Apteryx rowi and Ornithorhynchus_anatinus are the respective outgroups: their sequences are extracted with seqkit (Shen et al., 2016) and used to enrich the respective sister-groups alignments. Concatenated genes and used to recalculate the phylogeny with the outgroup, using the previously computed phylogeny to restrict the topology. shared_busco being the list of genes used for the phylogeny:
```
for gene in $(cat shared_busco); do
# align platypus sequence to the avian alignment
seqkit grep –n -p mammalia_”$gene”_final_align_NT.fasta -s Ornithorhynchus_anatinus | sed "s/-//g"  > Ornithorhynchus_anatinus.fasta # extract gapless outgroup gene
java -jar ~/bin/macse_v2.06.jar -prog enrichAlignment -align aves_"$gene"_final_mask_align_NT.aln -seq Ornithorhynchus_anatinus.fasta # align it to the avian alignment; output is aves_${gene}_final_mask_align_NT_AA.aln
sed -i 's/!/-/g' aves_"$gene"_final_mask_align_NT_AA.aln # Replace ! with - (symbol for frameshift mutations in MACSE)
sed -i 's/!/-/g' aves_"$gene"_final_mask_align_NT.aln

# align kiwi sequence to the mammalian alignment
seqkit grep –n -p aves_”$gene”_final_align_NT.fasta -s Apteryx_rowi | sed "s/-//g"  > Apteryx_rowi.fasta
java -jar ~/bin/macse_v2.06.jar -prog enrichAlignment -align mammalia_"$gene"_final_mask_align_NT.aln -seq Apteryx_rowi.fasta
sed -i 's/!/-/g' mammalia_"$gene"_final_mask_align_NT_AA.aln
sed -i 's/!/-/g' mammalia_"$gene"_final_mask_align_NT.aln
done

# concatenate
grep '>' mammalia*_final_mask_align_NT_AA.aln | cut -d: -f2 | sort -u | sed 's/>//' > list_esp1.txt # mammalia + outgroup species list
ls mammalia*_final_mask_align_NT_AA.aln > list_file1.txt # list of mammalia + outgroup alignment files
~/softwares/Concatenate list_file1.txt mammalia_outgroup.fst list_esp1.txt AA Fasta
grep '>' aves*_final_mask_align_NT_AA.aln | cut -d: -f2 | sort -u | sed 's/>//' > list_esp2.txt # aves + outgroup species list
ls $aves*_NT_AA.fasta > list_file2.txt # list of aves + outgroup alignment files
~/softwares/Concatenate list_file2.txt aves_outgroup.fst list_esp2.txt AA Fasta

# root trees with their outgroup
Rscript excludeSpeciesTree_edit.R mammalia_outgroup.fst mammalia_out_100.fst.treefile # outputs pruned new_mammalia_out_100.fst.treefile
iqtree -s mammalia_outgroup.fst -m LG+G -g new_mammalia_out_100.fst.treefile -pre mammalia_outgroup # use previously calculated mammalian topology to constrain the rooted tree
Rscript excludeSpeciesTree_edit.R aves_outgroup.fst aves_out_107.fst.treefile
iqtree -s aves_outgroup.fst -m LG+G -g new_aves_out_107.fst.treefile -pre aves_outgroup
```

The mammalian and avian phylogenies were joined into a tetrapod phylogeny with the tree editor Baobab (Dutheil & Galtier, 2002).
50 top-shared genes in Supplementary Table 5 were chosen to recompute the branch lengths. The 50 mammalian and avian nucleotide alignments were joined together into a tetrapod alignment:
```
for gene in $(cat 50sharedbuscogenes); do
java -jar ~/bin/macse_v2.06.jar -prog alignTwoProfiles -p1 mammalia_"$gene"_final_mask_align_NT.aln -p2 aves_"$gene"_final_mask_align_NT.aln -out_NT tetrapoda_"$gene"_NT.fasta -out_AA tetrapoda_"$gene"_AA.fasta
sed -i 's/!/-/g' tetrapoda_"$gene"_AA.fasta
done
```

The same workflow was repeated to join tetrapods to Actinopteri, Mollusca to Insecta, and finally vertebrates to protostomes.
With the full tree topology and the full gene alignments, 50 top-shared genes were chosen to recompute branch lengths and obtain `FULLTREE_brlengths_rooted.treefile`:
```
grep '>' fulltree*AA.fasta | cut -d: -f2 | sort -u | sed 's/>//' > list_esp.txt
ls fulltree*AA.fasta > list_file.txt
~/softwares/Concatenate list_file.txt fulltree_aligned.fasta list_esp.txt AA Fasta
# Estimate branch lengths on the full phylogeny
iqtree -s fulltree_aligned.fasta -st AA -m LG+G  -te $6 -pre fulltree_brle
```

## dN/dS estimation
Genes with more than 10% of missing sequence information are removed from the alignments. With preal being the unaligned fasta of a gene:
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
rename "s/^/${gene}_/" CM_counts_*

# Map substitutions by categories
sed -i 's/DnDs/Combination\(reg1=DnDs, reg2=SW\)/g' mapNH_currentgene.bpp
~/bin/MapNH_bppV3/mapnh param=mapNH_currentgene.bpp
rename "s/^/${gene}_/" counts_*

done
```

## TE annotation of genome assemblies
Earl Grey 1.3 (Baril et al., 2022) was used to annotate 29 dipteran assemblies with the command below:
```
earlGrey -g GCA_001542645.1.fa -s Anopheles_gambiae -r metazoa -o ./Anopheles_gambiae_earlgrey -t 8
```

## TE annotation from unassembled reads
The pipeline available at https://github.com/sigau/pipeline_dnapipe is based on two rounds of dnaPipeTE (Goubert et al., 2015) and was used to quantify TEs from short reads. The full process from reads download to second clean output can be automatized with
```
python3 sradownload_dnapipe.py table.tsv
```
where `table.tsv` can be any chunk of Supplementary Table 2 having an available reads experiment ID.
`dnaPT_summary_overall_recent.sh` was used to extract the overall and recent TE content. It leverages `dnapt_recentTEs.R` which adapts part of `dnaPT_landscapes.sh` from https://github.com/clemgoub/dnaPT_utils. The overall and recent TE content below 5% divergence were obtained with:
```
./dnaPT_summary_overall_recent.sh 5 recentTEs_5perc.tsv overallTEs.tsv
```

## Traits correlation under phylogenetic non-independence
Phylogenetic generalized Least Squares was performed with ape and nlme packages, using Supplementary Table 2 (`df`) and `FULLTREE_brlengths_rooted.treefile`:
```
library(nlme)
library(ape)
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
Bayesian inference of traits and substitution rates covariation was conducted with Coevol 1.6 (Lartillot & Poujol, 2011) on each clade separately.
**ADD WORKFLOW TO DEFINE GC-POOR AND GC-RICH GENESETS. Add script(s)?**

**ADD COMMAND LINE**


# References
Baril, T., Imrie, R. M., & Hayward, A. (2022). Earl Grey: a fully automated user-friendly transposable element annotation and analysis pipeline [Preprint]. In Review. doi: 10.21203/rs.3.rs-1812599/v1  

Dutheil, J., Gaillard, S., Bazin, E., Glémin, S., Ranwez, V., Galtier, N., & Belkhir, K. (2006). Bio++: a set of C++ libraries for sequence analysis, phylogenetics, molecular evolution and population genetics. BMC Bioinformatics, 7(1), 188. doi: 10.1186/1471-2105-7-188  

Goubert, C., Modolo, L., Vieira, C., ValienteMoro, C., Mavingui, P., & Boulesteix, M. (2015). De Novo Assembly and Annotation of the Asian Tiger Mosquito (Aedes albopictus) Repeatome with dnaPipeTE from Raw Genomic Reads and Comparative Analysis with the Yellow Fever Mosquito (Aedes aegypti). Genome Biology and Evolution, 7(4), 1192–1205. doi: 10.1093/gbe/evv050  

Guéguen, L., & Duret, L. (2018). Unbiased Estimate of Synonymous and Nonsynonymous Substitution Rates with Nonstationary Base Composition. Molecular Biology and Evolution, 35(3), 734–742. doi: 10.1093/molbev/msx308  

Guéguen, L., Gaillard, S., Boussau, B., Gouy, M., Groussin, M., Rochette, N. C., Bigot, T., Fournier, D., Pouyet, F., Cahais, V., Bernard, A., Scornavacca, C., Nabholz, B., Haudry, A., Dachary, L., Galtier, N., Belkhir, K., & Dutheil, J. Y. (2013). Bio++: Efficient Extensible Libraries and Tools for Computational Molecular Evolution. Molecular Biology and Evolution, 30(8), 1745–1750. doi: 10.1093/molbev/mst097  

Gurevich, A., Saveliev, V., Vyahhi, N., & Tesler, G. (2013). QUAST: quality assessment tool for genome assemblies. Bioinformatics, 29(8), 1072–1075. doi: 10.1093/bioinformatics/btt086  

Lartillot, N., & Poujol, R. (2011). A Phylogenetic Model for Investigating Correlated Evolution of Substitution Rates and Continuous Phenotypic Characters. Molecular Biology and Evolution, 28(1), 729–744. doi: 10.1093/molbev/msq244  

Manni, M., Berkeley, M. R., Seppey, M., Simão, F. A., & Zdobnov, E. M. (2021). BUSCO Update: Novel and Streamlined Workflows along with Broader and Deeper Phylogenetic Coverage for Scoring of Eukaryotic, Prokaryotic, and Viral Genomes. Molecular Biology and Evolution, 38(10), 4647–4654. doi: 10.1093/molbev/msab199  

Nguyen, L.-T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2015). IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. Molecular Biology and Evolution, 32(1), 268–274. doi: 10.1093/molbev/msu300  

Ranwez, V., Douzery, E. J. P., Cambon, C., Chantret, N., & Delsuc, F. (2018). MACSE v2: Toolkit for the Alignment of Coding Sequences Accounting for Frameshifts and Stop Codons. Molecular Biology and Evolution, 35(10), 2582–2584. doi: 10.1093/molbev/msy159  

Romiguier, J., Figuet, E., Galtier, N., Douzery, E. J. P., Boussau, B., Dutheil, J. Y., & Ranwez, V. (2012). Fast and Robust Characterization of Time-Heterogeneous Sequence Evolutionary Processes Using Substitution Mapping. PLOS ONE, 7(3), e33852. doi: 10.1371/journal.pone.0033852  

Scornavacca, C., Belkhir, K., Lopez, J., Dernat, R., Delsuc, F., Douzery, E. J. P., & Ranwez, V. (2019). OrthoMaM v10: Scaling-Up Orthologous Coding Sequence and Exon Alignments with More than One Hundred Mammalian Genomes. Molecular Biology and Evolution, 36(4), 861–862. doi: 10.1093/molbev/msz015  

Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE, 11(10), e0163962. doi: 10.1371/journal.pone.0163962  
