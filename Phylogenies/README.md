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
