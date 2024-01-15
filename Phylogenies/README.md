## Computing clade and full phylogeny

Scripts to compute clade phylogenies and merge them.

Clade phylogenies are computed over the same geneset (`concatenate.sh`, `iqtree.sh`).<br>
Briefly, a clade phylogeny is rooted with one outgroup species of the closest clade by enriching the alignments of the 107 genes with the homologous sequences of the outgroup. The same procedure is reiterated until the full topology is reached. Example between mammals and birds: `merge_phylogenies_example.sh`.
In parallel, alignment profiles of 50 top-shared busco genes are aligned to each other. The same procedure is reiterated until the full profile is reached for every gene (`align_profiles_example.sh`).
The concatenation of the 50 full-profile genes is used to compute branch lenghts for the full topology (`recompute_brlens_fullphylo.sh`).

Output full and clade rooted phylogenies are supplied (`*brlen_rooted.treefile`).

## References

Nguyen, L.-T., Schmidt, H. A., von Haeseler, A., & Minh, B. Q. (2015). IQ-TREE: A Fast and Effective Stochastic Algorithm for Estimating Maximum-Likelihood Phylogenies. Molecular Biology and Evolution, 32(1), 268–274. doi: 10.1093/molbev/msu300

Ranwez, V., Douzery, E. J. P., Cambon, C., Chantret, N., & Delsuc, F. (2018). MACSE v2: Toolkit for the Alignment of Coding Sequences Accounting for Frameshifts and Stop Codons. Molecular Biology and Evolution, 35(10), 2582–2584. doi: 10.1093/molbev/msy159

Scornavacca, C., Belkhir, K., Lopez, J., Dernat, R., Delsuc, F., Douzery, E. J. P., & Ranwez, V. (2019). OrthoMaM v10: Scaling-Up Orthologous Coding Sequence and Exon Alignments with More than One Hundred Mammalian Genomes. Molecular Biology and Evolution, 36(4), 861–862. doi: 10.1093/molbev/msz015

Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. PLOS ONE, 11(10), e0163962. doi: 10.1371/journal.pone.0163962
