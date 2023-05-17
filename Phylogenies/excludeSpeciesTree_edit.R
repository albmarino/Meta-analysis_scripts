#!/usr/bin/env Rscript

# Remove specific tip labels from a tree
# Example --> Rscript excludeSpeciesTree_edit.R mammalia_outgroup.fst mammalia_out_100.fst.treefile

args = commandArgs(trailingOnly=TRUE)

library("ape")

al<-read.dna(args[1], format="fasta", as.character=T, as.matrix=F)

tree<-read.tree(args[2])

sp2rm<-tree$tip.label[! tree$tip.label %in% names(al)]

print(sp2rm)
tree<-drop.tip(tree, sp2rm)
tree <- unroot(tree)

write.tree(tree, file=paste0("new_", args[2]))
