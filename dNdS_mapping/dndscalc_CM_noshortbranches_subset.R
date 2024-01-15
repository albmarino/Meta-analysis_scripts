#!/usr/bin/env Rscript

# Calculate dN, dS with Cumulated Mean approach: dN and dS for each species are obtained as sum(counts)/sum(opportunities/brlength) = sum(counts)*sum(brlength/opportunities), going over the given geneset

# Required input: clade nwk phylogeny as 1st argument, newline-seaprated list of gene prefixes that we want to include in the dN/dS computation as 2nd argument (e.g. remove outlier genes selected by PhylteR) (in the form falserm_<clade>_<genename>__final_mask_align_NT.aln), outfile name as 3rd argument
# Output files of run_bppml_mapnh_V3.sh must be in pwd
# Example command: Rscript dndscalc_CM_noshortbranches_subset.R actinopteri_out_107.fst.treefile actinopteri_genes2keep trim700_actinopteri_dNdSCM_phylter_noshortbranches.tsv

library(ape)
library(stringr)
args = commandArgs(trailingOnly=TRUE)

### list of all species in the treefile

mytree<-read.tree(args[1]) # all species tree
species<-mytree$tip.label
dfspecies<-data.frame(species)

## list of substitution trees (*dN.dnd files) for suffixes

genesprefix <- as.vector(read.table(args[2], header=F, sep="\t")$V1)

dfCM <- data.frame(matrix(ncol = 5, nrow = 0))

for (sp in dfspecies$species) {

	K_N = 0 # nonsyn counts
	K_S = 0 # syn counts
	S_N = 0 # sites subsceptible to nonsyn substitutions (Opport_N/brlen)
	S_S = 0 # sites subsceptible to syn substitutions (Opport_S/brlen)
	n_markers = 0

# For each gene of a species

	for (i in 1:length(genesprefix)) {

		treecounts_n <- read.tree(paste0(genesprefix[i],"_CM_splitcounts_dN.dnd"))
		treecounts_s <- read.tree(paste0(genesprefix[i],"_CM_splitcounts_dS.dnd"))
		treeopport_n <- read.tree(paste0(genesprefix[i], "_CM_splitcounts_dN_norm.dnd"))
		treeopport_s <- read.tree(paste0(genesprefix[i], "_CM_splitcounts_dS_norm.dnd"))
		treebrlens<- read.tree(paste0(genesprefix[i],"_final_mask_align_NT.aln.ml_h.dnd_1"))

		if (sp %in% treecounts_n$tip.label) {

			sp_position <- grep(sp, treecounts_n$tip.label) # get the position of that species in the current trees
			brlen <- treebrlens$edge.length[treebrlens$edge[,2]==sp_position]

# take into account the current gene only if it has a branch length longer than 0.001

			if (brlen > 0.001) {

# Sum the nonsyn and syn counts across all genes, respectively
				n_markers <- n_markers + 1
				K_N = K_N + treecounts_n$edge.length[treecounts_n$edge[,2]==sp_position]
				K_S = K_S + treecounts_s$edge.length[treecounts_s$edge[,2]==sp_position]

# Sum the sites subsceptible to nonsyn and syn substitutions, respectively

				S_N = S_N + treeopport_n$edge.length[treeopport_n$edge[,2]==sp_position]/treebrlens$edge.length[treebrlens$edge[,2]==sp_position]
				S_S = S_S + treeopport_s$edge.length[treeopport_s$edge[,2]==sp_position]/treebrlens$edge.length[treebrlens$edge[,2]==sp_position]
				}

			else {next}
			
			}
		}

# Once obtained sums of counts and of substitution sites, get cumulated mean dN/dS

	dN_CM = (K_N/S_N)
	dS_CM = (K_S/S_S)
	dNdS_CM = dN_CM/dS_CM

	dfCM <- rbind(dfCM, cbind(Species=sp, n_markers_phylter_noshortbr=n_markers, dN_CM_phylter_noshortbr=dN_CM, dS_CM_phylter_noshortbr=dS_CM, dNdS_CM_phylter_noshortbr=dNdS_CM))
	}

write.table(dfCM, file=args[3],quote=F,sep='\t',row.names=F,col.names=T)
