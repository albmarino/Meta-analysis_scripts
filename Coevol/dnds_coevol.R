#!/usr/bin/Rscript

# Script to extract terminal dS, dNdS and branch lenght estimated by coevol

# Required input: tab readcoevol outputs of all the clades (*postmeanbranchnonsynrate.tab and *postmeanbranchsynrate.tab) from both gcpoor and gcrich runs in pwd
# Example run: Rscript dnds_coevol.R dnds_coevol.tsv # 1st argument is the outfile name 

library(dplyr)
args = commandArgs(trailingOnly=TRUE)

bashcmd_poor_ds <- "cat *gcpoor_dsom.postmeanbranchsynrate.tab | grep -v '#'"
bashcmd_poor_dnds <- "cat *gcpoor_dsom.postmeanbranchnonsynrate.tab | grep -v '#'"

bashcmd_rich_ds <- "cat *gcrich_dsom.postmeanbranchsynrate.tab | grep -v '#'"
bashcmd_rich_dnds <- "cat *gcrich_dsom.postmeanbranchnonsynrate.tab | grep -v '#'"

ds_gcpoor <- read.table(header=F, sep = "\t", text=system(bashcmd_poor_ds, intern = T))
ds_gcrich <- read.table(header=F, sep = "\t", text=system(bashcmd_rich_ds, intern = T))
dnds_gcpoor <- read.table(header=F, sep = "\t", text=system(bashcmd_poor_dnds, intern = T))
dnds_gcrich <-read.table(header=F, sep = "\t", text=system(bashcmd_rich_dnds, intern = T))

species <- unique(c(ds_gcpoor$V1, ds_gcpoor$V2))

dnds <-  data.frame(matrix(nrow = 0, ncol = 7))

for (sp in species) {

	brlenpoor <- filter(ds_gcpoor, V1==sp & V2==sp)$V3
	dspoor <- filter(ds_gcpoor, V1==sp & V2==sp)$V4
	dndspoor <- filter(dnds_gcpoor, V1==sp & V2==sp)$V4
	brlenrich <- filter(ds_gcrich, V1==sp & V2==sp)$V3 
	dsrich <- filter(ds_gcrich, V1==sp & V2==sp)$V4
	dndsrich <- filter(dnds_gcrich, V1==sp & V2==sp)$V4

	dnds <- rbind(dnds, data.frame(sp, brlenpoor, dspoor, dndspoor, brlenrich, dsrich, dndsrich))
	}

colnames(dnds) <- c("Species", "brlen_coevol_gcpoor", "dS_coevol_gcpoor", "dNdS_coevol_gcpoor", "brlen_coevol_gcrich", "dS_coevol_gcrich", "dNdS_coevol_gcrich")

write.table(dnds, file=args[1], quote=F,sep="\t", row.names=F, col.names=T)
