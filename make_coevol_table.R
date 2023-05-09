#!/usr/bin/Rscript

#Rscript make_coevol_table.R Actinopteri

library(dplyr)
args = commandArgs(trailingOnly=TRUE)

df <- read.table(file="table_coevol.tsv", header=T,sep="\t")
sp2rm <- read.table("sp2rm_50percgenespersp", header=F)$V1
df <- filter(df, ! Species %in% sp2rm)
clade <- filter(df, Clade==args[1]) %>% select(Species, Unified_Cvalues_methods, BodyLength_cm, BasalMetRate_mLO2hr_pantheria, AgeatFirstBirth_d_pantheria, PopulationDensity_nperkm2_pantheria, MaxLongevity_y, Mass_g, Metabolic_rate_w, SexMaturity_d, DepthRange, Overall_repeat_bp, dnaPT_DNA_bp, dnaPT_RC_bp, dnaPT_LTR_bp, dnaPT_PLE_bp, dnaPT_LINE_bp, dnaPT_SINE_bp, dnaPT_Retroposon_bp, dnaPT_Unknown_bp, DNA_recent_bp, LINE_recent_bp, LTR_recent_bp, Overall_TEs_recent_bp, PLE_recent_bp, RC_recent_bp, Retroposon_recent_bp, SINE_recent_bp, Unknown_recent_bp)

clade$Unified_Cvalues_methods<-round(clade$Unified_Cvalues_methods, 0)
colnames(clade)[2]<-"GenomeSize"
colnames(clade)[12:20]<-c("OverallRepeat_bp","DNA_bp", "RC_bp","LTR_bp","PLE_bp","LINE_bp","SINE_bp","Retroposon_bp","Unknown_bp")
colnames(clade)[24]<-"OverallTEs_recent_bp"

for (c in colnames(clade)[-1]) {
	if ( unique(is.na(clade[[c]])) && length(unique(is.na(clade[[c]])))==1 ) {
		clade <- clade[, !colnames(clade) %in% c ]
	}
}
		
		
		
ntaxa<-nrow(clade)
ntraits<-ncol(clade)-1
secondline <- c(ntaxa, ntraits, colnames(clade)[-1])
write("#TRAITS", file="header")
write(secondline, file="header", append=T, sep=" ", ncolumns=length(secondline))
write.table(clade, file="traits", sep=" ", quote=F, na="-1", row.names=F, col.names=F)
		
bashcmd <- paste0("cat header traits > ", args[1], "_coevol.txt && rm header traits")
system(bashcmd)
