#!/usr/bin/Rscript

## GC3 content overview
# Produce boxplots and histograms to show the overall distribution of average GC3 content, average missing data in a sequence, and average gene sharing across species (this for each dataset). Supply three integer arguments to sample genes shared from at least the 90% of the reference dataset, from them sample the 300 top GC-rich genes, from which the 100 top shared genes are extracted (normally the quantity of missing data can be ignored). Three boxplots are produced for the extracted geneset to show (1) the distribution of average GC3 content, (2) the distribution of average missing data in a sequence, and (3) the distribution of gene sharing across species (in % relative to the reference dataset). The venn diagram shows the overlapping among the selected genesets of the several datasets.
## Number of genes per species overview
# Produce barplot with the gene proportion of each species compared to the total dataset (954 genes) and output a list of species under the given threshold (args[6]) to potentially exclude for downstream analyses (namely dN/dS calculation, coevol run, etc.). The sp2rm list can be used as input in gcstat.sh whose output can be used again in turn in this script. 
# Example --> Rscript ~/bin/gc3plots.R 150 100 90 gc3_report_new.tsv genespersp.tsv 70

.libPaths( c( "/home/alba.marino/R/x86_64-pc-linux-gnu-library/4.2/" , .libPaths() ) )
library(ggplot2)
library(ggpubr)
library(dplyr)
library(stringr)
library(ggVennDiagram)
library(RColorBrewer)
args = commandArgs(trailingOnly=TRUE)

df<- read.table(args[4], header=T,sep="\t")
df_species <- read.table(args[5],header=T,sep="\t")
appliedfilt <- str_remove(args[5], "genespersp_") %>% str_remove(".tsv")

# Show GC3 and species sharing stats in boxplots
pdf(paste0("genes_gc3_boxpl_report_", appliedfilt,".pdf"), width=9, height=15)
par(mfrow=c(2,1), cex.axis=0.7)
boxplot(df$GC3_average_total, df$GC3_average_actinopteri, df$GC3_average_aves, df$GC3_average_insecta, df$GC3_average_mammalia, df$GC3_average_mollusca, names=c("","","","","",""), col=c("gray34","tomato1","yellow4","forestgreen","dodgerblue3","magenta1"), notch=T, ylab="Average GC3 content")

boxplot(df$Missingdata_average_total, df$Missingdata_actinopteri, df$Missingdata_aves, df$Missingdata_insecta, df$Missingdata_mammalia, df$Missingdata_mollusca, names=c("","","","","",""), col=c("gray34","tomato1","yellow4","forestgreen","dodgerblue3","magenta1"), notch=T, ylab="Average missing data in a sequence")

boxplot(df$Gene_nsp_total/df$Sp_total*100, df$Gene_nsp_actinopteri/df$Sp_actinopteri*100, df$Gene_nsp_aves/df$Sp_aves*100, df$Gene_nsp_insecta/df$Sp_insecta*100, df$Gene_nsp_mammalia/df$Sp_mammalia*100, df$Gene_nsp_mollusca/df$Sp_mollusca*100, col=c("gray34","tomato1","yellow4","forestgreen","dodgerblue3","magenta1"),  names=c(paste0("Total (N=",nrow(df_species), ")"),paste0("Actinopteri (N=",nrow(df_species[df_species$Clade=="Actinopteri",]), ")"),paste0("Aves (N=",nrow(df_species[df_species$Clade=="Aves",]), ")"),paste0("Insecta (N=",nrow(df_species[df_species$Clade=="Insecta",]), ")"),paste0("Mammalia (N=",nrow(df_species[df_species$Clade=="Mammalia",]), ")"),paste0("Mollusca (N=",nrow(df_species[df_species$Clade=="Mollusca",]), ")")), notch=T, ylab="Shared species across a dataset (%)")
dev.off()


# Show GC3 and species sharing stats in histograms
total_gc3 <- ggplot(df, aes(x=GC3_average_total))+
geom_histogram(color="white", fill="gray34")+
labs(y="Gene count", x="Average GC3 content in all dataset", title=paste0("Total dataset N = ", nrow(df_species)))+
xlim(0.24, 0.97)

total_miss <- ggplot(df, aes(x=Missingdata_average_total))+
geom_histogram(color="white", fill="gray34")+
labs(y="Gene count", x="Average missing data in a sequence across all dataset", title=paste0("Total dataset N = ", nrow(df_species)))

total_sharedsp <- ggplot(df, aes(x=Gene_nsp_total/Sp_total*100))+
geom_histogram(color="white", fill="gray34")+
labs(y="Gene count", x="Shared species across all dataset (%)", title=paste0("Total dataset N = ", nrow(df_species)))

actinopteri_gc3 <- ggplot(df, aes(x=GC3_average_actinopteri))+
geom_histogram(color="white", fill="tomato1")+
labs(y="Gene count", x="Average GC3 content in Actinopteri", title=paste0("Actinopteri dataset N = ", nrow(df_species[df_species$Clade=="Actinopteri",])))+
xlim(0.24, 0.97)

actinopteri_miss <- ggplot(df, aes(x=Missingdata_actinopteri))+
geom_histogram(color="white", fill="tomato1")+
labs(y="Gene count", x="Average missing data in a sequence across actinopteri", title=paste0("Actinopteri dataset N = ", nrow(df_species[df_species$Clade=="Actinopteri",])))

actinopteri_sharedsp <- ggplot(df, aes(x=Gene_nsp_actinopteri/Sp_actinopteri*100))+
geom_histogram(color="white", fill="tomato1")+
labs(y="Gene count", x="Shared species across Actinopteri (%)", title=paste0("Actinopteri dataset N = ", nrow(df_species[df_species$Clade=="Actinopteri",])))

aves_gc3 <- ggplot(df, aes(x=GC3_average_aves))+
geom_histogram(color="white", fill="yellow4")+
labs(y="Gene count", x="Average GC3 content in Aves", title=paste0("Aves dataset N = ", nrow(df_species[df_species$Clade=="Aves",])))+
xlim(0.24, 0.97)

aves_miss <- ggplot(df, aes(x=Missingdata_aves))+
geom_histogram(color="white", fill="yellow4")+
labs(y="Gene count", x="Average missing data in a sequence across aves", title=paste0("Aves dataset N = ", nrow(df_species[df_species$Clade=="Aves",])))

aves_sharedsp <- ggplot(df, aes(x=Gene_nsp_aves/Sp_aves*100))+
geom_histogram(color="white", fill="yellow4")+
labs(y="Gene count", x="Shared species across Aves (%)", title=paste0("Aves dataset N = ", nrow(df_species[df_species$Clade=="Aves",])))

insecta_gc3 <- ggplot(df, aes(x=GC3_average_insecta))+
geom_histogram(color="white", fill="forestgreen")+
labs(y="Gene count", x="Average GC3 content in Insecta", title=paste0("Insecta dataset N = ", nrow(df_species[df_species$Clade=="Insecta",])))+
xlim(0.24, 0.97)

insecta_miss <- ggplot(df, aes(x=Missingdata_insecta))+
geom_histogram(color="white", fill="forestgreen")+
labs(y="Gene count", x="Average missing data in a sequence across insecta", title=paste0("Insecta dataset N = ", nrow(df_species[df_species$Clade=="Insecta",])))

insecta_sharedsp <- ggplot(df, aes(x=Gene_nsp_insecta/Sp_insecta*100))+
geom_histogram(color="white", fill="forestgreen")+
labs(y="Gene count", x="Shared species across Insecta (%)", title=paste0("Insecta dataset N = ", nrow(df_species[df_species$Clade=="Insecta",])))

mammalia_gc3 <- ggplot(df, aes(x=GC3_average_mammalia))+
geom_histogram(color="white", fill="dodgerblue3")+
labs(y="Gene count", x="Average GC3 content in Mammalia", title=paste0("Mammalia dataset N = ", nrow(df_species[df_species$Clade=="Mammalia",])))+
xlim(0.24, 0.97)

mammalia_miss <- ggplot(df, aes(x=Missingdata_mammalia))+
geom_histogram(color="white", fill="dodgerblue3")+
labs(y="Gene count", x="Average missing data in a sequence across mammalia", title=paste0("Mammalia dataset N = ", nrow(df_species[df_species$Clade=="Mammalia",])))

mammalia_sharedsp <- ggplot(df, aes(x=Gene_nsp_mammalia/Sp_mammalia*100))+
geom_histogram(color="white", fill="dodgerblue3")+
labs(y="Gene count", x="Shared species across Mammalia (%)",title=paste0("Mammalia dataset N = ", nrow(df_species[df_species$Clade=="Mammalia",])))

mollusca_gc3 <- ggplot(df, aes(x=GC3_average_mollusca))+
geom_histogram(color="white", fill="magenta1")+
labs(y="Gene count", x="Average GC3 content in Mollusca", title=paste0("Mollusca dataset N = ", nrow(df_species[df_species$Clade=="Mollusca",])))+
xlim(0.24, 0.97)

mollusca_miss <- ggplot(df, aes(x=Missingdata_mollusca))+
geom_histogram(color="white", fill="magenta1")+
labs(y="Gene count", x="Average missing data in a sequence across mollusca", title=paste0("Mollusca dataset N = ",nrow(df_species[df_species$Clade=="Mollusca",])))

mollusca_sharedsp <- ggplot(df, aes(x=Gene_nsp_mollusca/Sp_mollusca*100))+
geom_histogram(color="white", fill="magenta1")+
labs(y="Gene count", x="Shared species across Mollusca (%)",title=paste0("Mollusca dataset N = ", nrow(df_species[df_species$Clade=="Mollusca",])))

pdf(paste0("genes_gc3_histo_report_", appliedfilt, ".pdf"),width=11)
print(ggarrange(total_gc3, total_miss, total_sharedsp, nrow=1,ncol=2))
print(ggarrange(actinopteri_gc3, actinopteri_miss, actinopteri_sharedsp, nrow=1,ncol=2))
print(ggarrange(aves_gc3, aves_miss, aves_sharedsp, nrow=1,ncol=2))
print(ggarrange(insecta_gc3, insecta_miss, insecta_sharedsp, nrow=1,ncol=2))
print(ggarrange(mammalia_gc3, mammalia_miss, mammalia_sharedsp, nrow=1,ncol=2))
print(ggarrange(mollusca_gc3, mollusca_miss, mollusca_sharedsp, nrow=1,ncol=2))
dev.off()

# Show proportion of genes for every species relatively to the clade and to the total gene number
		
totalgenes <- filter(df_species, Clade=="Mollusca") %>% pull(ngenes_total) %>% unique
mol_totgenes <- filter(df_species, Clade=="Mollusca") %>% pull(ngenes_clade) %>% unique
mam_totgenes <- filter(df_species, Clade=="Mammalia") %>% pull(ngenes_clade) %>% unique
ins_totgenes <- filter(df_species, Clade=="Insecta") %>% pull(ngenes_clade) %>% unique
ave_totgenes <- filter(df_species, Clade=="Aves") %>% pull(ngenes_clade) %>% unique
act_totgenes <- filter(df_species, Clade=="Actinopteri") %>% pull(ngenes_clade) %>% unique
		
#molplot <- filter(df_species, Clade=="Mollusca") %>% ggplot(aes(x=Species, y=ngenes_sp/ngenes_clade*100))+
#eom_bar(stat="identity", fill="magenta1")+
#labs(y="genes % out of Mollusca", title=paste0("Mollusca total genes = ", mol_totgenes, "   (Total = ", totalgenes, ")"))+
#theme(axis.title.y=element_blank())+
#ylim(0, 100)+
#coord_flip()
molplot_tot <- filter(df_species, Clade=="Mollusca") %>% ggplot(aes(x=Species, y=ngenes_sp/ngenes_total*100))+
geom_bar(stat="identity", fill="magenta1")+
labs(y="genes % out of Total",  title=paste0("Mollusca total genes = ", mol_totgenes, "   (Total = ", totalgenes, ")"))+
theme(axis.title.y=element_blank())+
ylim(0, 100)+
coord_flip()

#mamplot <- filter(df_species, Clade=="Mammalia") %>% ggplot(aes(x=Species, y=ngenes_sp/ngenes_clade*100))+
#geom_bar(stat="identity", fill="dodgerblue3")+
#labs(y="genes % out of Mammalia", title=paste0("Mammalia total genes = ", mam_totgenes, "   (Total = ", totalgenes, ")"))+
#theme(axis.title.y=element_blank())+
#ylim(0, 100)+
#coord_flip()
mamplot_tot <- filter(df_species, Clade=="Mammalia") %>% ggplot(aes(x=Species, y=ngenes_sp/ngenes_total*100))+
geom_bar(stat="identity", fill="dodgerblue3")+
labs(y="genes % out of Total",  title=paste0("Mammalia total genes = ", mam_totgenes, "   (Total = ", totalgenes, ")"))+
theme(axis.title.y=element_blank())+
ylim(0, 100)+
coord_flip()
		
#insplot <- filter(df_species, Clade=="Insecta") %>% ggplot(aes(x=Species, y=ngenes_sp/ngenes_clade*100))+
#geom_bar(stat="identity", fill="forestgreen")+
#labs(y="genes % out of Insecta", title=paste0("Insecta total genes = ", ins_totgenes, "   (Total = ", totalgenes, ")"))+
#heme(axis.title.y=element_blank())+
#lim(0, 100)+
#coord_flip()
insplot_tot <- filter(df_species, Clade=="Insecta") %>% ggplot(aes(x=Species, y=ngenes_sp/ngenes_total*100))+
geom_bar(stat="identity", fill="forestgreen")+
labs(y="genes % out of Total",  title=paste0("Insecta total genes = ", ins_totgenes, "   (Total = ", totalgenes, ")"))+
theme(axis.title.y=element_blank())+
ylim(0, 100)+
coord_flip()
		
#aveplot <- filter(df_species, Clade=="Aves") %>% ggplot(aes(x=Species, y=ngenes_sp/ngenes_clade*100))+
#geom_bar(stat="identity", fill="yellow4")+
#labs(y="genes % out of Aves", title=paste0("Aves total genes = ", ave_totgenes, "   (Total = ", totalgenes, ")"))+
#heme(axis.title.y=element_blank())+
#ylim(0, 100)+
#coord_flip()
aveplot_tot <- filter(df_species, Clade=="Aves") %>% ggplot(aes(x=Species, y=ngenes_sp/ngenes_total*100))+
geom_bar(stat="identity", fill="yellow4")+
labs(y="genes % out of Total",  title=paste0("Aves total genes = ", ave_totgenes, "   (Total = ", totalgenes, ")"))+
theme(axis.title.y=element_blank())+
ylim(0, 100)+
coord_flip()
		
#actplot <- filter(df_species, Clade=="Actinopteri") %>% ggplot(aes(x=Species, y=ngenes_sp/ngenes_clade*100))+
#geom_bar(stat="identity", fill="tomato1")+
#labs(y="genes % out of Actinopteri", title=paste0("Actinopteri total genes = ", act_totgenes, "   (Total = ", totalgenes, ")"))+
#theme(axis.title.y=element_blank())+
#lim(0, 100)+
#coord_flip()
actplot_tot <- filter(df_species, Clade=="Actinopteri") %>% ggplot(aes(x=Species, y=ngenes_sp/ngenes_total*100))+
geom_bar(stat="identity", fill="tomato1")+
labs(y="genes % out of Total",  title=paste0("Actinopteri total genes = ", act_totgenes, "   (Total = ", totalgenes, ")"))+
theme(axis.title.y=element_blank())+
ylim(0, 100)+
coord_flip()
		
pdf(paste0("genespersp_barplot_", appliedfilt, ".pdf"),width=11, height=30)
print(actplot_tot)
print(aveplot_tot)
print(insplot_tot)
print(mamplot_tot)
print(molplot_tot)
dev.off()

# Produce list of species to remove upon a minimum threshold of e.g. 70% of genes each

if (is.na(args[6])) { print("Threshold for removing species not supplied, sp2rm list not produced")} else { sp2rm <- filter(df_species, ngenes_sp/ngenes_total*100<as.numeric(args[6])) %>% pull(Species)
write.table(sp2rm, file=paste0("sp2rm_", args[6], "percgenespersp"),row.names=F,col.names=F,quote=F,sep="\n")
}

# Show GC3 and species sharing in the GC-rich selected geneset for every dataset, filtering for minimum perc of gene sharing across species (Total and in each clade)

bashcmd <- paste0("sed '1d' ", args[4]," | awk '{FS=\"\\t\"} {if ($2/$3*100>=", args[3],") {print $0}}' | sort -r -k5| head -", as.numeric(args[1]), "| sort -r -k2 |head -", as.numeric(args[2]), "| cut -f1,2,3,5")
subgc_tot <- read.table(sep = "\t", text=system(bashcmd, intern = T))
subgc_tot <- mutate(subgc_tot, perc=V2/V3*100, Group=rep("Total",nrow(subgc_tot)))

bashcmd <- paste0("sed '1d' ",args[4]," | awk '{FS=\"\\t\"} {if ($7/$8*100>=", args[3],") {print $0}}' | sort -r -k10| head -", as.numeric(args[1]), "| sort -r -k7 |head -", as.numeric(args[2]), "| cut -f1,7,8,10")
subgc_act <- read.table(sep = "\t", text=system(bashcmd, intern = T))
subgc_act <- mutate(subgc_act, perc=V2/V3*100, Group=rep("Actinopteri",nrow(subgc_act)))

bashcmd <- paste0("sed '1d' ", args[4]," | awk '{FS=\"\\t\"} {if ($12/$13*100>=", args[3],") {print $0}}' | sort -r -k15| head -", as.numeric(args[1]), "| sort -r -k12 |head -", as.numeric(args[2]), "| cut -f1,12,13,15")
subgc_aves <- read.table(sep = "\t", text=system(bashcmd, intern = T))
subgc_aves <- mutate(subgc_aves, perc=V2/V3*100, Group=rep("Aves",nrow(subgc_aves)))

bashcmd <- paste0("sed '1d' ", args[4], " | awk '{FS=\"\\t\"} {if ($17/$18*100>=", args[3],") {print $0}}' | sort -r -k20| head -", as.numeric(args[1]), "| sort -r -k17 |head -", as.numeric(args[2]), "| cut -f1,17,18,20")
subgc_ins <- read.table(sep = "\t", text=system(bashcmd, intern = T))
subgc_ins <- mutate(subgc_ins, perc=V2/V3*100, Group=rep("Insecta",nrow(subgc_ins)))

bashcmd <- paste0("sed '1d' ", args[4]," | awk '{FS=\"\\t\"} {if ($22/$23*100>=", args[3],") {print $0}}' | sort -r -k25| head -", as.numeric(args[1]), "| sort -r -k22 |head -", as.numeric(args[2]), "| cut -f1,22,23,25")
subgc_mam <- read.table(sep = "\t", text=system(bashcmd, intern = T))
subgc_mam <- mutate(subgc_mam, perc=V2/V3*100, Group=rep("Mammalia",nrow(subgc_mam)))

bashcmd <- paste0("sed '1d' ", args[4]," | awk '{FS=\"\\t\"} {if ($27/$28*100>=", args[3],") {print $0}}' |sort -r -k30| head -", as.numeric(args[1]), "| sort -r -k27 |head -", as.numeric(args[2]), "| cut -f1,27,28,30")
subgc_mol <- read.table(sep = "\t", text=system(bashcmd, intern = T))
subgc_mol <- mutate(subgc_mol, perc=V2/V3*100, Group=rep("Mollusca",nrow(subgc_mol)))

subgc<-rbind(subgc_tot, subgc_act, subgc_aves, subgc_ins, subgc_mam, subgc_mol)

sharedsp <- ggplot(subgc, aes(x=Group, y=perc, fill=Group))+
geom_boxplot()+
labs(y="Shared species in ref dataset (%)", title=paste0("Min ", args[3], "% of species sharing per gene in a dataset --> ", args[1], " GC3-richest --> ", args[2], " top shared genes"), subtitle=paste0("Actinopteri=",nrow(df_species[df_species$Clade=="Actinopteri",]),"    Aves=", nrow(df_species[df_species$Clade=="Aves",]),   "    Insecta=", nrow(df_species[df_species$Clade=="Insecta",]), "    Mammalia=", nrow(df_species[df_species$Clade=="Mammalia",]),"    Mollusca=", nrow(df_species[df_species$Clade=="Mollusca",]),"    Total=", nrow(df_species)))+
scale_fill_manual(values=c("tomato1", "yellow4", "forestgreen", "dodgerblue3", "magenta1", "gray34"))

gc3cont <- ggplot(subgc, aes(x=Group, y=V4, fill=Group))+
geom_boxplot()+
labs(y="GC3 content")+
scale_fill_manual(values=c("tomato1", "yellow4", "forestgreen", "dodgerblue3", "magenta1", "gray34"))

vec <- list(Total=filter(subgc, Group=="Total") %>% pull(V1), Actinopteri=filter(subgc, Group=="Actinopteri") %>% pull(V1), Aves=filter(subgc, Group=="Aves") %>% pull(V1), Insecta=filter(subgc, Group=="Insecta") %>% pull(V1), Mammalia=filter(subgc, Group=="Mammalia") %>% pull(V1), Mollusca=filter(subgc, Group=="Mollusca") %>% pull(V1))

vendiag <- ggVennDiagram(vec, label = "count", label_alpha = 0, set_color = c("Total" = "gray34","Actinopteri" ="tomato1",'Aves' = 'yellow4', "Insecta"="forestgreen", "Mammalia" = "dodgerblue3", "Mollusca" = "magenta1"))+
scale_color_manual(values = c("gray34","tomato1","yellow4","forestgreen","dodgerblue3","magenta1"))+
scale_fill_distiller(palette = "Blues", direction = 1)+
labs(title=paste0("Set of ", args[2], " genes"))

#Write out lists of GC-rich selected genesets for each clade
write.table(subgc_tot$V1,file=paste0("gcrichset_tot_", appliedfilt, "_",  args[3],"_",args[1],"_",args[2],".tsv"),sep="\n",row.names=F,col.names=F,quote=F)
write.table(subgc_act$V1,file=paste0("gcrichset_actinopteri_",appliedfilt, "_",args[3],"_",args[1],"_",args[2],".tsv"),sep="\n",row.names=F,col.names=F,quote=F)
write.table(subgc_aves$V1,file=paste0("gcrichset_aves_",appliedfilt, "_",args[3],"_",args[1],"_",args[2],".tsv"),sep="\n",row.names=F,col.names=F,quote=F)
write.table(subgc_ins$V1,file=paste0("gcrichset_insecta_",appliedfilt, "_",args[3],"_",args[1],"_",args[2],".tsv"),sep="\n",row.names=F,col.names=F,quote=F)
write.table(subgc_mam$V1,file=paste0("gcrichset_mammalia_",appliedfilt, "_",args[3],"_",args[1],"_",args[2],".tsv"),sep="\n",row.names=F,col.names=F,quote=F)
write.table(subgc_mol$V1,file=paste0("gcrichset_mollusca_",appliedfilt, "_",args[3],"_",args[1],"_",args[2],".tsv"),sep="\n",row.names=F,col.names=F,quote=F)

pdf(paste0("GCrich_subsets_", appliedfilt, "_", args[3], "_", args[1], "_", args[2], ".pdf"),width=11)
print(ggarrange(sharedsp + rremove("x.text")+rremove("xlab"), gc3cont, ncol=1, nrow=2))
print(vendiag)
dev.off()
