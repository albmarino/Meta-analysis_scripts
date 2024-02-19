#!/usr/bin/Rscript

# Commands to calculate phylogenetic Independent Contrasts (PIC) among traits and dN/dS - example with Bio++ dN/dS (dNdS_CM_phylter_noshortbr)
# Required input: full phylogeny FULLTREE_brlengths_rooted.treefile and TableS2.tsv
# Rscript pic_Code.R TableS2.tsv FULLTREE_brlengths_rooted.treefile

library(ape)
library(caper)
library(nlme)
library(dplyr)
library(ggplot2)
library(ggpubr)
args = commandArgs(trailingOnly=TRUE)

df <- read.table(args[1],header=T,sep="\t")
fulltree <- read.tree(args[2])
sp2rm <- fulltree$tip.label[! fulltree$tip.label %in% rownames(df)]
prunedtree <- drop.tip(fulltree, sp2rm)
df <-  df[prunedtree$tip.label, ]
prunedtree$node.label = NULL
actinopteri <- filter(df, Clade=="Actinopteri")
aves <- filter(df, Clade=="Aves")
insecta <- filter(df, Clade=="Insecta")
mammalia <- filter(df, Clade=="Mammalia")
mollusca <- filter(df, Clade=="Mollusca")

# PIC data for overall dataset and each clade
df_ic <- comparative.data(prunedtree, df, Species, na.omit = F , vcv=TRUE )
actinopteri_ic <- comparative.data(prunedtree, actinopteri, Species, na.omit = F , vcv=TRUE )
aves_ic <- comparative.data(prunedtree, aves, Species, na.omit = F , vcv=TRUE )
insecta_ic <- comparative.data(prunedtree, insecta, Species, na.omit = F , vcv=TRUE )
mammalia_ic <- comparative.data(prunedtree, mammalia, Species, na.omit = F , vcv=TRUE )
mollusca_ic <- comparative.data(prunedtree, mollusca, Species, na.omit = F , vcv=TRUE )

# LHTs vs dN/dS independent contrasts

ic_logmass_logdnds <- crunch(Mass_g_log~log(dNdS_CM_phylter_noshortbr) , data=df_ic)
summary(ic_logmass_logdnds)
contrast_logmass_logdnds <- caic.table(ic_logmass_logdnds) # for plotting

ic_loglong_logdnds <- crunch(MaxLongevity_y_log~log(dNdS_CM_phylter_noshortbr) , data=df_ic)
summary(ic_loglong_logdnds)
contrast_loglong_logdnds <- caic.table(ic_loglong_logdnds) # for plotting

# Genome size vs total TE content

ic_loggs_logTE <- crunch(log(Unified_Cvalues_methods)~log(Overall_repeat_bp) , data=df_ic)
contrast_loggs_logTE <- caic.table(ic_loggs_logTE) # for plotting
ic_actinopteri_loggs_logTE <- crunch(log(Unified_Cvalues_methods)~log(Overall_repeat_bp) , data=actinopteri_ic)
ic_aves_loggs_logTE <- crunch(log(Unified_Cvalues_methods)~log(Overall_repeat_bp) , data=aves_ic)
ic_insecta_loggs_logTE <- crunch(log(Unified_Cvalues_methods)~log(Overall_repeat_bp) , data=insecta_ic)
ic_mammalia_loggs_logTE <- crunch(log(Unified_Cvalues_methods)~log(Overall_repeat_bp) , data=mammalia_ic)
ic_mollusca_loggs_logTE <- crunch(log(Unified_Cvalues_methods)~log(Overall_repeat_bp) , data=mollusca_ic)

summary(ic_loggs_logTE)
summary(ic_actinopteri_loggs_logTE)
summary(ic_aves_loggs_logTE)
summary(ic_insecta_loggs_logTE)
summary(ic_mammalia_loggs_logTE)
summary(ic_mollusca_loggs_logTE)

# Genome size vs dN/dS 
ic_loggs_logdnds <- crunch(log(Unified_Cvalues_methods)~log(dNdS_CM_phylter_noshortbr) , data=df_ic)
contrast_loggs_logdnds <- caic.table(ic_loggs_logdnds) # for plotting
ic_actinopteri_loggs_logdnds <- crunch(log(Unified_Cvalues_methods)~log(dNdS_CM_phylter_noshortbr) , data=actinopteri_ic)
ic_aves_loggs_logdnds <- crunch(log(Unified_Cvalues_methods)~log(dNdS_CM_phylter_noshortbr) , data=aves_ic)
ic_insecta_loggs_logdnds <- crunch(log(Unified_Cvalues_methods)~log(dNdS_CM_phylter_noshortbr) , data=insecta_ic)
ic_mammalia_loggs_logdnds <- crunch(log(Unified_Cvalues_methods)~log(dNdS_CM_phylter_noshortbr) , data=mammalia_ic)
ic_mollusca_loggs_logdnds <- crunch(log(Unified_Cvalues_methods)~log(dNdS_CM_phylter_noshortbr) , data=mollusca_ic)

summary(ic_loggs_logdnds)
summary(ic_actinopteri_loggs_logdnds)
summary(ic_aves_loggs_logdnds)
summary(ic_insecta_loggs_logdnds)
summary(ic_mammalia_loggs_logdnds)
summary(ic_mollusca_loggs_logdnds)

# Total TE content vs dN/dS
ic_logTE_logdnds <- crunch(log(Overall_repeat_bp)~log(dNdS_CM_phylter_noshortbr) , data=df_ic)
contrast_logTE_logdnds <- caic.table(ic_logTE_logdnds) # for plotting
ic_actinopteri_logTE_logdnds <- crunch(log(Overall_repeat_bp)~log(dNdS_CM_phylter_noshortbr) , data=actinopteri_ic)
ic_aves_logTE_logdnds <- crunch(log(Overall_repeat_bp)~log(dNdS_CM_phylter_noshortbr) , data=aves_ic)
ic_insecta_logTE_logdnds <- crunch(log(Overall_repeat_bp)~log(dNdS_CM_phylter_noshortbr) , data=insecta_ic)
ic_mammalia_logTE_logdnds <- crunch(log(Overall_repeat_bp)~log(dNdS_CM_phylter_noshortbr) , data=mammalia_ic)
ic_mollusca_logTE_logdnds <- crunch(log(Overall_repeat_bp)~log(dNdS_CM_phylter_noshortbr) , data=mollusca_ic)

summary(ic_logTE_logdnds)
summary(ic_actinopteri_logTE_logdnds)
summary(ic_aves_logTE_logdnds)
summary(ic_insecta_logTE_logdnds)
summary(ic_mammalia_logTE_logdnds)
summary(ic_mollusca_logTE_logdnds)

# Recent TE content vs dN/dS
ic_logrecentTE_logdnds <- crunch(log(Overall_TEs_recent_bp)~log(dNdS_CM_phylter_noshortbr) , data=df_ic)
contrast_logrecentTE_logdnds <- caic.table(ic_logrecentTE_logdnds) # for plotting
ic_actinopteri_logrecentTE_logdnds <- crunch(log(Overall_TEs_recent_bp)~log(dNdS_CM_phylter_noshortbr) , data=actinopteri_ic)
ic_aves_logrecentTE_logdnds <- crunch(log(Overall_TEs_recent_bp)~log(dNdS_CM_phylter_noshortbr) , data=aves_ic)
ic_insecta_logrecentTE_logdnds <- crunch(log(Overall_TEs_recent_bp)~log(dNdS_CM_phylter_noshortbr) , data=insecta_ic)
ic_mammalia_logrecentTE_logdnds <- crunch(log(Overall_TEs_recent_bp)~log(dNdS_CM_phylter_noshortbr) , data=mammalia_ic)
ic_mollusca_logrecentTE_logdnds <- crunch(log(Overall_TEs_recent_bp)~log(dNdS_CM_phylter_noshortbr) , data=mollusca_ic)

summary(ic_logrecentTE_logdnds)
summary(ic_actinopteri_logrecentTE_logdnds)
summary(ic_aves_logrecentTE_logdnds)
summary(ic_insecta_logrecentTE_logdnds)
summary(ic_mammalia_logrecentTE_logdnds)
summary(ic_mollusca_logrecentTE_logdnds)
