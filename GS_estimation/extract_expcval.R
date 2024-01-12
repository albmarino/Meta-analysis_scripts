#!usr/bin/env Rscript

# build the WLS model from the training dataset (Table S3), extract expected C-values for the full dataset (Table S1) and test reads effect - tab-separated files
# example run: Rscript extract_expcval.R train_dataset.tsv full_dataset.tsv

args <- commandArgs(trailingOnly = TRUE)
library(dplyr)

trainset <- read.table(args[1], header=T,sep="\t")

fullset <- read.table(args[2], header=T,sep="\t") %>% select(Species, Quast_ContigN50, Assembly_size, Reads_assembly)

# WLS model with residuals of a linear model as weights 
m1 <- lm(Cvalue~Assembly_size, data=trainset)
wt_m1 <- 1 / lm(abs(m1$residuals) ~ m1$fitted.values)$fitted.values^2
wls_m1 <- lm(Cvalue ~ Assembly_size, data = args[1], weights=wt_m1)

#Predict values from all assembly sizes with Quast_ContigN50 >= 50kb in the full dataset

n50set <- filter(fullset, Quast_ContigN50 >= 50000)
as <- data.frame(Assembly_size=n50set$Assembly_size)
corrected_as <- data.frame(Species=n50set$Species, Expected_Cvalue = predict(wls_m1,newdata = as))

fullset <- left_join(fullset, corrected_as, by="Species")
write.table(fullset, file="assemblysize_expectedcvalues.tsv", ,quote=F, sep="\t", row.names=F, col.names=T)

# test read type effect

newheader <- "Reads_category"

newdf <- data.frame(matrix(nrow=0, ncol=length(newheader)))
for (i in trainset$Reads_assembly) {
	if (i=="LR" || i == "LR-SR") {newdf <- rbind(newdf, "LR")}
	else if (i=="SR") {newdf <- rbind(newdf, "SR")}
	else {newdf <- rbind(newdf, "")}
	}
colnames(newdf) <- newheader

trainset <- cbind(trainset, newdf)
modreads <- lm(Cvalue ~ Assembly_size + Reads_category + Assembly_size:Reads_category, data=trainset)
wt_modreads <- 1 / lm(abs(modreads$residuals) ~ modreads$fitted.values)$fitted.values^2
wls_modreads <- lm(Cvalue ~ Assembly_size + Reads_category + Assembly_size:Reads_category, data = trainset, weights=wt_modreads)

summary(wls_modreads)

