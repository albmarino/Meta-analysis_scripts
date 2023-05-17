#!/usr/bin/Rscript

# use sorted_blast3 and: (1) weight reads counts with the total read count from renamed.blasting_reads.fasta; (2) use also blasted reads bp (to weight by genome size if we want). (2) probably not necessary unless we prefer a measure in bp instead of genome percentage
# Usage --> Rscript dnapt_recentTEs.R 5 $species $startpath $outfile # numeric argument with the maximum blastn divergence to sample as recent TEs, species, initial working directory, and output filename to append recent TEs info to
# output table with one species for each line, with overall and subclass by subclass TEs (bp and genome coverage?)  
args = commandArgs(trailingOnly=TRUE)
library(dplyr)
library(tidyr)
library(ggplot2)

# edit of dnaPT_utils script

bashcmd <- "join -a1 -12 -21 -o 1.3,1.4,2.4,2.5  sorted_blast3 one_RM_hit_per_Trinity_contigs | awk '{if (NF == 2) {print $1, \"\\t\", $2\"\\tUnknown_repeat\\tUnknown\"} else {print $0}}' | awk '/LINE/ { print $0 \"\\t\" $4; next} /LTR/ {print $0 \"\\t\" $4; next} /SINE/ {print $0 \"\\tSINE\"; next} /DNA/ {print $0 \"\\tDNA\"; next} /MITE/ {print $0 \"\\tMITE\";next} /RC/ {print $0 \"\\tRC\";next} /Unknown/ {print $0 \"\\tUnknown\";next} !/Simple_repeat|Low_complexity|Satellite|srpRNA|rRNA|tRNA|snRNA|ARTEFACT/ {if (NF == 4) {print $0\"\\tOthers\"} else {print $0\"\\tNA\\tUnknown\\tUnknown\"}}' | sed 's/ /\t/g;s/\t\t\t/\\t/g'"
df<-read.table(sep = "\t", text=system(bashcmd, intern = T))
df<- separate(df, V4, c("TE_subclass", "TE_superfamily"), sep = "/",fill = "right")
df$V1 <- 100-df$V1 # convert similarity to blastn distance
colnames(df)<-c("blastndiv", "Bp_length", "RM_annotation","TE_subclass", "TE_superfamily", "TE_family")
reads.c<-as.numeric(system("grep -c '>' ../renamed.blasting_reads.fasta", intern = T))
subclasslist <- sort(scan(file="subclasslist", what="character"))
outdf = data.frame(matrix(ncol = 1, nrow = length(subclasslist)))
colnames(outdf) <- "TE_subclass"
outdf$TE_subclass<- subclasslist

# Check that the landscape is the same as the one produced by dnaPT_utils 
#p_copycount<-ggplot(df, aes(x=blastndiv, fill=TE_subclass))+
#geom_histogram(aes(y=..count../reads.c*100), binwidth = 1)+
#labs(x="blastn divergence", y="Genome percentage (% reads count/total reads count)")

#pdf("tryplot_copycount.pdf", width=11)
#print(p_copycount)
#dev.off()

# Check landscape with bplength on y
bin <- seq(0,100,1)
df_discretebins <- data.frame()
for (i in 1:(length(bin)-1)) {
j = i-1
df_bin <- filter(df, blastndiv >= bin[i] & blastndiv < bin[i+1])
df_bin <- df_bin %>% mutate(blastndiv_bin = rep((i-1),nrow(df_bin)))
df_discretebins <- rbind(df_discretebins, df_bin)
}


#p_bpcount <- ggplot(df_discretebins_aggr, aes(x=blastndiv_bin, y=Length_bp, fill=TE_subclass))+
#geom_bar(stat="identity")+
#labs(x="blastn divergence", y="reads coverage (bp)")

#pdf("tryplot_bplength.pdf", width=11)
#print(p_bpcount)
#dev.off()

# Subset TE info for blastn divergence <5%: use df_discretebins for copy counts info

recent_df <- filter(df_discretebins, blastndiv<as.numeric(args[1]))
overall_te <- nrow(recent_df)/reads.c*100

vec <- c(overall_te)
names(vec) <- "Overall_TEs"
i = 2
for (subcl in unique(recent_df$TE_subclass)) {
	vec <- append(vec, filter(recent_df, TE_subclass==subcl) %>% nrow/reads.c*100)
	names(vec)[i]<-paste0(subcl)
	i = i + 1
	}

outdf$vec<-vec[match(outdf$TE_subclass,names(vec))]
#colnames(outdf)[colnames(outdf)=="vec"] <- args[2]
rownames(outdf)<- paste0(outdf$TE_subclass, "_recent")
#lineout <- as.data.frame(t(outdf)) %>% mutate(Species=args[2], .before=1) %>% filter(!row_number() %in% 1)
lineout <- c(args[2], outdf$vec)
setwd(args[3])
write(lineout, sep ="\t", append=TRUE, file=args[4], ncolumns=length(lineout))

#recent_bp_df <- filter(df_discretebins_aggr, blastndiv_bin<=as.numeric(args[1]))
#for (subcl in unique(recent_bp_df$TE_subclass)) {
	
#	assign(paste0(subcl, "_recent_bp"), filter(recent_bp_df, TE_subclass==subcl) %>% pull(Length_bp) %>% sum)
#	}
