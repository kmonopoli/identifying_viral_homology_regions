# install.packages("rJava")
# install.packages("xlsxjars")
# install.packages("xlsx")
# install.packages("motifStack")
# install.packages("rio")
# install.packages("dplyr")
# 
# library(readxl)
# library("rJava")
# library("xlsxjars")
# library("xlsx")
# library(motifStack)
# library("rio")
# library(dplyr)



setwd("/Users/kmonopoli/Desktop/covid19_predictions")
All=read_excel("COVID-19_Sequences.xlsx", sheet="human_COVID19__output", col_names=T)
All=as.data.frame(All)
All=as.data.frame(All[72:29909])
colnames(All) = c("Oligo_count", "Oligo_ID", "Sequence", "Gene_region", "GC_content", "Revcomp", "Oligo_ID_sense", "Oligo_ID_antisense", "Score", "Region", "Hits2microRNA", "refseq_seed", "hits2accessions", "hits2gIDs", "gIDs", "Accessions", "hits2mouse", "hits2nhp", "hits2human", "hits2rat","homology_score")#,

## Get rid of sequences with hts2gIDs>1
All = subset(All, hits2gIDs <= 1, select = c(Oligo_count, Oligo_ID, Sequence, Gene_region, GC_content, Revcomp, Oligo_ID_sense, Oligo_ID_antisense, Score, Region,Hits2microRNA, refseq_seed, hits2accessions, hits2gIDs, gIDs, Accessions, hits2mouse, hits2nhp, hits2human, hits2rat,homology_score))

## Get rid of sequences with miRNA hits > 1
All = subset(All, Hits2microRNA < 1 , select = c(Oligo_count, Oligo_ID, Sequence, Gene_region, GC_content, Revcomp, Oligo_ID_sense, Oligo_ID_antisense, Score, Region,Hits2microRNA, refseq_seed, hits2accessions, hits2gIDs, gIDs, Accessions, hits2mouse, hits2nhp, hits2human, hits2rat,homology_score))

## Get rid of sequences with GC content > 56%
All = subset(All, GC_content <= 56, select = c(Oligo_count, Oligo_ID, Sequence, Gene_region, GC_content, Revcomp, Oligo_ID_sense, Oligo_ID_antisense, Score, Region,Hits2microRNA, refseq_seed, hits2accessions, hits2gIDs, gIDs, Accessions, hits2mouse, hits2nhp, hits2human, hits2rat,homology_score))


#get rid of GGGG, CCCC, AAAA, or UUUU stretches
remove.list <- paste(c("GGGG", "CCCC", "AAAA", "UUUU"), collapse = '|')
All = All %>% filter(!grepl(remove.list, Sequence))


All$Oligo_count <- as.numeric(as.character(All$Oligo_count))
All$Score <- as.numeric(as.character(All$Score))
All$homology_score <- as.numeric(as.character(All$homology_score))
All$ScoCOVID_Genere <-""
All$score_plus <-All$Score+(0.5)*All$homology_score # combination of homology score and weight matrix score




#Select orf1a
orf1a = subset(All, Oligo_count >=266 & Oligo_count<= 13483, select = c(Oligo_count, Oligo_ID, Sequence, Gene_region, GC_content, Revcomp, Oligo_ID_sense, Oligo_ID_antisense, Score, Region,Hits2microRNA, refseq_seed, hits2accessions, hits2gIDs, gIDs, Accessions, hits2mouse, hits2nhp, hits2human, hits2rat,homology_score,score_plus))
#Sort by Score then homology_score
# orf1a= orf1a[order(orf1a$Score,orf1a$homology_score,decreasing=TRUE),]
orf1a= orf1a[order(orf1a$score_plus,decreasing=TRUE),]

orf1a=orf1a[1:12,]
# add column to indicate which COVID gene it targets
orf1a$COVID_Gene = "orf1a"
  
  
#Select orf1ab
orf1ab = subset(All, Oligo_count >=13483 & Oligo_count<= 21555, select = c(Oligo_count, Oligo_ID, Sequence, Gene_region, GC_content, Revcomp, Oligo_ID_sense, Oligo_ID_antisense, Score, Region,Hits2microRNA, refseq_seed, hits2accessions, hits2gIDs, gIDs, Accessions, hits2mouse, hits2nhp, hits2human, hits2rat,homology_score,score_plus))
# orf1ab= orf1ab[order(orf1ab$Score,orf1ab$homology_score,decreasing=TRUE),]
orf1ab= orf1ab[order(orf1ab$score_plus,decreasing=TRUE),]

orf1ab=orf1ab[1:12,]
# add column to indicate which COVID gene it targets
orf1ab$COVID_Gene = "orf1ab"

#Select S
S = subset(All, Oligo_count >=21563 & Oligo_count<= 25384, select = c(Oligo_count, Oligo_ID, Sequence, Gene_region, GC_content, Revcomp, Oligo_ID_sense, Oligo_ID_antisense, Score, Region,Hits2microRNA, refseq_seed, hits2accessions, hits2gIDs, gIDs, Accessions, hits2mouse, hits2nhp, hits2human, hits2rat,homology_score,score_plus))
# S= S[order(S$Score,S$homology_score,decreasing=TRUE),]
S= S[order(S$score_plus,decreasing=TRUE),]

S=S[1:12,]
# add column to indicate which COVID gene it targets
S$COVID_Gene = "S"

#Select 3a
ThreeA = subset(All, Oligo_count >=25393 & Oligo_count<= 26220, select = c(Oligo_count, Oligo_ID, Sequence, Gene_region, GC_content, Revcomp, Oligo_ID_sense, Oligo_ID_antisense, Score, Region,Hits2microRNA, refseq_seed, hits2accessions, hits2gIDs, gIDs, Accessions, hits2mouse, hits2nhp, hits2human, hits2rat,homology_score,score_plus))
# ThreeA= ThreeA[order(ThreeA$Score,ThreeA$homology_score,decreasing=TRUE),]
ThreeA= ThreeA[order(ThreeA$score_plus,decreasing=TRUE),]

ThreeA=ThreeA[1:12,]
# add column to indicate which COVID gene it targets
ThreeA$COVID_Gene = "3a"

#Select E
E = subset(All, Oligo_count >=26245 & Oligo_count<= 26472, select = c(Oligo_count, Oligo_ID, Sequence, Gene_region, GC_content, Revcomp, Oligo_ID_sense, Oligo_ID_antisense, Score, Region,Hits2microRNA, refseq_seed, hits2accessions, hits2gIDs, gIDs, Accessions, hits2mouse, hits2nhp, hits2human, hits2rat,homology_score,score_plus))
# E= E[order(E$Score,E$homology_score,decreasing=TRUE),]
E= E[order(E$score_plus,decreasing=TRUE),]


E=E[1:12,]
# add column to indicate which COVID gene it targets
E$COVID_Gene = "E"

#Select M
M = subset(All, Oligo_count >=26523 & Oligo_count<= 27191, select = c(Oligo_count, Oligo_ID, Sequence, Gene_region, GC_content, Revcomp, Oligo_ID_sense, Oligo_ID_antisense, Score, Region,Hits2microRNA, refseq_seed, hits2accessions, hits2gIDs, gIDs, Accessions, hits2mouse, hits2nhp, hits2human, hits2rat,homology_score,score_plus))
# M= M[order(M$Score,M$homology_score,decreasing=TRUE),]
M= M[order(M$score_plus,decreasing=TRUE),]

M=M[1:12,]
# add column to indicate which COVID gene it targets
M$COVID_Gene = "M"

#Select 7a
SevenA = subset(All, Oligo_count >=27394 & Oligo_count<= 27759, select = c(Oligo_count, Oligo_ID, Sequence, Gene_region, GC_content, Revcomp, Oligo_ID_sense, Oligo_ID_antisense, Score, Region,Hits2microRNA, refseq_seed, hits2accessions, hits2gIDs, gIDs, Accessions, hits2mouse, hits2nhp, hits2human, hits2rat,homology_score,score_plus))
# SevenA= SevenA[order(SevenA$Score,SevenA$homology_score,decreasing=TRUE),]
SevenA= SevenA[order(SevenA$score_plus,decreasing=TRUE),]

SevenA=SevenA[1:12,]
# add column to indicate which COVID gene it targets
SevenA$COVID_Gene = "7a"

#Select 8b
EightB = subset(All, Oligo_count >=27894 & Oligo_count<= 28259, select = c(Oligo_count, Oligo_ID, Sequence, Gene_region, GC_content, Revcomp, Oligo_ID_sense, Oligo_ID_antisense, Score, Region,Hits2microRNA, refseq_seed, hits2accessions, hits2gIDs, gIDs, Accessions, hits2mouse, hits2nhp, hits2human, hits2rat,homology_score,score_plus))
# EightB= EightB[order(EightB$Score,EightB$homology_score,decreasing=TRUE),]
EightB= EightB[order(EightB$score_plus,decreasing=TRUE),]

EightB=EightB[1:12,]
# add column to indicate which COVID gene it targets
EightB$COVID_Gene = "8b"




#Select N
N = subset(All, Oligo_count >=28274 & Oligo_count<= 29533, select = c(Oligo_count, Oligo_ID, Sequence, Gene_region, GC_content, Revcomp, Oligo_ID_sense, Oligo_ID_antisense, Score, Region,Hits2microRNA, refseq_seed, hits2accessions, hits2gIDs, gIDs, Accessions, hits2mouse, hits2nhp, hits2human, hits2rat,homology_score,score_plus))
# N= N[order(N$Score,N$homology_score,decreasing=TRUE),]
N= N[order(N$score_plus,decreasing=TRUE),]

N=N[1:12,]
# add column to indicate which COVID gene it targets
N$COVID_Gene = "N"



# 
# # ID overlaps in orf1a and orf1ab top hits
# Overlap_orf1a_orf1ab = merge(orf1a, orf1ab, by = "Oligo_count")
# Overlap_orf1a_orf1ab = Overlap_orf1a_orf1ab[,1:22]
# colnames(Overlap_orf1a_orf1ab) = c("Oligo_count", "Oligo_ID", "Sequence", "Gene_region", "GC_content", "Revcomp", "Oligo_ID_sense", "Oligo_ID_antisense", "Score", "Region", "Hits2microRNA", "refseq_seed", "hits2accessions", "hits2gIDs", "gIDs", "Accessions", "hits2mouse", "hits2nhp", "hits2human", "hits2rat","homology_score","COVID_Gene")
# h = matrix(nrow=nrow(Overlap_orf1a_orf1ab), ncol=1)
# h[,]= "orf1a_orf1ab"
# h = as.data.frame(h)
# colnames(h)="Gene_Target"
# 
# 
# #ID unique orf1a and orf1ab hits
# orf1a_only = anti_join(orf1a, Overlap_orf1a_orf1ab)
# i = matrix(nrow=nrow(orf1a_only), ncol=1)
# i[,]= "orf1a"
# i = as.data.frame(i)
# colnames(i)="Gene_Target"
# 
# orf1ab_only = anti_join(orf1ab, Overlap_orf1a_orf1ab)
# j = matrix(nrow=nrow(orf1ab_only), ncol=1)
# j[,]= "orf1ab"
# j = as.data.frame(j)
# colnames(j)="Gene_Target"

#join all gene target names
# Gene_Target = rbind(orf1a_only,orf1ab_only,Overlap_orf1a_orf1ab,S,ThreeA,E,M,SevenA,EightB,N)

#join all sequence info
All_Top = rbind(orf1a,orf1ab,S,ThreeA,E,M,SevenA,EightB,N)
# All_Top = cbind(Gene_Target, All_Top)

export(All_Top, "Human_COVID_19_Sequences.xlsx")


