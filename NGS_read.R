
library(tidyverse)
library(ggrepel)
library(dplyr)
library(DESeq2)
library(EnhancedVolcano)

#Specify pair of Conditions - should match row labels in csv

A_cond <- 'Spec0'
B_cond <- 'Spec800'
target <- 'TCPTPcat/TCPTPfull'

read_thresh <- 6

#path to spreadsheet
NGS_cts_df <- data.frame(read.csv("/Users/nickles/Desktop/TBS/NGS/NGS_cts_test_TC.csv", header = TRUE))


#NGS_cts_df <- data.frame(read.csv("/Users/nickles/Desktop/NGS/NGS_cts_test1.csv", header = TRUE, nrows = 5))

NGS_cts_df <- setNames(data.frame(t(NGS_cts_df[,-1])), paste0(NGS_cts_df[,1]))

rownames(NGS_cts_df) <- substr(rownames(NGS_cts_df), 1, nchar(rownames(NGS_cts_df))-1) # remove . in names

Areps <- c(paste0(A_cond, '_R1'), paste0(A_cond, '_R2'))
Breps <- c(paste0(B_cond, '_R1'), paste0(B_cond, '_R2'))

# removal of members prior to deseq?
NGS_cts_filt <- subset(NGS_cts_df, Areps[1] > read_thresh & Areps[2] > read_thresh & 
                             Breps[1] > read_thresh & Breps[2] > read_thresh) #filter plasmids with few reads

countdata <- as.matrix(bind_cols(NGS_cts_filt[Areps[1]], NGS_cts_filt[Areps[2]], NGS_cts_filt[Breps[1]], NGS_cts_filt[Breps[2]]))
countdata <- as.matrix(bind_cols(NGS_cts_df[Areps[1]], NGS_cts_filt[Areps[2]], NGS_cts_filt[Breps[1]], NGS_cts_filt[Breps[2]]))
