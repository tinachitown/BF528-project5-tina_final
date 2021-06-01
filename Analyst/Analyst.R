library(ggplot2)
library(dplyr)
library(plyr)

####    ANALYST ROLE 

setwd('/usr4/bf528/bcole99/BF528-project5-tina_final/project2_final/')

#6.1 Load gene_exp fil
dfx<-read.table('gene_exp.diff', header = TRUE)
df<-arrange(dfx, q_value, decreasing = TRUE)%>%select(gene, value_1, value_2, log2.fold_change., p_value, q_value)
DE_10<-head(df,10)

#Fixing column names
colnames(df)[colnames(df) == "log2.fold_change."] <- "log2fc"
colnames(df)



#Fixing margins
par(mar=c(1,1,1,1))

#6.2 produce histogram of logfc  for ALL genes
hist(df$log2.fold_change.,breaks = 50, col = "Blue", main = "Histogram of logfc on all genes")

#DE genes p < 0.01 for deliverables section 
df_pvalue01<-filter(df, p_value<0.01)
up_reg_p01<-filter(up_reg, p_value<0.01)
down_regp01<-filter(down_reg, p_value<0.01)

#There are 1999 genes at < 0.01 p-value


#6.3 Create new dataframe with ONLY genes, last column labeled significant (listing  ONLY "yes" results)
df.sub<-subset(df, df$significant == "yes")

##This results in 1,637 from the original data set of 36,329


#filter by important values: Subset by log fold change > 1.5, then select top 10 DE genes
df.logfc<-subset(df.sub, df.sub$log2.fold_change. > 1.5)



##This results in 527 from the subset of 1,637 

write.csv(DE_10<-select(df.logfc, gene, log2.fold_change., significant, p_value)%>%arrange(p_value)%>%head(10)%>%arrange(log2.fold_change.), "Top_10_DE_Genes.csv")


#6.4 Create another histogram with logfc ONLY for significant genes 
hist(df.sub$log2.fold_change.,breaks = 30, col = "Green", main = "Histogram of logfc on all Significant genes", freq = TRUE, axes =  TRUE)


#6.5 Subset significant gene dataframe into 2 dataframes, #1 UP #2 DOWN regulated genes using log2.foldchange column. (Use the # of UP / DOWN regulated genes in final report)
up_reg<-subset(df.sub, df.sub$log2.fold_change. >= 0)
down_reg<-subset(df.sub, df.sub$log2.fold_change. < 0)

#There are 995 genes down regulated (less than 0 fc) and there are 642 genes up-regulated = > 0

#Histogram of up/ down regulated genes . As described in paper, the reverse log fold change is down and up is positive 
hist(up_reg$log2.fold_change.,breaks = 30, col = "Yellow", main = "Histogram of logfc Up-Regulated", freq = TRUE, axes =  TRUE)
hist(down_reg$log2.fold_change.,breaks = 30, col = "Purple", main = "Histogram of logfc Down-Regulated", freq = TRUE, axes =  TRUE)



#6.8 Save the UP/DOWN regulated genes using write function 
write.csv(up_reg, "Up_Regulated.csv")
write.csv(down_reg, "Down_Regulated.csv")

#tina

#6.7 Use DAVID gene enrichment software to run analysis on both UP / DOWN regulated genes. Save the tables for the final report.

## As of up-regulated genes, of the 642 genes 601 were identified within David for Mus musculus species. Running under medium GOTERM_BP_FAT, GOTERM_CC_FAT, and 
##GOTERM_MF_FAT resulted in 294 clusters, with 13 clusters above 4.0 Enrichment score.  There were 492 terms not clustered. 

## As of down-regulated genes, of the 995 genes 950 were identified within David for Mus musculus species. Running under medium GOTERM_BP_FAT, GOTERM_CC_FAT, and 
##GOTERM_MF_FAT resulted in 390 clusters, with 25 above 5.0 Enrichment score.  725 terms were unclustered.  


