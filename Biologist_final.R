suppressPackageStartupMessages({
  install.packages('caret')
  install.packages('devtools')
  devtools::install_github("riatelab/cartography")
  install.packages("cartography")
  install.packages("eulerr")   
  install.packages("ggdendro")
  install.packages("ggsci")
  install.packages('hrbrthemes')
  hrbrthemes::import_roboto_condensed()
  install.packages('panelr')
  install.packages("patchwork")
  install.packages("reshape2")
  install.packages("R.filesets") 
  install.packages("Seurat")
  install.packages("shiny")  
  install_github('SingleCellAssay', 'RGLab', build_vignettes=FALSE)
  vignette('SingleCellAssay-intro')
  install.packages('UpSetR')
  install.packages("viridisLite")
 })


suppressPackageStartupMessages({
  library(caret)
  library(cowplot) 
  library(DescTools)
  library(devtools)
  library(dplyr) 
  library(glue) 
  library(GGally)
  library(gridExtra)
  library(ggplot2)
  library(ggdendro)
  library(ggrepel)
  library(grid)
  library(hrbrthemes) 
  library(magrittr)
  library(Matrix)
  library(panelr)
  library(patchwork)
  library(pheatmap)
  library(purrr)
  library(RColorBrewer)
  library(reshape2)
  library(Rtsne)
  library(R.filesets) 
  library(Seurat) 
  library(tidyverse)
  library(viridis)
  
})
#1. Plot fold change for genes from Figure 1D, compare/contrast to part
## Read in P0 v AD and P4 v P7 to acquire fold change values
P0_AD<-read.delim2('gene_P01.diff')
P4_P7<-read.delim2('gene_P4vP7.diff')
 

#Genes to select for sarcomere, mitochondria and cell cycle genes
s_genes = c("Tcap","Pygm", "Myoz2", "Des", "Csrp3", "Cryab")
m_genes = c("Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh")
ccx_genes = c("Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "E2f1", "Cdc27", "Cdc45", "Rad51","Aurkb", "Cdc23")

#STEP 1. Filter by genes from paper - Sarcomere
P0_ADs<-filter(P0_AD, gene %in% s_genes)%>%select(gene, value_1, value_2)
P4_P7s<-filter(P4_P7, gene %in% s_genes)%>%select(gene, value_1, value_2)
gene = P0_ADs$gene
b = as.numeric(P0_ADs$value_1)
c = as.numeric(P0_ADs$value_2)
d = as.numeric(P4_P7s$value_1)
e = as.numeric(P4_P7s$value_2)
P0 = log2(b)
AD= log2(c)
P4 = log2(d)
P7=  log2(e)
sar<-tbl_df(data.frame(gene,P0, P4, P7,AD))


#REPEAT for Mitochondria and Cell Cycle
##Mitochondria 
P0_ADm<-filter(P0_AD, gene %in% m_genes)%>%select(gene, value_1, value_2)
P4_P7m<-filter(P4_P7, gene %in% m_genes)%>%select(gene, value_1, value_2)

gene_m = P0_ADm$gene
j = as.numeric(P0_ADm$value_1)
k = as.numeric(P0_ADm$value_2)
l = as.numeric(P4_P7m$value_1)
m = as.numeric(P4_P7m$value_2)
P0m = log2(j)
ADm= log2(k)
P4m = log2(l)
P7m=  log2(m)
mit<-tbl_df(data.frame(gene_m,P0m, P4m, P7m,ADm))

##Cell cycle
P0_ADc<-filter(P0_AD, gene %in% ccx_genes)%>%select(gene, value_1, value_2)
P4_P7c<-filter(P4_P7, gene %in% ccx_genes)%>%select(gene, value_1, value_2)

gene_c = P0_ADc$gene
w = as.numeric(P0_ADc$value_1)
x = as.numeric(P0_ADc$value_2)
y = as.numeric(P4_P7c$value_1)
z = as.numeric(P4_P7c$value_2)
P0c = log2(w)
ADc= log2(x)
P4c = log2(y)
P7c=  log2(z)
cc<-tbl_df(data.frame(gene_c,P0c,P4c,P7c,ADc ))

##Plotting all Genes


#Fixing margins
par(mar=c(1,1,1,1))

#plot Sarcomere

ggparcoord(sar,
           columns = 2:5, groupColumn = 1, 
           scale = "std",
           showPoints = TRUE, 
           title = "Comparing Results: Sarcomere",
           alphaLines = 0.6
) + 
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum()+
 
  xlab("Sample") + ylab("log2") + theme_dark()


#plot Mitochondria

ggparcoord(mit,
           columns = 2:5, groupColumn = 1, 
           scale = "std",
           showPoints = TRUE, 
           title = "Comparing Results: Mitochondria",
           alphaLines = 0.6
) + 
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum()+
  
  xlab("Sample") + ylab("log2") + theme_light()

#plot Cell Cycle

ggparcoord(cc,
           columns = 2:5, groupColumn = 1, 
           scale = "std",
           showPoints = TRUE, 
           title = "Comparing Results: Cell Cycle",
           alphaLines = 0.6
) + 
  scale_color_viridis(discrete=TRUE) +
  theme_ipsum()+
  
  xlab("Sample") + ylab("log2") + theme_ft_rc()

#7.2 Expand table to include 6 enriched biological pathways that comparison with paper
 
##Read in DAVID enrichment genes for sarcomere, mitochoria & cell cycle genes
sar_up_genes<-read.csv("sarcomere.csv")
mito_up_genes<-read.csv("mitochondrian.csv")
cell_cycle_down_genes<-read.csv("cell_cycle.csv")

#Genes to select for sarcomere, mitochondria & cell cycle genes
sarcomere = c("Jph2", "Pdlim5", "ltgb1pb2", "Pgm5", "Synpo2", "Tcap","Pygm", "Myoz2", "Des", "Csrp3", "Cryab")
mitochondria = c("Uqcrfs1", "Ndufa8", "Dist", "Sdhc", "Sucla2", "Cs", "Prdx3", "Acat1", "Echs1", "Slc25a11", "Phyh")
cell_cycle = c("Cdc25a", "Cinp", "Cdc5l", "Mapk6", "Cdc123", "Mapk1", "Cdc7", "E2f8", "Cdk7", "Cdc26", "Cdc6", "E2f1", "Cdc27", "Bora", "Cdc45", "Rad51","Aurkb", "Cdc23")

#Filter DAVID results for sarcomere, mitochondria & cell cycle genes
sar<-filter(sar_up_genes, ID %in% sarcomere)
mit<-filter(mito_up_genes, ID %in% mitochondria)
cc<-filter(cell_cycle_down_genes, ID %in% cell_cycle)


#Read in top up and down regulated genes
up_reg_genes<-read.csv("Up-Regulated.csv")
down_genes<-read.csv("Down-Regulated.csv")

#Cleaning up for extra column & creating strings for matched results 
up_reg_genes<-up_reg_genes[,-9]
down_genes<-down_genes[-(10:14),]

m<-mit$ID 
sx<-sar$ID
s<-c("Pdlim5, Jph2, Mhyoz, Tcap")
cx<-cc$ID
c<-c("Cdc23,E2f1,E2f8,Rad51,Aurkb,Cdc7")


#Adding column with matching gene names from paper
for(i in 1:nrow(up_reg_genes)) {
  if(up_reg_genes$Enrichment.Term[i] == "mitochondrion") {
  up_reg_genes$Match[i]<-m 
  } else if (up_reg_genes$Enrichment.Term[i] == "sarcomere") {
  up_reg_genes$Match[i]<-s 
  } else up_reg_genes$Match[i]<- "N/A" 
}

#Down regulated matched genes from paper

for(i in 1:nrow(down_genes)) {
  if(down_genes$Enrichment.Term[i] == "cell cycle") {
    down_genes$Match[i]<-c 
   } else down_genes$Match[i]<- "N/A" 
}



#Save updated, matched top up_regulated and down_regulated DAVID results
write.csv(up_reg_genes, "Up_Regulated_Genes_Matched_David.csv")
write.csv(down_genes, "Down_Regulated_Genes_Matched_David.csv")



#3a. Create new FPKM matrix by combining the FPKM columns from each tracking table from 7.1
  ##Use multimerge to join greater than 3 dataframes
#Filter dataset for mitochondria genes

##Read in FPKM tables
P0_1_FPKM<-read.delim2("genes.fpkm_P0_1_tracking")
P0_2_FPKM<-read.delim2("genes.fpkm_P0_2_tracking")
AD_1_FPKM<-read.delim2("genes.fpkm_AD_1_tracking")
AD_2_FPKM<-read.delim2("genes.fpkm_AD_2_tracking")
P4_1_FPKM<-read.delim2("genes.fpkm_P4_1_tracking")
P4_2_FPKM<-read.delim2("genes.fpkm_P4_2_tracking")
P7_1_FPKM<-read.delim2("genes.fpkm_P7_1_tracking")
P7_2_FPKM<-read.delim2("genes.fpkm_P7_2_tracking")

#Read in gene_exp.diff, header = TRUE
DE_diff<-read.table('gene_P01.diff', header = TRUE)
DE_genes<-filter(DE_diff, p_value < 0.05)
 
###When selecting for DE genes at p_value < 0.05 the dataset went from 36,329 to 3,597

#Get DE gene names of DE_genes and filter all dataframes, then combine FPKM
DE_names<-DE_genes%>%select(gene)

#Filter dataframes by DE_genes - result was 4524 genes 
P01<-filter(P0_1_FPKM, gene_short_name %in% DE_names$gene)%>%select(gene_short_name, FPKM)
P02<-filter(P0_2_FPKM, gene_short_name %in% DE_names$gene)%>%select(gene_short_name, FPKM)
P4_1<-filter(P4_1_FPKM, gene_short_name %in% DE_names$gene)%>%select(gene_short_name, FPKM)
P4_2<-filter(P4_2_FPKM, gene_short_name %in% DE_names$gene)%>%select(gene_short_name, FPKM)
P7_1<-filter(P7_1_FPKM, gene_short_name %in% DE_names$gene)%>%select(gene_short_name, FPKM)
P7_2<-filter(P7_2_FPKM, gene_short_name %in% DE_names$gene)%>%select(gene_short_name, FPKM)
AD_1<-filter(AD_1_FPKM, gene_short_name %in% DE_names$gene)%>%select(gene_short_name, FPKM)
AD_2<-filter(AD_2_FPKM, gene_short_name %in% DE_names$gene)%>%select(gene_short_name, FPKM)

## remove duplicates               
P01x<-P01[!duplicated(P01$gene_short_name), ]
P02x<-P02[!duplicated(P02$gene_short_name), ]
P4_1x<-P4_1[!duplicated(P4_1$gene_short_name), ]
P4_2x<-P4_2[!duplicated(P4_2$gene_short_name), ]
P7_1x<-P7_1[!duplicated(P7_1$gene_short_name), ]
P7_2x<-P7_2[!duplicated(P7_2$gene_short_name), ]
AD1x<-AD_1[!duplicated(AD_1$gene_short_name), ]
AD2x<-AD_2[!duplicated(AD_2$gene_short_name), ]

#### After removing duplicates, we end up with 3519 genes per data set. 


# change column names
colnames(P01x)<-c("gene", "P01")
colnames(P02x)<-c("gene", "P02")
colnames(P4_1x)<-c("gene","P4_1")
colnames(P4_2x)<-c("gene","P4_2")
colnames(P7_1x)<-c("gene", "P7_1")
colnames(P7_2x)<-c("gene", "P7_2")
colnames(AD1x)<-c("gene", "AD_1")
colnames(AD2x)<-c("gene", "AD_2")

#Merge dataframes


FPKM_all<-MultMerge(P01x, P02x, P4_1x,P4_2x, P7_1x, P7_2x, AD1x, AD2x, by = "gene")

#3b. Subset top 1000 DE genes (5.4) Reduce list by p_value < 0.01. 

DE_genes01<-filter(DE_genes, p_value < 0.01, DE_genes$significant == "yes")
##Noted: All the genes denoted as "yes" for significance have p_value <0.01.  Reduces list to 1,637 genes

#Get top 1000 DE genes, running a little more to get top 1000 for heat map 
top1k_DE<-top_n(DE_genes01, wt= p_value, -1015)  
#get names 
top1k_names<-top1k_DE$gene


#3c. Clustered heatmap of FPKM values, top 1000 DE genes found in P0 v Ad.
#Genes on rows and samples on columns w/dendrogram

###Filter FPKM_all by top 1000 genes
top_heat<-filter(FPKM_all, gene %in% top1k_names)

# Run clustering
top.matrix <- data.matrix(top_heat)
top.matrix<-top.matrix[,-1]
rownames(top.matrix)<-top_heat$gene

#Scale matrix
top.matrix.scaled<-scale(top.matrix[,2:9])
rownames(top.matrix.scaled)<-top_heat$gene

#Fixing margins - ALWAYS RUN BEFORE PLOTTING
par(mar=c(1,1,1,1))
dev.off()

#Make heatmap
heatmap.plot<-pheatmap(top.matrix.scaled, color = colorRampPalette(c("navy", "white", "firebrick3"))(50), cexRow=0.5, cexCol=0.1, margins = c(5,5),
          scale="row",fontsize_row = 5, cluster_rows = TRUE, main = "Top 1000 Genes")
# Create dendro
plot.dendro <- as.dendrogram(hclust(d = dist(x = top.matrix.scaled)))

# Create dendro
dendro.plot <- ggdendrogram(data = plot.dendro, rotate = TRUE)
dendro.plot <- dendro.plot + theme(axis.text.y = element_text(size = 4))
# Preview the plot
print(dendro.plot)

grid.newpage()
print(heatmap.plot, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro.plot, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))


