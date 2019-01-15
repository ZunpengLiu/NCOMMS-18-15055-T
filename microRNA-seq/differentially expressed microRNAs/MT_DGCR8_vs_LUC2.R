library(DESeq2)
library(pheatmap)
library(readr)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)

merge_mature <- read.csv("/Volumes/Zunpeng/08_Result/04_DGCR8/microRNA_batch2_20181120/merge.mature.csv",sep='\t')


dataFrame<-merge_mature[,1:7]

DESeq2_result_prefix<-"MT_DGCR8_vs_LUC"

countMatrix <- as.matrix(dataFrame[2:7])
# countMatrix <- as.matrix(dataFrame[2:5])
rownames(countMatrix) <- dataFrame$miRNA_ID
condition <- factor(c(rep("LUC",3), rep("MT_DGCR8",3)), levels = c("LUC","MT_DGCR8"))

# Build the Sample Information Table. First method, sample information table, also named as colData, or design matrix
# sample names
sampleNames <- colnames(countMatrix)

# Build the Object colData
# Build the Sample Information table. Second method, sample information table, also named as colData, or design matrix
colData <- data.frame(Samplename= sampleNames, condition)
# add the rowname for colData using the vector sampleNames
rownames(colData) <- sampleNames

######################### Build the Object dds (DESeqDataSet)
dds <- DESeqDataSetFromMatrix(countMatrix, colData, design= ~ condition)

dds <- dds[rowSums(counts(dds)) > 1,]

######################### normalize data
dds2 <- DESeq(dds)

# acquire the results using function results(), and assign to res
res <- results(dds2)

res <- res[order(res$padj),]
resdata <- merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)
names(resdata)[1]<-"MiRBase_ID"

diff_microRNA <-subset(resdata, padj < 0.05 & (log2FoldChange >= 1 | log2FoldChange <= -1))
write.csv(resdata, paste(DESeq2_result_prefix,"_all_microRNA.csv",sep=""),row.names = F)
write.csv(diff_microRNA, paste(DESeq2_result_prefix,"_diff_microRNA.csv",sep=""),row.names = F)
up_microRNA<-subset(diff_microRNA,diff_microRNA$log2FoldChange >=1)
write.csv(up_microRNA,paste(DESeq2_result_prefix,"_up_microRNA.csv",sep=""),row.names = F)
down_microRNA<-subset(diff_microRNA,diff_microRNA$log2FoldChange <= -1)
write.csv(down_microRNA,paste(DESeq2_result_prefix,"_down_microRNA.csv",sep=""),row.names = F)


########################################
miR_Family_Info <- read_csv("/Volumes/Zunpeng/06_Database/07_TargetScanhuman/miR_Family_Info.csv")

resdata_anno<-merge(resdata,miR_Family_Info,by="MiRBase_ID")
write.csv(resdata_anno, paste(DESeq2_result_prefix,"_all_microRNA_anno.csv",sep=""),row.names = F)
diff_microRNA_anno<-merge(miR_Family_Info,diff_microRNA,by="MiRBase_ID")
write.csv(diff_microRNA_anno, paste(DESeq2_result_prefix,"_diff_microRNA_anno.csv",sep=""),row.names = F)
up_anno<-merge(miR_Family_Info,up_microRNA,by="MiRBase_ID")
write.csv(up_anno,paste(DESeq2_result_prefix,"_up_microRNA_anno.csv",sep=""),row.names = F)
down_anno<-merge(miR_Family_Info,down_microRNA,by="MiRBase_ID")
write.csv(down_anno,paste(DESeq2_result_prefix,"_down_microRNA_anno.csv",sep=""),row.names = F)

########################################
rld <- rlog(dds2)

Plot_prefix<-"MT_DGCR8_vs_LUC"

pdf(paste(Plot_prefix,"_log_scatter_diagram_WT.pdf",sep=""),width=6, height=6)
plot( log2( 1+counts(dds2, normalized=TRUE)[, 1:2] ), col="black", pch=20, cex=0.3 )
dev.off()
pdf(paste(Plot_prefix,"_log_scatter_diagram_KO.pdf",sep=""),width=6, height=6)
plot( log2( 1+counts(dds2, normalized=TRUE)[, 3:4] ), col="black", pch=20, cex=0.3 )
dev.off()

pdf(paste(Plot_prefix,"_rlog_scatter_diagram_WT.pdf",sep=""),width=6, height=6)
plot( assay(rld)[, 1:2], col="black", pch=20, cex=0.3 )
dev.off()
pdf(paste(Plot_prefix,"_rlog_scatter_diagram_KO.pdf",sep=""),width=6, height=6)
plot( assay(rld)[, 3:4], col="black", pch=20, cex=0.3 )
dev.off()

########################################

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, seq(2),sep="-")
colnames(sampleDistMatrix) <- NULL
sampleDist_colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

pdf(paste(Plot_prefix,"_Euclidean_distances_heatmap.pdf",sep=""),width=8, height=8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=sampleDist_colors)
dev.off()


########################################
pca_data <- plotPCA(rld, intgroup=c("condition"), returnData=TRUE)
write.csv(paste(Plot_prefix,"_PCA_data.csv",sep=""),row.names = F)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pdf(paste(Plot_prefix,"_PCA_plot.pdf",sep=""),width=6,height=6)
ggplot(pca_data, aes(PC1, PC2, color=condition, shape=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +ylim(-8,8)+
  theme_few()
dev.off()

########################################

pdf(paste(Plot_prefix,"_Expression_value_boxplot.pdf",sep=""),width=8, height=8)
expression_value_data<-assay(rld)
expression_value_color<-rainbow(7)
dev.off()
########################################
pdf(paste(Plot_prefix,"_MA_Plot.pdf",sep=""),width=8, height=8)
plotMA(res, main=paste(Plot_prefix,"_MA_Plot.",sep=""), ylim=c(-5,5))
dev.off()


########################################
pdf(paste(Plot_prefix,"_Histogram_p-value.pdf",sep=""),width=8, height=8)
hist(resdata$pvalue, breaks=20,col="blue",xlab="p-value of the result data",main="Histogram of p-value")
dev.off()


########################################
Volcano_data <-resdata

Volcano_data$threshold <- as.factor(ifelse(Volcano_data$padj < 0.05 & 
                                             abs(Volcano_data$log2FoldChange) >=1,
                                           ifelse(Volcano_data$log2FoldChange >= 1 ,'Up regulated','Down regulated'),'Unchanged'))

Volcano_data$lg10<--log10(Volcano_data$padj)


Volcano_tile <- paste0(paste(Plot_prefix,"_volcano plot",sep=""),
                       '\nCutoff for padj is 0.05',
                       '\nThe number of up regulated gene is ',nrow(Volcano_data[Volcano_data$threshold =='Up regulated',]) ,
                       '\nThe number of down regulated gene is ',nrow(Volcano_data[Volcano_data$threshold =='Down regulated',]),
                       '\nThe number of unchanged gene is ',nrow(Volcano_data[Volcano_data$threshold =='Unchanged',]) )

volcano_sum<-summary(Volcano_data$threshold)
write.csv(volcano_sum,paste(DESeq2_result_prefix,"_sunmmary_result.csv",sep=""))

pdf(paste(Plot_prefix,"_volcano plot_8.pdf",sep=""), width = 6, height = 7)
ggplot( data=Volcano_data,
        aes(x=log2FoldChange, y =-log10(padj),colour=threshold,fill=threshold)) +
  scale_color_manual(values=c("green", "blue","red"))+
  geom_point(alpha=1, size=1) +
  xlim(c(-4,4)) +ylim(c(0,5))+
  ggtitle(Volcano_tile) +
  theme_few() +
  geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
  theme(plot.title = element_text(hjust = 0.5,size=12),legend.title = element_blank())+
  labs(x="log2 (fold change)",y="-log10 (padj)")
dev.off()



Volcano_data$lg10[Volcano_data$lg10>=5]<-5
Volcano_data$log2FoldChange[Volcano_data$log2FoldChange<=-3]<--3
Volcano_data$log2FoldChange[Volcano_data$log2FoldChange>=3]<-3

volcano_sum<-summary(Volcano_data$threshold)
write.csv(volcano_sum,paste(DESeq2_result_prefix,"_sunmmary_result.csv",sep=""))

pdf(paste(Plot_prefix,"_volcano plot_10.pdf",sep=""), width = 6, height = 7)
ggplot( data=Volcano_data,
        aes(x=log2FoldChange, y =lg10,colour=threshold,fill=threshold)) +
  scale_color_manual(values=c("green", "blue","red"))+
  geom_point(alpha=1, size=1.2) +
  xlim(c(-3,3)) +ylim(c(0,5))+
  ggtitle(Volcano_tile) +
  theme_few() +
  geom_vline(xintercept=c(-1,1),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
  theme(plot.title = element_text(hjust = 0.5,size=12),legend.title = element_blank())+
  labs(x="log2 (fold change)",y="-log10 (padj)")
dev.off()






########################################
NAME<-Plot_prefix
heatmap_data<-diff_microRNA
row.names(heatmap_data)<-heatmap_data$MiRBase_ID
heatmap_data<-heatmap_data[,c(8,9,10,11,12,13)]
heatmap_data<-data.matrix(heatmap_data)



pheatmap(heatmap_data,
         main = paste(NAME,"_heatmap",sep=""),
         cluster_rows = TRUE,
         cluster_cols =FALSE,
         scale = "row",#"none", "row", "column"
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = FALSE,
         show_colnames = TRUE,
         filename = paste(Plot_prefix,"_heatmap.pdf",sep=""),
         width =6,         
         height = 7)

pheatmap(heatmap_data,
         main = paste(NAME,"_heatmap",sep=""),
         cluster_rows = TRUE,
         cluster_cols =FALSE,
         clustering_distance_rows = "correlation",#euclidean
         clustering_distance_cols = "correlation",#euclidean
         clustering_method = "complete",
         scale = "row",#"none", "row", "column"
         #color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         #c("skyblue","yellow","tomato")
         #color=colorRampPalette(c('green','yellow','red'), bias=1)(50),
         border_color = NA, #"grey60",
         fontsize = 14,
         fontsize_row = 10, #fontsize,
         fontsize_col = 10, #fontsize,
         fontsize_number = 0.8* fontsize,
         #margins=c(5,10)
         #clustering_callback = identity2,
         kmeans_k = NA,
         breaks = NA,
         cutree_rows = NA,
         cutree_cols = NA,
         #treeheight_row = ifelse(cluster_rows, 50, 0),
         #treeheight_col = ifelse(cluster_cols, 50, 0),
         show_rownames = FALSE,
         show_colnames = TRUE,
         legend = TRUE,
         legend_breaks = NA,
         legend_labels = NA,
         annotation = NA,
         annotation_col= NA,
         #annotation_row = annotation_row,#NA,
         #annotation_colors = ann_colors,
         annotation_legend = TRUE,
         drop_levels = TRUE,
         #display_numbers = F,
         #display_numbers = TRUE, number_format = "%.2f", number_color="purple")
         #number_format = "%.2f",
         #number_color = "grey30",
         gaps_row = NULL,
         gaps_col = NULL,
         labels_row = NULL,
         labels_col = NULL,
         cellwidth = NA,
         cellheight = NA,
         #revC=TRUE,
         filename = paste(Plot_prefix,"_heatmap.pdf",sep=""),
         width =6,         
         height = 7)


#########################
wordcloud2
pdf("DGCR8_KO_vs_WT_miRNA_up_wordcloud", width = 8, height = 6)
wordcloud_up<-data.frame(word=up_microRNA$MiRBase_ID,
                         freq=round(up_microRNA$log2FoldChange,0))
wordcloud2(wordcloud_up,size = 0.3,shape = "circle")
dev.off()

DGCR8_KO_vs_WT_miRNA_down_wordcloud
wordcloud_down<-data.frame(word=down_microRNA$MiRBase_ID,
                           freq=round(down_microRNA$log2FoldChange,0))
wordcloud2(wordcloud_down,size = 0.1,shape = "circle")



################################################################
library(readr)
hsa_MTI <- read_csv("~/03_database/06_miRTarBase/hsa_MTI.csv")

up_MTI<-subset(hsa_MTI,hsa_MTI$miRNA%in%up_microRNA$MiRBase_ID)
write.csv(up_MTI,"up_regulated_miRNA_hsa_MTI.csv")

down_MTI<-subset(hsa_MTI,hsa_MTI$miRNA%in%down_microRNA$MiRBase_ID)
write.csv(down_MTI,"down_regulated_miRNA_hsa_MTI.csv")

################################################################
Predicted_Targets_Info_default_predictions <- read_csv("~/03_database/07_TargetScanhuman/Predicted_Targets_Info.default_predictions.csv")

up_Predicted_Targets<-subset(Predicted_Targets_Info_default_predictions,
                             Predicted_Targets_Info_default_predictions$`miR Family`%in%up_anno$`miR family`)
write.csv(up_Predicted_Targets,"up_regulated_miRNA_TargetScan.csv",row.names = F)

down_Predicted_Targets<-subset(Predicted_Targets_Info_default_predictions,
                               Predicted_Targets_Info_default_predictions$`miR Family`%in%down_anno$`miR family`)
write.csv(down_Predicted_Targets,"down_regulated_miRNA_TargetScan.csv",row.names = F)


source("https://bioconductor.org/biocLite.R")
biocLite("topGO")
################## Heatmap #######################
res<-DGCR8_KO_vs_WT_miRNA_all_microRNA
res$sum<-res$DN1+res$DN2+res$DT1+res$DT2
res2<-data.matrix(res[order(res$sum,decreasing = TRUE),])

heatmap_data<-res[,c(8,9,10,11)]

heatmap_data2<-log2(heatmap_data+1)
pheatmap(heatmap_data2,
         main = "all_microRNA__heatmap",
         cluster_rows = FALSE,
         cluster_cols =FALSE,
         clustering_distance_rows = "correlation",#euclidean
         clustering_distance_cols = "correlation",#euclidean
         clustering_method = "complete",
         scale = "row",#"none", "row", "column"
         #color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         #c("skyblue","yellow","tomato")
         #color=colorRampPalette(c('green','yellow','red'), bias=1)(50),
         border_color = NA, #"grey60",
         fontsize = 14,
         fontsize_row = 10, #fontsize,
         fontsize_col = 10, #fontsize,
         fontsize_number = 0.8* fontsize,
         #margins=c(5,10)
         #clustering_callback = identity2,
         kmeans_k = NA,
         breaks = NA,
         cutree_rows = NA,
         cutree_cols = NA,
         #treeheight_row = ifelse(cluster_rows, 50, 0),
         #treeheight_col = ifelse(cluster_cols, 50, 0),
         show_rownames = FALSE,
         show_colnames = TRUE,
         legend = TRUE,
         legend_breaks = NA,
         legend_labels = NA,
         annotation = NA,
         annotation_col= NA,
         #annotation_row = annotation_row,#NA,
         #annotation_colors = ann_colors,
         annotation_legend = TRUE,
         drop_levels = TRUE,
         #display_numbers = F,
         #display_numbers = TRUE, number_format = "%.2f", number_color="purple")
         #number_format = "%.2f",
         #number_color = "grey30",
         gaps_row = NULL,
         gaps_col = NULL,
         labels_row = NULL,
         labels_col = NULL,
         cellwidth = NA,
         cellheight = NA,
         #revC=TRUE,
         filename = "all_microRNA_log2_deseq2_scale_row_heatmap_1227.pdf",
         width =5,         
         height = 7)


library(pheatmap)
