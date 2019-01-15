library(GenomicFeatures)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ReactomePA)
library(ChIPseeker)
library(ggplot2)
library(ggthemes)

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
chr<-c(paste("chr",seq(22),sep = "") ,"chrX")

all_bedPeaksFile = file.path("/Volumes/Zunpeng/08_Result/04_DGCR8/01_ChIP_seq_H3K9me3/07_mountain/WT_merge_d2000_up10000_merge_100000.final_2.bed")
all_peak<-readPeakFile(all_bedPeaksFile)

loss_bedPeaksFile = file.path("/Volumes/Zunpeng/08_Result/04_DGCR8/01_ChIP_seq_H3K9me3/07_mountain/WT_merge_only_merge_d2000_up10000_merge_100000.final_2.bed")
loss_peak<-readPeakFile(loss_bedPeaksFile)

genome_bedPeaksFile = file.path("/Volumes/Zunpeng/08_Result/04_DGCR8/01_ChIP_seq_H3K9me3/07_mountain/")
genome_peak<-readPeakFile(genome_bedPeaksFile)

peaks <- list(b_mountain=all_peak,c_Loss_mountain=loss_peak,a_genome=genome_peak)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

H9_anno<-peakAnnoList$H9
write.csv(H9_anno,"H9.anno.csv",row.names = FALSE)
KO_anno<-peakAnnoList$KO
write.csv(KO_anno,"KO.anno.csv",row.names = FALSE)

pdf('p1_H3K9me3 ChIP peaks over Chromosomes.pdf',width=20, height=11)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX"), weightCol="V5") + facet_grid(chr ~ .id)
dev.off()

pdf('p1_H3K9me3 ChIP peaks over Chromosomes 2.pdf',width=20,height=11)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX"), weightCol="V5") + facet_grid(chr ~ .id)+
  scale_color_manual(values=c("lightgrey","mediumturquoise","red"))+  scale_fill_manual(values=c("lightgrey","mediumturquoise","red"))
dev.off()


pdf('p1_H3K9me3 ChIP peaks over Chromosomes 3.pdf',width=20,height=11)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX"), weightCol="V5") + facet_grid(chr ~ .id)+
  scale_color_manual(values=c("lightgrey","mediumturquoise","red"))+  scale_fill_manual(values=c("lightgrey","mediumturquoise","red"))
dev.off()


pdf('p1_H3K9me3 ChIP peaks over Chromosomes 4.pdf',width=20,height=11)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX"), weightCol="V5") +  facet_grid(chr ~ .id)+
  scale_color_manual(values=c("lightgrey","mediumturquoise","lightcoral"))+ 
  scale_fill_manual(values=c("lightgrey","mediumturquoise","lightcoral"))
dev.off()


pdf('p1_H3K9me3 ChIP peaks over Chromosomes 5.pdf',width=8,height=11)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX"), weightCol="V5") +
  scale_fill_manual(values=c("gainsboro","mediumturquoise","red"))
dev.off()



pdf('p1_H3K9me3 ChIP peaks over Chromosomes 5.pdf',width=11,height=11)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX"), weightCol="V5") + 
  scale_color_manual(values=c("mediumturquoise","lightcoral"))+  scale_fill_manual(values=c("mediumturquoise","lightcoral"))
dev.off()

pdf('p1_H3K9me3 ChIP peaks over Chromosomes 6.pdf',width=8,height=11)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX"), weightCol="V5") + 
  scale_color_manual(values=c("grey","lightcoral"))+  scale_fill_manual(values=c("grey","lightcoral"))
dev.off()



pdf('p1_H3K9me3 ChIP peaks over Chromosomes 7.pdf',width=8,height=11)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX"), weightCol="V5") + 
  scale_color_manual(values=c("grey","red"))+  scale_fill_manual(values=c("grey","red"))
dev.off()

pdf('p1_H3K9me3 ChIP peaks over Chromosomes 8.pdf',width=12,height=11)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX"), weightCol="V5") + 
  facet_grid(chr ~ .id)
dev.off()

pdf('p1_H3K9me3 ChIP peaks over Chromosomes 9.pdf',width=8,height=11)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX"), weightCol="V5")
dev.off()


pdf('p1_H3K9me3 ChIP peaks over Chromosomes 10.pdf',width=8,height=11)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX"), weightCol="V5")+scale_color_manual(values=c("grey","green3"))
dev.off()

pdf('p1_H3K9me3 ChIP peaks over Chromosomes 12.pdf',width=7,height=11)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX"), weightCol="V5")+scale_color_manual(values=c("grey","lightseagreen"))
dev.off()

pdf('p1_H3K9me3 ChIP peaks over Chromosomes 13.pdf',width=7,height=11)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX"), weightCol="V5")+
  scale_color_manual(values=c("grey","cyan3"))
dev.off()


















pdf('p2_TagHeatmap of H3K9me3 ChIP peaks binding to TSS region_all.pdf',width=8,height=12)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=c("lightcoral","grey"),xlab="Enrichment of H3K9me3 at TSS region")
dev.off()

png('p2_TagHeatmap of H3K9me3 ChIP peaks binding to TSS region_all.tiff',width=600,height=1200)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=c("lightcoral","mediumturquoise"),xlab="Enrichment of H3K9me3 at TSS region")
dev.off()

png('p2_TagHeatmap of H3K9me3 ChIP peaks binding to TSS region_all_2.tiff',width=800,height=1200)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=c("lightcoral","mediumturquoise"),xlab="Enrichment of H3K9me3 at TSS region")
dev.off()


png('p2_TagHeatmap of H3K9me3 ChIP peaks binding to TSS region_all_3.tiff',width=700,height=1200)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=c("lightcoral","mediumturquoise"),xlab="Enrichment of H3K9me3 at TSS region")
dev.off()


pdf('p3_Average Profile of H3K9me3 ChIP peaks binding to TSS region.pdf',width=10,height=6)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

pdf('p3_Average Profile of H3K9me3 ChIP peaks binding to TSS region 2.pdf',width=10,height=6)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), resample=500, facet="row",xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

pdf('p3_Average Profile of H3K9me3 ChIP peaks binding to TSS region_all_3.pdf',width=10,height=10)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row",xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

pdf('p3_Average Profile of H3K9me3 ChIP peaks binding to TSS region_all_4.pdf',width=10,height=10)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.99,resample=500, facet="row",xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

pdf('p5_Bar-summarize the distribution of peaks over different type of features_all.pdf',width=10,height=4)
plotAnnoBar(peakAnnoList)
dev.off()

pdf('p9_Distribution of transcription factor-binding loci relative to TSS_all.pdf',width=10,height=4)
plotDistToTSS(peakAnnoList, title="Distribution of transcription factor-binding loci relative to TSS")
dev.off()


