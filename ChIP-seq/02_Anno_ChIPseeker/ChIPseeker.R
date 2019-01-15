setwd("/Volumes/Zunpeng/08_Result/04_DGCR8/ChIP_seq_H3K9me3/02_Anno_ChIPseeker")

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
chr<-c(paste("chr",seq(22),sep = "") ,"chrX","chrY")

KO_bedPeaksFile = file.path("/Volumes/Zunpeng/08_Result/04_DGCR8/ChIP_seq_H3K9me3/02_Anno_ChIPseeker/DGCR8_KO_K9.macs2.peaks.bed")
KO_peak<-readPeakFile(KO_bedPeaksFile)

H9_bedPeaksFile = file.path("/Volumes/Zunpeng/08_Result/04_DGCR8/ChIP_seq_H3K9me3/02_Anno_ChIPseeker/DGCR8_WT_K9.macs2.peaks.bed")
H9_peak<-readPeakFile(H9_bedPeaksFile)

peaks <- list(H9=H9_peak,KO=KO_peak)
tagMatrixList <- lapply(peaks, getTagMatrix, windows=promoter)
peakAnnoList <- lapply(peaks, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)
H9_anno<-peakAnnoList$H9
write.csv(H9_anno,"H9.anno.csv",row.names = FALSE)
KO_anno<-peakAnnoList$KO
write.csv(KO_anno,"KO.anno.csv",row.names = FALSE)

pdf('p1_H3K9me3 ChIP peaks over Chromosomes .pdf',width=20,height=12)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX","chrY"), weightCol="V5") + facet_grid(chr ~ .id)
dev.off()

pdf('p1_H3K9me3 ChIP peaks over Chromosomes 2.pdf',width=14,height=12)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX","chrY"), weightCol="V5") + facet_grid(chr ~ .id)
dev.off()

pdf('p1_H3K9me3 ChIP peaks over Chromosomes 3.pdf',width=12,height=10)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX","chrY"), weightCol="V5") + facet_grid(chr ~ .id)+scale_color_manual(values=c("red","blue"))
dev.off()

pdf('p1_H3K9me3 ChIP peaks over Chromosomes 5.pdf',width=11,height=10)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX","chrY"), weightCol="V5") + facet_grid(chr ~ .id)+scale_color_manual(values=c("red","blue"))
dev.off()

pdf('p1_H3K9me3 ChIP peaks over Chromosomes 6.pdf',width=11,height=10)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX","chrY"), weightCol="V5") + facet_grid(chr ~ .id)+scale_color_manual(values=c("blue","red"))
dev.off()


pdf('p1_H3K9me3 ChIP peaks over Chromosomes 4.pdf',width=11,height=10)
covplot(peaks,chrs = c(paste("chr",seq(22),sep = "") ,"chrX","chrY"), weightCol="V5") + facet_grid(chr ~ .id)+scale_color_manual(values=c("grey","green3"))
dev.off()


pdf('p2_TagHeatmap of H3K9me3 ChIP peaks binding to TSS region_all.pdf',width=12,height=12)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL,xlab="Enrichment of H3K9me3 at TSS region")
dev.off()

png('p2_TagHeatmap of H3K9me3 ChIP peaks binding to TSS region_all.tiff',width=600,height=1200)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL,xlab="Enrichment of H3K9me3 at TSS region")
dev.off()

png('p2_TagHeatmap of H3K9me3 ChIP peaks binding to TSS region_all_2.tiff',width=800,height=1200)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL,xlab="Enrichment of H3K9me3 at TSS region")
dev.off()


png('p2_TagHeatmap of H3K9me3 ChIP peaks binding to TSS region_all_3.tiff',width=1000,height=1200)
tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL,xlab="Enrichment of H3K9me3 at TSS region")
dev.off()


pdf('p3_Average Profile of H3K9me3 ChIP peaks binding to TSS region.pdf',width=6,height=6)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")
dev.off()

pdf('p3_Average Profile of H3K9me3 ChIP peaks binding to TSS region 2.pdf',width=6,height=6)
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


