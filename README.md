# NCOMMS-18-15055-T

Welcome to the GitHub repository associated with our recent publication：

NCOMMS-18-15055-T （in revision）

## Abstract

DiGeorge syndrome critical region 8 (DGCR8) is a critical component of the canonical microprocessor complex for microRNA biogenesis. However, the non-canonical functions of DGCR8 have not been studied. Here, we demonstrate that DGCR8 plays an important role in maintaining heterochromatin organization and preventing aging. An N-terminal-truncated version of DGCR8 (DR8dex2) accelerated senescence in human mesenchymal stem cells (hMSCs) independent of its miRNA-processing activity, which is mediated by its C-terminal domains. Further studies revealed that DGCR8 maintained heterochromatin organization by interacting with the nuclear envelope protein Lamin B1, and heterochromatin-associated proteins, KAP1 and HP1r Overexpression of any of these proteins, including DGCR8, reversed premature senescent phenotypes in DR8dex2 hMSCs. Finally, DGCR8 was downregulated in pathologically and naturally aged hMSCs, whereas DGCR8 overexpression alleviated hMSC aging and osteoarthritis in mice. Taken together, these analyses uncovered a novel, miRNA-independent role for DGCR8 in maintaining heterochromatin organization and preventing senescence. DGCR8 may therefore represent a new therapeutic target for alleviating human aging-related disorders.


## What you'll find here
The file folders of each sequencing include pipelines (software/code) and some small (simulated or real) dataset to demo the software/code.

The Instruction for each pipeline has been included in the main fuction. 

1. The pipeline and NA-seq data processing. A pipeline used in which could processing from raw fastq reads to transcripts quantitation and differentially expressed genes (DEGs) analysis.
2. The pipeline of ChIP-seq data processing. A pipeline used in which could processing from raw fastq reads to peak calling analysis.
3. The pipeline of microRNA-seq data processing. A pipeline used in which could processing from raw fastq reads to peak calling analysis.


## Required softwares
* fastqc=/FastQC/fastqc
* fqvalue=/fqvalue_v2.4/fqvalue
* FASTX_fastq_quality_filter=/FASTX-Toolkit/fastq_quality_filter
* FASTX_clipper=/FASTX-Toolkit/fastx_clipper
* FASTX_trimmer=/FASTX-Toolkit/fastx_trimmer
* samtools=/samtools-1.6/samtools
* bedtools=/bedtools2/bin
* fastq2fasta=/mirdeep2-2.0.0.8-pl5.22.0_5/bin/fastq2fasta.pl
* miRDeep2=/mirdeep2-2.0.0.8-pl5.22.0_5/bin/miRDeep2.pl
* mapper=/mirdeep2-2.0.0.8-pl5.22.0_5/bin/mapper.pl
* quantifier=/mirdeep2-2.0.0.8-pl5.22.0_5/bin/quantifier.pl
* trim_galore=/TrimGalore-0.4.5/trim_galore
* tophat=/tophat-2.1.1.Linux_x86_64/tophat2
* hisat2=/hisat2-2.0.4/hisat2
* bowtie2=/bowtie2-2.2.9/bowtie2
* HTseq=/home/liuzunpeng/.local/bin/htseq-count
* stringTie=/stringtie-1.2.3.Linux_x86_64/stringtie
* bedtools=/bedtools2/bin
* bamToBed=/bedtools2/bin/bamToBed
* genomeCoverageBed=/bedtools2/bin/genomeCoverageBed
* bedGraphToBigWig=/bedGraphToBigWig/bedGraphToBigWig
* bamCoverage=/deepTools/local2/bin/bamCoverage
* macs2=/MACS2-2.1.1.20160309/bin/macs2



## Required R libraries
* library(DESeq2) 
* library(geneplotter)
* library(ggthemes)
* library(pheatmap)
* library(GO.db)
* library(org.Hs.eg.db) 
* library(topGO) 
* library(GSEABase)
* library(clusterProfiler) 
* library(biomaRt) 
* library(ggplot2)
* library(RColorBrewer)
* library(ggpubr)
* library(scatterplot3d)
* library(GenomicFeatures)
* library(GenomicRanges)
* library(TxDb.Hsapiens.UCSC.hg19.knownGene)
* library(ReactomePA)
* library(ChIPseeker)

## Questions/Issues

If you have any questions or issues, please feel free to contact with me zunpengAT163.com.


