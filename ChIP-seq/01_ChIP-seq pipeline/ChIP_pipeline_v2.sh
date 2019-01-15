#!usr/bin/sh
if [ $# -ne 1 ]
then
echo "sh `basename $0`  <IN_dir>  <sample name>   example sh $0 /public/agis/gaofei_group/liuzunpeng/02_Data/02_ChIP-seq/H3K9me3 KO_hMSC_H3K9me3_rep1"
exit
fi

indir=/pfs1/liuguanghui/liuzunpeng/02_Data/02_ChIP-seq/DGCR8_KAT7_H3K9me3/raw_data/DGCR8_H3K9me3
sample=$1
f1=$indir/$sample"_"1.fq.gz
f2=$indir/$sample"_"2.fq.gz

out=/pfs1/liuguanghui/liuzunpeng/05_Results/15_DGCR8_ChIP_seq
mkdir -p $out

qsu=/pfs1/liuguanghui/liuzunpeng/06_Qsub/15_DGCR8_ChIP_seq
mkdir -p $qsu

log_dir=/home/liuzunpeng/../../pfs1/liuguanghui/liuzunpeng/06_Qsub/15_DGCR8_ChIP_seq

################ Reference ################
ref=/pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/hg19.fa
gtf=/pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/Homo_sapiens.GRCh37.87_chr.chr.gtf
bowtie2_index=/pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/02_bowtie2_index/hg19

################    Src    ################
ATAC_seq_src=/pfs1/liuguanghui/liuzunpeng/04_Softwares/ATAC_seq_src/
macs_stat=$ATAC_seq_src/macs_stat.pl

################ Softwares ################
fastqc=/pfs1/liuguanghui/liuzunpeng/04_Softwares/fastqc/FastQC/fastqc
fqvalue=/pfs1/liuguanghui/liuzunpeng/04_Softwares/fqvalue/fqvalue_v2.4/fqvalue
samtools=/pfs1/liuguanghui/liuzunpeng/04_Softwares/samtools/samtools-1.6/samtools
bowtie2=/pfs1/liuguanghui/liuzunpeng/04_Softwares/bowtie/bowtie2-2.2.9/bowtie2
bedtools=/pfs1/liuguanghui/liuzunpeng/04_Softwares/bedtools2/bin
macs2=/pfs1/liuguanghui/liuzunpeng/04_Softwares/macs2/MACS2-2.1.1.20160309/bin/macs2

################ fastQCraw ################
mkdir -p $out/$sample/01_fastQC_raw
mkdir -p $qsu/01_fastQC_raw

fqv1=$out/$sample/01_fastQC_raw/$sample.R1.basic.xls
fqv2=$out/$sample/01_fastQC_raw/$sample.R2.basic.xls

echo -ne "######
#PBS -S /bin/bash
#PBS -N FastQC.$sample
#PBS -q XD64
#PBS -k eo
#PBS -o $log_dir/01_fastQC_raw/fastqc.$sample.log
#PBS -e $log_dir/01_fastQC_raw/fastqc.$sample.$sample.err
#PBS -l nodes=2:ppn=10

$fastqc --outdir $out/$sample/01_fastQC_raw -f fastq $f1 $f2 &&
$fqvalue -q 33 $f1 > $fqv1 &&
$fqvalue -q 33 $f2 > $fqv2 &&
echo \"Raw data FastQC&fqvalue of $sample has been done! \"
"> $qsu/01_fastQC_raw/fastqc.$sample.sh


################ TrimGalore clean ################
mkdir -p $out/$sample/02_TrimGalore
mkdir -p $qsu/02_TrimGalore

trim_galore=/pfs1/liuguanghui/liuzunpeng/04_Softwares/trim_galore/TrimGalore-0.4.5/trim_galore
trim_galore_log=$out/$sample/02_TrimGalore/$sample.trim_galore.log

echo -ne "######
#PBS -S /bin/bash
#PBS -N trim_galore.$sample
#PBS -q XD64
#PBS -k eo
#PBS -o $log_dir/02_TrimGalore/trim_galore.$sample.log
#PBS -e $log_dir/02_TrimGalore/trim_galore.$sample.err
#PBS -l nodes=2:ppn=10

$trim_galore --fastqc --path_to_cutadapt /home/liuzunpeng/.local/bin/cutadapt --stringency 3 --paired --output_dir $out/$sample/02_TrimGalore $f1 $f2 2>trim_galore_log &&
echo \"Clean data and FastQC done! \n All done! \"
" > $qsu/02_TrimGalore/trim_galore.$sample.sh

################ Alignment-to-reference-by-bowtie2 ################
cf1=$out/$sample/02_TrimGalore/$sample"_1_val_1.fq.gz"
cf2=$out/$sample/02_TrimGalore/$sample"_2_val_2.fq.gz"

mkdir -p $out/$sample/03_bowtie2
mkdir -p $qsu/03_bowtie2
bowtie2_Align_log=$out/$sample/03_bowtie2/$sample.bowtie2_Align.log
bowtie2_Align_result=$out/$sample/03_bowtie2/$sample.bowtie2_Align_result.xls

echo -ne "######
#PBS -S /bin/bash
#PBS -N Bowtie2.$sample
#PBS -q XD64
#PBS -k eo
#PBS -o $log_dir/03_bowtie2/bowtie2.$sample.log
#PBS -e $log_dir/03_bowtie2/bowtie2.$sample.err
#PBS -l nodes=2:ppn=10

$bowtie2 -p 10 -x $bowtie2_index -1 $cf1 -2 $cf2 -S $out/$sample/03_bowtie2/$sample.sam 2>$bowtie2_Align_log &&
echo \"$sample align done! \" &&
echo \"Transfer $sample sam to bam! \" &&
$samtools view -S -b $out/$sample/03_bowtie2/$sample.sam >$out/$sample/03_bowtie2/$sample.bam &&
echo \"Transfer $sample sam to bam done! \" &&
echo \"Sort $sample bam file! \" &&
$samtools sort $out/$sample/03_bowtie2/$sample.bam -o $out/$sample/03_bowtie2/$sample.sort.bam &&
echo \"Sort $sample bam file done! \" &&
echo \"Samtools index $sample sort.bam! \" &&
$samtools index $out/$sample/03_bowtie2/$sample.sort.bam &&
rm $out/$sample/03_bowtie2/$sample.sam &&
echo \"$sample sam files removed. \nall done! \n\"
" >$qsu/03_bowtie2/bowtie2.$sample.sh

################ Alignment-to-reference-by-bowtie2 ################

mkdir -p $qsu/03_rmdup
rmdup_log=$out/$sample/03_bowtie2/$sample.Remove_duplicates.log

echo -ne "######
#PBS -S /bin/bash
#PBS -N Remove_duplicates.$sample
#PBS -q XD64
#PBS -k eo
#PBS -o $log_dir/03_bowtie2/Remove_duplicates.$sample.log
#PBS -e $log_dir/03_bowtie2/Remove_duplicates.$sample.err
#PBS -l nodes=2:ppn=10

echo \"Samtools rmdup $sample sort.bam! \" &&
$samtools rmdup $out/$sample/03_bowtie2/$sample.sort.bam $out/$sample/03_bowtie2/$sample.sort.rmdup.bam 2>$rmdup_log &&
echo \"$sample rmdup done !\" &&
echo \"Samtools index $sample sort.rumdup.bam! \" &&
$samtools index $out/$sample/03_bowtie2/$sample.sort.rmdup.bam &&
echo \"Samtools index $sample sort.rumdup.bam done! \"
" >$qsu/03_rmdup/rmdup.$sample.sh


################ Calling-peak-by-MACS2 ################
mkdir -p $out/$sample/04_MACS2
mkdir -p $qsu/04_MACS2

MACS2_log=$out/$sample/04_MACS2/$sample.MACS2.log

Control_bowtie2_bam=/pfs1/liuguanghui/liuzunpeng/05_Results/15_DGCR8_ChIP_seq/Input/03_bowtie2/Input.sort.rmdup.bam

peak_stat=$out/$sample/04_MACS2/$sample.bowtie2_state.txt

echo -ne "######
#PBS -S /bin/bash
#PBS -N MACS2.$sample
#PBS -q XD64
#PBS -k eo
#PBS -o $log_dir/04_MACS2/MACS2.$sample.log
#PBS -e $log_dir/04_MACS2/MACS2.$sample.err
#PBS -l nodes=2:ppn=10

$macs2 callpeak -t $out/$sample/03_bowtie2/$sample.sort.rmdup.bam -c $Control_bowtie2_bam -f BAM -g hs -B -q 0.01 -n $sample --outdir $out/$sample/04_MACS2 2>$MACS2_log &&
echo \"$sample calling-peak-by-MACS2 done! \" &&
grep '^chr\S' $out/$sample/04_MACS2/$sample"_"peaks.xls |perl -alne '{print \"\$F[0]\\\t\$F[1]\\\t\$F[2]\\\t\$F[9]\\\t\$F[7]\\\t+\"}' >$out/$sample/04_MACS2/$sample.macs2.peaks.bed &&
echo \"Get $sample peak_bed done! \" &&
perl $macs_stat $out/$sample/04_MACS2/$sample\_summits.bed >$peak_stat &&
echo \"$sample peak stat done! \nAll done! \"
">$qsu/04_MACS2/MACS2.$sample.sh


################ Bamcovrage ################
mkdir -p $out/$sample/05_Bamcovrage
mkdir -p $qsu/05_Bamcovrage

Deeptools=/pfs1/liuguanghui/liuzunpeng/04_Softwares/deepTools/deepTools/local2/bin
bamCoverage=$Deeptools/bamCoverage
plotCoverage=$Deeptools/plotCoverage

computeMatrix=$Deeptools/computeMatrix
plotProfile=$Deeptools/plotProfile
plotHeatmap=$Deeptools/plotHeatmap
plotEnrichment=$Deeptools/plotEnrichment

Bamcovrage_log=$out/$sample/05_Bamcovrage/$sample.Bamcovrage.log

echo -ne "######
#PBS -S /bin/bash
#PBS -N Bamcovrage.$sample
#PBS -q XD64
#PBS -k eo
#PBS -o $log_dir/05_Bamcovrage/Bamcovrage.$sample.log
#PBS -e $log_dir/05_Bamcovrage/Bamcovrage.$sample.err
#PBS -l nodes=2:ppn=10
echo \"Start deeptools bamCoverage of $sample binSize 2000 bigwig! \" &&
$bamCoverage -b $out/$sample/03_bowtie2/$sample.sort.rmdup.bam \
--ignoreDuplicates --normalizeUsingRPKM --binSize 2000 \
-o $out/$sample/05_Bamcovrage/$sample.coverage_bin2000.bw  2>$Bamcovrage_log &&
echo \"Start deeptools bamCoverage of $sample binSize 2000 bedgraph! \" &&
$bamCoverage -b $out/$sample/03_bowtie2/$sample.sort.rmdup.bam \
--outFileFormat bedgraph --ignoreDuplicates --normalizeUsingRPKM --binSize 2000 \
-o $out/$sample/05_Bamcovrage/$sample.coverage_bin2000.bedgraph  2>>$Bamcovrage_log &&
echo \"Finish deeptools bamCoverage of $sample binSize 2000! \"
">$qsu/05_Bamcovrage/Bamcovrage.$sample.sh

echo -ne "######
#PBS -S /bin/bash
#PBS -N Bamcovrage.$sample
#PBS -q XD64
#PBS -k eo
#PBS -o $log_dir/05_Bamcovrage/Bamcovrage.$sample.log
#PBS -e $log_dir/05_Bamcovrage/Bamcovrage.$sample.err
#PBS -l nodes=2:ppn=10
echo \"Start deeptools bamCoverage of $sample binSize 100 bigwig! \" &&
$bamCoverage -b $out/$sample/03_bowtie2/$sample.sort.rmdup.bam \
--ignoreDuplicates --normalizeUsingRPKM --binSize 100 \
-o $out/$sample/05_Bamcovrage/$sample.coverage_bin100.bw  2>$Bamcovrage_log &&
echo \"Start deeptools bamCoverage of $sample binSize 100 bedgraph! \" &&
$bamCoverage -b $out/$sample/03_bowtie2/$sample.sort.rmdup.bam \
--outFileFormat bedgraph --ignoreDuplicates --normalizeUsingRPKM --binSize 100 \
-o $out/$sample/05_Bamcovrage/$sample.coverage_bin100.bedgraph  2>>$Bamcovrage_log &&
echo \"Finish deeptools bamCoverage of $sample binSize 100! \"
">$qsu/05_Bamcovrage/Bamcovrage.bin100.$sample.sh

