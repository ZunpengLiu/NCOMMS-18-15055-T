#!usr/bin/sh
if [ $# -ne 1 ]
then
echo "sh `basename $0`  <sample name> "
exit
fi

sample=$1

# gunzip /pfs1/liuguanghui/liuzunpeng/02_Data/05_microRNA_seq/MF_Heart_Liver/new/Rawdata/$sample/$sample"_"R1.fq.gz

fastq=/pfs1/liuguanghui/liuzunpeng/05_Results/58_micro_RNA_DGCR8/raw_data/$sample".fq"

out=/pfs1/liuguanghui/liuzunpeng/05_Results/58_micro_RNA_DGCR8/results

qsu=/pfs1/liuguanghui/liuzunpeng/05_Results/58_micro_RNA_DGCR8/qsub

log_dir=/home/liuzunpeng/../../pfs1/liuguanghui/liuzunpeng/05_Results/58_micro_RNA_DGCR8/qsub


################ Softwares ################
fastqc=/pfs1/liuguanghui/liuzunpeng/04_Softwares/fastqc/FastQC/fastqc
fqvalue=/pfs1/liuguanghui/liuzunpeng/04_Softwares/fqvalue/fqvalue_v2.4/fqvalue

FASTX_fastq_quality_filter=/pfs1/liuguanghui/liuzunpeng/04_Softwares/FASTX-Toolkit/fastq_quality_filter
FASTX_clipper=/pfs1/liuguanghui/liuzunpeng/04_Softwares/FASTX-Toolkit/fastx_clipper
FASTX_trimmer=/pfs1/liuguanghui/liuzunpeng/04_Softwares/FASTX-Toolkit/fastx_trimmer

samtools=/pfs1/liuguanghui/liuzunpeng/04_Softwares/samtools/samtools-1.6/samtools
bedtools=/pfs1/liuguanghui/liuzunpeng/04_Softwares/bedtools2/bin

fastq2fasta=/home/liuzunpeng/anaconda3/pkgs/mirdeep2-2.0.0.8-pl5.22.0_5/bin/fastq2fasta.pl
miRDeep2=/home/liuzunpeng/anaconda3/pkgs/mirdeep2-2.0.0.8-pl5.22.0_5/bin/miRDeep2.pl
mapper=/home/liuzunpeng/anaconda3/pkgs/mirdeep2-2.0.0.8-pl5.22.0_5/bin/mapper.pl
quantifier=/home/liuzunpeng/anaconda3/pkgs/mirdeep2-2.0.0.8-pl5.22.0_5/bin/quantifier.pl

################ Reference ################
genome_index=/pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/04_bowtie1.index/hg19
genome_fa=/pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/hg19.fa
gtf=/pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/Homo_sapiens.GRCh37.87.final.gtf

hairpin_fa=/pfs1/liuguanghui/liuzunpeng/03_Database/03_miRBase/hairpin.fa
mature_fa=/pfs1/liuguanghui/liuzunpeng/03_Database/03_miRBase/mature.fa
hairpin_hsa_fa=/pfs1/liuguanghui/liuzunpeng/03_Database/03_miRBase/hairpin.hsa.fa
mature_hsa_fa=/pfs1/liuguanghui/liuzunpeng/03_Database/03_miRBase/mature.hsa.fa
#mature_other_fa=/public/agis/gaofei_group/liuzunpeng/03_database/03_miRBase/mature.other.fa
#mature_cel_fa=/public/agis/gaofei_group/liuzunpeng/03_database/03_miRBase/mature.cel.fa

########################## step1 fastQCraw ##########################
mkdir -p $out/$sample/01_fastQC_raw
mkdir -p $qsu/01_fastQC_raw
fqv_raw=$out/$sample/01_fastQC_raw/$sample.raw.basic.xls

echo -ne "######
#PBS -S /bin/bash
#PBS -N FastQC.$sample
#PBS -q XD64
#PBS -o $log_dir/01_fastQC_raw/fastqc.$sample.log
#PBS -e $log_dir/01_fastQC_raw/fastqc.$sample.$sample.err
#PBS -l nodes=2:ppn=10

echo \"Start FastQC & fqvalue raw data  $sample ! \" &&
$fastqc --outdir $out/$sample/01_fastQC_raw $fastq &&
$fqvalue -q 33 $fastq > $fqv_raw &&
echo \"Raw data FastQC & fqvalue of $sample has been done! \"
######
\n"> $qsu/01_fastQC_raw/fastqc.$sample.sh

##########################           step2 Cut adapter & Clean              ##########################

mkdir -p $out/$sample/02_clean
mkdir -p $qsu/02_clean

fq_filter=$out/$sample/02_clean/$sample.filtered.fq
fq_clipper=$out/$sample/02_clean/$sample.clippered.fq
fq_clean=$out/$sample/02_clean/$sample.clean.fq

echo -ne "######
#PBS -S /bin/bash
#PBS -N clean.$sample
#PBS -q XD64
#PBS -o $log_dir/02_clean/clean.$sample.log
#PBS -e $log_dir/02_clean/clean.$sample.err
#PBS -l nodes=2:ppn=10

echo \"Start Clean raw data of $sample! \" &&
$FASTX_fastq_quality_filter -v -q 20 -p 80 -Q 33 -i $fastq -o $fq_filter &&
$FASTX_clipper -a TGGAATTCTCGGGTGCCAAGG -l 16 -v -Q 33  -i $fq_filter -o $fq_clipper &&
$FASTX_trimmer -v -f 1 -l 26 -Q 33 -i $fq_clipper -o $fq_clean &&
rm $fq_filter $fq_clipper &&
echo \"Clean raw data $sample has been done! \"
######
\n"> $qsu/02_clean/clean.$sample.sh


########################## step3 fastQCclean ##########################
mkdir -p $out/$sample/03_fastQC_clean
mkdir -p $qsu/03_fastQC_clean

fqv_clean=$out/$sample/03_fastQC_clean/$sample.clean.basic.xls

echo -ne "######
#PBS -S /bin/bash
#PBS -N fastQC_clean.$sample
#PBS -q XD64
#PBS -o $log_dir/03_fastQC_clean/fastQC_clean.$sample.log
#PBS -e $log_dir/03_fastQC_clean/fastQC_clean.$sample.err
#PBS -l nodes=2:ppn=10

echo \"Start FastQC clean data of $sample! \" &&
$fastqc --outdir $out/$sample/03_fastQC_clean $fq_clean &&
$fqvalue -q 33 $fq_clean > $fqv_clean &&
echo \"FastQC clean data has been done! \"
######
\n"> $qsu/03_fastQC_clean/fastQCclean.$sample.sh


########################## Step 4  mirdeep2_mapping ##########################
mkdir -p $out/$sample/04_mirdeep2_mapping/
mkdir -p $qsu/04_mirdeep2_mapping

fasta=$out/$sample/04_mirdeep2_mapping/$sample.fasta

fa_collapsed=$out/$sample/04_mirdeep2_mapping/$sample.reads_collapsed.fa
mapping_result=$out/$sample/04_mirdeep2_mapping/$sample.reads_collapsed_vs_genome.arf

mapping_log=$out/$sample/04_mirdeep2_mapping/$sample.mapping_log.txt

echo -ne "######
#PBS -S /bin/bash
#PBS -N mirdeep2_mapping.$sample
#PBS -q XD64
#PBS -o $log_dir/04_mirdeep2_mapping/mirdeep2_mapping.$sample.log
#PBS -e $log_dir/04_mirdeep2_mapping/mirdeep2_mapping.$sample.err
#PBS -l nodes=2:ppn=10
echo \"Start transform fastq to fasta of $sample! \" &&
$fastq2fasta $fq_clean >$fasta &&
echo \"transform fastq to fasta of $sample has been done\" &&
echo \"Start mirdeep2_mapping $sample! \" &&
$mapper $fasta -c -j -k TGGAATTCTCGGGTGCCAAGG -l 18 -m -p $genome_index -s $fa_collapsed -t $mapping_result -u -v 2>$mapping_log &&
echo \"mirdeep2_mapping $sample has been done ! \"
######
\n" >$qsu/04_mirdeep2_mapping/mapper.$sample.sh



########################## Step 5  fast quantitation of reads mapping to known miRBase precursors.##########################
mkdir -p $out/$sample/05_quantifier
mkdir -p $qsu/05_quantifier

quantifier_log=$out/$sample/05_quantifier/$sample.quantifier_log.txt

echo -ne "######
echo \"Start quantifier data of $sample! \" &&
$quantifier \
-p $hairpin_fa -m $mature_fa -r $fa_collapsed -t Human -y 20181121 2>$quantifier_log &&
echo \"Quantifier data of $sample has been done ! \"
######
\n" >$out/$sample/05_quantifier/Quantifier.$sample.sh



########################## Step 6 Indentification microRNA from deep sequencing data ##########################
mkdir -p $out/$sample/06_miRDeep2
mkdir -p $qsu/06_miRDeep2

miRDeep2_log=$out/$sample/06_miRDeep2/$sample.miRDeep2_log.txt

echo -ne "######
#PBS -S /bin/bash
#PBS -N mirdeep2.$sample
#PBS -q batch
#PBS -o $log_dir/04_mirdeep2_mapping/mirdeep2.$sample.log
#PBS -e $log_dir/04_mirdeep2_mapping/mirdeep2.$sample.err
#PBS -l nodes=2:ppn=10

echo \"$sample Start miRDeep2 data! \" &&
$miRDeep2 $fa_collapsed $genome_fa $mapping_result none none none 2>$miRDeep2_log &&
echo \"miRDeep2 data of $sample has been done ! \"
######
\n" >$qsu/06_miRDeep2/miRDeep2.$sample.sh



