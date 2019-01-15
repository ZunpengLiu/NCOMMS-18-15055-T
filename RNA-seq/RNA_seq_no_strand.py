#!/usr/bin/env python3
# coding: utf-8

from optparse import OptionParser
import os
import sys
import re

'''
Title:        RNA-seq pipeline
Description:  This pipeline is used for ATAC-seq analysis and dowanstream analysis.

Author:       Zunpeng Liu
Email:        zunpengliuAT163.com
Date:         12/6/2018
Version:      1.1
Python:       based on python 3.6
Citation:     please
'''

MY_USAGE = '''
python3  RNA-seq_v1.py  [options]

Example: python3 RNA-seq_v1.py -i /pfs1/liuguanghui/liuzunpeng/02_Data/03_ATAC-seq/04_Lamin-KI/raw_data -p M91_HHVK7CCXY_L1 \
-o /pfs1/liuguanghui/liuzunpeng/05_Results/20_Lamin_KI
-p /pfs1/liuguanghui/liuzunpeng/06_Qsub/20_Lamin_KI
-c /pfs1/liuguanghui/liuzunpeng/01_Script/03_ATAC-seq/07_Lamin_KI/config.txt
-r /pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/hg19.fa
-g /pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/Homo_sapiens.GRCh37.87_chr.chr.gtf
-x /pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/02_bowtie2_index/hg19
-b /pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/04_bed/hg19_RefSeq.bed

Author:       Zunpeng Liu
Email:        zunpengliu@163.com
Date:         16/5/2018

Blessings:    If you have any question or query of this pipeline, please feel free to contact with me.
'''

parser = OptionParser(MY_USAGE, version='Version 1.1')
parser.add_option('-i', '--in_dir', dest='in_dir', type='string',
                  help='You should sepecify the input file path, wich stored all raw data.')
parser.add_option('-p', '--prefix', dest='prefix', type='string',
                  help='You should specifty the prefix of samples, which will be used to label all processed results. For example, if the prefix was defined as "M91_HHVK7CCXY_L1", the prefix of file names will be precessed as "M91".')
parser.add_option('-o', '--out_dir', dest='out_dir', type='string',
                  help='You should specifty the output directory, which will store all processed results.')
parser.add_option('-q', '--qsub_dir', dest='qsub_dir', type='string',
                  help='You should specifty the qub directory, which will store all shell script derived from this pipeline.')
parser.add_option('-c', '--configure_file', dest='configure_file', type='string',
                  default='/pfs1/liuguanghui/liuzunpeng/01_Script/03_ATAC-seq/07_ATAC_seq_pipeline/config.txt',
                  help='The configure lib file should include the absolute path of softwares used in this pipeline')
parser.add_option('-g', '--gtf', dest='gtf', type='string',
                  default='/pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/Homo_sapiens.GRCh37.87_chr.chr.gtf',
                  help='Qub directory of all processing data.')
parser.add_option('-r', '--reference', dest='reference', type='string',
                  default='/pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/hg19.fa',
                  help='Selectable opotion. If giving, the reference will be used for building the index.')
parser.add_option('-x', '--index', dest='index', type='string',
                  default='/pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/02_bowtie2_index/hg19',
                  help='The index file of genome, which will be used for aligement by bowtie2 should be build by bowtie2-build tool before using.')
parser.add_option('-b', '--RefSeq_bed', dest='RefSeq_bed', type='string',
                  default='/pfs1/liuguanghui/liuzunpeng/03_Database/00_genome/01_hg19/03_Ensembl/01_fa_gtf_bed_chr/02_bowtie2_index/hg19',
                  help='The refseq bed used for genome annotation')
(options, args) = parser.parse_args()

''' Predefine variables '''
in_dir = str(options.in_dir)
prefix = str(options.prefix)
sample  = prefix
out_dir = str(options.out_dir)
qsu_dir = str(options.qsub_dir)
log_dir = '/home/liuzunpeng/../..' + qsu_dir
ref = str(options.reference)
gtf = str(options.gtf)
bowtie2_index = str(options.index)

''' ensure all directory has been created, if not exist,then creat the new ones. '''


def ensure_dir(directory):
    # directory = os.path.dirname(f)
    if not os.path.exists(directory):
        os.makedirs(directory)


''' make new directions '''
ensure_dir(out_dir)
ensure_dir(qsu_dir)

''' Configure all softwares and library '''
# configure all softwares and sources from configure.txt
sft = {}
config = options.configure_file
with open(config, 'r') as config:
    for line in config:
        #        if len(line)!=0 and not line.startswith('#'):
        if len(re.findall('=', line)) == 1:
            sft[re.findall(r'(\w+?)=.*', line)[0]] = re.findall(r'\w+=(.*)', line)[0]
globals().update(sft)

''' Raw data files '''
fq1 = in_dir + '/' + prefix + '_1.fq.gz'
fq2 = in_dir + '/' + prefix + '_2.fq.gz'

''' 1. FastQC rawdata '''
ensure_dir(str(out_dir + '/' + sample + '/01_fastQC_raw'))
ensure_dir(str(qsu_dir + '/01_fastQC_raw'))

fqv1 = out_dir + '/' + sample + '/01_fastQC_raw/' + sample + '.R1.basic.xls'
fqv2 = out_dir + '/' + sample + '/01_fastQC_raw/' + sample + '.R2.basic.xls'

FastQC_str = '\
#PBS -S /bin/bash\n\
#PBS -N fastQC_raw.' + sample + '\n\
#PBS -q batch\n\
#PBS -o ' + log_dir + '/01_fastQC_raw/fastQC_raw.' + sample + '.log\n\
#PBS -e ' + log_dir + '/01_fastQC_raw/fastQC_raw.' + sample + '.err\n\
#PBS -l nodes=1:batch:ppn=1\n\n\
\
echo "Start fastQC and fqvalue raw data of ' + sample + '" && \n' \
             + fastqc + ' --outdir ' + out_dir + '/' + sample + '/01_fastQC_raw -f fastq ' + fq1 + ' ' + fq2 + ' && \n' \
             + fqvalue + ' -q 33 ' + fq1 + ' > ' + fqv1 + ' && \n' \
             + fqvalue + ' -q 33 ' + fq2 + ' > ' + fqv2 + ' && \n' \
             + 'echo "FastQC and fqvalue raw data of ' + sample + ' has been done! " '

FastQC_shell_file = open(qsu_dir + '/01_fastQC_raw/fastQC_raw.' + sample + '.sh', 'w')
FastQC_shell_file.write(FastQC_str)
FastQC_shell_file.close()

''' 2. TrimGalore clean '''
ensure_dir(out_dir + '/' + sample + '/02_TrimGalore')
ensure_dir(qsu_dir + '/02_TrimGalore')

clean_fq1 = out_dir + '/' + sample + '/02_TrimGalore/' + prefix + '_1_val_1.fq.gz'
clean_fq2 = out_dir + '/' + sample + '/02_TrimGalore/' + prefix + '_2_val_2.fq.gz'

trim_galore_log = out_dir + '/' + sample + '/02_TrimGalore/' + sample + '.trim_galore.log'

trim_galore_str = '\
#PBS -S /bin/bash\n\
#PBS -N trim_galore.' + sample + '\n\
#PBS -q batch\n\
#PBS -o ' + log_dir + '/02_TrimGalore/trim_galore.' + sample + '.log\n\
#PBS -e ' + log_dir + '/02_TrimGalore/trim_galore.' + sample + '.err\n\
#PBS -l nodes=1:batch:ppn=12\n\n\
\
echo "Start trim raw data of ' + sample + '" && \n' \
                  + trim_galore + ' --fastqc --path_to_cutadapt ' + cutadapt + ' --stringency 3 --paired --output_dir ' + out_dir + '/' + sample + '/02_TrimGalore ' \
                  + fq1 + ' ' + fq2 + ' 2>' + trim_galore_log + ' && \n\
echo "Clean and FastQC of ' + sample + ' done! All done! \"'

trim_galore_shell_file = open(qsu_dir + '/02_TrimGalore/trim_galore.' + sample + '.sh', 'w')
trim_galore_shell_file.write(trim_galore_str)
trim_galore_shell_file.close()

''' 3. FastQC clean '''
ensure_dir(out_dir + '/' + sample + '/03_fastQCclean')
ensure_dir(qsu_dir + '/03_fastQCclean')

fqvC1 = out_dir + '/' + sample + '/03_fastQCclean/' + sample + '.R1.clean.basic.xls'
fqvC2 = out_dir + '/' + sample + '/03_fastQCclean/' + sample + '.R2.clean.basic.xls'

Basicresult = out_dir + '/' + sample + '/03_fastQCclean/' + sample + '.BasicInfo.xls'

FastQC_clean_str = '\
#PBS -S /bin/bash\n\
#PBS -N fastQCclean.' + sample + '\n\
#PBS -q batch\n\
#PBS -o ' + log_dir + '/03_fastQCclean/fastQCclean.' + sample + '.log\n\
#PBS -e ' + log_dir + '/03_fastQCclean/fastQCclean.' + sample + '.err\n\
#PBS -l nodes=1:batch:ppn=1\n\n\
\
echo "Start FastQC Clean data of ' + sample + '" && \n' \
                   + fqvalue + ' -q 33 ' + clean_fq1 + ' > ' + fqvC1 + ' && \n' \
                   + fqvalue + ' -q 33 ' + clean_fq2 + ' > ' + fqvC2 + ' && \n\
echo "Fqvalue raw data of ' + sample + ' has been done! \"'

FastQC_clean_shell_file = open(qsu_dir + '/03_fastQCclean/fastQCclean.' + sample + '.sh', 'w')
FastQC_clean_shell_file.write(FastQC_clean_str)
FastQC_clean_shell_file.close()

''' 4. Mapping by hisat2 '''
ensure_dir(out_dir + '/' + sample + '/04_hisat2')
ensure_dir(qsu_dir + '/04_hisat2')

hisat2_alignLog = out_dir + '/' + sample + '/04_hisat2/' + sample + '.hisat2_alignLog.txt'
hisat2_Alignresult = out_dir + '/' + sample + '/04_hisat2/' + sample + '.hisat2_AlignRate.xls'

hisat2_Alignresult_file = open(hisat2_Alignresult, 'w')
hisat2_Alignresult_str = 'Sample_ID\tTotalReadsPair\tUniqPairAlignConcordantly\tMultiPairAlignConcordantly\tuniqPairAlignDisconcordantly\tUniqSingleAlign\tMultiSingleAlign\tTotalMapped'
hisat2_Alignresult_file.write(hisat2_Alignresult_str)
FastQC_clean_shell_file.close()

hisat_sam = out_dir + '/' + sample + '/04_hisat2/' + sample + '.sam'
hisat_bam = out_dir + '/' + sample + '/04_hisat2/' + sample + '.bam'
hisat_sorted_bam = out_dir + '/' + sample + '/04_hisat2/' + sample + '.sort.bam'

hisat_str = '\
#PBS -S /bin/bash\n\
#PBS -N hisat2.' + sample + '\n\
#PBS -q XD64\n\
#PBS -o ' + log_dir + '/04_hisat2/hisat2.' + sample + '.log\n\
#PBS -e ' + log_dir + '/04_hisat2/hisat2.' + sample + '.err\n\
#PBS -l nodes=1:XD64:ppn=12\n\n\
\
echo "Start hisat2 mapping of ' + sample + '" && \n' \
            + hisat2 + ' -p 12 -x ' + hisat2_index + ' -1 ' + clean_fq1 + ' -2 ' + clean_fq2 + ' --dta  -S ' + hisat_sam + ' --un-gz ' + out_dir + '/' + sample + '/04_hisat2 2>' + hisat2_alignLog + ' && \n' \
            + 'echo "Alignment to reference of ' + sample + ' by hisat2 has been done! \" && \n' \
            + samtools + ' view -@ 12 -S -b -q 10 ' + hisat_sam + ' > ' + hisat_bam + ' && \n' \
            + 'echo "Samtools transfer sam to bam of ' + sample + ' has been done! \" && \n' \
            + samtools + ' sort -@ 12 ' + hisat_bam + ' -o ' + hisat_sorted_bam + ' && \n' \
            + 'echo "Samtools sort bam files of ' + sample + ' has been done! \" && \n' \
            + samtools + ' index -@ 12 ' + hisat_sorted_bam + ' && \n' \
            + samtools + ' stats -@ 12 ' + hisat_sorted_bam + ' && \n' \
            + 'echo "Samtools index and stats of ' + sample + ' has been done! \" && \n' \
            + 'rm ' + hisat_sam + ' && \n' \
            + 'echo "All done! " '

hisat_shell_file = open(qsu_dir + '/04_hisat2/hisat2.' + sample + '.sh', 'w')
hisat_shell_file.write(hisat_str)
hisat_shell_file.close()

''' 5. Mapping by tophat '''
ensure_dir(out_dir + '/' + sample + '/05_tophat')
ensure_dir(qsu_dir + '/05_tophat')

tophat_alignLog = out_dir + '/' + sample + '/05_tophat/' + sample + '.tophat_alignLog.txt'

tophat_sam = out_dir + '/' + sample + '/05_tophat/' + sample + '.sam'
tophat_bam = out_dir + '/' + sample + '/05_tophat/' + sample + '.bam'

tophat_str = '\
#PBS -S /bin/bash\n\
#PBS -N tophat.' + sample + '\n\
#PBS -q XD64\n\
#PBS -o ' + log_dir + '/05_tophat/tophat.' + sample + '.log\n\
#PBS -e ' + log_dir + '/05_tophat/tophat.' + sample + '.err\n\
#PBS -l nodes=1:XD64:ppn=12\n\n\
\
echo "Start tophat mapping of ' + sample + '" && \n' \
             + tophat + ' -p 12 --library-type fr-unstranded -o ' + out_dir + '/' + sample + '/05_tophat -G ' + gtf + ' ' + bowtie2_index + ' ' + clean_fq1 + '  ' + clean_fq2 + ' 2>' + tophat_alignLog + ' && \n' \
             + 'echo "Alignment to reference of ' + sample + ' by tophat has been done! \" && \n' \
             + samtools + ' view -@ 12 -S -b -q 10 ' + tophat_sam + ' > ' + tophat_bam + ' && \n' \
             + 'echo "Samtools transfer sam to bam of ' + sample + ' has been done! \" && \n' \
             + samtools + ' index -@ 12 ' + tophat_bam + ' && \n' \
             + samtools + ' stats -@ 12 ' + tophat_bam + ' && \n' \
             + 'echo "Samtools index and stats of ' + sample + ' has been done! \" && \n' \
             + 'echo "All done! " '

tophat_shell_file = open(qsu_dir + '/05_tophat/tophat.' + sample + '.sh', 'w')
tophat_shell_file.write(tophat_str)
tophat_shell_file.close()

''' 5. Get_bw file '''
ensure_dir(out_dir + '/' + sample + '/06_Bamcovrage')
ensure_dir(qsu_dir + '/06_Bamcovrage')

bin10_bw = out_dir + '/' + sample + '/06_Bamcovrage/' + sample + '.normalizeUsingRPKM_bin10.bw'
bin10_bed = out_dir + '/' + sample + '/06_Bamcovrage/' + sample + '.normalizeUsingRPKM_bin10.bed'
bedGraph_file = out_dir + '/' + sample + '/06_Bamcovrage/' + sample + '.sort_rmdup.bedGraph'
bw_file = out_dir + '/' + sample + '/06_Bamcovrage/' + sample + '.sort_rmdup.bw'
Bamcovrage_log = out_dir + '/' + sample + '/06_Bamcovrage/' + sample + '.Bamcovrage.log'

Bamcovrage_str = '\
#PBS -S /bin/bash\n\
#PBS -N Bamcovrage.' + sample + '\n\
#PBS -q batch\n\
#PBS -o ' + log_dir + '/06_Bamcovrage/Bamcovrage.' + sample + '.log\n\
#PBS -e ' + log_dir + '/06_Bamcovrage/Bamcovrage.' + sample + '.err\n\
#PBS -l nodes=1:batch:ppn=1\n\n\
\
echo "Start deeptools bamCoverage of ' + sample + ' at binSize 10 bp!" && \n' \
                 + bamCoverage + ' -b ' + hisat_sorted_bam + ' --normalizeUsingRPKM --binSize 10 --ignoreForNormalization chrX chrM chrY -o ' + bin10_bw + ' 2>' + Bamcovrage_log + ' && \n' \
                 + bamCoverage + ' -b ' + hisat_sorted_bam + ' --outFileFormat bedgraph --normalizeUsingRPKM --binSize 10 --ignoreForNormalization chrX chrM chrY -o ' + bin10_bed + ' 2>> ' + Bamcovrage_log + ' && \n' \
                 + 'echo "Finish deeptools bamCoverage of ' + sample + ' at binSize 10 bp!\" && \n' \
                 + genomeCoverageBed + ' -bg -split -g ' + chromosomesize + ' -ibam ' + hisat_sorted_bam + ' >' + bedGraph_file + ' 2>>' + Bamcovrage_log + ' && \n' \
                 + bedGraphToBigWig + ' ' + bedGraph_file + ' ' + chromosomesize + ' ' + bw_file + ' 2>>' + Bamcovrage_log + ' && \n' \
                 + 'echo "bedGraphToBigWig make bedGraph of ' + sample + ' has been done!\" && \n' \
                 + 'echo "All done! " '

Bamcovrage_shell_file = open(qsu_dir + '/06_Bamcovrage/Bamcovrage.' + sample + '.sh', 'w')
Bamcovrage_shell_file.write(Bamcovrage_str)
Bamcovrage_shell_file.close()

''' 6. Count reads by HTseq '''
ensure_dir(out_dir + '/' + sample + '/07_HTseq_count')
ensure_dir(qsu_dir + '/07_HTseq_count')

count_res = out_dir + '/' + sample + '/07_HTseq_count/' + sample + '.count'
count_txt = out_dir + '/' + sample + '/07_HTseq_count/' + sample + '.count.txt'
HTseq_log = out_dir + '/' + sample + '/07_HTseq_count/' + sample + '.HTseq.log.txt'

HTseq_str = '\
#PBS -S /bin/bash\n\
#PBS -N Count.' + sample + '\n\
#PBS -q XD64\n\
#PBS -o ' + log_dir + '/07_HTseq_count/Count.' + sample + '.log\n\
#PBS -e ' + log_dir + '/07_HTseq_count/Count.' + sample + '.err\n\
#PBS -l nodes=1:XD64:ppn=12\n\n\
\
echo "Start HTseq count of ' + sample + ' !" && \n' \
            + HTseq + ' -r pos -s no -a 20 -f bam ' + hisat_sorted_bam + ' ' + gtf + ' > ' + count_res + ' 2>' + HTseq_log + ' && \n' \
            + 'sed \'1i GeneID\\t' + sample + '\' ' + count_res + ' >' + count_txt + ' && \n' \
            + 'echo "Finish HTseq count of ' + sample + '!\" && \n' \
            + 'echo "All done! " '

HTseq_shell_file = open(qsu_dir + '/07_HTseq_count/Count.' + sample + '.sh', 'w')
HTseq_shell_file.write(HTseq_str)
HTseq_shell_file.close()

''' 8. Count reads by StringTie '''

ensure_dir(out_dir + '/' + sample + '/08_stirngTie')
ensure_dir(qsu_dir + '/08_stirngTie')

stirngTie_res = out_dir + '/' + sample + '/08_stirngTie/' + sample + '.gene_abund.tab'
stirngTie_txt = out_dir + '/' + sample + '/08_stirngTie/' + sample + 'gene_abund.txt'
stirngTie_FPKM = out_dir + '/' + sample + '/08_stirngTie/' + sample + 'gene_FPKM.txt'
stirngTie_TPM = out_dir + '/' + sample + '/08_stirngTie/' + sample + 'gene_TPM.txt'

stirngTie_gtf = out_dir + '/' + sample + '/08_stirngTie/' + sample + '.gtf'

stirngTie_log = out_dir + '/' + sample + '/08_stirngTie/' + sample + '.stirngTie.log.txt'

stirngTie_str = '\
#PBS -S /bin/bash\n\
#PBS -N stirngTie.' + sample + '\n\
#PBS -q XD64\n\
#PBS -o ' + log_dir + '/08_stirngTie/stirngTie.' + sample + '.log\n\
#PBS -e ' + log_dir + '/08_stirngTie/stirngTie.' + sample + '.err\n\
#PBS -l nodes=1:XD64:ppn=12\n\n\
\
echo "Start stirngTie count of ' + sample + ' !" && \n' \
                + stringTie + ' -A ' + stirngTie_res + ' -e -B -p 12 -G ' + gtf + ' -o ' + stirngTie_gtf + ' ' + hisat_sorted_bam + ' 2>' + stirngTie_log + ' && \n' \
                + 'cat ' + stirngTie_res + ' | awk \'{print $1"\\t"$2"\\t"$8"\\t"$9}\' |sed \'1c GeneID\\tGene_Name\\t' + sample + '_FPKM\\t' + sample + '_TPM\' >' + stirngTie_txt + ' && \n' \
                + 'cat ' + stirngTie_res + ' | awk \'{print $1"\\t"$8} \'|sed \'1c GeneID\\t' + sample + '_FPKM\' >' + stirngTie_FPKM + ' && \n' \
                + 'cat ' + stirngTie_res + ' | awk \'{print $1"\\t"$9}\' |sed \'1c GeneID\\t' + sample + '_TPM\' >' + stirngTie_TPM + ' && \n' \
                + 'echo "Finish StirngTie count of ' + sample + '!\" && \n' \
                + 'echo "All done! " '

stirngTie_shell_file = open(qsu_dir + '/08_stirngTie/stirngTie.' + sample + '.sh', 'w')
stirngTie_shell_file.write(stirngTie_str)
stirngTie_shell_file.close()

''' 9. Count reads by Salmon '''
ensure_dir(out_dir + '/' + sample + '/09_Salmon')
ensure_dir(qsu_dir + '/09_Salmon')

Salmon_out = out_dir + '/' + sample + '/09_Salmon/' + sample + '.quant'
Salmon_res = out_dir + '/' + sample + '/09_Salmon/' + sample + '.quant/quant.genes.sf'
Salmon_txt = out_dir + '/' + sample + '/09_Salmon/' + sample + 'genes.TPM.txt'
Salmon_log = out_dir + '/' + sample + '/09_Salmon/' + sample + '.Salmon.log.txt'

Salmon_str = '\
#PBS -S /bin/bash\n\
#PBS -N Salmon.' + sample + '\n\
#PBS -q batch\n\
#PBS -o ' + log_dir + '/09_Salmon/Salmon.' + sample + '.log\n\
#PBS -e ' + log_dir + '/09_Salmon/Salmon.' + sample + '.err\n\
#PBS -l nodes=1:batch:ppn=12\n\n\
\
echo "Start stirngTie count of ' + sample + ' !" && \n' \
             + Salmon + ' quant -i ' + salmon_index + ' -l A -g ' + gtf + ' -1 ' + clean_fq1 + ' -2 ' + clean_fq2 + ' -p 12 -o ' + Salmon_out + ' 2> ' + Salmon_log + ' && \n' \
             + 'awk \'{print $1"\\t"$4}\' ' + Salmon_res + ' |sed \'1c GeneID\\t' + sample + '_TPM\' |grep -v "ENST" >' + Salmon_txt + ' 2> ' + Salmon_log + ' && \n' \
             + 'echo "Finish Salmon count of ' + sample + '!\" && \n' \
             + 'echo "All done! " '

Salmon_shell_file = open(qsu_dir + '/09_Salmon/Salmon.' + sample + '.sh', 'w')
Salmon_shell_file.write(Salmon_str)
Salmon_shell_file.close()

''' 10. Count reads by cufflinks '''
ensure_dir(out_dir + '/' + sample + '/10_cufflinks')
ensure_dir(qsu_dir + '/10_cufflinks')

cufflinks_str = '\
#PBS -S /bin/bash\n\
#PBS -N cufflinks.' + sample + '\n\
#PBS -q batch\n\
#PBS -o ' + log_dir + '/10_cufflinks/cufflinks.' + sample + '.log\n\
#PBS -e ' + log_dir + '/10_cufflinks/cufflinks.' + sample + '.err\n\
#PBS -l nodes=1:batch:ppn=12\n\n\
\
echo "Start cufflinks count of ' + sample + ' !" && \n' \
                + cufflinks + ' -o ' + out_dir + '/' + sample + '/10_cufflinks -p 12 --library-type fr-firststrand -u -g ' + gtf + ' ' + tophat_bam + ' && \n' \
                + 'echo "Finish cufflinks count of ' + sample + '!\" && \n' \
                + 'echo "All done! " '

cufflinks_shell_file = open(qsu_dir + '/10_cufflinks/cufflinks.' + sample + '.sh', 'w')
cufflinks_shell_file.write(cufflinks_str)
cufflinks_shell_file.close()

''' 11. Repenrich_1_bowti2_mapping '''
ensure_dir(out_dir + '/' + sample + '/11_Repenrich_step1_bowti2_mapping')
ensure_dir(qsu_dir + '/11_Repenrich_step1_bowti2_mapping')

Repenrich_step1_bowtie2_Align_log = out_dir + '/' + sample + '/11_Repenrich_step1_bowti2_mapping/' + sample + '.bowtie2_Align_log.txt'

rep_sam = out_dir + '/' + sample + '/11_Repenrich_step1_bowti2_mapping/' + sample + '.sam'
rep_bam = out_dir + '/' + sample + '/11_Repenrich_step1_bowti2_mapping/' + sample + '.bam'

Repenrich_step1_bowti2_str = '\
#PBS -S /bin/bash\n\
#PBS -N Repenrich_1_bowtie2.' + sample + '\n\
#PBS -q batch\n\
#PBS -o ' + log_dir + '/11_Repenrich_step1_bowti2_mapping/Repenrich_1_bowti2.' + sample + '.log\n\
#PBS -e ' + log_dir + '/11_Repenrich_step1_bowti2_mapping/Repenrich_1_bowti2.' + sample + '.err\n\
#PBS -l nodes=1:batch:ppn=12\n\n\
\
echo "Start the alignment of ' + sample + ' to reference by bowtie2!" && \n' \
+ bowtie2 + ' -q -p 12 -x ' + bowtie2_index + ' -t -1 ' + clean_fq1 + ' -2 ' + clean_fq2 + ' -S ' + rep_sam + ' 2>' + Repenrich_step1_bowtie2_Align_log + ' && \n' \
+ 'echo "Alignment to reference of ' + sample + ' by bowtie2 has been done! \" && \n' \
+ 'echo "Start transfering the sam file of ' + sample + ' to bam file and remove chrM and chrY DNA!" && \n' \
 + 'awk \'$3!="chrM" && $3!="chrY" \' ' + rep_sam + ' | ' + samtools + ' view -S -b > ' + rep_bam + ' && \n' \
+ 'echo "Repenrich_step1_bowti2_mapping of ' + sample + ' done!" '

Repenrich_1_bowti2_mapping_shell_file = open(
    qsu_dir + '/11_Repenrich_step1_bowti2_mapping/Repenrich_step1_bowti2_mapping.' + sample + '.sh', 'w')
Repenrich_1_bowti2_mapping_shell_file.write(Repenrich_step1_bowti2_str)
Repenrich_1_bowti2_mapping_shell_file.close()

''' 12. Repenrich_step2_subset '''

ensure_dir(out_dir + '/' + sample + '/12_Repenrich_step2_subset_count')
ensure_dir(qsu_dir + '/12_Repenrich_step2_subset_count')

rep_r1fastq = out_dir + '/' + sample + '/12_Repenrich_step2_subset_count/' + sample + '_multimap_R1.fastq'
rep_r2fastq = out_dir + '/' + sample + '/12_Repenrich_step2_subset_count/' + sample + '_multimap_R2.fastq'
rep_uniqueBAM = out_dir + '/' + sample + '/12_Repenrich_step2_subset_count/' + sample + '_unique.bam'

Repenrich_step2_subset_log = out_dir + '/' + sample + '/12_Repenrich_step2_subset_count/' + sample + '.Repenrich_step2_subset_log.txt'

Repenrich_step2_str = '\
#PBS -S /bin/bash\n\
#PBS -N Repenrich_step2.' + sample + '\n\
#PBS -q XD128\n\
#PBS -o ' + log_dir + '/12_Repenrich_step2_subset_count/Repenrich_step2_subset.' + sample + '.log\n\
#PBS -e ' + log_dir + '/12_Repenrich_step2_subset_count/Repenrich_step2_subset.' + sample + '.err\n\
#PBS -l nodes=1:XD128:ppn=12\n\n\
\
echo "Start unique and multi-mapped reads of ' + sample + ' !" && \n' \
                      + python + ' ' + RepEnrich2_subset + ' ' + rep_bam + ' 30 ' + out_dir + '/' + sample + '/12_Repenrich_step2_subset_count/' + sample + ' --pairedend TRUE  2>' + Repenrich_step2_subset_log + ' && \n' \
                      + 'echo "Subset unique and multi-mapped reads of ' + sample + ' has been done! \" && \n' \
                      + python + ' ' + RepEnrich2 + ' ' + repeatmasker + ' ' + out_dir + '/' + sample + '/12_Repenrich_step2_subset_count ' + sample + ' ' + hg19setup + ' ' + rep_r1fastq + ' --fastqfile2 ' + rep_r2fastq + ' ' + rep_uniqueBAM + ' --cpus 12 --pairedend TRUE  2>' + Repenrich_step2_subset_log + ' && \n' \
                      + 'echo "Repenrich of ' + sample + ' done!" '

Repenrich_step2_shell_file = open(qsu_dir + '/12_Repenrich_step2_subset_count/Repenrich_step2.' + sample + '.sh', 'w')
Repenrich_step2_shell_file.write(Repenrich_step2_str)
Repenrich_step2_shell_file.close()






''' 13. Mapping by bowtie2 '''
ensure_dir(out_dir + '/' + sample + '/13_bowtie2')
ensure_dir(qsu_dir + '/13_bowtie2')

bowtie2_alignLog = out_dir + '/' + sample + '/13_bowtie2/' + sample + '.bowtie2_alignLog.txt'

bowtie2_sam = out_dir + '/' + sample + '/13_bowtie2/' + sample + '.sam'
bowtie2_bam = out_dir + '/' + sample + '/13_bowtie2/' + sample + '.bam'
bowtie2_sorted_bam = out_dir + '/' + sample + '/13_bowtie2/' + sample + '.sort.bam'

bowtie2_str = '\
#PBS -S /bin/bash\n\
#PBS -N bowtie2.' + sample + '\n\
#PBS -q XD64\n\
#PBS -o ' + log_dir + '/13_bowtie2/bowtie2.' + sample + '.log\n\
#PBS -e ' + log_dir + '/13_bowtie2/bowtie2.' + sample + '.err\n\
#PBS -l nodes=1:XD64:ppn=12\n\n\
\
echo "Start bowtie2 mapping of ' + sample + '" && \n' \
+ bowtie2 + ' -q -p 12 -x ' + bowtie2_rep_index + ' -t -1 ' + clean_fq1 + ' -2 ' + clean_fq2 + ' -S ' + bowtie2_sam + ' 2>' + bowtie2_alignLog + ' && \n' \
+ 'echo "Alignment to reference of ' + sample + ' by bowtie2 has been done! \" && \n' \
+ samtools + ' view -@ 12 -S -b -q 10 ' + bowtie2_sam + ' > ' + bowtie2_bam + ' && \n' \
+ samtools + ' sort -@ 12 -l 9 ' + bowtie2_bam + ' -o ' + bowtie2_sorted_bam + ' && \n' \
+ samtools + ' index -@ 12 ' + bowtie2_sorted_bam + ' && \n' \
+ samtools +' idxstats ' + bowtie2_sorted_bam + ' > ' + out_dir + '/'+sample+'/13_bowtie2/bow_rep.' + sample +'.idxstats.txt && \n'\
+ samtools +' flagstat ' + bowtie2_sorted_bam + ' > ' + out_dir + '/'+sample+'/13_bowtie2/bow_rep.' + sample +'.flagstat.txt && \n'\
+ samtools +' stats ' + bowtie2_sorted_bam + ' > ' + out_dir + '/'+sample+'/13_bowtie2/bow_rep.' + sample +'_.all_stats.txt && \n'\
+ 'echo "All done! " '

bowtie2_shell_file = open(qsu_dir + '/13_bowtie2/bowtie2.' + sample + '.sh', 'w')
bowtie2_shell_file.write(bowtie2_str)
bowtie2_shell_file.close()

