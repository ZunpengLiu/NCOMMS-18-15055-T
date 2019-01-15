#!/usr/bin/python
#coding:utf-8

import os
import re
import sys

input=sys.argv[1]
lsd = os.listdir(input)
lsd.sort()
lst=[]

py='/pfs1/liuguanghui/liuzunpeng/01_Script/01_RNA-seq/01_RNA_seq_pipelines/00_Hisat2-HTseq-DESeq2/55_DGCR8/RNA_seq_no_strand.py'
ii=input
out=re.findall(r'(.*/)',py)[0]+'step1.qsub.sh'
oo='/pfs1/liuguanghui/liuzunpeng/05_Results/57_RNA_seq_DGCR8'
cc='/pfs1/liuguanghui/liuzunpeng/01_Script/01_RNA-seq/01_RNA_seq_pipelines/00_Hisat2-HTseq-DESeq2/55_DGCR8/config.txt'
qq='/pfs1/liuguanghui/liuzunpeng/06_Qsub/57_RNA_seq_DGCR8'

for i in lsd:
    if i.endswith('_1.fq.gz'):
        pre=re.findall(r'(.*)_1.fq.gz',i)
        lst.append(pre[0])
uni_lst=list(set(lst))

f1=open(out,'w')
for x in uni_lst:
    tempstr='python3 '+py+' -i '+input+' -o '+oo+' -q '+qq+' -c '+cc+' -p '+x+'\n'
    f1.write(tempstr)

f1.close()
