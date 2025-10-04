#!/bin/bash

GATK=path/to/your/GATK3.8
java=path/to/your/java1.8
ref=path/to/your/haplotype1
genome=path/to/your/haplotype1
bwa=path/to/your/bwa
picard=path/to/your/picard
samtools=path/to/your/samtools

fq1=$1
fq2=$2
sample=$3
echo ${sample}
${bwa} mem -M -t 16 -R "@RG\tID:${sample}\tSM:${sample}\tLB:WGS\tPL:Illumina" ${genome} ${fq1} ${fq2} | ${samtools} sort -@ 16 -o ${sample}.M.bam -
${samtools} index ${sample}.M.bam

bam=${sample}.M.bam
${java} -Xmx4g -jar ${GATK} -T RealignerTargetCreator -nt 8 -R ${ref} -I ${bam} -o ${bam}.intervals
${java} -Xmx4g -jar ${GATK} -T IndelRealigner -R ${ref} -I ${bam} -targetIntervals ${bam}.intervals -o ${bam}.realigned.bam
${java} -jar ${picard} MarkDuplicates -I ${bam}.realigned.bam -M ${bam}_markdup_metrics.txt -O ${bam}.realigned.markdup.bam
${samtools} index ${bam}.realigned.markdup.bam
rm ${bam}.realigned.bam ${bam}.realigned.bam.bai ${bam} ${bam}.bai
