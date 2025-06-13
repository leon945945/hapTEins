#!/bin/bash
file=$1
cat ${file} | grep -v "#" | while read line
do
	fields=(${line}); bam=${fields[0]}; readsTag=${fields[1]}; fq1=${fields[2]}; fq2=${fields[3]}
	tmp=$(basename ${bam}); sample=${tmp%%_*}
	echo ${sample} ${bam} ${readsTag} ${fq1} ${fq2}
	#1.处理比对到TE的reads，并检测插入位点
	echo "1.detect reads aligned to TE sequence"
	./getPos.py ${bam} ${readsTag} ${sample}
	#2.处理没有比对到TE的discordant reads，并检测TE插入位点
	echo "2.detect discordant reads which not aligned to TE sequence"
	samtools view -F 1294 ${bam} | grep -v "TE" | cut -f 1-9 > ${sample}.sam
	samtools depth -g 256 -a ${bam} > ${sample}.depth
	SAM=${sample}.sam
	DEP=${sample}.depth
	./fltDiscordant.py ${SAM} ${bam} ${DEP} ${fq1} ${fq2} ${sample}
	#3.清理无用文件
	echo "3.remove sam and depth files"
	rm ${SAM} ${DEP}
done
