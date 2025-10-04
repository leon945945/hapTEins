#!/bin/bash
#!/bin/bash

show_help() {
    cat << EOF
usage: ${0##*/} [-h] [-f FILE]

parameters:
  -h            help message
  -f FILE       fastq path file; file format: fq1 fq2 samplename
EOF
}

# 默认值
file=""

# 解析参数
while getopts "hf:" opt; do
    case $opt in
        h)
            show_help
            exit 0
            ;;
        f)
            file=$OPTARG
            ;;
        ?)
            show_help >&2
            exit 1
            ;;
    esac
done

cat ${file} | grep -v "#" | while read line
do
	fqs=(${line}); fq1=${fqs[0]}; fq2=${fqs[1]}; sample=${fqs[2]}
	echo $fq1, $fq2, $sample
	#1.将fastq比对到M/P单倍型
	echo "##################################################"
	echo "1.align fastq to genome"
	./align-M.sh ${fq1} ${fq2} ${sample}
	./align-P.sh ${fq1} ${fq2} ${sample}
	echo "done"
	echo "##################################################"
	
	#2.提取reads并比较
	echo "##################################################"
	echo "2.extract reads alignment"
	~/learning/TEdev/getReadsMateTag ${sample}.M.bam.realigned.markdup.bam . > ${sample}.M.readsTag.txt
	~/learning/TEdev/getReadsMateTag ${sample}.P.bam.realigned.markdup.bam . > ${sample}.P.readsTag.txt
	echo "done"
	echo "##################################################"
	
	#3.比较reads
	echo "##################################################"
	echo "3.compare reads alignment"
	./readsCompare.py ${sample}.M.readsTag.txt ${sample}.P.readsTag.txt ${sample}
	echo "done"
	echo "##################################################"
	
	#4.分配reads
	echo "##################################################"
	echo "4.assign reads"
	seqkit grep -f ${sample}.M.reads ${fq1} -o ${sample}.M.reads.fq1.gz
	seqkit grep -f ${sample}.M.reads ${fq2} -o ${sample}.M.reads.fq2.gz
	seqkit grep -f ${sample}.P.reads ${fq1} -o ${sample}.P.reads.fq1.gz
	seqkit grep -f ${sample}.P.reads ${fq2} -o ${sample}.P.reads.fq2.gz
	echo "done"
	echo "##################################################"
	
	#5.将分配到M/P单倍型的fastq分别比对到M/P单倍型
	echo "##################################################"
	echo "5.align assigned reads"
	./align-M_TE.sh ${sample}.M.reads.fq1.gz ${sample}.M.reads.fq2.gz ${sample}
	./align-P_TE.sh ${sample}.P.reads.fq1.gz ${sample}.P.reads.fq2.gz ${sample}
	echo "done"
	echo "##################################################"
	
	#6.比对后的bam提取reads比对情况
	~/learning/TEdev/getReadsMateTag ${sample}.M_TE.bam.realigned.markdup.bam . > ${sample}.M_TE.readsTag.txt
	~/learning/TEdev/getReadsMateTag ${sample}.P_TE.bam.realigned.markdup.bam . > ${sample}.P_TE.readsTag.txt
	
	#7.清理文件
	rm ${sample}.M.bam.realigned.markdup.bam* ${sample}.P.bam.realigned.markdup.bam* ${sample}.M.readsTag.txt ${sample}.P.readsTag.txt ${sample}.M.reads ${sample}.P.reads ${sample}.M.reads.fq1.gz ${sample}.M.reads.fq2.gz ${sample}.P.reads.fq1.gz ${sample}.P.reads.fq2.gz *intervals *_metrics.txt *realigned.bai
done
