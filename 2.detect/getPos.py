#!/home/songlizhi/software/miniconda3/bin/python
from collections import Counter
import pandas as pd
import subprocess
import os
import sys
bam = sys.argv[1]
readsTag = sys.argv[2]
sample = sys.argv[3]
TEs = ["TE3002","TE4474","TE5564","TE5813","TE6353","TE6936","TE12547"]
for TE in TEs:
	mL = int(TE.replace("TE",""))
	#1.提取比对到TE的reads
	os.system("grep {TE} {RT} | sort -k 12n > {sample}.txt".format(TE=TE,RT=readsTag,sample=sample))
	dat = pd.read_csv("{sample}.txt".format(sample=sample),sep="\t",names=["readName","chr","begin","AS","NM","alignLength","cigarNum","readLength","ratio","type","Mchr","Mbegin"])
	#2.提取discordant reads
	dat = dat.loc[~dat["Mchr"].str.match("TE"),:].reset_index(drop=True)
	if dat.empty:
		continue
	#3.根据reads长度计算插入区间
	beg = 0
	binL = []
	begP = dat.loc[0,"Mbegin"]
	for idx in dat.index:
		if dat.loc[idx,"Mbegin"] - begP > mL:
			binL.append(dat.iloc[beg:idx,])
			beg = idx
			begP = dat.loc[idx,"Mbegin"]
		else:
			continue
	binL.append(dat.iloc[beg:,])
	#4.对每个插入区间计算reads的平均插入位置
	avePos = []
	for df in binL:
		pos = df["Mbegin"].apply("mean")
		avePos.append(int(pos))
	posChr = [i.reset_index().loc[0,"Mchr"] for i in binL]
	#5.在插入平均位置提取soft clip位置
	softC = []
	for chrom,pos in zip(posChr,avePos):
		TEregion = chrom+":"+str(pos-500)+"-"+str(pos+500) if pos> 500 else chrom+":"+"1"+"-"+str(pos+500)
		os.system("/home/songlizhi/learning/TEdev/findSplit {bam} {TEregion} > {sample}.softClip".format(bam=bam,TEregion=TEregion,sample=sample))
		softC.append(pd.read_csv("{sample}.softClip".format(sample=sample),sep="\t",names=["reads","clipPos"]))
	#6.对reads的soft clip位置统计数量，数量最多的前两个的距离小于20bp，认为是TE插入的候选区间
	softcN = []; softClipN = []; regionIs = []; regionLs = []; regionRs = []
	for i in softC:
		counter = Counter(i["clipPos"])
		softcN.append(counter)
	candTEins = [] 
	for chrom,clip in zip(posChr,softcN):
		top2 = clip.most_common(2)
		if len(top2)<2:
			continue
		if abs(top2[0][0]-top2[1][0]) < 20:
			softClipN.append(top2[0][1]+top2[1][1])
			Len = abs(top2[0][0]-top2[1][0])
			S   = top2[0][0] if top2[0][0] < top2[1][0] else top2[1][0]
			ins = chrom+"-"+str(top2[0][0])
			#使用softlip reads数量评分
			regionI = chrom+":"+str(S+1)+"-"+str(S+1+Len)
			regionL = chrom+":"+str(S-Len)+"-"+str(S)
			regionR = chrom+":"+str(S+1+Len)+"-"+str(S+1+Len+Len)
			depI   = int(subprocess.getoutput("echo $(samtools depth -a -g 256 -r {region} {bam} | cut -f 3 | tr '\n' '+' | sed 's/+$//') | bc".format(region=regionI,bam=bam)))/Len
			depL   = int(subprocess.getoutput("echo $(samtools depth -a -g 256 -r {region} {bam} | cut -f 3 | tr '\n' '+' | sed 's/+$//') | bc".format(region=regionL,bam=bam)))/Len
			depR   = int(subprocess.getoutput("echo $(samtools depth -a -g 256 -r {region} {bam} | cut -f 3 | tr '\n' '+' | sed 's/+$//') | bc".format(region=regionR,bam=bam)))/Len
			#根据softclip的深度比例调整得分
			candTEins.append(ins); regionIs.append(depI); regionLs.append(depL); regionRs.append(depR)
	#7.对候选插入位置的深度进行过滤，如果超过100则过滤掉
	candTEIns = []; softClipN1 = []; regionIs1 = []; regionLs1 = []; regionRs1 = []
	for idx,reg in enumerate(candTEins):
		chrom = "-".join(reg.split("-")[0:2])
		pos   = reg.split("-")[-1]
		region = chrom+":"+pos+"-"+pos
		dep   = subprocess.getoutput("samtools depth -a -g 256 -r {region} {bam}".format(region=region,bam=bam)).split("\t")[-1]
		print(region,dep)
		if int(dep)<100 and int(dep)>10:
			candTEIns.append(candTEins[idx]); softClipN1.append(softClipN[idx]); regionIs1.append(regionIs[idx]); regionLs1.append(regionLs[idx]); regionRs1.append(regionRs[idx])

	df = pd.DataFrame({"pos":candTEIns,"softclipN":softClipN1,"regionIdep":regionIs1,"regionLdep":regionLs1,"regionRdep":regionRs1})
	df.to_csv(sample+"-"+TE+".TEalign.csv")
	#print(*candTEIns,sep="\n",file=open(sample+"-"+TE+".TEalign.txt","w"))
