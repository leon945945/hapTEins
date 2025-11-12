#!/home/songlizhi/software/miniconda3/bin/python
from collections import Counter
import pandas as pd
import subprocess
import os
import sys
SAM = sys.argv[1]
bam = sys.argv[2]
Dep = sys.argv[3]
fq1 = sys.argv[4]
fq2 = sys.argv[5]
sample = sys.argv[6]
#1.读入sam文件，和每个位点的深度
print("2.1 read in sam file and get the depth of each position")
sam = pd.read_csv(SAM,sep="\t",names=["qName","FLAG","RNAME","POS","MAPQ","CIGAR","MRNM","MPOS","ISIZE"])
sam["ID"] = sam["RNAME"]+"-"+sam["POS"].astype(str)
DEP = pd.read_csv(Dep,sep="\t",names=["chrom","pos","depth"],on_bad_lines="warn")
DEP["ID"] = DEP["chrom"]+"-"+DEP["pos"].astype(str)
dep = pd.merge(sam,DEP,on="ID",how="left")
del DEP
#2.对深度进行过滤，保留深度10-50之间的位点
print("2.2 filter depth, keep position whose depth is 10-50")
dep = dep.loc[(dep["depth"]>10) & (dep["depth"]<50),:]
sam = sam.loc[dep.index,:]
#3.以1Kb对插入区间分区
print("2.3 binning insertion region with 1Kb window")
beg = dep.index[0]
binL = []
begP = dep.loc[beg,"pos"]
for idx in dep.index:
	if abs(dep.loc[idx,"pos"] - begP) > 1000:
		#print(beg,idx,begP,dep.loc[idx,"pos"])
		binL.append(dep.loc[beg:idx,])
		beg = idx
		begP = dep.loc[idx,"pos"]	
	else:
		continue
binL.append(dep.loc[beg:,])
#4.对每个插入区间计算reads的平均插入位置
print("2.4 get average insertion position of each insertion region")
avePos = []
for df in binL:
	pos = df["pos"].apply("mean")
	avePos.append(int(pos))
posChr = [i.reset_index().loc[0,"chrom"] for i in binL]
#5.在插入平均位置提取soft clip位置
print("2.5 get soft clip position of each average position")
softC = []
for chrom,pos in zip(posChr,avePos):
	TEregion = chrom+":"+str(pos-500)+"-"+str(pos+500) if pos > 500 else chrom+":"+"1"+"-"+str(pos+500)
	#print(bam,TEregion)
	os.system("/home/songlizhi/learning/TEdev/findSplit {bam} {TEregion} > {sample}.softClip".format(bam=bam,TEregion=TEregion,sample=sample))
	softC.append(pd.read_csv("{sample}.softClip".format(sample=sample),sep="\t",names=["reads","clipPos"]))
#6.对reads的soft clip位置统计数量，数量最多的前两个的距离小于20bp，认为是TE插入的候选区间
print("2.6 stats the num of soft clip, the distance of top2 soft clip position less than 20bp was considered as candidate region")
softcN = []
for i in softC:
	counter = Counter(i["clipPos"])
	softcN.append(counter)
candTEins = []; softClipN = []; regionIs = []; regionLs = []; regionRs = []
for chrom,clip in zip(posChr,softcN):
	top2 = clip.most_common(2)
	if len(top2)<2:
		continue
	if abs(top2[0][0]-top2[1][0]) < 20 and top2[0][1]+top2[1][1] > 10:
		ins = chrom+"-"+str(top2[0][0])
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
		candTEins.append(ins); regionIs.append(depI); regionLs.append(depL); regionRs.append(depR)
#7.再对候选位点的深度进行过滤
print("2.7 filter the depth of candidate position")
DEP = pd.read_csv(Dep,sep="\t",names=["chrom","pos","depth"],on_bad_lines="warn")
DEP["ID"] = DEP["chrom"]+"-"+DEP["pos"].astype(str)
candTEins_d = {"ID":candTEins,"softclipN":softClipN,"regionIdep":regionIs,"regionLdep":regionLs,"regionRdep":regionRs}
candTEins_df = pd.DataFrame(candTEins_d)
candTEins_df = pd.merge(candTEins_df,DEP,on="ID",how="left")
candTEins_df = candTEins_df.loc[(candTEins_df["depth"]>10) & (candTEins_df["depth"]<50),:]
candTEIns = list(candTEins_df["ID"]); softClipN1 = list(candTEins_df["softclipN"]); regionIs1 = list(candTEins_df["regionIdep"]); regionLs1 = list(candTEins_df["regionLdep"]); regionRs1 = list(candTEins_df["regionRdep"])
del DEP
#8.再对候选位点左右100bp是否有比对到其他染色体的reads进行过滤，少于5条比对到其他染色体的reads则过滤掉
print("2.8 detect whether existing reads mapped to other chromosome among flanking 100bp of insertion position")
allRT = []; candTEins = []; softClipN2 = []; regionIs2 = []; regionLs2 = []; regionRs2 = []
for idx,reg in enumerate(candTEIns):
	print(reg)
	chrom = "-".join(reg.split("-")[0:2])
	pos   = int(reg.split("-")[-1])
	region = chrom+":"+str(pos-50)+"-"+str(pos+50) if pos > 50 else chrom+":"+"1"+"-"+str(pos+50)
	os.system("/home/songlizhi/learning/TEdev/getReadsMateTag {bam} {region} > {sample}.readsTag".format(bam=bam,region=region,sample=sample))
	try:
		RT = pd.read_csv("{sample}.readsTag".format(sample=sample),sep="\t",header=0)
	except:
		print("warning: empty dataframe!")
		continue
	RT["ID"] = reg
	notChrN = sum(RT["Mchr"]!=chrom)
	if notChrN>5:
		allRT.append(RT)
		candTEins.append(reg)
		softClipN2.append(softClipN1[idx]); regionIs2.append(regionIs1[idx]); regionLs2.append(regionLs1[idx]); regionRs2.append(regionRs1[idx])
allRT_df = pd.concat(allRT,axis=0)
#9.将候选位点左右50bp区间内的reads比对到TE，过滤掉比对不上的位点
print("2.9 map reads of flanking 50bp of candidate position, filter out position that not mapped to TE")
print(*list(allRT_df["readName"]),sep="\n",file=open("{sample}.all.reads".format(sample=sample),"w"))
os.system("seqkit grep -f {sample}.all.reads {fq1} > {sample}.fq1".format(fq1=fq1,sample=sample))
os.system("seqkit grep -f {sample}.all.reads {fq2} > {sample}.fq2".format(fq2=fq2,sample=sample))
TEtype = []; candTE = []; softClipN3 = []; regionIs3 = []; regionLs3 = []; regionRs3 = []
for reg,dat in allRT_df.groupby("ID"):
	print(reg)
	print(*list(dat["readName"]),sep="\n",file=open("{sample}.reads".format(sample=sample),"w"))
	os.system("seqkit grep -f {sample}.reads {sample}.fq1 > tmp.{sample}.fq1".format(sample=sample))
	os.system("seqkit grep -f {sample}.reads {sample}.fq2 > tmp.{sample}.fq2".format(sample=sample))
	os.system("bwa mem /home/songlizhi/genome/swoT2T/T2T_TE/TEs.fa tmp.{sample}.fq1 tmp.{sample}.fq2 | grep -v @ | cut -f 1-9 > {sample}.sam".format(sample=sample))
	df = pd.read_csv("{sample}.sam".format(sample=sample),sep="\t",names=["qName","FLAG","RNAME","POS","MAPQ","CIGAR","MRNM","MPOS","ISIZE"])
	print(df)
	if sum(df["RNAME"].str.match("TE"))==0:
		continue
		#candTEins.pop(candTEins.index(reg))
	else:
		TEtype.append(Counter([i for i in df["RNAME"] if i!="*"]))
		candTE.append(reg)
		softClipN3.append(softClipN2[candTEins.index(reg)]); regionIs3.append(regionIs2[candTEins.index(reg)]); regionLs3.append(regionLs2[candTEins.index(reg)])
		regionRs3.append(regionRs2[candTEins.index(reg)])
TEpos = []
for t,p in zip(TEtype,candTE):
	tmp = t.most_common(1)
	print(tmp)
	TEpos.append(tmp[0][0]+"-"+p)

df = pd.DataFrame({"pos":TEpos,"softclipN":softClipN3,"regionIdep":regionIs3,"regionLdep":regionLs3,"regionRdep":regionRs3})
#print(*TEpos,sep="\n",file=open(sample+".disAlign.txt","w"))
df.to_csv(sample+".disAlign.csv")
os.system("rm {sample}.fq1 {sample}.fq2 {sample}.reads {sample}.readsTag tmp.{sample}.fq1 tmp.{sample}.fq2".format(sample=sample))
