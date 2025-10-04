#!/home/songlizhi/software/miniconda3/bin/python
import pandas as pd
import glob
import os
import re
import sys
#modify point
igvScript = "TE.P.igv.batch"
sample = sys.argv[1]
bamDir = "/home/songlizhi/sports/swo/align/"
bam    = sample + ".P_TE.bam.realigned.markdup.bam"
disAln = pd.read_csv("/home/songlizhi/sports/swo/align/P_TE_detect/"+sample+".P.disAlign.txt",header=None,names=["id"])
TEAlns = glob.glob("/home/songlizhi/sports/swo/align/P_TE_detect/"+sample+".P-TE*.TEalign.txt")
dfs = []
if not os.path.exists("{sample}.P".format(sample=sample)):
        os.system("mkdir {sample}.P".format(sample=sample))
for TEAln in TEAlns:
        tp = re.findall("TE[0-9]+",TEAln)[0]
        df = pd.read_csv(TEAln,header=None,names=["id"])
        df["id"] = tp+"-"+df["id"]
        dfs.append(df)
dfs.append(disAln)
Adf = pd.concat(dfs,axis=0).sort_values("id").drop_duplicates().reset_index(drop=True)
def splt(x):
        return pd.Series([x.split("-")[0], "-".join(x.split("-")[1:3]), x.split("-")[-1]],index=["type","chrom","pos"])
Odf = Adf["id"].apply(splt).sort_values(["type","chrom","pos"],axis="index")
Odf.to_csv("/home/songlizhi/sports/swo/align/P_TE_detect/{sample}.P.TEins.csv".format(sample=sample))
for index in Adf.index:
        item = Adf.loc[index,"id"]
        chr = "-".join(item.split("-")[1:3])
        pos = item.split("-")[-1]
        start = int(pos) - 300 if int(pos)>300 else 1
        end = int(pos) + 300
        name = item + ".png"
        with open(igvScript, "w") as script:
                print("new", file = script)
                print("genome /home/songlizhi/genome/swoT2T/T2T_TE/swoP_TE.fa", file = script)
                print("load " + bamDir + bam, file = script)
                print("snapshotDirectory {sample}.P".format(sample=sample), file = script)
                print("goto " + chr + ":" + str(start) + "-" + str(end), file = script)
                print("setSleepInterval 800", file = script)
                print("snapshot " + item + ".png", file = script)
                print("exit", file = script)
        #modify point
        os.system("/home/songlizhi/software/IGV_Linux_2.14.1/igv.sh --batch TE.P.igv.batch")