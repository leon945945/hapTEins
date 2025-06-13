#!/home/songlizhi/software/miniconda3/bin/python
import pandas as pd
import sys
Mf = sys.argv[1]
Pf = sys.argv[2]
sample = sys.argv[3]
Mout = open(sample+".M.reads","w")
Pout = open(sample+".P.reads","w")
stats = open(sample+".reads.stats","w")
Mreads = pd.read_csv(Mf,sep="\t",header=0)
Preads = pd.read_csv(Pf,sep="\t",header=0)
Mreads_pri = Mreads.loc[Mreads["type"]=="pri",["readName","AS","NM"]]
Preads_pri = Preads.loc[Preads["type"]=="pri",["readName","AS","NM"]]
Mreads_priSum = Mreads_pri.groupby("readName")["AS","NM"].agg("sum")
Preads_priSum = Preads_pri.groupby("readName")["AS","NM"].agg("sum")
AllReads = pd.merge(Mreads_priSum,Preads_priSum,left_index=True,right_index=True,how="outer",suffixes=["_M","_P"])
AllReads = AllReads.fillna(0)
AS_M = AllReads.loc[AllReads["AS_M"]>AllReads["AS_P"],:]
AS_P = AllReads.loc[AllReads["AS_P"]>AllReads["AS_M"],:]
NM_M = AllReads.loc[(AllReads["AS_M"]==AllReads["AS_P"]) & (AllReads["NM_M"]<AllReads["NM_P"]),:]
NM_P = AllReads.loc[(AllReads["AS_M"]==AllReads["AS_P"]) & (AllReads["NM_P"]<AllReads["NM_M"]),:]
MP   = AllReads.loc[(AllReads["AS_M"]==AllReads["AS_P"]) & (AllReads["NM_P"]==AllReads["NM_M"]),:]
M_readsN = AS_M.shape[0] + NM_M.shape[0]
P_readsN = AS_P.shape[0] + NM_P.shape[0]
MP_readsN = MP.shape[0]
M_reads = list(AS_M.index) + list(NM_M.index) + list(MP.index)
P_reads = list(AS_P.index) + list(NM_P.index) + list(MP.index)
print("No.M reads: {M}\nNo.P reads: {P}\nNo.MP reads: {MP}".format(M=M_readsN,P=P_readsN,MP=MP_readsN),file=stats)
print(*M_reads,sep="\n",file=Mout)
print(*P_reads,sep="\n",file=Pout)
