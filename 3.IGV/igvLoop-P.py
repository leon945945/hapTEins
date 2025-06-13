#!/bin/env python3
import pandas as pd
import openpyxl
import os
import re
#modify point
igvScript = "sv.P.igv.batch"
sheets = openpyxl.load_workbook("/home02/songlizhi/sports/hifi/calling/sv/SWO.P.sv.noXG.filter.csq.xlsx").sheetnames
sheets = sheets[0:len(sheets)-2]
sheets = ["HQ_TL"]

bamDir = "/home02/songlizhi/sports/hifi/align/"
#files  = os.listdir("/home02/songlizhi/sports/hifi/calling/sv/igv/sv.P")
for sheet in sheets:
	sv  = pd.read_excel("/home02/songlizhi/sports/hifi/calling/sv/SWO.P.sv.noXG.filter.csq.xlsx", sheet_name = sheet)
	for index in sv.index:
		chr = sv.loc[index,"chr"]
		start = int(sv.loc[index,"pos"]) - 200 if int(sv.loc[index,"pos"])>200 else 1
		end = int(sv.loc[index,"pos"]) + abs(sv.loc[index,"svlen"]) + 200 if sv.loc[index,"svtype"]=="DEL" else int(sv.loc[index,"pos"]) + 200
		name = sheet + "-" + chr + "-" + str(start) + "-" + str(end) + ".png"
		#if not chr=="V-P9":
			#continue

		with open(igvScript, "w") as script:
			print("new", file = script)
			print("genome /home/songlizhi/genome/swoT2T/V_T2T-P.fa", file = script)
			Bams = os.listdir(bamDir)
			Bams.sort()
			for bam in Bams:
				#modify point
				if re.findall(".*P.reads.bam$", bam):
					print("load " + bamDir + bam, file = script)
			#modify point
			print("snapshotDirectory sv.P", file = script)
			print("goto " + chr + ":" + str(start) + "-" + str(end), file = script)
			print("setSleepInterval 800", file = script)
			print("snapshot " + sheet + "-" + chr + "-" + str(start) + "-" + str(end) + ".png", file = script)
			print("exit", file = script)

		#modify point
		os.system("/home/songlizhi/software/IGV_Linux_2.14.1/igv.sh --batch sv.P.igv.batch")
