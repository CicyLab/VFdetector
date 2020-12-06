from __future__ import division
import argparse
import os
import sys
import logging
from . import Alignment
from . import Utils

def FormatData(args, Total_reads, REFCOUNT_UNIQ_DICT, REFCOUNT_MULT_DICT, OUTFILE):	
	###
	DICT_GENELENGHT={}
	DICT_GENE={}
	DICT_VF={}
	DICT_SUBTYPE={}
	
	if args.type == "prot":
		glenfile = args.data_path+"/ANNO/VFDB_setA_pro.length.xls"
	elif args.type == "nucl":
		glenfile = args.data_path+"/ANNO/VFDB_setA_nt.length.xls"
	
	F=open(glenfile,"r")
	for line in F.readlines():
		array=line.strip().split("\t")
		DICT_GENELENGHT.update({array[0] : float(array[1])})

	F.close()

	vfgfile = args.data_path+"/ANNO/VFG2VF.xls"
	F=open(vfgfile,"r")
	colnames=[]
	for line in F.readlines():
		if line.startswith("VFGID"):
			colnames=line.strip().split("\t")
		else:
			line=line.strip().split("\t")
			DICT_GENE.update({line[0]:0})
			if line[0] in REFCOUNT_UNIQ_DICT.keys():
				GeneLength = DICT_GENELENGHT[line[0]]
				DICT_GENE[line[0]] = DICT_GENE[line[0]] + REFCOUNT_UNIQ_DICT[line[0]] * 1000000 /(GeneLength * Total_reads)
			
			if line[0] in REFCOUNT_MULT_DICT.keys():
				GeneLength = DICT_GENELENGHT[line[0]]
				DICT_GENE[line[0]] = DICT_GENE[line[0]] + REFCOUNT_MULT_DICT[line[0]] * 1000000 /(GeneLength * Total_reads) 

			if line[1] not in DICT_VF.keys():
				DICT_VF.update({line[1]:0})
			else:
				DICT_VF[line[1]] = DICT_VF[line[1]] + DICT_GENE[line[0]]
	F.close()
	

	###DATABASE
	DB_VF2ANNO = {}
	vffile = args.data_path+"/ANNO/VF.annotation.xls"
	names=[]
	F=open(vffile,"r")
	for line in F.readlines():
		if line.startswith("VFID"):
        		names=line.strip().split("\t")
		else:
			line=line.strip().split("\t")
#			print(line[0])
			DB_VF2ANNO.update({line[0]:{"VF_Name":line[1]}})
			DB_VF2ANNO[line[0]].update({"Bacteria":line[2]})
			DB_VF2ANNO[line[0]].update({"Characteristics":line[3]})
			DB_VF2ANNO[line[0]].update({"Function":line[5]})
			DB_VF2ANNO[line[0]].update({"Mechanism":line[6]})
	F.close()
	
	DB_VF2SUBTYPE = {}
	DB_SUBTYPE2TYPE = {}
	DB_SUBTYPE2VF = {}

	dbfile = args.data_path+"/ANNO/VF2VFCLASS.xls"
	names=[]
	with open(dbfile) as f :
		for line in f :
			if line.startswith("types"):
				names=line.strip().split("\t")
			else:
				line=line.strip().split("\t")
				tmp = line[2]
				DB_SUBTYPE2VF.update({line[1]:0})
				for i in range(len(line)):
					if i > 2:
						tmp = tmp + line[i]
						if line[i] not in DICT_VF.keys():
							DICT_VF.update({line[i]:0})
						DB_SUBTYPE2VF[line[1]]=DB_SUBTYPE2VF[line[1]]+DICT_VF[line[i]]

					if line[i] in DB_VF2SUBTYPE.keys():
						DB_VF2SUBTYPE.update({line[i]:DB_VF2SUBTYPE[line[i]]+"|"+line[1]})
					else:
						DB_VF2SUBTYPE.update({line[i]:line[1]})
				DB_SUBTYPE2TYPE.update({line[1]:line[0]})

	f.close()

	#To store the raw information about genes
	FILE1=OUTFILE+".VFG.XLS"
	FW1=open(FILE1,"w")
	FW1.write("%s\t%s\t%s\t%s\t%s\n" % ("GID","GeneLength","UniqCount","MultiCount","Abundance"))
	for gene,value in DICT_GENE.items():
		GeneLength=DICT_GENELENGHT.setdefault(gene, 0)
		UniqCount=REFCOUNT_UNIQ_DICT.setdefault(gene, 0)
		MultiCount=REFCOUNT_MULT_DICT.setdefault(gene, 0)
		FW1.write("%s\t%3f\t%3f\t%3f\t%3f\n" % (gene,GeneLength,UniqCount,MultiCount,value))
	FW1.close()

	#To store the user-friendly result about virulence factors and subtypes
	FILE2=OUTFILE+".VF.XLS"
	FW2=open(FILE2,"w")
	FW2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("Virulence Factors", "Value", "VF Name", "Subtype", "Associated Bacteria", "Mechanism", "Function"))
	for VF,value in DICT_VF.items():
		FW2.write("%s\t%3f\t%s\t%s\t%s\t%s\t%s\n" % (VF, value, DB_VF2ANNO[VF]['VF_Name'], DB_VF2SUBTYPE[VF], DB_VF2ANNO[VF]["Bacteria"], DB_VF2ANNO[VF]["Mechanism"], DB_VF2ANNO[VF]["Function"]))
	FW2.close()
	
