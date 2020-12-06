import os
from . import Utils

def MakeXY(args, samfile, tsvfile):
	
	if args.type == "prot":
		dbfile = args.data_path+"/PROTEIN/VFDB_setA_pro.fas.short.filt.fa"
	elif args.type == "nucl":
		dbfile = args.data_path+"/DNA/VFDB_setA_nt.fas.short.filt.fa"
		
	DICT_refgene = {}
	FI=open(dbfile, "r")
	for line in FI.readlines():
		if line.startswith(">"):
			vfg=line.strip().replace(">","")
			DICT_refgene.update({vfg: 1})
	FI.close()

	DICT_score = {}
	FI=open(samfile,"r")
	for line in FI.readlines():
		if line.startswith("@"):
			continue
		else:
			array=line.strip().split("\t")
			if(array[2] != "*"):
				read_id = array[0]
				refgene = array[2]
				tmp   = array[11].split(":")
				score = tmp[2]
				Utils.TwoDimDict(DICT_score, read_id, refgene, score)

	FI.close()

	FO=open(tsvfile,"w")
	FO.write("ID")
	for refgene in DICT_refgene.keys():
		FO.write("\t")
		FO.write(refgene)
	FO.write("\n")

	for read_id, value in DICT_score.items():
		FO.write(read_id)
		for refgene in DICT_refgene.keys():
			if refgene in value.keys():
				score=DICT_score[read_id][refgene]
			else:
				score=0
			FO.write("\t")
			FO.write(str(score))
	
		FO.write("\n")

	FO.close()
