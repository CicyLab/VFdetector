from __future__ import division
import argparse
import os
import sys
import logging

def align(args):
	cmd=""
	if args.type == "prot":
		args.aligner = "DIAMOND"
		cmd = " ".join(['diamond ', 'blastx',
			'-q', args.input_file,
			'-d', args.data_path+"/PROTEIN/VFDB.setA.pro.db.dmnd",
			'-p', args.threads,
			'-k', str(args.arg_num_alignments_per_entry),
			'--id', str(args.arg_alignment_identity),
			'--sensitive',
			'-e', str(args.arg_alignment_evalue),
			'-o', args.output_file + ".sam",
			'-f', '101' 
			])
	elif args.type == "nucl":
		args.aligner = "BOWTIE2"
		cmd = " ".join(['bowtie2 ',
			'-q', args.input_file,
			'-x', args.data_path + "/DNA/VFDB",
			'-p', args.threads,
			'-S', args.output_file + ".sam",
			'--quiet',
			'--no-unal',
		#	'| /mnt/osf1/soft/samtools-1.10/samtools view -b -o', args.output_file+".bam -"	
			])
	return(cmd)

def align_unique(args):
	cmd=""
	if args.type == "prot":
		args.aligner = "DIAMOND"
		cmd = " ".join(['diamond ', 'blastx',
			'-q', args.input_file,
			'-d', args.data_path+"/PROTEIN/VFDB.setA.pro.db.dmnd",
			'-p', args.threads,
			'-k', str(1),
			'--id', str(args.arg_alignment_identity),
			'--sensitive',
			'-e', str(args.arg_alignment_evalue),
			'-o', args.output_file + ".sam",
			'-f', '101'
			])
	elif args.type == "nucl":
		args.aligner = "BOWTIE2"
		cmd = " ".join(['bowtie2 ',
			'-q', args.input_file,
			'-x', args.data_path + "/DNA/VFDB",
			'-p', args.threads,
			'-S', args.output_file + ".sam"
			])

	return(cmd)

def process_alignment(args, samfile, predfile, tsvfile):
	FI=open(samfile,"r")
	FP=open(predfile,"r")
	FO=open(tsvfile,"w")

	FO.write("\t".join(["ReadId","Alignment", "VirulenceType","Probability","PredictedType","BestHit","DiamondIdentity","DiamondLength","DiamondBitscore","DiamondEvalue","Bowtie2MappingQual","Bowtie2MappingLength"]))
	FO.write("\n")	
	dict_prob = {}	
	dict_pred = {}
	for line in FP.readlines():
		array= line.strip().split("\t")
		dict_pred.update({array[0]:array[2]})
		dict_prob.update({array[0]:array[1]})
	FP.close()
	
	if args.type == "prot":
		for line in FI.readlines():
			if not line.startswith("@"):
				array= line.strip().split("\t")
				if not array[2] == "*":
					ReadId= array[0]
					Alignment = array[2]
					VirulenceType = "unknown"
					Probability   = dict_prob[array[0]]
					PredictedType = dict_pred[array[0]]
					bitscore=array[11].split(":")
					DiamondBitscore=bitscore[2]
					BestHit= array[2]
					evalue=array[15].split(":")
					DiamondEvalue=evalue[2]
					id=array[16].split(":")
					DiamondIdentity=id[2]
					len=array[19].split(":")
					DiamondLength=len[2]
					Bowtie2MappingQual="unknown"
					Bowtie2MappingLength="unknown"
					a="\t"
					FO.write(a.join([ReadId, Alignment, VirulenceType, Probability, PredictedType, BestHit, DiamondIdentity, DiamondLength, DiamondBitscore, DiamondEvalue, Bowtie2MappingQual, Bowtie2MappingLength]))
					FO.write("\n")
	elif args.type == "nucl":
		for line in FI.readlines():
			if not line.startswith("@"):
				array= line.strip().split("\t")
				if not array[2] == "*":
					ReadId= array[0]
					Alignment = array[2]
					VirulenceType = "unknown"
					Probability   = dict_prob[array[0]]
					PredictedType = dict_pred[array[0]]
					DiamondBitscore="unknown"
					BestHit= array[2]
					DiamondEvalue="unknown"
					DiamondIdentity="unknown"
					DiamondLength="unknown"
					Bowtie2MappingQual= array[4]
					Bowtie2MappingLength= array[5]
					a="\t"
					FO.write(a.join([ReadId, Alignment, VirulenceType, Probability, PredictedType, BestHit, DiamondIdentity, DiamondLength, DiamondBitscore, DiamondEvalue, Bowtie2MappingQual, Bowtie2MappingLength]))
					FO.write("\n")
	
	FI.close()
	FO.close()

def process_alignment_unique(samfile, tsvfile):
	FI=open(samfile,"r")
	FO=open(tsvfile,"w")
	FO.write("\t".join(["ReadId","Alignment", "VirulenceType","Probability","PredictedType","BestHit","DiamondIdentity","DiamondLength","DiamondBitscore","DiamondEvalue","Bowtie2MappingQual","Bowtie2MappingLength"]))
	FO.write("\n")
	for line in FI.readlines():
		if not line.startswith("@"):
			array= line.strip().split("\t")
			if not array[2] == "*":
				ReadId= array[0]
				Alignment = array[2]
				VirulenceType = "unknown"
				Probability   = "unknown"
				PredictedType = "unknown"
				bitscore=array[11].split(":")
				DiamondBitscore=bitscore[2]
				BestHit= array[2]
				evalue=array[15].split(":")
				DiamondEvalue=evalue[2]
				id=array[16].split(":")
				DiamondIdentity=id[2]
				len=array[19].split(":")
				DiamondLength=len[2]
				Bowtie2MappingQual="unknown"
				Bowtie2MappingLength="unknown"
				a="\t"
				FO.write(a.join([ReadId, Alignment, VirulenceType, Probability, PredictedType, BestHit, DiamondIdentity, DiamondLength, DiamondBitscore, DiamondEvalue, Bowtie2MappingQual, Bowtie2MappingLength]))
				FO.write("\n")
	FI.close()
	FO.close()
