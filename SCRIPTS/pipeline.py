from __future__ import division
import argparse
import os
import sys
import logging

def mkdir(path):
	try:
		os.makedirs(path)
	except:
		logger.info("Directory {} already exists or couldn't create new directory".format(path))

def predict(args):
	cmd=""
	if args.type == "prot":
		args.aligner = "DIAMOND"
		cmd = " ".join(['diamond ', 'blastx',
			'-q', args.input_file,
			'-d', args.data_path+"DATABASE/DIAMOND/VFDB.setA.pro.db.dmnd",
			'-p', args.threads,
			'-k', str(args.arg_num_alignments_per_entry),
			'--id', str(args.arg_alignment_identity),
			'--sensitive',
			'-e', str(args.arg_alignment_evalue),
			'-o', args.output_file + ".sam",
			'-f', '101' 
			])
	elif args.type == "nulc":
		args.aligner = "BOWTIE2"
		cmd = " ".join(['bowtie2 ',
			'-q', args.input_file,
			'-x', args.data_path + "DATABASE/BOWTIE2/VFDB",
			'-p', args.threads,
			'-S', args.output_file + ".sam"
			])
	
	print(cmd)
	os.system(cmd)
	ofile=args.output_file + ".sam"
	
	(align_u, align_m, total) = count_reads(args.output_file + ".sam")
	FormatData(total, align_u, align_m, args.output_file)


def count_reads(samfile):
	Total_reads = 0
	READ1_DICT={}
	READ2_DICT={}
	COUNT1_DICT={}
	COUNT2_DICT={}
	F=open(samfile,"r")
	for line in F:
		if line.startswith("@"):
			continue
		else:
			array=line.strip().split("\t")
			newid=array[0].split("/")[0]
			if array[2] == "*":
				continue
			if array[0].endswith("/1"):
				COUNT1_DICT[newid] = COUNT1_DICT.setdefault(newid, 0) + 1
				tmp = READ1_DICT.setdefault(newid, array[2])
				if tmp != array[2]:
					READ1_DICT[newid] = READ1_DICT[newid] + ":" + array[2]
			if array[0].endswith("/2"):
				COUNT2_DICT[newid] = COUNT2_DICT.setdefault(newid, 0) + 1
				tmp = READ2_DICT.setdefault(newid, array[2])
				if tmp != array[2]:
					READ2_DICT[newid] = READ2_DICT[newid] + ":" + array[2]

	F.close()

	REFCOUNT_UNIQ_DICT = {}
	REFCOUNT_MULT_DICT = {}

	#Uniq
	for newid,value1 in COUNT1_DICT.items():
		value2 = COUNT2_DICT.setdefault(newid, 0)
		if (value1 == 1) & (value2 == 1):
			if (READ1_DICT[newid] == READ2_DICT[newid]):
				Total_reads = Total_reads + 1
				REFCOUNT_UNIQ_DICT[READ1_DICT[newid]]=REFCOUNT_UNIQ_DICT.setdefault(READ1_DICT[newid], 0) + 1
		elif ((value1 > 1) & (value2 > 0)) | ((value1 > 0) & (value2 > 1)):
			gene_from_dict1 = READ1_DICT[newid].split(":")
			gene_from_dict2 = READ2_DICT[newid].split(":")
			count = 0	
			for i in range(len(gene_from_dict1)):
				if gene_from_dict1[i] in gene_from_dict2:
					count = count + 1

			for i in range(len(gene_from_dict1)):
				if gene_from_dict1[i] in gene_from_dict2:
					REFCOUNT_MULT_DICT[gene_from_dict1[i]] = REFCOUNT_MULT_DICT.setdefault(gene_from_dict1[i], 0) +  1 / count

			Total_reads = Total_reads + 1

	return(REFCOUNT_UNIQ_DICT, REFCOUNT_MULT_DICT, Total_reads)

def download_data(args):
	print("XXX")

def FormatData(Total_reads, REFCOUNT_UNIQ_DICT, REFCOUNT_MULT_DICT, OUTFILE):	
	###
	DICT_GENELENGHT={}
	DICT_GENE={}
	DICT_VF={}
	DICT_SUBTYPE={}
	
	dbfile = "/root/SOFTWARE/CNNVF/DATABASE/ANNO/VFDB_setA_pro.length.xls"
	F=open(dbfile,"r")
	for line in F.readlines():
		array=line.strip().split("\t")
		DICT_GENELENGHT.update({array[0] : float(array[1])})

	F.close()

	dbfile = "/root/SOFTWARE/CNNVF/DATABASE/ANNO/VFG2VF.xls"
	F=open(dbfile,"r")
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
	dbfile = "/root/SOFTWARE/CNNVF/DATABASE/ANNO/VF.annotation.xls"
	names=[]
	F=open(dbfile,"r")
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

	dbfile = "/root/SOFTWARE/CNNVF/DATABASE/ANNO/VF2VFCLASS.xls"
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
	for gene,value in DICT_GENE.items():
		FW1.write("%s\t%3f\n" % (gene,value))
	FW1.close()

	#To store the user-friendly result about virulence factors and subtypes
	FILE2=OUTFILE+".VF.XLS"
	FW2=open(FILE2,"w")
	FW2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ("Virulence Factors", "Value", "VF Name", "Subtype", "Associated Bacteria", "Mechanism", "Function"))
	for VF,value in DICT_VF.items():
		FW2.write("%s\t%3f\t%s\t%s\t%s\t%s\t%s\n" % (VF, value, DB_VF2ANNO[VF]['VF_Name'], DB_VF2SUBTYPE[VF], DB_VF2ANNO[VF]["Bacteria"], DB_VF2ANNO[VF]["Mechanism"], DB_VF2ANNO[VF]["Function"]))
	FW2.close()
	
	#To save the result from 
	FILE3=OUTFILE+".SUBTYPE.XLS"
	FW3=open(FILE3,"w")
	FW3.write("%s\t%s\n" % ("Subtype", "Value"))
	for VF, value in DB_SUBTYPE2VF.items():
		FW3.write("%s\t%3f\n" % (VF, value))
	FW3.close()

def TwoDimDict(thedict, key_a, key_b, val):
	if key_a in thedict.keys():
		thedict[key_a].update({key_b: val})
	else:
		thedict.update({key_a:{key_b: val}})

def Normalization(DICT_GENE, DICT_GENE_LENGTH, total):
	for k in DICT_GENE.keys():
		DICT_GENE[k] = DICT_GENE[k] * 1000/ (DICT_GENE_LENGTH[k] * total) 
		
	return DICT_GENE

def plot_data(args):
	print("XXX")

def main():
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers()

    	# use deeparg section
	reads = subparsers.add_parser("PREDICT", help="Predict virulence factors from next generation sequencing reads")
	reads.add_argument('-i', '--input-file', required=False,
                       help='Input file (Fastq input file)')
	reads.add_argument('-o', '--output-file', required=True,
                       help='Output file where to store results')
	reads.add_argument('-d', '--data-path', required=False,
                       help="Path where data was downloaded [see deeparg download-data --help for details]")
	reads.add_argument('--type', default='nucl',
                       help='Molecular data type prot/nucl [Default: nucl]')
	reads.add_argument('--method', default='VFdetector',
                       help='Method adopted to detect virulent factors [VFdetector, BestHit]')
	reads.add_argument('-t','--threads',required=False,
                       help='number of CPUs to allocate [Default: 1]')
	reads.add_argument('--min-prob', default=0.8,
                       help='Minimum probability cutoff [Default: 0.8]')
	reads.add_argument('--arg-alignment-identity', default=50,
                       help='Identity cutoff for sequence alignment [Default: 50]')
	reads.add_argument('--arg-alignment-evalue', default=1e-10,
                       help='Evalue cutoff [Default: 1e-10]')
	reads.add_argument('--arg-alignment-overlap', default=0.8,
                       help='Alignment overlap cutoff [Default: 1e-10]')
	reads.add_argument('--arg-num-alignments-per-entry', default=1000,
                       help='Diamond, minimum number of alignments per entry [Default: 1000]')
	
	reads.set_defaults(func=predict)


   	# Download section
	download = subparsers.add_parser("DOWNLOAD", 
		       help="Download the databases and models used in CNNVF")
	download.add_argument('-o','--output_path',required=False, 
		       help='Output Directory to store database')
	download.set_defaults(func=download_data)
	
	# Draw
	draw = subparsers.add_parser("DRAW", 
		       help="Draw the coverage information for each virulence factors")
	draw.add_argument("--TSV","-T", 
		       help = "The aligned file from DIAMOND")
	draw.add_argument("--BAM","-B", 
		       help = "The aligned file from BOWTIE2")
	draw.add_argument("--GENE","-G", 
		       help = "The Gene ID to extract and draw")
	draw.add_argument("--PLOT","-P",
		       help = "The file to draw")
	draw.set_defaults(func=plot_data)
	
	args = parser.parse_args()
	parser.parse_args()
	args.func(args)

	pass

if __name__ == '__main__':
    main()
