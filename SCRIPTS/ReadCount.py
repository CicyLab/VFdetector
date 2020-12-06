from __future__ import division
import argparse
import os
import sys
import logging
from . import Utils
from . import Store

def count_reads(args, tabfile):
	Total_reads = 0
	READ1_DICT={}
	READ2_DICT={}
	COUNT1_DICT={}
	COUNT2_DICT={}

	F=open(tabfile,"r")
	for line in F:
		if line.startswith("@"):
			continue
		else:
			array=line.strip().split("\t")
			newid=array[0].split("/")[0]
			if array[1] == "*":
				continue
			if array[0].endswith("/1"):
				COUNT1_DICT[newid] = COUNT1_DICT.setdefault(newid, 0) + 1
				tmp = READ1_DICT.setdefault(newid, array[1])
				if tmp != array[1]:
					READ1_DICT[newid] = READ1_DICT[newid] + ":" + array[1]
			if array[0].endswith("/2"):
				COUNT2_DICT[newid] = COUNT2_DICT.setdefault(newid, 0) + 1
				tmp = READ2_DICT.setdefault(newid, array[1])
				if tmp != array[1]:
					READ2_DICT[newid] = READ2_DICT[newid] + ":" + array[1]

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
					ucount = 1
					count = count + ucount

			for i in range(len(gene_from_dict1)):
				if gene_from_dict1[i] in gene_from_dict2:
					ucount = 1
					REFCOUNT_MULT_DICT[gene_from_dict1[i]] = REFCOUNT_MULT_DICT.setdefault(gene_from_dict1[i], 0) +  ucount / count

			Total_reads = Total_reads + 1

	return(REFCOUNT_UNIQ_DICT, REFCOUNT_MULT_DICT, Total_reads)

