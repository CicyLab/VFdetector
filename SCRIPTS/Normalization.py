import argparse
import os
import sys
import logging

def FRPM(DICT_GENE, DICT_GENE_LENGTH, total):
	for k in DICT_GENE.keys():
		DICT_GENE[k] = DICT_GENE[k] * 1000/ (DICT_GENE_LENGTH[k] * total) 
		
	return DICT_GENE
