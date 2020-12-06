from __future__ import division

import argparse
import csv
import logging
import math
import numpy as np
import os
import random
import sys
import shutil 
import time
from SCRIPTS import pipeline

def CreateResultDirectory(path):
	""" 
	Return the location of the directory for saving results
	"""
	if not os.path.isdir(path):
		try:
			print("Creating output directory: " + path)
			os.mkdir(path)
		except:
			logging.info("Directory {} already exists or couldn't create new directory".format(path))

def check_requirements():
	###Checking the directory to write 
	if not os.path.isdir(path):
		try:
			print("Creating output directory: " + path)
			os.mkdir(path)
		except:
			logging.info("Directory {} already exists or couldn't create new directory".format(path))

def main():
	parser = argparse.ArgumentParser()
	subparsers = parser.add_subparsers()

    	# use machine learning method
	reads = subparsers.add_parser("PREDICT", help="Predict virulence factors from next generation sequencing reads")
	reads.add_argument('-i', '--input-file', required=True,
                       help='Input file (Fastq input file)')
	reads.add_argument('-o', '--output-file', required=True,
                       help='Output file where to store results')
	reads.add_argument('-d', '--data-path', required=True,
                       help="Directory to store database")
	reads.add_argument('--type', default='nucl',
                       help='Molecular data type prot/nucl [Default: nucl]')
	reads.add_argument('-t','--threads',required=False,
                       help='number of CPUs to allocate [Default: 1]')
	reads.add_argument('--min-prob', default=0.5,
                       help='Minimum probability cutoff [Default: 0.5]')
	reads.add_argument('--arg-alignment-identity', default=30,
                       help='Identity cutoff for sequence alignment [Default: 30]')
	reads.add_argument('--arg-alignment-evalue', default=1e-3,
                       help='Evalue cutoff [Default: 1e-3]')
	reads.add_argument('--arg-alignment-overlap', default=0.8,
                       help='Alignment overlap cutoff [Default: 0.8]')
	reads.add_argument('--arg-num-alignments-per-entry', default=1000,
                       help='Diamond, minimum number of alignments per entry [Default: 1000]')
	
	reads.set_defaults(func=pipeline.predict)
	
	# use best-hit method
	reads = subparsers.add_parser("BESTHIT", help="Predict virulence factors from next generation sequencing reads")
	reads.add_argument('-i', '--input-file', required=True,
			help='Input file (Fastq input file)')
	reads.add_argument('-o', '--output-file', required=True,
			help='Output file where to store results')
	reads.add_argument('-d', '--data-path', required=True,
			help="Directory to store database")
	reads.add_argument('--type', default='nucl',
			help='Molecular data type prot/nucl [Default: nucl]')
	reads.add_argument('-t','--threads',required=False,
			help='number of CPUs to allocate [Default: 1]')
	reads.add_argument('--arg-alignment-identity', default=20,
			help='Identity cutoff for sequence alignment [Default: 20]')
	reads.add_argument('--arg-alignment-evalue', default=1e-3,
			help='Evalue cutoff [Default: 1e-3]')
	reads.add_argument('--arg-alignment-overlap', default=0.8,
			help='Alignment overlap cutoff [Default: 0.8]')
	reads.add_argument('--arg-num-alignments-per-entry', default=1000,
			help='Diamond, minimum number of alignments per entry [Default: 1000]')
	
	reads.set_defaults(func=pipeline.best_hit)
	
	# Draw
	draw = subparsers.add_parser("DRAW", 
		       help="Draw the coverage information for each virulence gene")
	draw.add_argument("-g","--gid", required=True,
		       help = "The gene id in database")
	draw.add_argument("-s","--sam", required=True,
		       help = "The alignment file from BOWTIE2")
	draw.add_argument("-o","--outfile", required=True,
		       help = "The file to draw")
	draw.add_argument("-d","--database", required=True,
		       help = "The gene length database")
	draw.set_defaults(func=pipeline.draw)
	
	args = parser.parse_args()
	parser.parse_args()
	args.func(args)

	pass

if __name__ == '__main__':
    main()
