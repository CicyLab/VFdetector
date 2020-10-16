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

	try:
		os.makedirs(path)
	except:
		logger.info("Directory {} already exists or couldn't create new directory".format(path))

def check_requirements():
	


def download_data(args):
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
	
	reads.set_defaults(func=pipeline.predict)


   	# Download section
	download = subparsers.add_parser("DOWNLOAD", 
		       help="Download the databases and models used in CNNVF")
	download.add_argument('-o','--output_path',required=False, 
		       help='Output Directory to store database')
	download.set_defaults(func=pipeline.download_data)
	
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
	draw.set_defaults(func=pipeline.plot_data)
	
	args = parser.parse_args()
	parser.parse_args()
	args.func(args)

	pass

if __name__ == '__main__':
    main()
