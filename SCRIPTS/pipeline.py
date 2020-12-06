from __future__ import division
import argparse
import os
import sys
import logging
from . import Alignment
from . import ReadCount
from . import Store
from . import MakeXY
from . import Model
from . import Draw

def mkdir(path):
	try:
		os.makedirs(path)
	except:
		logging.info("Directory {} already exists or couldn't create new directory".format(path))

def predict(args):
	cmd=Alignment.align(args)
	print(cmd)
	os.system(cmd)
	ofile=args.output_file + ".sam"
	tmpfile=args.output_file + ".tsv"
	tsvfile=args.output_file + ".data.tsv"
	predfile=args.output_file + ".pred.tsv"
	tabfile=args.output_file + ".tabular.tsv"
	outfile=args.output_file + ".final"

	MakeXY.MakeXY(args, ofile, tmpfile)
	Model.ModelRidge(args, tmpfile, predfile)
	Alignment.process_alignment(args, ofile, predfile, tabfile)
	(align_u, align_m, total) = ReadCount.count_reads(args, tabfile)
	Store.FormatData(args, total, align_u, align_m, args.output_file)

def best_hit(args):
	cmd=Alignment.align_unique(args)
	print(cmd)
	os.system(cmd)
	ofile=args.output_file + ".sam"
	tmpfile=args.output_file + ".tsv"
	tsvfile=args.output_file + ".data.tsv"
	predfile=args.output_file + ".pred.tsv"
	tabfile=args.output_file + ".tabular.tsv"
	outfile=args.output_file + ".final"
	Alignment.process_alignment_unique(ofile, tabfile)
	(align_u, align_m, total) = ReadCount.count_reads(args, tabfile)
	Store.FormatData(args, total, align_u, align_m, args.output_file)
	
	
def draw(args):
	cmd=Draw.genecov(args)
	
