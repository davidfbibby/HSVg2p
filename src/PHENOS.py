#!/usr/bin/env python

import argparse
import os
import sys

import pandas as pd

from utils import g2pU, gU, pU
from data_init.g2pTables import *

################################################################################
def parse_arguments():

	ap = argparse.ArgumentParser()

	ap.add_argument(
		"-p", dest="previous",
		help="Report all phenotypic data from the same sample(s)"
	)
	ap.add_argument(
		"-s", "--statistics", action="store_true",
		help="Give description of databases (overrides other options)"
	)

	ap.add_argument(
		"-f", dest="single_file", help="path/to/pheno/file. OVERRIDES -d"
	)
	ap.add_argument(
		"-d", dest="directory", default=os.getcwd(),
		help="path/to/dir/containing/pheno/file(s) [.]"
	)
	ap.add_argument(
		"--rundate", help="Run date (YYMMDD, overrides automatic parsing)"
	)
	ap.add_argument(
		"--force", action="store_true",
		help="Force overwrite of existing records with identical data"
	)
	ap.add_argument(
		"--recursive", action="store_true",
		help="Set to look in all subdirectories of <directory>"
	)

	args = ap.parse_args()

	return args

################################################################################
def find_PRA_files():

	log.info("*-- Getting PRA file(s) --*")

	PRAs = g2pU.find_input_files(args, r"^[\w\.-]+\.xlsx?$")

	df = pd.DataFrame(map(parse_input_file, PRAs), columns=phe.cols)
	df = phe.append(df)[1]
	phe.write()

	return df

#-------------------------------------------------------------------------------
def parse_input_file(PRA):
	path, fname = os.path.split(PRA)
	return [path, fname, gU.MOLIS_name(fname, spaces=True)]

################################################################################
def data_import(df, table):
	"""
	If <args.force>, deletes entries matching <df> from <table> and iteratively
		removes corresponding entries from child tables.
	"""

	if table is g2pU.phe and args.force:
		print("Forcibly overwriting PHENOS. Removing PRA files and their "
			  "SIR data from the PHENOs and SIRs tables.")

		indexes = table.delete(df)
		print("Removing old PRA(s) from archives")
		gU.remove([*map(archive, table.df.iloc[indexes].itertuples())])

	print(f"Importing any new data into {gU.filestem(table.fname)}")
	return table.append(df, write=False)

################################################################################
def main(args):

	global log

	log = g2pU.getLog("phenos")
	g2pU.log = gU.log = log

	if args.statistics: return pU.statistics()

	df = find_PRA_files()
	g2pU.analyse_data(df, ec50, pU.parse_PRA, "PRA")

################################################################################
if __name__ == "__main__":

	args = parse_arguments()
	sys.exit(main(args))

################################################################################