#!/usr/bin/env python
"""
Module for handling sequence data.

Takes FASTAs (or a directory thereof) and processes them for variants, then
processes the variants for resistance association (literature) and phenotypic
resistance (in-house data from the PHENOS module).
If importation is flagged, the sequence file(s) are copied to the database and
dynamic tables are appended ("FASTAS" and "VARIANTS").
"""

import argparse
import os
import re
import sys

import itertools as it
import numpy as np
import pandas as  pd

from datetime import datetime as dt
from glob import glob

from utils import g2pU, gU, sU
from components import rP
from data_init.g2pTables import *

################################################################################
def parse_arguments():
	"""
	Parse command-line arguments.
	"""

	ap = argparse.ArgumentParser(prog="HSVgeno2pheno.py sequences")

	ap.add_argument(
		"-i", "--import_data", action="store_true",
		help="Import sequences to FASTAS table and copy file to ./data/archive"
	)
	ap.add_argument(
		"-r", "--report_data", action="store_true",
		help="Report interpretation (default report settings)"
	)
	ap.add_argument(
		"-s", "--statistics", action="store_true",
		help="Give description of databases (overrides other options)"
	)

	ap.add_argument(
		"--report_format", default=0, type=int,
		help="Report settings (see docs for details) [0]"
	)

	ap.add_argument(
		"-f", dest="single_file", help="path/to/fasta/file. *OVERRIDES -d*"
	)
	ap.add_argument(
		"-d", dest="directory", default=os.getcwd(),
		help="path/to/dir/containing/fasta/file(s) [.]"
	)
	ap.add_argument(
		"--recursive", action="store_true",
		help="Set to look in all subdirectories of <directory>"
	)

	ap.add_argument(
		"--runid", help="Run ID (overrides automatic parsing)"
	)
	ap.add_argument(
		"--rundate", help="Run date (YYMMDD, overrides automatic parsing)"
	)
	ap.add_argument(
		"--force", action="store_true",
		help="Force overwrite of existing records with identical data"
	)

	args = ap.parse_args()

	if not (args.statistics or args.import_data or args.report_data):
		ap.print_help()
		log.info("At least one of -i, -rm or -s must be passed")
		return 65

	rP.format_code = args.report_data

	return args

################################################################################
def find_FASTA_files():

	log.info("*-- Getting FASTA file(s) --*")

	FASTAs = g2pU.find_input_files(args, r"^[\w\.-]+\.fas?(?:ta)?$")
	index = []
	df = pd.DataFrame(map(parse_input_file, FASTAs), columns=fil.cols)
	if not df.empty:
		index, df = fil.append(df)
		fil.write()

	return index, df

#-------------------------------------------------------------------------------
def parse_input_file(fas):
	path, fname = os.path.split(fas)

	regex = re.compile(r"^(\d{6})_(.*?)(_(TK|POL).*)?$")
	for target in map(gU.filestem, (fas, dirname := os.path.dirname(fas))):
		if match := regex.match(target.upper()):
			date, runid = match.groups()[:2]
			break
	else:
		date = dt.fromtimestamp(os.path.getmtime(fas))
		runid = gU.filestem(dirname)

	# Passed rundate and runid values override derived values
	date = args.rundate if args.rundate else pd.to_datetime(date, format="%y%m%d")
	runid = args.runid if args.runid else runid

	return [path, fname, date, runid]

################################################################################
def data_import(df, table):
	"""
	If not <args.import_data>, returns.
	If <args.force>, deletes entries matching <df> from <table> and iteratively
		removes corresponding entries from child tables.
	"""

	if not args.import_data: return df.index, df
	if table is g2pU.fil and args.force:
		log.info("Forcibly overwriting FASTA. Removing sequences and their "
			  "variants from the FASTAs and VARIANTs tables.")

		indexes = table.delete(df)
		log.info("Removing old FASTA(s) from archives")
		gU.remove([*map(archive, table.df.iloc[indexes].itertuples())])

		# tmp_df = pd.concat((table.df, df))
		# dup = tmp_df.duplicated(keep="last")
		#
		# print("Removing old FASTA(s) from archives")
		# gU.remove([*map(archive, table.df.iloc[dup].itertuples())])
		#
		# print("Removing old entries from FILES, FASTAS, and VARIANTS")
		# indexes = tmp_df.index.values[dup]
		# table.delete(indexes)
		# next_table = table.child
		# while next_table is not None:
		# 	indexes = next_table.filter(("PARENT_ID", indexes), index=True)
		# 	next_table.delete(indexes)
		# 	next_table = next_table.child

	log.info(f"Importing any new data into {gU.filestem(table.fname)}")
	return table.append(df, write=False)

################################################################################
def main(args):

	global log

	log = g2pU.getLog("sequences")
	g2pU.log = gU.log = rP.log = log

	if args.statistics: return sU.statistics()

	index, fil_df = find_FASTA_files()
	fas_df = g2pU.analyse_data(fil_df, fas, sU.parse_FASTA, "FASTA")
	print(fas_df)
	input()

	var_df = g2pU.analyse_data(fas_df, var, sU.parse_variants, "SEQ")
	print(var_df)
	input()

	if not args.import_data: fil.delete(index)	# Remove "new" data
	if not args.report_data: return

	log.debug(f"REPORT format: {args.report_data}")
	log.info("*-- Merging FASTA and VARIANT data for reporting --*")

	df = pd.merge(
		fil_df.drop("LOCATION", axis=1), fas_df.reset_index(),
		how="outer", left_index=True, right_on="PARENT_ID"
	).set_index("index")

	df = pd.merge(
		df.drop(["PARENT_ID", "SEQ"], axis=1), var_df,
		how="outer", left_index=True, right_on="PARENT_ID"
	)

	# Send sub-DFs on a per FASTA basis
	gU.threaded(rP.generate_FASTA_report, df.groupby("PARENT_ID"), 16, iterate=True)

################################################################################

if __name__ == "__main__":

	args = parse_arguments()

	returncode = args if type(args) is int else main(args)
	sys.exit(returncode)