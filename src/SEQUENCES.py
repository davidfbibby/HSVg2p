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

import pandas as  pd

from glob import glob

from utils import g2pU, gU
from components import parse_variants as pV

################################################################################
id_cols = ["FILENAME", "DATE", "RUNID", "LOCATION"]

################################################################################
def parse_arguments():
	"""
	Parse command-line arguments.
	"""

	ap = argparse.ArgumentParser()

	ap.add_argument(
		"-i", "--import_data", action="store_true",
		help="Import sequences to FASTAS table and copy file to ./data/archive"
	)
	ap.add_argument(
		"-r", "--report_data",
		help="Control reporting of interpretation (see docs)"
	)
	ap.add_argument(
		"-f", dest="fasta", help="path/to/fasta/file. OVERRIDES -d"
	)
	ap.add_argument(
		"-d", dest="directory", default=os.getcwd(),
		help="path/to/dir/containing/fasta/file(s) [.]"
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

	return ap.parse_args()

################################################################################
def get_fastas():
	"""
	Parses the -f argument if present, else the -d argument. Returns a list of
	FASTA file(s).
	"""

	if args.fasta is not None: return [os.path.abspath(args.fasta)]

	map_func = lambda f: re.search(r"^\w+\.fas?(?:ta)?$", os.path.basename(f))

	return [*filter(map_func, glob(f"{os.path.abspath(args.directory)}/*"))]

################################################################################
def get_df(fas):
	"""
	Takes an input FASTA file, imports the FASTAs into a DataFrame, and derives
	metadata from the filename and/or passed arguments. These are used to
	compare new sequences to existing data.
	"""

	location, filename = os.path.split(fas)
	print(f"Processing {filename} in {location}")

	df =  gU.fas2df(fas)

	df["SEQ"] = df.SEQ.str.replace("-", "")
	df["MOLIS"] = df.NAME.apply(g2pU.MOLIS_name)

	if date_str := re.match(r"^\d{6}", filename):
		date = pd.to_datetime(date_str.group(0), format="%y%m%d")
		runid = get_runid(gU.filestem(basename).split("_" , 1)[1])
	else:
		date = pd.NaT
		runid = gU.filestem(location)

	date = args.rundate if args.rundate else date
	runid = args.runid if args.runid else runid

	metadata = (filename, date, runid, location)
	df[id_cols] = metadata

	return df

################################################################################
def import_sequences(df):
	"""
	If --force, any existing rows in the "FASTAS" table are deleted, and the
	corresponding entries in the "VARIANTS" table are removed
	In all cases, the new FASTA file is copied into the archive (overwriting an
	existing one), the resulting "FASTAS" table is appended with the new data,
	and NAME/id_col duplicates are removed (i.e. if the data already exists and
	not --force, it will be deduplicated immediately after being appended,
	leaving the table unchanged.)
	"""
	#
	if args.force:
		indexes = get_indexes(df[id_cols[:-1]].iloc[0].values)
		g2pU.fas.df = g2pU.fas.filter(
			filters=("index", indexes), inverse=True, copy=True
		)
		g2pU.fas.df = g2pU.fas.filter(
			filters=("FAS_ID", indexes), inverse=True, copy=True
		)

	gU.shell("cp", fas, f"{g2pU.data_dir}/archive")
	g2pU.fas.append(df=df, dup_cols=["NAME"] + id_cols[:-1])

################################################################################
def report_sequences(df):
	pass


################################################################################
def main(args):

	# Test args
	if not (args.import_data or args.report_data):
		print("Either -i must be set and/or -r must have an argument")
		return 65

	if (fastas := get_fastas()) == []: return 67 	# No FASTAs to process

	for fas in fastas:
		df = get_df(fas)
		if args.import_data: import_sequences(df)
		if args.report_data: report_sequences(df)

################################################################################
def get_indexes(values):
	"""Finds indexes of a FASTA file in the "FASTAS" table"""

	return g2pU.fas.filter(filters=[*zip(id_cols, values)], index=True)


################################################################################
if __name__ == "__main__":

	args = parse_arguments()
	sys.exit(main(args))