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

	ap = argparse.ArgumentParser(prog="HSVgeno2pheno.py sequences")

	ap.add_argument(
		"-i", "--import_data", action="store_true",
		help="Import sequences to FASTAS table and copy file to ./data/archive"
	)
	ap.add_argument(
		"-r", "--report_data", nargs="?", const=1,
		help="Control reporting of interpretation (see docs for format codes) [1]"
	)

	ap.add_argument(
		"-f", dest="fasta", help="path/to/fasta/file. *OVERRIDES -d*"
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

	args = ap.parse_args()

	if not (args.import_data or args.report_data):
		ap.print_help()
		print("At least one of -i or -r must be passed")
		return 65

	return args

################################################################################
def get_fastas():
	"""
	Parses the -f argument if present, else the -d argument. Returns a list of
	FASTA file(s).
	"""

	map_func = lambda f: re.search(r"^[\w\.]+\.fas?(?:ta)?$", os.path.basename(f))

	if args.fasta is not None:
		return [os.path.abspath(args.fasta)]
	else:
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
	df.index += 1

	df["SEQ"] = df.SEQ.str.replace("-", "")
	df["MOLIS"] = df.NAME.apply(g2pU.MOLIS_name)
	df["FAS_ID"] = int(0)

	date = pd.NaT
	if date_str := re.match(r"^\d{6}", runid := gU.filestem(location)):
		date = pd.to_datetime(date_str.group(0), format="%y%m%d")
		runid = runid.split("_" , 1)[1]

	# Passed rundate and runid values override derived values
	date = args.rundate if args.rundate else date
	runid = args.runid if args.runid else runid

	# Associate all metadata with each sample
	metadata = (filename, date, runid, location)
	df[id_cols] = metadata

	return df

################################################################################
def import_sequences(df, fas):
	"""
	If --force, any existing rows in the "FASTAS" table are deleted, and the
	corresponding entries in the "VARIANTS" table are removed
	In all cases, the new FASTA file is copied into the archive (overwriting an
	existing one), the resulting "FASTAS" table is appended with the new data,
	and NAME/id_col duplicates are removed (i.e. if the data already exists and
	not --force, it will be deduplicated immediately after being appended,
	leaving the table unchanged.)
	"""

	if args.force:
		print("Forcing overwrite")

		indexes = get_indexes(df[id_cols[:-1]].iloc[0].values)
		print("Previous indexes of input FASTA file sequences:",
			  ", ".join(map(str, indexes))
		)

		print("Deleting from FASTAS table")
		g2pU.fas.df = g2pU.fas.filter(
			filters=("index", indexes), inverse=True
		)

		print("Deleting associated variants from VARIANTS table")
		g2pU.var.df = g2pU.fas.filter(
			filters=("FAS_ID", indexes), inverse=True
		)

	print("Appending new data to FASTAS table, and deleting if duplicated.")
	g2pU.fas.append(df=df.drop("FAS_ID", axis=1), dup_cols=["NAME"] + id_cols[:-1])

	indexes = get_indexes(df[id_cols[:-1]].iloc[0].values)
	df.loc[:, "FAS_ID"] = indexes

	archive = f"{g2pU.data_dir}/archive/{indexes[0]}_{os.path.basename(fas)}"
	print(f"Copying FASTA files to {archive}")
	gU.shell("cp", fas, archive)

	return df

################################################################################
def process_variants(df):

	print("Processing variants in FASTA(s)")

	tested = df.FAS_ID.isin(g2pU.var.df.FAS_ID)
	var_df = pd.concat((
		pV.get_vars(pV.map_seqs(df[~tested])),
		g2pU.var.filter(filters=("FAS_ID", df.FAS_ID.values))
	))

	if args.import_data:
		g2pU.var.append(df=var_df, dup_cols=["HGVS", "FAS_ID"])

	print(df)
	print(var_df)
	
################################################################################
def report_sequences(df):
	pass


################################################################################
def main(args):

	# Test args
	if (fastas := get_fastas()) == []: return 67 	# No FASTAs to process

	for fas in fastas:
		df = get_df(fas)
		if args.import_data: df = import_sequences(df, fas)
		var_df = process_variants(df)
		if args.report_data: report_sequences(df)

################################################################################
def get_indexes(values):
	"""Finds indexes of a FASTA file in the "FASTAS" table"""

	return g2pU.fas.filter(filters=[*zip(id_cols, values)], index=True)


################################################################################
if __name__ == "__main__":

	args = parse_arguments()
	returncode = args if type(args) is int else main(args)
	sys.exit(returncode)