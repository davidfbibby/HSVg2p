#!/usr/bin/env python

import argparse
import os
import re
import sys

from utils import g2pU, gU

################################################################################
"Parse command-line arguments"

def parse_arguments():

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

def main(args):

	# Test args
	if not (args.import_data or args.report_data):
		print("Either -i must be set and/or -r must have an argument")
		return 65

	# Establish FASTA file "list"
	regex = re.compile(r"^\w+\.fas?(?:ta)?$")
	if args.fasta is not None: args.fasta=[args.fasta]
	else: args.fasta = [*filter(regex.search, os.listdir(args.directory))]

################################################################################
if __name__ == "__main__":

	sys.exit(main(parse_arguments()))