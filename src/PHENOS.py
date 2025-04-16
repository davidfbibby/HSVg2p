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
		"-p", dest="previous",
		help="Report all phenotypic data from the same sample(s)"
	)
	ap.add_argument(
		"-f", dest="pheno", help="path/to/pheno/file. OVERRIDES -d"
	)
	ap.add_argument(
		"-d", dest="directory", default=os.getcwd(),
		help="path/to/dir/containing/pheno/file(s) [.]"
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


################################################################################
def main(args):

	# Establish PHENO file "list"
	regex = re.compile(r"^\w+\.xlsx?$")
	if args.pheno is not None: args.pheno=[args.pheno]
	else: args.pheno = [*filter(regex.search, os.listdir(args.directory))]

################################################################################
if __name__ == "__main__":

	sys.exit(main(parse_arguments()))