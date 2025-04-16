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
		"-l", "--literature", help="Report information from literature"
	)
	ap.add_argument(
		"-p", "--phenos", help="Report information from in-house phenotyping"
	)
	ap.add_argument(
		"--recent", action="store_true",
		help="Report only the most recent phenotyping"
	)


	return ap.parse_args()


################################################################################


################################################################################
def main(args):

	if args.table not in ("FASTAS", "PHENOS"):
		print("Table argument must be either FASTAS or PHENOS")
		sys.exit(65)

	pass

################################################################################
if __name__ == "__main__":

	sys.exit(main(parse_arguments()))