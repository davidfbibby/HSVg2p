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
		"--table", required=True, help="Table to amend"
	)
	ap.add_argument(
		"--amendment", required=True, help="Amendment (see docs)"
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