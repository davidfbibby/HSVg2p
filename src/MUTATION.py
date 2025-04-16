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
		"--range", default=0, type=int, help="Explore N flanking loci [0]"
	)
	ap.add_argument(
		"--homology", action="store_true",
		help="Also interrogate the homologous loci from the other HSV type"
	)
	ap.add_argument(
		"mutation", help="Mutation/locus in HGVS format"
	)

	return ap.parse_args()


################################################################################


################################################################################
def main(args):

	if len(elements := args.mutation.split(".")) != 4:
		print("Mutation not in valid HGVS format")
		return 65

	pass

################################################################################
if __name__ == "__main__":

	sys.exit(main(parse_arguments()))