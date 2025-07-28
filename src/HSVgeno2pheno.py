#!/usr/bin/env python
"""
Very lightweight top-level script for managing HSV geno2pheno.

Handles the first passed argument and distributes the remainder to the
appropriate module. Also receives the error codes from said module and reports
relevant messages.
"""

import argparse
import os
import sys

import subprocess as sp

from collections import defaultdict
from glob import glob

from utils import gU

################################################################################
code_dict = defaultdict(
	lambda: "Unspecified error",
	{
		 2: "Invalid arguments",
		65: "Invalid arguments",
		66: "Tool name not valid",
		67: "No FASTAs to process",
})

################################################################################
def parse_arguments():
	"""
	Parses the first argument only. Returns the rest as <others>.
	"""

	ap = argparse.ArgumentParser(add_help=False)
	ap.add_argument(
		"tool", nargs="?", default="",
		help="Choose from 'sequences', 'phenos', 'mutation', 'molis', 'modify' "
		     "and 'suppress'"
	)
	ap.add_argument(
		"-h", "--help", action="store_true",
		help="show this help message and exit"
	)

	args, others = ap.parse_known_args()

	if (tool := args.tool.upper()): return tool, others

	ap.print_help()
	if not args.help: error_report(65)

	sys.exit()

################################################################################
def error_report(code):
	"""
	Prints relevant error messages based upon <code>.
	"""

	print(f"Failed with exit code {code}: {code_dict[code]}")
	sys.exit(code)

################################################################################
def HSVgeno2pheno(tool, others):
	"""
	Takes the first passed argument and uses it to call the appropriate module
	with the remaining arguments (<others>).
	"""

	module = f"{os.path.dirname(__file__)}/{tool.upper()}.py"
	try:
		print(f"Running {gU.filestem(module)}")
		returncode = sp.run([module, *others]).returncode
	except FileNotFoundError:
		returncode = 66
	return returncode

################################################################################
if __name__ == "__main__":
	"""
	Calls the argument parser, passes <tool> and <others> to _HSVgeno2pheno_, and
	exits on a zero return code. Reports the error otherwise.
	"""

	if (returncode := HSVgeno2pheno(*parse_arguments())) == 0: exit(0)
	error_report(returncode)

################################################################################