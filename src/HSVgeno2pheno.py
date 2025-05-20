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

from utils import gU

################################################################################
code_dict = defaultdict(
	lambda: "Unspecified error",
	{
		65: "Invalid arguments",
		66: "Tool name not valid",
		67: "No FASTAs to process",
})

################################################################################
def parse_arguments():
	"""
	Parses the first argument only. Returns the rest as <others>.
	"""

	ap = argparse.ArgumentParser()
	ap.add_argument("tool")

	args, others = ap.parse_known_args()
	tool = args.tool.upper()

	return tool, others

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

	if not gU.exists(module := f"{os.path.dirname(__file__)}/{tool}.py"):
		print("Tool must be one of SEQUENCES, PHENOS, MOLIS, MUTATIONS, "
			 f"SUPPRESS, or MODIFY ({tool} entered)")
		return 66

	cp = sp.run([module, *others])
	return cp.returncode

################################################################################
if __name__ == "__main__":
	"""
	Calls the argument parser, passes <tool> and <others> to _HSVgeno2pheno_, and
	exits on a zero return code. Reports the error otherwise.
	"""

	if (returncode := HSVgeno2pheno(*parse_arguments())) == 0: exit(0)
	error_report(returncode)

################################################################################