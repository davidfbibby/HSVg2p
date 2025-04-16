#!/usr/bin/env python

import argparse
import sys

import subprocess as sp

from collections import defaultdict

from utils import gU

################################################################################
"Parse command-line arguments"

def parse_arguments():

	ap = argparse.ArgumentParser()
	ap.add_argument("tool")

	args, others = ap.parse_known_args()
	tool = args.tool.upper()

	return tool, others

################################################################################
def error_report(code):
	code_dict = {
		65: "Invalid arguments",
		66: "Tool name not valid"
	}
	code_dict = defaultdict(lambda: "Unspecified error", code_dict)

	print(f"HSVgeno2pheno.py failed with exit code {code}: {code_dict[code]}")
	sys.exit(code)

################################################################################
def HSVgeno2pheno(tool, others):

	if not gU.exists(module := f"{g2pU.src}/{tool}.py"):
		print("Tool must be one of SEQUENCES, PHENOS, MOLIS, MUTATIONS, "
			 f"SUPPRESS, or MODIFY ({tool} entered)")
		return 66

	cp = sp.run([module, *others])
	return cp.returncode

################################################################################
if __name__ == "__main__":

	if (returncode := HSVgeno2pheno(*parse_arguments())) == 0: exit()
	error_report(returncode)
