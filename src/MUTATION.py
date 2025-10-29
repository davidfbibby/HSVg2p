#!/usr/bin/env python

import argparse
import os
import re
import sys

from functools import partial

import itertools as it
import numpy as np

from utils import g2pU, gU
from components import rP
from data_init.g2pTables import *

################################################################################
"Parse command-line arguments"

def parse_arguments():

	ap = argparse.ArgumentParser(add_help=False)

	"""
	1.pol.p.S724N is a *mutation*
	1.pol.p.724 is a *locus*
	"""

	ap.add_argument(
		"-m", "--mutation", type=str,
		help="Query a mutation in HGVS format - if specifics are omitted, the "
			 "locus is queried."
	)

	ap.add_argument(
		"-e", "--extend", default=0, type=int, choices=[0, 2, 5, 10],
		help="Explore N loci flanking a locus (N=2, 5, 10) [0]."
	)
	ap.add_argument(
		"-a", "--homology", action="store_true",
		help="Also interrogate homologous locus/loci from the other HSV type(s)"
	)

	args = ap.parse_args()

	return args



def report_mut():

	raw, lit, var = (table.filter(("HGVS", args.mutation))
						for table in (g2pU.raw, g2pU.lit, g2pU.var))

	if raw.empty:
		print("Unreported variant")
	else:
		print(raw.sort_values(["DRUG", "SUSCEPTIBILITY"]))
		print(lit.sort_values(["DRUG", "SUSC"]))
	if var.empty:
		print("Variant not in sequence database")
	else:
		@gU.memo
		def drug2drugs(drug):
			df = g2pU.targets.loc[:, ~np.isnan(g2pU.targets.loc[drug])]
			return df.dropna(how="all").index.values


		def func(molis, drug):
			indices = g2pU.fas.filter(("MOLIS", molis), index=True)
			df = g2pU.var.filter(("PARENT_ID", indices))
			# print(df)
			df = g2pU.lit.filter((
				("HGVS", df[df.HGVS!=args.mutation].HGVS),
				("DRUG", drug2drugs(drug)),
				("SUSC", ["R", "R*", "?"])
			))
			if df.empty: return ""
			return ", ".join(df.HGVS.unique())

		molis = g2pU.fas.filter(("index", var.PARENT_ID)).MOLIS
		ec50 = g2pU.phe.filter(("MOLIS", molis)).merge(
			g2pU.ec50.df, left_index=True, right_on="PARENT_ID", how="inner"
		)
		for drug, df in ec50.groupby("DRUG"):
			df["OTHER_MUTS"] = df.MOLIS.apply(partial(func, drug=drug))
			df = rP.evaluate_ec50s(df)
			df = df[df.VALID=="pass"]
			print(df.drop(["FILENAME", "LOCATION", "PARENT_ID", "VALID", "HSV"], axis=1))


def main(args):

	mut = g2pU.Mutation(args.mutation, args.homology, args.extend)


	variants, literature = mut.get_mutations()

	print(variants)
	[print(v.literature(), "\n") for v in variants]
	[print(v.phenotypes(), "\n") for v in variants]

	print(literature)
	[print(v.literature(), "\n") for v in literature]
	[print(v.phenotypes(), "\n") for v in literature]

	exit()


################################################################################
if __name__ == "__main__":

	args = parse_arguments()
	main(args)
