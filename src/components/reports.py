"""
Writes a report for HSV geno2pheno.

input
-----

df		DataFrame in g2pU.var format
format	Integer detailing extent of reported material

Depending upon the format value passed, more or less information is generated
from the report, ranging from basic - literature reports on the detected
variants for the three basic drugs - ACV, FOS & CDV, PRI depending upon
domain(s) - to an all-inclusive report based upon phenotypic data and literature
information on all reportable drugs (i.e. including BVDU, GCV, AMN etc.)

Upstream of this module, it is expected that report generation will be requested
in the following circumstances:
- Analysis of a single domain sequence
- Analysis of a single MOLIS no. (with multiple domains)
- Analysis of a single mutation
- Analysis of a range of mutations


The <format> integer is built up from bits, a la SAM FLAGs. In most cases, the
flag will overwrite defaults:

1		EXPAND drug selection to include all possible drugs for the domain(s)
2		EXCLUDE UKHSA phenotypic data in the analysis
4		DO NOT merge reports with the same MOLIS no.

"""
import os
import re

import itertools as it
import numpy as np
import pandas as pd

from collections import defaultdict
from types import SimpleNamespace

from utils import g2pU, gU
from data_init.g2pTables import *

format_code = 0
log = gU.dummy_log()
args = SimpleNamespace()

################################################################################
class Report():
	"""
	Top-level class for reports. Is subclassed for different reporting needs
	(e.g. FASTA, MUTATION, MOLIS, etc.). All instances have parameters that
	govern output:

	drugs		If 1, include all drugs for a report. If 2, then limit to the
				headline drug(s).
	susc		If False, include only non-"S" mutations when expanding a locus.
	molis		If True, collates data from all reports matching the MOLIS ID.
	expand		The number of mutations to look either side of a locus when
				reporting upon novel mutations.
	homology	If True, also look at homologous loci in other HSVs when
				reporting upon novel mutations.
	phenotypes	If True, report local phenotypes of samples with ambiguous or
				novel mutations.

	Because several parameters govern reporting across multiple report types,
	this top-level class has some class methods:

	get_lit		Returns the literature of a mutation (passed as <hgvs>), with
				filters according to drugs and the class susceptibility
				threshold. Adds the hgvs to a new column (for missing data), and
				creates dummy rows for novel mutations (with 0 CITATIONS).
	"""

	def __init__(self, drugs=2, molis=True, expand=10, homology=True, phenotypes=True):
		self.drugs = drugs
		self.molis = molis
		self.expand = expand
		self.homology = homology
		self.phenotypes = phenotypes

	def get_lit(self, hgvs):
		mut = g2pU.Mutation(hgvs, self.homology, self.expand)
		muts = mut()
		drugs = get_drugs_by_code(mut.d, self.drugs)
		lit_df = lit.filter(filters=(("HGVS", muts), ("DRUG", drugs)))

		if lit_df.empty: lit_df = pd.DataFrame(
			data=[[hgvs, drug, 0, "?"] for drug in drugs],
			columns=lit.cols
		)

		lit_df["VARIANT"] = hgvs

		return lit_df

class FASTA_report(Report):
	"""
	Subclass of *Report*, designed to report upon FASTA requests made through
	SEQUENCES.py (HSVgeno2pheno.py sequences -r ....). In addition to the

	"""
	def __init__(self, df, **kwargs):
		super().__init__(**kwargs)
		self.df = df

	def interpret(self):

		self.lit_df = pd.concat(map(self.get_lit, self.df.HGVS))

@gU.memo
def get_drugs_by_code(domain, drug_code=2):
	return tgt.filter((domain, drug_code), index=True)


################################################################################
"""Takes a DF merged from <g2pU.fil>-, <g2pU.fas>- & <g2pU.var>-like DFs and
generates a report. Uses the <PARENT_ID> field (of <g2pU.var>) to separate
individual reports. Varying amounts of historic, literature, and phenotyping
data can be included in the report. For phenotyping, this can include (a) any
phenotypes corresponding to the MOLIS ID, and/or (b) phenotypes from variants
with either ambiguous or no literature reports. This latter can lead to
iterative interrogation of resistant samples, to determine if other, known
genotypic markers can account for the phenotype, rather than the variant under
scrutiny.
"""

def generate_FASTA_report(index, df):
	"Generates a report for a single sample/domain combination"

	log.info(f"Generating report for FASTA {index} (MOLIS={df.iloc[0].MOLIS})")
	report = FASTA_report(df)
	report.interpret()
	print(report.lit_df)
	print()


def get_phenos(hgvs):

	var_df = g2pU.var.filter(("HGVS", hgvs))

	molis = g2pU.fas.filter(("index", var_df.PARENT_ID)).MOLIS.unique()
	phe_df = g2pU.phe.filter(("MOLIS", molis))


	return var_df, molis, phe_df


################################################################################
"Takes a DF subsetted from <g2pU.phe> and generates a phenotyping report"
def phenotyping_reports(phe_df):
	pass

################################################################################
"Takes a DF like <g2pU.ec50> and evaluates against <g2p2U.thr> (valid / SIR)"
def evaluate_ec50s(df):

	df[["SIR", "VALID"]] = df.apply(evaluate_ec50, axis=1, result_type="expand")
	return df

#-------------------------------------------------------------------------------
def evaluate_ec50(row):
	"Determines susceptibility (S/I/R) and validity for a single drug"

	SI, IR = get_thresholds(str(row.HSV)[0], row.DRUG)
	if SI == 0:
		susc, valid = ("-", "inv")
	else:
		ec50, ctrl = map(lambda x: float(str(x).strip("<>")), (row.EC50, row.CTRL))
		susc = "S" if ec50 < SI else "R" if ec50 >= IR else "I"
		valid = "pass" if ctrl < SI else "fail"

	return (susc, valid)

#-------------------------------------------------------------------------------
@gU.memo
def get_thresholds(HSV, drug):
	thr = g2pU.thr.filter((("HSV", HSV), ("DRUG", drug)))
	if thr.empty: return (0, 0)
	return tuple(thr.iloc[0][["S_I", "I_R"]].astype(float))
################################################################################