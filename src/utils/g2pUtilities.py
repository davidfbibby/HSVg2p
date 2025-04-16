import collections
import os
import re
import sys
import yaml

import itertools as it
import pandas as pd

from functools import partial, reduce
from glob import glob

 # To get generalUtilities and taskHandler
sys.path.append("/home/dbibby/DevDirectory/Genomancer/src")

from . import gU

"CONSTANTS"

src = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
data_dir = os.path.join(src, os.pardir, "data")
src, data_dir = map(os.path.abspath, (src, data_dir))

"SHORT FUNCTIONS"
def read_tsv(tsv, **kwargs):
	"Wraps gU.read_tsv (specifying index_col), which wraps pd.read_csv"

	return gU.read_tsv(tsv, index_col=0, **kwargs)

def MOLIS_name(name):
	return gU.MOLIS_name(name, spaces=True, suffix=False)

def log_list(data):
	if len(data) > 10: data = data[:5] + ["..."] + data[-5:]
	return data

"CLASSES"

class Table(object):

	def __init__(self, name_or_df, cols=[]):
		self.cols = cols
		t = type(name_or_df)
		if t == str: self.df = self.read(name_or_df)
		elif t == type(pd.DataFrame()): self.df = name_or_df
		else: raise

	def __repr__(self):
		return self.df.to_string()

	def read(self, name):
		self.fname = f"{data_dir}/{name}.tsv"
		df = read_tsv(self.fname, cols=self.cols, comment="\"").fillna("").astype(str)
		self.cols = df.columns
		if "DATE" in df.columns:
			df.DATE = pd.to_datetime(df.DATE, format="%Y-%m-%d")
		return df

	def append(self, df, *,
			   fill="", dup_cols=[], sort_cols=[], reset=True, write=True):
		"""
		Appends <df> data to a Table.df by simple pd.concat. The following
		actions are performed in order (i.e. reset will follow a sort).
		- index		Adds max(Table.df.index) + 1 to the <df> index
		- fill		The fill value for NA cells
		- dup_cols	Removes duplicates of <dup_cols> columns (keeps first)
		- sort		Sorts the Table.df by the columns in <sort>
		- reset		Resets the index after concatenation
		- write		Finally, the table is saved to file by default
		"""

		self.df = pd.concat((self.df, df)).fillna(fill)
		if dup_cols:
			self.df = self.df.drop_duplicates(dup_cols)
		if sort_cols:
			self.df = self.df.sort_values(sort_cols)
		if reset:
			self.df = self.df.reset_index(drop=True)
		if write:
			self.write()

	def filter(self, filters, inverse=False, index=False, copy=False, setop=set.intersection):
		"""
		Returns the sub_df where all <filters> are satisfied. <filters> is a
		list of (col, values) pairs (<col> can be "index"), returning indices
		where col == values. The <inverse> argument defaults to all False, but
		can be a tuple of the same length as <filters>, specifying whether a
		filter returns the inverse or not.
		The intersection of all indices
		"""
		def subdf(cv, i):
			c, v = cv
			v = gU.parse_input_data(v)
			if c == "index": boolean = self.df.index.isin(v)
			else: boolean = self.df[c].isin(v)
			if i: boolean = ~boolean
			return set(self.df[boolean].index)

		if type(filters[0]) != tuple: filters = (filters, )
		if type(inverse) == bool: inverse = it.repeat(inverse, len(filters))

		indices = sorted(
			setop(*it.starmap(subdf, zip(filters, inverse)))
		)
		if index: return indices
		df = self.df.loc[indices]
		if copy: return Table(df.copy())
		return df

	def write(self):
		gU.write_tsv(self.df, self.fname, index=True)

################################################################################

"SET UP TABLES FROM DATA TSVs"

def setup_tables():

	# Variable tables - these get updated with new data
	#	fas		FASTA sequences, with source RUN ID, date of run etc.
	#	var		Variants detected in FASTAs, with an ID link to <fas>
	#	phe		PRA files
	#	sir		Susceptibilities from PRAs, with an ID link to <phe>

	cfg.fas = Table("FASTAS", ["RUNID", "DATE", "NAME", "MOLIS", "SEQ"])
	cfg.var = Table("VARIANTS", ["HGVS", "FAS_ID"])
	cfg.var.df = cfg.var.df.astype({"FAS_ID": int})

	cfg.phe = Table("PHENOS", ["NAME", "MOLIS"])
	cfg.sir = Table("SUSCEPTIBILITIES", ["DATE", "DRUG",  "CONTROL", "EC50", "HSVTYPE", "SIR", "VALID", "PHE_ID"])
	cfg.sir.df = cfg.sir.df.astype({"PHE_ID": int})

	# Static tables - these contain fixed reference information
	#	hsv		Maps older MOLIS IDs to HSV type, based upon an MMD download
	#	thr		Configures the Susceptible/Intermediate/Resistant PRA cut-offs
	#	lit		Table of reported variants and their impact upon drug
	#			susceptibilities
	#	ref		Table of the references used in <lit>
	#	tgt		Table of drugs and the domains targeted (for mutation queries)

	cfg.hsv, cfg.thr, cfg.lit, cfg.ref = map(
		Table, ("HSVTYPES", "THRESHOLDS", "LITERATURE", "REFERENCES")
	)

	cfg.targets = yaml.safe_load(open(f"{data_dir}/external/targets.yaml"))


"TABLE QUERY FUNCTIONS"

def reindex(df, ref_df, on=None, how="inner"):
	df = ref_df.reset_index().merge(df, on=on, how=how).set_index("index")
	return df

"PARSING VARIANTS"

def get_drugs(domain, Core=False):
	drugs = cfg.targets[domain]
	if Core: return drugs["Core"]
	return [*gU.chain(drugs.values())]

def lit_report(hgvs):
	"Reports on a single variant (in HGVS format). See README for rules"

	elements = hgvs.split(".")
	drugs = get_drugs(elements[1])

	#	If HGVS encodes missing data, None is returned. This should be handled upstream
	if elements[2] == "m": return

	lit = hgvs2lit(hgvs)

	# Rule 1 - no reports at all, i.e. completely novel
	if lit.df.empty: return {(hgvs, drug): ("NOVEL", None) for drug in drugs}

	for drug, df in lit.df.groupby("DRUG"):
		groups = get_groups(df.SOURCE)



def get_groups(ref_ids):

	regex = re.compile(r"(?:^|;) ?(.*?),")

	print(ref_ids)
	authors = cfg.ref.filter(("index", ref_ids)).AUTHORS.values
	print(len(authors))
	authors = gU.consolidate([set(regex.findall(a)) for a in authors])

	print(len(authors))

	return authors

def lit_eval_func(sir, counts):
	"""
	Evaluates the collected literature reports in a table.
	If no reports, report '-' (i.e. novel)
	If mixed SEN & RES, report 'S/R?'
	If >=1 report of SEN, report 'S*'
	Else report "R*" (i.e. likely resistant)
	"""

	SR_arr = np.isin(["S", "R"], sir)
	if (~SR_arr).all(): return "-"
	if SR_arr.all(): return "S/R?"
	if SR_arr[0]: return "S*"
	return "R*"

def phe_eval_func(sir, counts):
	"""
	Evaluates the overall resistance levels from a table of PRA results.
	If any are SEN, report 'S'
	If only one is RES, report 'R?', i.e. possible resistance - can occur when
	a novel mutation is seen in conjunction with a resistance mutation.
	"""
	sir = dict(zip(sir, counts))
	if sir["S"] > 0: return "S"
	if sir["R"] == 1: return "R?"
	return "R"

def missing_hgvs(hgvs):
	"Reports the collated resistance interpretation of all missed loci"
	drugs = cfg.targets[mut[1]]
	lit = pd.Series(data=get_missing_lit(*mut), index=drugs, name=hgvs)
	phe = pd.Series(data="?", index=drugs, name=hgvs)

def get_missing_lit(hsvtype, domain, _, loci):
	"""
	Generates a report from all the mutations at the missing loci as per a
	set of FASTA variants. Creates a summary row.
	"""
	loci = np.array(loci.split("_"), dtype=np.int32)[[0, -1]] + [0, 1]
	filters = (
		("AA_POS", map(str, np.arange(*loci))),
		("HSVTYPE", hsvtype),
		("DOMAIN", domain)
	)
	hgvs_list = cfg.lit.filter(filters).HGVS.unique()
	df = pd.concat(map(get_lit, hgvs_list), axis=1).T
	return get_summary_row(df)

def get_summary_row(df):
	"""
	Takes a table of S-R-*-? values and creates a summary row. Used to report
	overall resistance from a FASTA and collated resistance from a missing
	region of a domain.
	"""

	def get_report(ser):
		"Reports the first level observed from descending levels"
		levels = [
			("R", "RES"), ("R*", "LIKELY RES"), ("R?", "POSS RES"), ("-", "?"),
			("?", "?"), ("S/R?", "AMBIGUOUS"), ("S*", "LIKELY SEN"), ("S", "SEN")
		]
		for char, report in levels:
			if (ser == char).any(): return report
		return "?"

	return df.apply(get_report, axis=0)



"SPECIFIC TABLE FUNCTIONS"

hgvs2lit = lambda x: Table(cfg.lit.filter(("HGVS", x)))

hgvs2fasids = lambda x: cfg.var.filter(("HGVS", x)).FAS_ID
fasids2molis = lambda x: cfg.fas.filter(("index", x)).MOLIS
molis2phe = lambda x: cfg.phe.filter(("MOLIS", x))
pheids2sir = lambda x: cfg.sir.filter(("PHE_ID", x))

def hgvs2sir(hgvs):
	fas_ids = hgvs2fasids(hgvs)
	molis = fasids2molis(fas_ids)
	phe_df = molis2phe(molis)
	sir_df = pheids2sir(phe_df.index)
	return Table(pd.merge(sir_df, phe_df, left_on="PHE_ID", right_index=True))