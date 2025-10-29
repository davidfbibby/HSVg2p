"""
Sub-component of HSVg2p that looks in the new/PRAs directory for new
HSV phenotyping files. Thus far, each such file deals with a single sample,
with PRA data against one or more drugs, often three or four. After cross-
referencing the new data against the existing data in processed/PHENO.tsv, new
files are parsed for their EC50 data. Control and sample values are obtained,
and by reference to config-determined thresholds, run validity and
susceptibility is assigned.
These come in one of several predictable formats, thus requiring individual
parsing functions.

*** REQUIRES xlrd AND openpyxl TO PARSE THE PRA SPREADSHEETS IN PANDAS ***
"""

import os
import xlrd

import itertools as it
import numpy as np
import pandas as pd

from datetime import datetime as dt
from functools import partial

from g2pUtils import g2pU, gU

################################################################################
"Takes a PRA file and returns the MOLIS name, run date, and drug data"

def import_file(PRA, *args, **kwargs):

	basename = os.path.basename(PRA)

	return pd.DataFrame(
		data=[[basename, g2pU.MOLIS_name(basename)]],
		columns=cfg.phe.cols
	)

################################################################################
"Parses all new PRA files according to their type. Calculates EC50s"

def process(df):
	"Finds new PRAs"

	cfg.log_prefix = "parse_PRAs"
	cfg.log.info(f"Parsing {len(df)} new phenotype files")
	if (df := parse_PRAs(df)) is None: return

	cfg.log.info("Assigning susceptibilities and validities")
	get_susc(df)

################################################################################
"Obtains EC50s (and HSV type and run date) from PRAs"
def parse_PRAs(df):


	if df.empty: return cfg.log.info("\tNo new phenotypes to calculate")

	df = g2pU.reindex(cfg.hsv.df, df, on="MOLIS", how="left")
	dfs = [*filter(
		lambda x: x is not None,
		gU.threaded(func=parse_PRA, data=[*df.iterrows()], procs=16)
	)]

	if len(dfs) == 0: return cfg.log.error("\tNo parseable data")
	df = pd.concat(dfs).reset_index(drop=True)

	return df

#-------------------------------------------------------------------------------

def parse_PRA(index, row):
	"Parses a single PRA, returning drug-EC50-control EC50 data, plus run date"

	PRA = f"{cfg.import_folder}/{row.NAME}"
	dfs = pd.read_excel(PRA, sheet_name=None)
	dfs.pop("Results", None)

	if "ACV" in dfs: data = (it.starmap(parse_old, dfs.items()), row.HSVTYPE)
	elif (df := dfs.get("HSV PRA")) is not None: data = parse_new(make_arr(df))
	else: data = None
	if data is None:
		return cfg.log.error(f"{row.NAME} cannot be parsed.")

	# Populate DataFrane, update fields and format the Date
	data, HSVTYPE = data
	df = pd.DataFrame(filter(None, data), columns=cfg.sir.cols[:4])
	df[["PHE_ID", "HSVTYPE"]] = index, HSVTYPE
	if (rundate := min(df.pop("DATE"))) is pd.NaT:
		rundate = dt.fromtimestamp(os.path.getmtime(PRA))
	df["DATE"] = pd.to_datetime(rundate, format="%Y-%m-%d")

	return df

#-------------------------------------------------------------------------------
def parse_old(drug, df):
	"Parses the old-style PRA worksheets"

	drug = drug.replace("PEN", "PCV")
	if cfg.thr.filter(("DRUG", drug)).empty: return

	arr = make_arr(df)
	rundate = parse_timestr(find_cells(arr, "Date", 1, dtype=str)[0][0])

	# x_arr = No drug plaque counts
	x_arr = find_cells(arr, "No Drug", 1, 0, 3, 3)[[0, 2]]

	# p_arr = test plaque counts, then divided by x_arr
	p_arr = find_cells(arr, "Plaque", 0, 2, 3, 5)[[0, 2]] / get_xDrug(x_arr)

	# c_arr = concentrations for each member of p_arr
	c_arr = find_cells(arr, "Plaque", -1, 2, 1, 5)[0]

	func = partial(get_EC50, c_arr=c_arr)

	return [rundate, drug, *map(func, p_arr)]

#-------------------------------------------------------------------------------
def parse_new(arr):
	"Parses the new-style PRA worksheets"

	def func(text):
		cells = find_cells(arr, text, 0, 2, 1, 3, dtype=str)[0]
		nz = cells.nonzero()[0]
		return "" if nz.size == 0 else cells[nz[-1]]

	rundate, HSVTYPE = map(func, ("ixed", "type"))
	rundate = parse_timestr(rundate)

	# x_arr = No drug plaque counts
	rows, cols = cell_coords(arr, "xDrug")
	x_arr = np.stack([arr[r + 1:, c][:3] for r, c in zip(rows, cols)])
	if np.any(~np.char.isnumeric(x_arr)): return

	arr = arr[rows[0]:, :cols[0]]

	def parse_drug(col):
		rows, _ = cell_coords(arr, drug := str(arr[0, col]))
		if cfg.thr.filter(("DRUG", drug)).empty: return

		try:
			pc = np.stack([arr[r + 5:r:-1, col] for r in rows]).astype(float)
		except Exception:
			return

		# p_arr = test plaque counts, then divided by x_arr
		p_arr = pc[:2] / get_xDrug(x_arr)

		# c_arr = concentrations for each member of p_arr
		func = partial(get_EC50, c_arr=pc[2])

		return [rundate, drug, *map(func, p_arr)]

	return (map(parse_drug, arr[0].nonzero()[0]), HSVTYPE)

################################################################################
"Assigns valid/invalid S/I/R"

def get_susc(df):

	sir_df = pd.DataFrame(
		gU.threaded(func=assign_susc, data=[*df.iterrows()], procs=16),#, iterate=True),
		columns=["index", "SIR", "VALID"]
	).set_index("index")

	cfg.log.info(f"Evaluated {len(sir_df)} susceptibilities and validities")

	cfg.sir.append(df.join(sir_df), dup_cols=["PHE_ID", "DRUG"], sort_cols=["PHE_ID"])

################################################################################

def make_arr(df):
	"Converts a DF into a np. array"

	return df.fillna("").to_numpy().astype(str)

def parse_timestr(timestr):
	"Parses Excel and strings into datetimes"

	if timestr in("0", "nan"): return pd.NaT
	if timestr.isnumeric(): return xlrd.xldate_as_datetime(int(timestr), 0)
	try:
		return pd.to_datetime(timestr, dayfirst=True)
	except:
		return pd.NaT

def cell_coords(arr, text):
	"Returns coordinates of cells in an np array matching a passed string pattern"

	return np.where(np.char.find(arr, text) != -1)

def find_cells(arr, text, row_offset=0, col_offset=0, n_rows=1, n_cols=1, dtype=float):
	"Finds the first cell from cell_coords and returns an offset range of cells"

	cells = cell_coords(arr, text)
	if cells[0].size == 0: return []
	r, c = cells
	r = r[0] + row_offset
	c = c[0] + col_offset

	cells = arr[r:, c:][:n_rows, :n_cols].astype(dtype)

	return cells

def get_xDrug(x_arr):
	"Converts string plaque counts from no-drug controls into a divideable array"

	return np.expand_dims(np.mean(x_arr.astype(int), axis=1), 1)

@gU.memo
def concs(arr):
	return np.log(arr)

def get_EC50(p_arr, c_arr):
	"Calculates the EC50 from matching plaque% & conc arrays (ascending)"

	index = (p_arr >= 0.5).nonzero()[0]

	if index.size == 0: EC50 = f"<{c_arr[-1]}"
	elif (index := index[0]) == 0: EC50 = f">{c_arr[0]}"
	else:
		c_arr = concs(c_arr)
		proportion = (p_arr[index] - 0.5) / (p_arr[index] - p_arr[index - 1])
		conc = c_arr[index] - ((c_arr[index] - c_arr[index - 1]) * proportion)
		EC50 = np.format_float_positional(np.exp(conc), 4, fractional=False)

	return EC50

################################################################################
def assign_susc(index, row):
	"Determines susceptibility (S/I/R) and validity for a single drug"

	thr = cfg.thr.filter((("HSVTYPE", row.HSVTYPE), ("DRUG", row.DRUG)))
	if thr.empty: return [index, "-", "-"]
	SI, IR = thr.iloc[0][["S_I", "I_R"]].astype(float)

	func = lambda x: float(x.strip("<>"))
	ctrl, test = map(func, (row.CONTROL, row.EC50))

	susc = "S" if test < SI else "R" if test >= IR else "I"
	valid = "valid" if ctrl < SI else "inv"

	return [index, susc, valid]

################################################################################