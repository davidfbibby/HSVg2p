"""Functions that do the heavy lifting for the PHENOS module"""

import os
import xlrd

import itertools as it
import numpy as np
import pandas as pd

from datetime import datetime as dt
from functools import partial, reduce

from . import gU, g2pU
from data_init.g2pTables import ec50, tgt

################################################################################
def parse_PRA(index, row):

	path = os.path.join(row.LOCATION, row.FILENAME)

	sheets = pd.read_excel(path, sheet_name=None)
	sheets.pop("Results", None)

	func = parse_old_PRA if "ACV" in sheets else parse_new_PRA

	if not (df := func(sheets)).empty:
		df["PARENT_ID"] = index
		if (date := min(df.pop("DATE"))) is pd.NaT:
			date = dt.fromtimestamp(os.path.getmtime(path))
		df["DATE"] = pd.to_datetime(date, format="%Y-%m-%d")

	return df

################################################################################
def parse_old_PRA(sheets, hsvtype=1):
	def parse_drug_sheet(drug, sheet):
		"Parses the old-style PRA worksheets - each sheet is one drug"

		drug = drug.replace("PEN", "PCV")
		index = ec50.cols[:4]

		if drug not in tgt.df.index: return pd.DataFrame(index=index)

		arr = g2pU.make_arr(sheet)

		date_cells = find_cells(arr, "Date", 1, dtype=str)
		date = pd.NaT if not date_cells else parse_timestr(date_cells[0][0])

		x_arr = find_cells(arr, "No Drug", 1, 0, 3, 3)[[0, 2]]
		p_arr = find_cells(arr, "Plaque", 0, 2, 3, 5)[[0, 2]] / get_xDrug(x_arr)
		c_arr = find_cells(arr, "Plaque", -1, 2, 1, 5)[0]

		ctrl, test = map(partial(calculate_EC50, c_arr), p_arr)
		ser = pd.Series(data=[date, drug, ctrl, test], index=index)

		return ser

	df = pd.concat(it.starmap(parse_drug_sheet, sheets.items()), axis=1).T

	return df

################################################################################
def parse_new_PRA(sheets, hsvtype=1):
	"Parses the new-style PRA worksheets"

	def parse_drug(col):
		indexes = find_txt(arr, drug := str(arr[0, col])).T[0]

		try: p_arr = np.stack([arr[i + 5:i:-1, col] for i in indexes]).astype(np.float64)
		except Exception: return

		return [drug, *map(partial(calculate_EC50, p_arr[2]), p_arr[:2] / xDrug)]

	if (sheet := sheets.get("HSV PRA")) is None: return pd.DataFrame()

	arr = g2pU.make_arr(sheet)

	date_cells = find_cells(arr, "ixed", 0, 2, 1, 4, dtype=str)[0]
	date = "" if (x := date_cells.nonzero()[0]).size == 0 else date_cells[x[-1]]

	xDrug_indexes = find_txt(arr, "xDrug")
	x_arr = np.stack([arr[r + 1:, c][:3] for r, c in xDrug_indexes])
	if np.any(~np.char.isnumeric(x_arr)): return pd.DataFrame()

	xDrug = get_xDrug(x_arr)
	r, c = xDrug_indexes[0]
	arr = arr[r:, :c]

	data = filter(lambda x: x is not None, map(parse_drug, arr[0].nonzero()[0]))
	df = pd.DataFrame(data=data, columns=ec50.cols[1:4])
	df["DATE"] = parse_timestr(date)

	return df

################################################################################
find_txt = lambda a, t: np.stack(np.where(np.char.find(a, t) != -1)).T

#-------------------------------------------------------------------------------
def parse_timestr(timestr):
	"Parses Excel and strings into datetimes"

	if timestr in("0", "nan"): return pd.NaT
	if timestr.isnumeric(): return xlrd.xldate_as_datetime(int(timestr), 0)
	try:
		return pd.to_datetime(timestr, dayfirst=True)
	except:
		return pd.NaT

#-------------------------------------------------------------------------------
def find_cells(arr, text, row_os=0, col_os=0, row_n=1, col_n=1, dtype=float):
	"Finds the first cell from find_txt and returns an offset range of cells"

	if (index := find_txt(arr, text)).size == 0: return []
	row, col = index[0] + [row_os, col_os]
	return arr[row:, col:][:row_n, :col_n].astype(dtype)

#-------------------------------------------------------------------------------
def get_xDrug(x_arr):
	"Converts string plaque counts from no-drug controls into a divideable array"

	return np.expand_dims(np.mean(x_arr.astype(int), axis=1), 1)

#-------------------------------------------------------------------------------
def calculate_EC50(c_arr, p_arr):
	"Calculates the EC50 from matching plaque% & conc arrays (ascending)"

	index = (p_arr >= 0.5).nonzero()[0]
	if index.size == 0: return f"<{c_arr[-1]}"
	if (index := index[0]) == 0: return f">{c_arr[0]}"
	c_arr = np.log(c_arr)
	ratio = (p_arr[index] - 0.5) / (p_arr[index] - p_arr[index - 1])
	conc = c_arr[index] - ((c_arr[index] - c_arr[index - 1]) * ratio)
	return np.format_float_positional(np.exp(conc), 4, fractional=False)

################################################################################
def statistics():

	pass

################################################################################