import collections
import logging
import os
import re
import sys
import xlrd
import yaml

import itertools as it
import pandas as pd

from collections import defaultdict
from functools import partial
from glob import glob

from utils import gU

"CONSTANTS"

src = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir)
data_dir = os.path.join(src, os.pardir, "data")
archive_dir = os.path.join(data_dir, "archive")
src, data_dir = map(os.path.abspath, (src, data_dir))
res_dict = {
	"R": "RESISTANT", "R*": "LIKELY RESISTANT", "R?": "POSSIBLY RESISTANT",
	"S": "SENSITIVE", "S*": "LIKELY SENSITIVE", "?": "AMBIGUOUS",
	"H": "HYPER-SENSITIVE", "H*": "LIKELY HYPER-SENSITIVE"
}

################################################################################
"PHENOS and SEQUENCES process"

def analyse_data(df, table, func, datatype, inner_func=None):

	df = df[~df.index.isin(table.df.PARENT_ID)]
	if df.empty: return log.info(f"No new {datatype} analysis required")

	log.info(f"*-- Analysing {len(df)} {datatype}s --*")

	if inner_func: df = inner_func(df)
	df = pd.concat(gU.threaded(func, df.iterrows(), 16))
	index, df = table.append(df)
	table.write()

	log.info(f"*-- Analysis of {datatype}s complete --*")

	return df


################################################################################
class Table(object):
	"""
	The Table class wraps a pandas DataFrame in its <df> attribute. In addition
	to being able to support standard operations through <df>, there are class
	properties and methods. Table is subclassed into DynamicTable and
	StaticTable to reflect the data and reference Tables, respectively.

	<name> and <location> combine to give <fname>. <location> specifies the sub-
	folder within /HSVg2p/data in which the tab-separated flatfile containing
	the Table data is stored. <cols> stores the column names and is useful for a
	number of operations. **kwargs are passed to the <read_tsv> module function
	to allow subclasses to specify behaviours (useful when the file doesn't yet
	exist).

	methods
	-------
	filter	Takes a list of column-value (or index-value) pairs and gets the
			sub-DataFrame where all conditions are true, i.e. those rows where
			each column contains a member of the list of corresponding values.
			Returnvalue can be the sub-DF or its indexes. The filtering
			operation can be inverted (i.e. where none of the values are in the
			columns), and the condition of all having to be true can be changed
			by passing a function. The default is set.intersection. For any of
			the filters to be true, then pass set.union.

	"""
	def __init__(self, name, location, **kwargs):

		self.name = name
		self.location = location
		self.fname = f"{data_dir}/{self.location}/{name}.tsv"

		kwargs.setdefault("index_col", 0)
		self.df = gU.read_tsv(self.fname, **kwargs)
		self.cols = self.df.columns

	def __repr__(self):
		return self.df.to_string()

	def filter(self, filters, inverse=False, index=False, setop=set.intersection):
		"""
		Returns the sub_df where all <filters> are satisfied. <filters> is a
		list of (col, values) pairs (<col> can be "index"), returning indices
		where col == values. The <inverse> argument defaults to all False, but
		can be a tuple of the same length as <filters>, specifying whether a
		filter returns the inverse or not.
		(The intersection of all indices...)
		"""
		def single_filter(col_val, inv):
			col, val = col_val
			val = gU.parse_input_data(val)

			if col == "index":
				if isinstance(self.df.index, pd.MultiIndex):
					val = [*gU.iter_zip(self.df.index.nlevels, val)]
				boolean = self.df.index.isin(val)
			else: boolean = self.df[col].isin(val)
			if inv: boolean = ~boolean
			return set(self.df[boolean].index)

		if self.df.empty: return self.df
		if type(filters[0]) != tuple: filters = (filters, )
		if type(inverse) == bool: inverse = it.repeat(inverse, len(filters))

		indices = sorted(
			setop(*it.starmap(single_filter, zip(filters, inverse)))
		)
		if index: return indices
		df = self.df.loc[indices]
		return df

#-------------------------------------------------------------------------------
class DynamicTable(Table):
	"""
	A subclass of Table, with additional __init__ actions and new functions:
	<delete>	allows propagation into "child" Tables so that all data relating
				to an entry is removed from all files.
	<append>	merges a new DataFrame with the existing Table, thus precluding
				duplicates. Returns the new DF but with revised indexes, plus
				the indexes of new entries only.
	<write>		saves the Table "as-is" to file. Not automatic, so needs to be
				explicit in the script.
	"""

	def __init__(self, name, cols, child=[]):
		self.child = child
		super().__init__(name, location="dynamic", cols=cols)

		if "DATE" in self.df.columns:
			self.df.DATE = pd.to_datetime(self.df.DATE, format="%Y-%m-%d")
		if "PARENT_ID" in self.df.columns:
			self.df = self.df.astype({"PARENT_ID": int})

	def append(self, df, *,
			   fill="", sort_cols=[], reset=False):
		"""
		Appends <df> data to a Table.df, checking (and removing) duplicates via
		pd.merge, and filling spaces with <fill>. The following actions are
		performed in order (i.e. reset will follow a sort):
		- sort		Sorts the Table.df by the columns in <sort>
		- reset		Resets the index after concatenation
		- write		Finally, the table is saved to file by default

		returns:
		- index

		The <result> variable can take one of two values, determining how much
		of the input DataFrame is returned. Both returned DFs are re-indexed
		with the corresponding indexes of the the updated Table.:
		["all"]		Returns the entire input DF.
		"new"		Returns only that part of the input DF that was not already
					present in the Table prior to the update.
		"""

		log.debug("Appending new data to existing Table")
		index = self.df.index
		self.df = self.df.merge(df, how="outer").fillna(fill)

		if sort_cols: self.df = self.df.sort_values(sort_cols)
		if reset: self.df = self.df.reset_index(drop=True)

		dups = pd.concat((self.df, df)).duplicated(keep="last")
		df = self.df.loc[dups[dups].index]
		index = df.index.difference(index)

		return index, df

	def delete(self, indexes):

		if (child := self.child) is not None:
			indexes = child.df[child.df.PARENT_ID.isin(indexes)].index.values
			child.delete(indexes)
		self.df = self.df[~self.df.index.isin(indexes)]
		self.write()

	def write(self):
		gU.write_tsv(self.df, self.fname, index=True)

#-------------------------------------------------------------------------------
class StaticTable(Table):
	"""
	A subclass of Table, with additional __init__ actions, namely the ability to
	make an index from selected columns (i.e. not the default [0]), and "" marks
	can be added to the read_csv <comment> parameter for the Citation table.
	"""

	def __init__(self, name, **kwargs):
		super().__init__(name, location="static", **kwargs)

#-------------------------------------------------------------------------------
class LiteratureTable(StaticTable):
	"""
	A subclass of StaticTable, specifically for Tables with geno-to-pheno info
	gleaned from the literature. Expands the <__init__> function to parse the
	HGVS.
	"""

	def __init__(self, name, **kwargs):
		super().__init__(name, **kwargs)
		print(self.name)
		print(self.df)
		if not self.df.empty:
			df = pd.DataFrame.from_records(
				self.df.index.map(parse_HGVS),
				index=self.df.index,
				columns=list("HDPVL")
			)
			self.df = pd.concat((self.df, df), axis=1)
			input(self.df)

#-------------------------------------------------------------------------------
def parse_HGVS(hgvs):

	parsed = hgvs.split(".")
	parsed.append(int(re.search(r"\d{1,4}", parsed[-1]).group(0)))

	return parsed

"SET UP TABLES FROM DATA TSVs"

# DynamicTables - these get updated with new data
#	mol		MOLIS IDs, with HSV types - True for present in sample.
#			Also maps older MOLIS IDs to HSV type, based upon an MMD download
# 	fil		FILES, with source RUN ID, date of run etc.
#	fas		FASTA sequences, with an ID link (child relationship) to <fil>
#	var		Variants detected in FASTAs, with a child relationship to <fas>
#	phe		PRA files
#	sir		Susceptibilities from PRAs, with a child relationship to <phe>

var = DynamicTable(
	"VARIANTS", cols=["HGVS", "PARENT_ID"]
)

fas = DynamicTable(
	"FASTAS", child=var, cols=["NAME", "MOLIS", "SEQ", "PARENT_ID"]
)

fil = DynamicTable(
	"FILES", child=fas, cols=["LOCATION", "FILENAME", "DATE", "RUNID"]
)

ec50 = DynamicTable(
	"EC50S",  cols=["DATE", "DRUG",  "CTRL", "EC50", "PARENT_ID"]
)

phe = DynamicTable(
	"PHENOS", child=ec50, cols=["LOCATION", "FILENAME", "MOLIS"]
)

mol = DynamicTable(
	"MOLIS", cols=["MOLIS", "1", "2", "2v"]
)

# StaticTables - these contain fixed reference information
#	thr		Configures the Susceptible/Intermediate/Resistant PRA cut-offs
#   raw		Table in which each row is a variant/drug/citation combination
#	cit		Table of the citations used in <raw>
#	lit		Table of reported variants and their impact upon drug
#			susceptibilities. An automated summary of <raw> generated using
#			<summarise_raw>
#   res		Table of master interpretations, curated and amended in the light
#			of data arising during routine work. Is the ultimate arbiter of
#			g2p, combining literature, phenotyping and other sources. Overrules
#			both <raw> and <lit>.

thr  = StaticTable("THRESHOLDS", index_col=[0, 1])
raw = StaticTable("LITERATURE", index_col=[0, 1, 2])
cit = StaticTable("CITATIONS", index_col=0, comment="\"")
tgt = StaticTable("TARGETS", index_col=0)
lit, res = map(LiteratureTable, ("INTERPRETATIONS", "RESOLVED"))
