import logging
import os
import re

import itertools as it
import numpy as np
import pandas as pd

from . import gU
from data_init.g2pConstants import *
from data_init.g2pTables import *

################################################################################
"SHORT FUNCTIONS"

hgvs_regex = re.compile(r"([12]v?)\.(\w+)\.([pcm])\.([A-Z]?)" \
						r"(\d{1,4}(?:-\d{1,4})?).*?([A-Z\*]|ins[A-Z]+|del)?")

def group_format(loci):
	"String formatting for start/end loci"
	start, end = loci
	return "_".join(map(str, sorted(set((start + 1, end)))))

def find_input_files(args, regex):
	"Returns a list of either FASTA or PRA file(s) to process"

	if args.single_file: return [os.path.abspath(args.single_file)]

	regex = re.compile(regex)
	paths = [os.path.abspath(args.directory)]

	if args.recursive:
		isdir = lambda d: os.path.isdir(os.path.join(paths[0], d))
		paths.extend([
			os.path.join(paths[0], subdirectory)
			for subdirectory in filter(isdir, os.listdir(paths[0]))
		])

	return sorted(
		os.path.join(d, f)
		for d in paths
		for f in filter(regex.search, os.listdir(d))
	)

def move_to_archive(df, prefix):
	"Moves new files to the archive data store (<g2pU.data_dir>/archive)"
	def archive(row):
		return f"{data_dir}/archive/{prefix}{str(row.Index).zfill(6)}_{row.FILENAME}"

	df["ARCHIVE"] = [*map(archive, df.itertuples())]
	df = df.loc[~df.ARCHIVE.apply(gU.exists)]

	if df.empty: return log.debug("No new files to archive")
	df.loc[:, "SRC"] = df.LOCATION.str.cat(df.FILENAME, sep="/")

	log.info("Archiving input files")
	for src, dst in zip(df.LOCATION.str.cat(df.FILENAME, sep="/"), df.ARCHIVE):
		log.debug(f"\t{src} -> <archive>/{gU.filestem(dst)}")
		gU.shell("cp", src, dst)

def make_arr(df):
	"Converts a DF into a np. array"
	return df.fillna("").to_numpy().astype(str)

################################################################################
"LOGGING"

def getLog(name, scrLevel=logging.INFO, filLevel=logging.DEBUG, fileLog=True, logDir=data_dir):

	log = logging.getLogger(name)
	log.setLevel(filLevel)
	formatter = logging.Formatter(
		"%(asctime)s %(levelname)s: %(message)s", datefmt="%H:%M:%S"
	)
	scr = logging.StreamHandler()
	scr.setLevel(scrLevel)
	scr.setFormatter(formatter)
	log.addHandler(scr)

	if fileLog:
		fil = logging.FileHandler(f"{logDir}/{name}.log")
		fil.setLevel(filLevel)
		fil.setFormatter(formatter)
		log.addHandler(fil)

	return log

log = getLog("g2pU", fileLog=False)

################################################################################
"PHENOS and SEQUENCES process"

def analyse_data(src_df, table, func, datatype):
	"""
	Takes data prepared from one Table and processes it for import into a child
	Table. Returns that part of the child Table containing data from the input.
	"""
	print(src_df)
	if (tmp_df := src_df[~src_df.index.isin(table.df.PARENT_ID)]).empty:
		log.info(f"No new {datatype} analysis required")

	else:
		input(tmp_df)
		log.info(f"*-- Analysing {len(tmp_df)} {datatype}s --*")

		if table is var: tmp_df = sU.map_fasta_seqs(tmp_df)
		tmp_df = pd.concat(gU.threaded(func, tmp_df.iterrows(), 16))
		table.append(tmp_df)
		table.write()

		log.info(f"*-- Analysis of {datatype}s complete --*")

	return table.filter(("PARENT_ID", src_df.index))

################################################################################
"REPORTS"

################################################################################
"CLASSES"

#@gU.memo
class HGVS(object):
	"""
	This the object for a single variant, i.e. not those that have ambiguity in
	their format - for those, an object of the subclass Mutation is required,
	which will return a list of HGVS objects where necessary.
	"""

	def __init__(self, hgvs):
		self.hgvs = hgvs
		if not (match := hgvs_regex.search(self.hgvs)):
			raise ValueError("Not in proper HGVS format (see docs)")
		self.h, self.d, self.p, self.r, self.l, self.a = match.groups()

	def __repr__(self):
		return self.hgvs

	def literature(self):
		return lit.filter(("HGVS", self.hgvs))

	def phenotypes(self):
		fas_index = var.filter(("HGVS", self.hgvs)).PARENT_ID
		molis = fas.filter(("index", fas_index)).MOLIS
		phe_index = phe.filter(("MOLIS", molis)).index.values
		phenos = ec50.filter(("PARENT_ID", phe_index))
		return phenos


class Mutation(HGVS):
	"""
	A Mutation's only method is <get_mutations>, which returns all variants
	within its specification. In its simplest form, this list contains a single
	variant in HGVS form, but may contain all the variants within locus ranges,
	missing data and homologous sites in other HSV types, depending upon the
	submitted <hgvs> and <homology> arguments.
	For a specific mutation, i.e. one with a complete HGVS format, then it will
	be returned "as is". Where the "ref" and/or the "alt" amino acid(s) are
	missing, or the indel definition of a nucleotide mutation is missing, then
	all valid mutations within the <lookaround> and <homology> extended ranges
	are returned, i.e. only those with an entry in <lit> and/or <VARIANT>
	tables.
	"""

	def __init__(self, hgvs, homology=False, extend=0):
		super().__init__(hgvs)
		self.homology = homology
		self.extend = extend

	def __repr__(self):
		return f"{self.hgvs}, homology={self.homology}, range={self.extend}"

	def __call__(self):
		"Takes the raw HGVS and processes for range, missing data, and homology"

		self.muts = []
		self.pc = "c" if self.p == "c" else "p"
		if not self.a is None: return [self.hgvs]



		self.l2 = self.l.split("-")[-1]
		self.l, self.l2 = map(int, (self.l, self.l2))

		params = [(self.h, self.l, self.l2)]
		if self.homology:
			params.append((3 - int(self.h[0]), *self.get_homologous_loci()))

		print("params:", params)

		hits = list(gU.chain(it.starmap(self.get_potential_muts, params)))
		print("HITS:", hits)
		input(f"self.a is None - HGVS={self.hgvs}")
		return hits

	def get_potential_muts(self, h, l, l2):

		loci = "|".join(map(str, range(l-self.extend, l2+self.extend+1)))
		hits_regex = f"{h}v?\.{self.d}\.{self.pc}\.[A-Z]?({loci})[idA-Z\*]"

		print(f"h: {h}, l: {l}, l2: {l2}")
		print("hits_regex:",hits_regex)

		def filter_table(table):
			ser = table.df.HGVS.str.extract(hits_regex).dropna()
			df = table.df.join(ser, how="right")
			return df

		df = pd.concat(map(filter_table, (var, lit)))
		print(df)
		hits = pd.unique(df.sort_values(0).HGVS)

		input(hits)

		return hits

	def get_homologous_loci(self):

		aln = f"{data_dir}/static/{self.d}.{self.pc}.aln"
		fastas = sorted([*gU.fasta_parser(aln)], key=lambda x: x[0] == self.h)
		is_base = np.stack([
			np.cumsum([*map(lambda x: x != "-", fasta[1])])
			for fasta in fastas
		])

		func = lambda x: np.where(is_base[-1]==x)[0][0]

		return is_base[0][[*map(func, (self.l, self.l2))]]

################################################################################
