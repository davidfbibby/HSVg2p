"""Functions that do the heavy lifting for the SEQUENCES module"""

import os

import itertools as it
import numpy as np
import operator as op
import pandas as pd

from functools import partial, reduce

from . import gU, g2pU

################################################################################
def parse_FASTA(index, row):

	path = os.path.join(row.LOCATION, row.FILENAME)
	if (df := gU.fas2df(path)).empty: return pd.DataFrame()

	df["SEQ"] = df.SEQ.str.replace("-", "")
	df["MOLIS"] = df.NAME.apply(gU.MOLIS_name, spaces=True)
	df["PARENT_ID"] = index
	df = df[df.MOLIS.str.match(r"([Hh]\d{9}|RS\d{8})")]

	return df

################################################################################
def map_fasta_seqs(df):

	"BWA maps sequences against the reference index and returns SAM as a <df>"

	if df.empty: return pd.DataFrame()
	tmp = gU.df2fas(df.reset_index(), "tmp.fas", cols=["index", "SEQ"])

	pc = gU.bwa_mem("-t1", "-O12", target=f"{g2pU.data_dir}/static/bwa", fastq=tmp)
	pc = gU.samtools("view", F=8191, stdin=pc.stdout).communicate()[0]
	gU.remove(tmp)

	df = pd.DataFrame(map(gU.get_sam_line, pc.strip().split("\n")))
	if df.empty: return df

	func = lambda x: str(x.RNAME).split("-")[1].split("_")
	df[["HSV", "DOMAIN"]] = df.apply(func, axis=1, result_type="expand")
	df = df.drop(["FLAG", "MAPQ", "RNEXT", "PNEXT", "TLEN", "QUAL"], axis=1)
	df = df.set_index("QNAME")
	df.index = df.index.astype(int)

	return df

################################################################################
def parse_variants(index, row):

	r_arr = gU.get_seq(row.RNAME, f"{g2pU.data_dir}/static/ref_seqs.fas")
	c_arr = gU.cigar2arr(row.CIGAR)
	p_arr = np.full_like(c_arr, int(row.POS) - 1, dtype=np.int32)
	p_arr += np.cumsum((c_arr != 8))
	p_arr -= (c_arr != 4).nonzero()[0][0]
	s_arr = np.zeros_like(c_arr)
	s_arr[~(c_arr == 1)] = gU.seq2arr(row.SEQ)

	############

	arr = np.stack((r_arr, np.full_like(r_arr, 15)))
	arr[1, p_arr[c_arr < 4] - 1] = s_arr[c_arr < 4]
	arr = arr[:, :-3]

	############

	arr, missing = get_missing(arr)
	arr, indels = get_indels(arr, c_arr, s_arr)
	SNPs = get_SNPs(arr)

	variants = [*gU.chain((missing, indels, SNPs))]
	if not variants: variants = [f"{row.HSV}.None"]

	############

	df = pd.DataFrame(data=variants, columns=["HGVS"])
	df["HGVS"] = f"{row.HSV}.{row.DOMAIN}." + df.HGVS
	df["PARENT_ID"] = index

	return df

#-------------------------------------------------------------------------------
def get_missing(arr):
	"""
	Finds missing amino acid loci, i.e. where Ns cause a >3-residue translation.
	Returns amino acid positions where at least one of the bases is an N and
	the translation has > 3 possible amino acids (i.e. it doesn't flag where
	a 3rd-base N leads to only 1 or 2 translated amino acids)
	"""

	def parse_missing(loci):
		arr[1, slice(*loci * 3)] = arr[0, slice(*loci * 3)]
		return f"m.{g2pU.group_format(loci)}"

	############

	@gU.memo
	def check_codon(start):
		return len(gU.translate_mixed(arr[1, start * 3:][:3])) > 3

	############

	missing_nt = (arr[1] == 15).nonzero()[0]
	missing_loci = [*filter(check_codon, np.unique(missing_nt // 3))]
	arr[1, missing_nt] = arr[0, missing_nt]
	missing = [*map(parse_missing, gU.groups(missing_loci))]

	return arr, missing

#-------------------------------------------------------------------------------
def get_indels(arr, c_arr, s_arr):
	"Finds nt and (frame-aware) aa indels, pushing them to their furthest 3' loci."

	def push_indel(loci, t_arr):
		while op.eq(*t_arr[loci]):
			loci += 1
		return loci

	############

	def parse_del(loci):
		"Parses loci where sample arr == 0 (i.e. base in reference)"

		# Overwrite deleted position with reference base to allow <get_SNPs>
		arr[1, slice(*loci)] = arr[0, slice(*loci)]
		loci = push_indel(loci, arr[0])
		pc = "c"

		if (op.sub(*loci) % 3) == 0:
			pc = "p"
			loci = (loci - loci[0] % 3) // 3

		return f"{pc}.{g2pU.group_format(loci)}del"

	############

	def parse_ins(loci):
		"Parses loci where c_arr == 8 (i.e. insertion in the CIGAR string)"

		get_ins = lambda seq: "".join(map(str, seq))

		loci = push_indel(loci, s_arr)

		fs = ((op.sub(*loci) % 3) == 0)
		pc, func = (("c", gU.arr2seq), ("p", gU.translate_seq))[1 * fs]
		divisor = (2 * fs + 1)
		loci -= (loci[0] % divisor)
		locus = loci // divisor
		ins = [*map(get_ins, it.product(*func(s_arr[slice(*loci)])))]

		return [f"{pc}.{locus[0]}ins{_ins}" for _ins in ins]

	############

	indels = [
		*map(parse_del, gU.groups((arr[1] == 0).nonzero()[0])),
		*gU.chain(map(parse_ins, gU.groups((c_arr == 8).nonzero()[0])))
	]

	return arr, indels

#-------------------------------------------------------------------------------
def get_SNPs(arr):
	"Returns non-wt amino acids at each variant locus - freq used for mixtures"

	SNPs = [
		f"p.{list(ref)[0]}{locus}{aa}"
		for locus, (ref, alt) in enumerate(zip(*map(gU.translate_seq, arr)), 1)
		for aa in alt - ref
	]

	return SNPs

################################################################################
def statistics():
	"Returns statistics of the three sequence tables - FILES, FASTAS & VARIANTS"

	log.info("Table statistics")
	ser = g2pU.fil.df.DATE.dt.year
	df = pd.DataFrame.from_dict(
		{year: [len(df), len(g2pU.fas.filter(("PARENT_ID", df.index)))]
		 for year, df in g2pU.fil.df.groupby(ser)},
		"index", columns=["FILES", "FASTAS"]
	)
	df.loc["Total"] = df.sum()
	log.info(df)
	log.info(f"{g2pU.fas.df.MOLIS.nunique()} MOLIS IDs")
	log.info("Most recent addition(s)")
	[*map(log.info, g2pU.fil.filter(("DATE", max(g2pU.fil.df.DATE))).FILENAME.values)]