
import numpy as np
import operator as op
import pandas as pd

from utils import gU, g2pU

refseqs = dict(gU.fasta_parser(f"{g2pU.data_dir}/static/ref_seqs.fas"))
################################################################################

def map_seqs(df):

	if df.empty: return pd.DataFrame()

	print("\tMapping with BWA")

	fas = gU.df2fas(df.reset_index(), "tmp.fas", cols=["index", "SEQ"])

	pc = gU.bwa_mem("-t1", "-O12", target=f"{g2pU.data_dir}/static/bwa", fastq=fas)
	pc = gU.samtools("view", F=8191, stdin=pc.stdout).communicate()[0]
	gU.remove(fas)

	df = pd.DataFrame(map(gU.get_sam_line, pc.strip().split("\n")))

	print(df)

	func = lambda x: str(x.RNAME).split("-")[1].split("_")
	df[["HSVTYPE", "DOMAIN"]] = df.apply(func, axis=1, result_type="expand")

	print("\tMapping complete")

	return df.set_index("QNAME")

################################################################################

def get_vars(df):

	if df.empty: return pd.DataFrame()

	print("\tFinding variants in the new FASTAs")

	data = [*df.iterrows()]
	df = pd.concat(gU.threaded(func=parse_fas, data=data, procs=16))
	print(f"\t{len(df)} variants found in {len(data)} sequences")

	return df.reset_index(drop=True)

#-------------------------------------------------------------------------------

"Takes a SAM line from a FASTA -> targets mapping, and derives variants"

def parse_fas(index, line):

	r_arr = gU.seq2arr(refseqs[line.RNAME])
	c_arr = gU.cigar2arr(line.CIGAR)
	p_arr = np.full_like(c_arr, int(line.POS) - 1, dtype=np.int32)
	p_arr += np.cumsum((c_arr != 8))
	p_arr -= (c_arr != 4).nonzero()[0][0]
	s_arr = np.zeros_like(c_arr)
	s_arr[~(c_arr == 1)] = gU.seq2arr(line.SEQ)

	############

	arr = np.stack((r_arr, np.full_like(r_arr, 15)))
	arr[1, p_arr[c_arr < 4] - 1] = s_arr[c_arr < 4]
	arr = arr[:, :-3]

	############

	arr, missing = get_missing(arr)
	arr, indels = get_indels(arr, c_arr, s_arr)
	SNPs = get_SNPs(arr)

	variants = [*gU.chain((missing, indels, SNPs))]
	if not variants: variants = ["None"]

	############

	df = pd.DataFrame(data=variants, columns=["HGVS"])
	df["HGVS"] = f"{line.HSVTYPE}.{line.DOMAIN}." + df.HGVS
	df["FAS_ID"] = index

	return df

#-------------------------------------------------------------------------------

def group_format(loci):
	"String formatting for start/end loci"

	start, end = loci
	return "_".join(map(str, sorted(set((start + 1, end)))))

#-------------------------------------------------------------------------------
"Finds missing amino acid loci, i.e. where Ns cause a >3-residue translation"

def get_missing(arr):
	"""
	Returns amino acid positions where at least one of the bases is an N and
	the translation has > 3 possible amino acids (i.e. it doesn't flag where
	a 3rd-base N leads to only 1 or 2 translated amino acids)
	"""

	def parse_missing(loci):
		arr[1, slice(*loci * 3)] = arr[0, slice(*loci * 3)]
		return f"m.{group_format(loci)}"

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
"Finds nt and (frame-aware) aa indels, pushing them to their furthest 3' loci."

def get_indels(arr, c_arr, s_arr):

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

		return f"{pc}.{group_format(loci)}del"

	############

	def parse_ins(loci):
		"Parses loci where c_arr == 8 (i.e. insertion in the CIGAR string)"

		get_ins = lambda seq: "".join(map(str, seq))

		loci = push_indel(loci, s_arr)
		# log.info(f"Insertion at {loci[0] + 1}")

		fs = ((op.sub(*loci) % 3) == 0)
		pc, func = (("c", gU.arr2seq), ("p", gU.translate_seq))[1 * fs]
		divisor = (2 * fs + 1)
		loci -= (loci[0] % divisor)
		# locus = group_format((loci // divisor)[[0, 0]])
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
"Returns all amino acid substitutions"

def get_SNPs(arr):
	"Returns non-wt amino acids at each variant locus - freq used for mixtures"

	SNPs = [
		f"p.{list(ref)[0]}{locus}{aa}"#, round(100 / len(a), 2))
		for locus, (ref, alt) in enumerate(zip(*map(gU.translate_seq, arr)), 1)
		for aa in alt - ref
	]

	return SNPs

################################################################################
