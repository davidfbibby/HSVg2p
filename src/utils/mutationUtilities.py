"""Functions that do the heavy lifting for the PHENOS module"""

import os

import itertools as it
import numpy as np
import pandas as pd

from functools import partial, reduce

from . import gU, g2pU

################################################################################
class Mutation():
	"""
	A Mutation's only method is <get_mutations>, which returns all variants
	within its specification. In its simplest form, this list contains a single
	variant in HGVS form, but may contain all the variants within locus ranges,
	missing data and homologous sites in other HSV types, depending upon the
	submitted <hgvs> and <homology> arguments.
	"""

	regex = r"([12]v?)\.(\w+)\.([pcm])\.([A-Z]?)" \
			r"(\d{1,4}(?:-\d{1,4})?).*?([A-Z]|ins|del)?"

	def __init__(self, hgvs, homology, lookaround):
		self.hgvs = hgvs
		self.homology = homology
		self.lookaround = lookaround

	def __repr__(self):
		return f"{self.hgvs}, homology={self.homology}, range={self.lookaround}"

	def get_mutations(self):
		"Takes the raw HGVS and processes for range, missing data, and homology"

		if not (match := re.search(Mutation.regex, self.hgvs)): return []
		self.h, self.d, self.p, self.r, self.l, self.a = match.groups()
		self.l2 = self.l.split("-")[-1]

		if self.p == "m": self.missing()

	def missing(self):
		for locus in range(self.l, self.l2 + 1):
			df = lit.filter((""))



	def add_homology(self):
		aln = f"{data_dir}/static/{self.domain}.{self.pc}.aln"

		fastas = sorted([*gU.fasta_parser(aln)], key=lambda x: x[0] == self.hsvtype)
		ref_seq = fastas.pop()[1]
		pos = np.where([*map(lambda x: x != "-", ref_seq)])[0][args.locus]
		self.homolog = get_homolog(pos)
		def get_homolog(pos):
			if self.ref and ref_seq[:pos][-1] != self.ref: return
			hsvtype = str(3 - int(self.hsvtype[0]))
			mut = "".join((self.ref, str(pos), self.alt))
			return Mutation(".".join((hsvtype, self.domain, self.pcm, mut)))

		def func(n, s): return (n, (s := s[:pos])[-1], len(s.replace("-", "")))
		return [*it.starmap(func, fastas)]

################################################################################

################################################################################