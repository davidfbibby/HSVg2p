# HSV genotype and phenotype database and query tool  

Author: David F Bibby, UKHSA  

## Overview  

The Antiviral Unit (AVU) at UKHSA routinely tests samples from patients infected with HSV, with both genotypic and phenotypic assays being available. The former makes use of a body of published literature to determine the presence or absence of resistance-associated variants in target gene sequences. The latter compares experimentally determined EC50 values to tables of thresholds to determine susceptibility to a range of drugs. Notably, upon completion of a sample's testing requirements and all analysis has been reported on the LIMS, the data is not processed further. It is stored on remote drives, specifically either the AVU space on FILECOL15 or on the HPC cluster (mainly FASTQ files and the outputs of their bioinformatic analysis).  
This package intends to create a database of both sequence and phenotype data, with on-demand reporting, querying, and analysis. For the database component, the top level tables store information about the raw data emitted from service assays, parsing essential identifiers distinguishing it from similar outputs (e.g. MOLIS no, date, run ID, etc.). Although it sits above the other layers, analysis without importation is possible. A second layer analyses the raw data entering the primary layer, returning lists of variants from sequences and/or EC50s from phenotyping. A reporting layer extracts information from the second layer, interprets it through a series of static reference tables, and returns collated information regarding the properties of the queried unit. A unit in this sense can be a FASTA sequence, a phenotyping dataset, a MOLIS ID, or a specific mutation. Depending on the type of information and passed arguments, this can range from succinct to highly verbose. Formats suitable for clinical reporting can be configured for routine use.  

---

## How to use  

All functions are accessed via a single top-level module - `HSVgeno2pheno.py`. Six subparsers are used to direct the analysis to first-level modules. Each has a range of parameters that can be either declared on the command line or passed via a configuration file containnig workflow-specific arguments. The last two commands modify the tables in-place, and should be restricted to authorised users only.  

| ID | Module | Description | Main configurables |
| :--- | :--- | :--- | :--- |
| 1 | `SEQUENCES` | Processes SEQUENCE data | ?Import, ?Report, ?NGS |
| 2 | `PHENOS` | Processes PHENOTYPIC data | ?Report all MOLIS |
| 3 | `MOLIS ID` | Reports on a MOLIS ID | ?Variants, ?Susceptibilities |
| 4 | `MUTATION` | Reports on a mutation (HGVS format) | |
| 5 | `MODIFY` | Modifies a STATIC table | ?check, ?password |
| 6 | `SUPPRESS` | Flags DYNAMIC table data as invalid | ?check, ?password |

---

## Database Table structure

There are two top layer tables - FASTAS and PHENOS. Each contains broad metadata about either a single FASTA sequence or a single phenotypic assay. Each entry points to primary data held in a *processed* folder; neither the primary sequences nor the phenotypic assay data are held in the database. It is not necessary for a sequence to be imported into the database/processed folder - variant analysis can be performed independently.  
There are two second layer tables - VARIANTS and EC50S. The former links sequence variants to specific FASTAs. For Sanger sequences, the variant information derives from the outputs of a module that takes a single sequence as input, aligns it to references and returns the difference(s) in a defined format based upon the HGVS notation. Each variant has four components, separated by periods - the type of HSV (1, 2 or 2v), the domain (TK, pol, UL5, UL52), the type of variant ("p" for amino acid change, "c" for nucleotide frameshift indel, or "m" for missing amino acid loci), and the specifics of the variant. For example, a common variant in HSV-1 TK that confers resistance to ACV is an insertion of a G in the homopolymeric tract at nucleotide positions 430-436. This would be encoded with 1.TK.c.436insG. A deletion at the same position - 1.TK.c.436del - also confers ACV resistance. In HSV-2 pol, the substitution of serine for asparagine at amino acid position 729 is encoded thus: 2.pol.p.S729N.  
Whilst most FASTQ analysis will result in a single FASTA file that can be analysed as per Sanger inputs, occasionally, sequence emerges with multiple frameshift indels at various frequencies, often precluding the generation of a single FASTA that can correctly identify them all. Here, the VARIANTS table can be directly populated from the QuasiBAM outputs, ensuring all variants are represented. In both Sanger and NGS analyses, mixtures can be seen - a frequency column in the table captures this information.  
The EC50S table is derived from the phenotypic assay files. These typically comprise a single file per sample per test, and contain data about the EC50s for one or more drugs. Because the validity of an assay output is dependent upon the data for contemporary control assays, these outputs are included in this table. Each row thus contains control and sample EC50 data calculated from a single assay file and for a single drug (multi-drug assays therefore have a line for each drug). As with the VARIANTS, susceptibility is not determined at the point of this table being populated. Rather the interpretations are generated at the time of reporting, using the most up-to-date reference information.  

## Static Tables  

A number of reference datasets support the analysis and reporting of HSV data:  
1.	Reference sequence location: For processing sequences into variants, a reference sequence for each HSV type/domain combination is needed. This table 