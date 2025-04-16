# HSV genotype and phenotype database and query tool  

Author: David F Bibby, UKHSA  

## Overview  

The Antiviral Unit (AVU) at UKHSA routinely tests samples from patients infected with HSV, with both genotypic and phenotypic assays being available. The former makes use of a body of published literature to determine the presence or absence of resistance-associated variants in target gene sequences. The latter compares experimentally determined EC50 values to tables of thresholds to determine susceptibility to a range of drugs. Notably, upon completion of a sample's testing requirements and all analysis has been reported on the LIMS, the data is not processed further. It is stored on remote drives, specifically either the AVU space on FILECOL15 or on the HPC cluster (mainly FASTQ files and the outputs of their bioinformatic analysis). 
 
This package intends to create a database of both sequence and phenotype data, with on-demand reporting, querying, and analysis. For the database component, the top level tables store information about the raw data emitted from service assays, parsing essential identifiers distinguishing it from similar outputs (e.g. MOLIS no, date, run ID, etc.). Although it sits above the other layers, analysis without importation is possible. A second layer analyses the raw data entering the primary layer, returning lists of variants from sequences and/or EC50s from phenotyping. A reporting layer extracts information from the second layer, interprets it through a series of static reference tables, and returns collated information regarding the properties of the queried unit. A unit in this sense can be a FASTA sequence, a phenotyping dataset, a MOLIS ID, or a specific mutation. Depending on the type of information and passed arguments, this can range from succinct to highly verbose. Formats suitable for clinical reporting can be configured for routine use.  

---

## How to use  

All functions are accessed via a single top-level module - `HSVgeno2pheno.py`. Six subparsers are used to direct the analysis to first-level modules. Each has a range of parameters that can be either declared on the command line or passed via a configuration file containnig workflow-specific arguments. The last two commands modify the tables in-place, and should be restricted to authorised users only.  

| ID | Subparser (and module) | Description | Main configurables |
| :--- | :--- | :--- | :--- |
| 1 | [`sequences`](#sequences) | Processes SEQUENCE data | ?Import, ?Report, ?NGS |
| 2 | [`phenos`](#phenos) | Processes PHENOTYPIC data | ?Report all MOLIS |
| 3 | [`molis`](#molis) | Reports on a MOLIS ID | ?Variants, ?Susceptibilities |
| 4 | [`mutation`](#mutation) | Reports on a mutation (HGVS format) | ?loci range |
| 5 | `modify` | Modifies a STATIC table | ?check, ?password |
| 6 | `suppress` | Flags DYNAMIC table data as invalid | ?check, ?password |

<a name="workflows">All subparsers can be configured with the `-w` argument, specifying a "workflow" configuration. within `config/` are `yaml` files for each subparser specifying combinations of arguments for analysis shortcuts. Workflow arguments override defaults, and command line arguments override workflow arguments.

---

## <a name="DTs">Database Table structure

There are two top layer tables - FASTAS and PHENOS. Each contains broad metadata about either a single FASTA sequence or a single phenotypic assay. Each entry points to primary data held in `data/archive`; neither the primary sequences nor the phenotypic assay data are held in the database. It is not necessary for a sequence to be imported into the `data/archive` folder - variant analysis can be performed independently.  

There are two second layer tables - VARIANTS and EC50S. The former links sequence variants to specific FASTAs. For Sanger sequences, the variant information derives from the outputs of a module that takes a single sequence as input, aligns it to references and returns the difference(s) in a defined format based upon the HGVS notation. Each variant has four components, separated by periods - the type of HSV (1, 2 or 2v), the domain (TK, pol, UL5, UL52), the type of variant ("p" for amino acid change, "c" for nucleotide frameshift indel, or "m" for missing amino acid loci), and the specifics of the variant. For example, a common variant in HSV-1 TK that confers resistance to ACV is an insertion of a G in the homopolymeric tract at nucleotide positions 430-436. This would be encoded with 1.TK.c.436insG. A deletion at the same position - 1.TK.c.436del - also confers ACV resistance. In HSV-2 pol, the substitution of serine for asparagine at amino acid position 729 is encoded thus: 2.pol.p.S729N.  

Whilst most FASTQ analysis will result in a single FASTA file that can be analysed as per Sanger inputs, occasionally, sequence emerges with multiple frameshift indels at various frequencies, often precluding the generation of a single FASTA that can correctly identify them all. Here, the VARIANTS table can be directly populated from the QuasiBAM outputs, ensuring all variants are represented. In both Sanger and NGS analyses, mixtures can be seen - a frequency column in the table captures this information.  
The EC50S table is derived from the phenotypic assay files. These typically comprise a single file per sample per test, and contain data about the EC50s for one or more drugs. Because the validity of an assay output is dependent upon the data for contemporary control assays, these outputs are included in this table. Each row thus contains control and sample EC50 data calculated from a single assay file and for a single drug (multi-drug assays therefore have a line for each drug). As with the VARIANTS, susceptibility is not determined at the point of this table being populated. Rather the interpretations are generated at the time of reporting, using the most up-to-date reference information.  

## Static Tables  

A number of reference datasets support the analysis and reporting of HSV data:  
| ID | Name | Description |
| :--- | :--- | :--- |
| 1 | Literature | A table of genotype-to-phenotype correlations derived from literature sources. Each row is a single source-variant-drug combination, together with its documented phenotype |
| 2 | References | A table containing the sources from which the Literature table (#1) was obtained |
| 3 | Thresholds | A table of EC50 cutoffs for each drug-HSV type combination |
| 4 | Targets | A table with domains as rows and drugs as columns and the relation between the two. 0: Drug does not target domain; 1: (Secondary) drug targets domain; 2: (Primary) drug targets domain |

---

## <a name="sequences">sequences

To process sequence data, use the `sequences` subparser. Major options are as follows:

| Argument | Description | Default |
| :--- | :--- | :--- |
| [`-i` / `--import_data`](#sequences_import) | Set to import FASTAs into the FASTAS top-layer table | |
| [`-r` / `--report_data`](#sequences_report) | Set to report the interpretation of each FASTA | |
| `-f` / `--fasta` | `path/to/fasta/file` | |
| `-d` / `--directoy` | `path/to/directory/containing/fastas` | Current directory |
| `--runid` | name of the run (e.g. LIMS tracking ID) | Inferred from directory name |
| `--rundate` | date of the run | Inferred from directory name |
| `--force` | Set to overwrite existing data in the FASTAS table | |

`-i` and `-r` apply to all sequences processed by a single `HSVgeno2pheno.py sequences` command. At least one must be present. If `-r` is passed, its argument will control the report format (can be part of a [configured workflow](#workflows)).

The `-f` and `-d` flags are mutually exclusive. If `-d` is selected, each FASTA file in the directory is processed sequentially, each in the same way as if it had been selected singly using the `-f` flag. Each sequence within a FASTA file is processed in turn.

 It is assumed that most sequences will originate from a directory named according to a loose convention that implies parseability by a regex. If this is not the case, and the `-i` flag is set, then the `--runid` and `--rundate` arguments can be used to manually specify these two labels for the sequence data. 

`--force` processes an imported sequence, even if it has been imported previously. This allows a re-analysis generating a new "truth" to be established.

### <a name="sequences_import">`-i` / `--import_data`

If import is selected, metadata about the sequence is extracted, namely the filename containing the sequence, the run ID, and the run date. (The latter two labels are either taken from the passed arguments or inferred from the directory name using regex matching.) If an entry corresponding to these labels is already present in the database and the `--force` flag is not set, the import process halts, and the index of the pre-existing record is passed to the reporting module. Otherwise, the file is copied into the `data/archive` directory, and in FASTAS, either the pre-existing record or a new record is updated with the metadata. The index of the imported record (new or pre-existing) is passed to the next stage - parsing variants.

### Parsing variants

After any import step and before any reporting step, a *new sequence* is analysed by the `parse_variants` module in `src/components`. The first step is for all entries in the VARIANTS second-level table whose FASTA ID matches the passed index to be deleted. For new (previously unused index) or non-imported sequences (index = 0), this will be redundant; it is only for the overwriting of sequence data that step is important. Next, the sequence is mapped to the set of references and the best matches for each domain are established. Each domain is analysed for the three types of variant - missing data, indels, and substitutions - and those detected are encoded with the modifed HGVS notation described above in the [Database Table structure](#DTs) section. The output is returned in the same format as the VARIANTS table and concatenated to VARIANTS. 

The VARIANTS table is populated with the outputs of `parse_variants` - one record per variant. For non-zero indexes, successful population is followed by setting a flag in the FASTAS table indicating the sequence has been parsed in this way. (Any interruption in this process could lead to incomplete and corrupt data).

The index is returned, allowing the reporting stage to access the variant and file metadata, if necessary.

### <a name="sequences_report">`-r` / `--report_data`

If report is selected, then the `report` module in `src/components` is used to generate a report. The extent of such reports can vary from a simple literature search giving susceptibilities to the "primary" drugs, to inclusion of secondary drugs and even cross-referencing phenotypic data. Where there is discordance between the literature and the phenotypic data, any genotypes associated with the phenotypically discordant samples can be interrogated. Alternatively, the MUTATION module can be used to look at these instances of confounding data. The `report` module is also used to generate phenotypic assay reports, and other combined reports from the MOLIS_ID and MUTATION modules, hence the list of passed arguments can be long. See the `report` section for more information.

### Post-processing

The final step is to remove entries from VARIANTS with a FASTA index of 0 (i.e. non-imported sequences). This is achieved with a simple `sed` script.

## <a name="phenos">phenos

To process phenotypic data, use the `phenos` subparser. Importing and reporting are set to True for all phenotypic assay files. Other options are as follows:

| Argument | Description | Default |
| :--- | :--- | :--- |
| `-p` / `--previous` | Set to report all phenotypic data sharing the same MOLIS ID | |
| `-f` / `--pheno` | `path/to/pheno/file` | |
| `-d` / `--dir` | `path/to/directory/containing/pheno file` | Current directory |
| `--runid` | name of the run (e.g. LIMS tracking ID) | Inferred from directory name |
| `--rundate` | date of the run | Inferred from directory name |
| `--force` | Set to overwrite existing data in the PHENOS table | |

The process map for this subparser is simpler than for `sequences`. Import is identical, save for the name of the top-level table. Instead of parsing variants with `parse_variants`, phenotypic data is parsed from the data file using the `parse_pheno` module in `src/components`. Unlike sequence data, where the input data is a series of FASTA formatted sequences, the nature of the phenotypic assays (currently exclusively Plaque Reduction Assays) is such that the format of the Microsoft Excel data files has varied over time, and a function to identify the data file version, and a set of separate functions for each specific file version are necessary to derive and calculate the EC50s for each drug within each file. The output of `parse_pheno` is a table of control EC50 and sample EC50 rows, one for each drug with a result. Again, the PHENOS table has a flag set once the EC50S table is successfully updated.

For reporting, the `-p` flag governs whether the single phenotypic datafile is reported or whether all data from all phenotypic assays for the same MOLIS ID are reported. Whether these are displayed separately or in a single table, and sorted by date and drug or not at all, is governed by parameters passed to the `report` module. As noted above, the versatility of the `report` module favours configuration. Some common configurations are preset in `config/report.yaml`.

## <a name="molis">molis

To obtain a report on all data for a given MOLIS ID, ise the `molis` subparser. Bypassing all primary analysis, this option directly engages the `report` module of `src/components`, returning a set of extended reports by default. Other options are as follows:

| Argument | Description | Default |
| :--- | :--- | :--- |
| `-l` / `--literature` | Return information from the literature | |
| `-p` / `--phenos` | Return data from phenotyping assays | |
| `--recent` | Return only the most recent result(s) | |

## <a name="mutation">mutation

Use the `mutation` subparser to interrogate the data on a specific mutation. 

| Argument | Description | Default |
| :--- | :--- | :--- |
| `--range` | Examine <INT> loci either side of the target mutation | -1 |
| `--homology` | Examine the corresponding locus/loci in the *other* HSV type | |
| `mut` | Mutation (usually in HGVS format) | *positional*, *required* |


The `mut` argument must be present as a positional argument. For most queries, this will be in the modified HGVS format (see [here](#DTs)) but for queries involving novel mutations, this can have an amino acid position as the final element rather than a fully encoded mutation. This is used in conjunction with the `--range` and `--homology` arguments, so as to enable scanning of the neighbourhood both in the same HSV type and in the other type (i.e. HSV-1 if the original query is for HSV-2).

At least one of `-l` and/or `-p` *must* be selected (or else there would be nothing to report!). The `-l` option performs a literature search either for the specific variant, or all mutations within the scope of the `--range` and `--homology` arguments. The default for `--range` is set to -1 so as to allow the option to pass 0, if all that is sought are other amino acid variants at that locus only. 

If `-p` is set, the `--range` argument will be ignored; `--homology` can be set to report any effects of the homologous mutation in the other HSV type. 

For all queries for a specific variant, the VARIANTS table is interrogated to find all FASTA IDs whose sequence contains the variant. From the FASTAS table, their MOLIS IDs áre obtained, and these are used to interrogate the PHENOS and EC50S table to obtain susceptibiilty outputs. In cases where resistance is detected, all *other* variants in that sample are interrogated for their resistance associations.
