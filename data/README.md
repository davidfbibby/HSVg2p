# Organisation of `HSVg2p/data`

Three categories of data are found within this project, and each are located in their own subdirectory:  
1. `archive`  
These are the raw datafiles to which the dynamic tables point when referencing genotypic or phenotypic data. Essentially, they are FASTAs or PRA in one of several Excel formats. An `--import` argument to the `sequences` subparser will copy a new file to this location.  
2. `static`  
The reference tables for analysis are kept here. These include the literature sources for interpreting variants, the interpretations themselves, EC50 cut-offs for each HSV type-target-drug combination, and drug-domain target information.  
3. `dynamic`  
The main database is held here in two pairs of tables - one for genotypic and one for phenotypic information. Each pair contains a top-level table containing metadata for all `--import`ed files stored in the `archive`, and a secondary table derived from the data pointed to by the first. For genotypic data, this is the parsed variants coded by HSV type, domain, locus, mutation type, and specific detail. For phenotypic data, this is the EC50s for each drug, calculated from the primary data (e.g. plaque counts and drug concentrations). These four tables are stored in `feather` format for rapid read-write and efficient (ZSTD).

A final table is held in the top directory - `HSVTYPES.tsv`. This is a simple mapping of MOLIS ID to HSV type for those samples for which this was determined by immunofluorescence prior to conducting a titration-PRA. More recently, the genotypic data (at least for TK) is generally output well before the need to know the HSV type for a phenotypic test, so the type is known ahead of time. Hence the genotypic data informs the choice of phenotypic cut-offs.

## Data flows

Genotypes can be analysed without importation. In these cases, the variants are parsed, and a report is generated with reference to the literature and, where stored samples have been phenotyped, the phenotypes (and associated genotypes). Population of the `dynamic` tables requires the `--import` flag to be set for genotypes. All phenotypes are automatically imported. The first step of importation is to extract the metadata and compare it to existing records. Duplicates are ignored unless the `--force` flag is set in the `sequences` or `phenos` subparser. For new data, the file is copied to the `archive` and the metadata and filename is entered into the top-level table. Each record is initialised with two Boolean fields set to False - *analysed* & *suppress* - and an empty *HSV_type* field. The *analysed* and *HSV_type* fields are only populated once the immediate downstream analysis is complete. When the *suppress* field is set to True (by using the `suppress` subparser), the record is rendered invisible to many analysis tools. 
