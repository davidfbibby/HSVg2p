# Component modules in HSVgeno2pheno

There are 3 modules in `src/componenets`, one for processing sequence data, one for processing phenotypic data, and a large one for generating a wide variety of reports.

## `parse_variants`

Takes a FASTA sequence and returns a `pandas` DataFrame of the variants. Uses a defined set of references, located in the /HSVamp directory within $GENOMANCER_REFS. 