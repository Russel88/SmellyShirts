# SmellyShirts

#### Scripts and data to reproduce results for the paper:
The T-shirt microbiome is individual and shaped by washing and fabric type (Submitted)

## How to reproduce results
### Preparation
* Get sequences from http://www.ncbi.nlm.nih.gov/bioproject/594290
* Run data/Cutadapt_seq.txt to remove primers from sequences
* Get reference databases from doi:10.5281/zenodo.801827
* Run data/Cutadapt_ref.txt to remove primers and everything outside V3V4 region from sequences

### R analysis
* Run data/load_data.R script (R 3.4.3)
* Run the remaning .R scripts (R 3.5.2)

## Software
* Cutadapt 1.14 (doi:10.14806/ej.17.1.200)
* mafft v7.376 (doi:10.1093/molbev/mst010)
* FastTree 2.1.10 (doi:10.1371/journal.pone.0009490)
* GNU parallel (O. Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine)
* Bowtie v2.2.3 (doi:10.1038/nmeth.1923)
* Samtools v1.5 (doi:10.1093/bioinformatics/btp352)
* Gephi 0.9.2
* R 3.4.3 (no link means they are from CRAN)
	* dada2 1.7.7 (doi:10.1038/nmeth.3869)
	* phyloseq 1.22.3 (doi:10.1371/journal.pone.0061217)
	* foreach 1.4.4
	* ips 0.0-7
	* ape 5.0
* R 3.5.2 (no link means they are from CRAN)
	* DAtest 2.7.15 (doi:10.1101/241802)
	* ggplot2 3.1.0
	* vegan 2.5-4
	* nlme 3.1-137
	* multcomp 1.4-1
	* foreach 1.4.4
	* phyloseq 1.26.1
	* COEF 0.0.0.9 (https://github.com/Russel88/COEF)
	* eulerr 5.1.0		
		
## Databases
* Taxonomy: RDP 16 release 11.5 train set and species set (doi:10.1093/nar/gkt1244, doi:10.5281/zenodo.801828)
* 16S rRNA gene copy numbers: rrnDB-5.4 RDP (doi:10.1093/nar/gku1201)