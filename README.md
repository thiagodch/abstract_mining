# abstract_mining
Repository for scripts used for the abstract mining project for the Data Science 2016 Course at Poli USP

## download_abs.sh
Bash script for downloading abstracts. The query term of interest should be placed in the "queryterm" placeholder. The output is a list with all relevant abstracts concatenated into a single file. In the script, the placeholder is "juvenile idiopathic arthritis", the query used in our project, for reproducibility.

`esearch -db pubmed -query <queryterm> | efetch -format medline` > abstracts.txt

## get_abs.pl
Perl script for separating the result of `download_abs.sh`. The input is the abstracts file, and the output is a folder containing all separated abstracts.

`perl get_abs.pl <abstractsfile>`

## processing.R
R functions for interactive importing and analysis of the `get_abs.pl` script. The outputs are several PDF files containing the results of analyses

## Rebember to:
Create a folder called abstracts
Name file containing the query results as abstracts.txt
