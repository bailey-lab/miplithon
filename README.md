# miplithon

This is a snakemake pipeline for analyzing drug resistance prevalences. The
pipeline expects 3 tables of read counts similar to the ones output by the
variant calling step of miptools. Each table consists of columns that are
mutations and rows that are samples, and a given (column, row) value is the
count of reads associated with a mutation in a given sample.

The tables are:
  - coverage_AA_table.csv: the total number of reads that span the genomic
  region.
  - reference_AA_table.csv: the count of reads that did not support
  the mutation.
  - alternate_AA_table.csv: the count of reads that supported the mutation.

The user also provides a list of mutations of interest (in a yaml file) along
with coverage and alternate thresholds to apply to each mutation, and a metadata
file that assigns categorical values to each sample in a column (or columns) of
interest.

The program outputs a table of categorical values for the column of interest,
as well as how many samples passed the thresholds for each category.

For example, a user might have counts for samples collected from several
geographic regions, and metadata that links individual samples to individual
regions. The user might be interested in 10 mutations at thresholds of coverage 3
and alternate count of 1. The program would output how many samples in each
geographic region had counts that passed these thresholds for each of the 10
mutations.

## Installation

To install this program, first download snakemake (e.g. pip install snakemake).

Next, cd to a folder of interest and git clone the pipeline, e.g.
```bash
git clone https://github.com/bailey-lab/miplithon.git
```

## Usage

To use this pipeline, edit the yaml file using the comments in the yaml file.
Some example values have been filled in.

## Results

alternate_freqs.csv gives the unfiltered within sample allele frequencies of
each mutation in each sample.

This is graphed as a plotly heatmap in frequency_heatmap.html

The main results can be found in the
prevalences/{summarize_by}:all_{cov}_{alt}_summary.tsv file, where {summarize_by}
is the column of interest, {cov} is the threshold coverage, and {alt} is the
threshold alternate count (number of times the mutant has been seen. Mutations
in this table will match those provided by the user in the yaml file.

prevalences_by_threshold/{cov}_{alt}_{count}_summary.tsv gives all mutations
(even those not explicitly asked for) that passed coverage, alternate count, and
minimum sample count thresholds from the yaml file.

background_mutations/561_on_Asian_backgrounds.tsv is currently hardcoded to give
the number of samples that have 561H and Asian mutations. This will be
generalized in the future.

prevalences/mdr1-NFD_{cov}_{alt}_cov.txt and
prevalences/mdr1-NFD_{cov}_{alt}_alt.txt give tables that show the number of
samples that contain the mdr1-Asn86, mdr1-Tyr184Phe, and mdr1-Asp1246 mutations
(number of samples that pass coverage thresholds and number of samples that pass
coverage and alternate thresholds, respectively).
