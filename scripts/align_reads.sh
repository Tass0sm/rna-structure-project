#!/usr/bin/env bash

# RNA_section__454_9627.fasta
# ./dna/072517_HIVDNA_filtered.fastq
# ./dna/alignment.sam

minimap2 -ax map-ont $1 $2 > $3
