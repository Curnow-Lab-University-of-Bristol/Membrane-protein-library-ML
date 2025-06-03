# CytbX-library-ML

This repository contains datafiles and additional code for the CytbX ML paper.

Branch 'Nanopore-sequencing' contains code for NGS basecalling and sequence filtering.

### Files available here
|Filename|Description|
|--------|-----------|
|original_designs.fasta|12,248 unique protein sequences from the original Rosetta design run, as described in Hardy et al (2023) PNAS 120 (16) e2300137120|
|oligo_pool.fasta|DNA sequences backtranslated to nucleotide from ***original_designs.fasta***, using the _E. coli_ codon usage table for high-expressing sequences|
|curnow_bright_seqs.fasta|FACS-selected sequences - high-GFP phenotype|
|curnow_bright_seqs_protein.fasta|Translation of FACS-selected sequences - high-GFP phenotype|
|curnow_dim_seqs.fasta|FACS-selected sequences - low GFP phenotype|
|curnow_dim_seqs_protein.fasta|Translation of FACS-selected sequences - low GFP phenotype|
|no_selected_seqs_DNA.fasta|Oligo pool nucelotide sequences with FACS-selected sequences	removed. Unlabelled, used for inference|
|no_selected_seqs_protein.fasta|Traslated oligo pool sequences with FACS-selected sequences removed. Unlabelled, used for inference|
|protein_predicted.csv|Model inference of expression level of each sequence	from ***no_selected_seqs_protein.fasta***|


