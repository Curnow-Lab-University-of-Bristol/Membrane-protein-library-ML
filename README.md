# CytbX-library-ML

This repository contains datafiles and additional code for the CytbX ML paper.

### Files available here
|Filename|Description|
|--------|-----------|
|original_designs.fasta|12,248 unique protein sequences from the original Rosetta design run, as described in Hardy et al (2023) PNAS 120 (16) e2300137120|
|resfile.txt|The original design resfile dictating sequence positions and the design alphabet|
|oligo_pool.fasta|DNA sequences backtranslated to nucleotide from ***original_designs.fasta***, using the _E. coli_ codon usage table for high-expressing sequences|
|curnow_bright_seqs.fasta|FACS-selected sequences - high-GFP phenotype|
|curnow_bright_seqs_protein.fasta|Translation of FACS-selected sequences - high-GFP phenotype|
|curnow_dim_seqs.fasta|FACS-selected sequences - low GFP phenotype|
|curnow_dim_seqs_protein.fasta|Translation of FACS-selected sequences - low GFP phenotype|
|no_selected_seqs_DNA.fasta|Oligo pool nucelotide sequences with FACS-selected sequences	removed. Unlabelled, used for inference|
|no_selected_seqs_protein.fasta|Translated oligo pool sequences with FACS-selected sequences removed. Unlabelled, used for inference|
|protein_predicted.csv|Model inference of expression level of each sequence	from ***no_selected_seqs_protein.fasta***|

# Extracting design sequences
The Rosetta design run in Hardy et al (2023) PNAS 120 (16) e2300137120 generated 18,339 independent computational design decoys. These design decoys were extracted from 50 Rosetta silent files using the ***extract.pdbs*** command. Silent file #22 was prematurely terminated because of a time limit imposed during the design run and so was treated separately from the other 49 files, which were extracted in a single batch. The relevant flags were:
```
-in:file:silent <filename1> <filename2> <filename3>…
-linmem_ig 20
-ignore_unrecognized_res true
-extra_res_fa HEMred.params
-out:path:pdb output
-out:path:all output
-missing_density_to_jump true
```
Using Rosetta 3.14 and the Slurm command:
```
srun extract_pdbs.serialization.linuxgccrelease
```
