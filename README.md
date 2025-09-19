# Membrane-protein-library-ML

This repository contains datafiles and additional code for the (unpublished) membrane protein ML paper, using computational libraries to train predeictive models of membrane protein expresion in E. coli.

### Files available here
|Filename|Description|
|--------|-----------|
|original_designs.fasta|12,248 unique protein sequences from the original Rosetta design run, as described in Hardy et al (2023) PNAS 120 (16) e2300137120|
|resfile.txt|The original design resfile dictating sequence positions and the design alphabet. Note the original model lacks a start Methionine, add +1 for equivalent residues in expressed proteins.|
|oligo_pool.fasta|DNA sequences backtranslated to nucleotide from ***original_designs.fasta***, using the _E. coli_ codon usage table for high-expressing sequences|
|curnow_bright_seqs.fasta|FACS-selected sequences - high-GFP phenotype|
|curnow_bright_seqs_protein.fasta|Translation of FACS-selected sequences - high-GFP phenotype|
|curnow_dim_seqs.fasta|FACS-selected sequences - low GFP phenotype|
|curnow_dim_seqs_protein.fasta|Translation of FACS-selected sequences - low GFP phenotype|
|no_selected_seqs_DNA.fasta|Oligo pool nucelotide sequences with FACS-selected sequences removed. Unlabelled|
|no_selected_seqs_protein.fasta|Translated oligo pool sequences with FACS-selected sequences removed. Unlabelled|
|no_selected_seqs_protein.csv|Translated oligo pool sequences with FACS-selected sequences removed. Unlabelled, used for inference|
|labelled_fasta.csv|Protein sequences from FACS-sorted populations labelled bright=1, dim=0|
|predicted_protein.csv|Model inference of expression level of each sequence	from ***no_selected_seqs_protein.fasta***|
|membrane_protein_ML.ipynb|Jupyter notebook for all machine learning|
|labelled_fasta_DNA.csv|DNA sequences from FACS-sorted populations labelled bright=1, dim=0. Enables DNA level analysis, not used here|
|no_selected_seqs_DNA.csv|Nucelotide sequences from the original oligo pool with FACS-selected sequences removed. Not used here|

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
Using Rosetta 3.14 and the Slurm command ```srun extract_pdbs.serialization.linuxgccrelease```

FASTA sequences were obtained from each extracted decoy with ```$ python pdb2fasta.py *.pdb > all_seqs.fasta```

Seqkit was installed via homebrew ```% brew install seqkit```

Sequences were unduplicated with ```% seqkit rmdup -s all_seqs.fasta > unique_seq.fasta```

Leaving 12,248 unique protein sequences. _EMBOSS_ tools were then used to backtranslate each amino acid sequence to a cognate nucleotide sequence based upon the codon-usage table for high-expressing _Escherichia coli_ proteins:
```
% brew install emboss
% embossdata -fetch -file Eecoli_high.cut
% backtranseq -cfile Eecoli.cut unique_seq.fasta
```
Universal adaptors were introduced at either end of all library sequences with the seqkit mutate function:
```
% cat unique_DNA_high.fasta | seqkit mutate -i 0:ACTTTAAGAAGGAGATATACCATG > 5_prime_extend_high.fasta     

% cat 5_prime_extend_high.fasta | seqkit mutate -i -1:GCGGCCGCACTCGAGCTGGTGCCGCGCGGCAGCA > end_ext_with_Ala_spacer.fasta
```
# Nanopore sequencing
Nanopore basecalling used Oxford Nanopore Dorado 0.7.1 in duplex mode on the University of Bristol HPC cluster BluePebble.

Raw sequencing data (POD5) was split by channel information on the local computer.

### Splitting POD5 files by channel:
Install POD5 python tools:
```
% pip install pod5
```
Then:
```
% pod5 view /path/to/pod5/folder/ --include "read_id, channel" --output summary.tsv
```
This creates file 'summary.tsv' with channel information. 
```
% pod5 subset /path/to/pod5/folder/ --summary summary.tsv --columns channel --output split_by_channel
```
### Basecalling
Dorado 0.7.1 for linux64 was downloaded from https://github.com/nanoporetech/dorado. To install on HPC:
```
$ tar -xzvf dorado-0.7.1-linux-x64.tar.gz
```
Run using Slurm submission on a single GPU node, specifying FASTQ output. The working directory requires the split pod5 files in subdirectory 'split_pod5s'. 
```
#!/bin/bash
#SBATCH --job-name=dorado_basecalling
#SBATCH --nodes=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=20GB
#SBATCH --output=dorado_basecalling_%j.log

module add cuda/12.4.1
export PATH="path/to/directory/nanopore/dorado-0.7.1-linux-x64/bin:$PATH"
dorado duplex --emit-fastq sup split_pod5s/ > seq_file.fastq
```
This takes ~5h to complete at 8e+04 bases/s, duplex rate 16%.

# Sequence analysis
Install EMBOSS and seqkit via homebrew and fastq-filter (https://github.com/LUMC/fastq-filter)
```
% brew install seqkit
% brew install EMBOSS
% pip install fastq-filter
```
Remove first 150 nt of each nanopore run, since these are of lower quality
```
% seqkit seq_file.fastq -r 150:-1 > seq_file_trim.fastq
```
Filter on length to remove short and very long reads
```
% seqkit seq seq_file_trim.fastq -m 1000 -M 14000 >  seq_file_short.fastq
```
Filter the remaining reads by quality:
```
% fastq-filter -e 0.001 -o seq_file_Q30.fastq seq_file_short.fastq.fastq
```
To extract the open reading frame for library-gfp fusions:
```
% getorf -sequence seq_file_Q30.fastq -minsize 1119 -maxsize 1119 -reverse Y -find 3 -outseq seq_file_Q30_orfs_nt.fasta
```
Flag ‘-find 3’ locates ORFS between start and stop codon, and retains the output as nucleotide. 
To isolate the library sequences:
```
% seqkit subseq seq_file_Q30_orfs_nt.fasta -r 1:339 > seq_file_Q30_orfs_library_nt.fasta
```
Finally, remove duplicate sequences and print one instance of each unique sequence to file
```
% seqkit rmdup seq_file_Q30_orfs_library_nt.fasta -s -d Q30_nt_dups.fasta -o Q30_nt_unique.fasta
```    
# Machine Learning
Machine learning was implemented in a Jupyter Notebook (ver. 1.1.1). The input files are ***labelled_fasta.csv*** which is the FACS-selected sequences with appropriate labels and ***no_selected_seqs_protein.csv*** which is the remainder of the oligo library used for inference. The analysis presented is at protein level, equivalent datafiles are provided here to allow analysis at DNA level if required. Packages were installed via Anaconda navigator and versions were, to best knowledge: pandas 2.2.3, seaborn 0.13.2, umap-learn 0.5.4, shap 0.47.2, matplotlib 3.10.0.

# Raw sequencing files
The raw DNA sequencing data is in the Nanopore .pod5 format. Raw files are deposited at the EMBL-EBI BioStudies repository as project S-BSST2184.
