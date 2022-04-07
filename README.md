# Scripts from Massively parallel, computationally-guided design of a pro-enzyme

This repository contains scripts used in our manuscript, "Massively parallel, computationally guided design of a proenzyme," which describes a computational approach and screening method to design pro-enzymes of an FDA-approved therapeutic enzyme, carboxypeptidase G2.

XML scripts used to generate are maintained in the [rosetta_scripts_scripts repository](https://github.com/RosettaCommons/rosetta_scripts_scripts) and are available with public releases of Rosetta.

## generate_sequences.py

The generate_sequences.py script loads Rosetta poses from various locations on disk, extracts the protein sequence from the appropriate chain, and reverse translates it according to a custom codon table.  Sequences are output as protein sequences, DNA sequences, and DNA sequences with Gibson assembly adapters.

## deep_seq_magicblast.R

The deep_seq_magicblast.R script analyses next-generation sequencing data.  The FASTQ files from various NGS samples were previously aligned to a set of "reference" sequences using Magic BLAST.  These alignments are used as the input to this script, which will determine the number of reads for each of the reference sequences.  Across biological replicates of similar samples, the average reads and error are determined.

# Publication Information

The manuscript that uses these scripts is published.  Please cite the following article if you use these scripts:
Massively parallel, computationally guided design of a proenzyme
Brahm J. Yachnin, Laura R. Azouz, Ralph E. White III, Conceição A. S. A. Minetti, David P. Remeta, Victor M. Tan, Justin M. Drake, Sagar D. Khare.
*Proc Natl Acad Sci U S A*. **119** (15): e2116097119; doi: https://doi.org/10.1073/pnas.2116097119
