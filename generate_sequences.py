#!/usr/bin/python3

# This script loads Rosetta poses from various locations on disk, extracts the protein sequence from the 
# appropriate chain, and reverse translates it according to a custom codon table.  Sequences are output 
# as protein sequences, DNA sequences, and DNA sequences with Gibson assembly adapters.

## Based on a file listing designs from various paths, generate a list of sequences.

import os, csv

infile = "poselist.csv"  # csv file with header fname, path, and class.  fname is the filename, path is its location, and class is the design set (given in outputprefixes below)
outputprefixes = ["reg", "hbnet", "hotspot"]  # The three design classes: reg (regular fixed bb design), hbnet (hbnet design), and hotspot (the original hotspot set)
# "Gibson assembly" adaptor sequences:
dna5prime = "CTTTATTTTCAAGGGGGGTCT"
dna3prime = "TAAAAGCTTAATTAGCTGAG"

#Functions to open and close files.
def init_outfiles(prefixes):
	proteinfile = []
	dnafile = []
	dnanohififile = []
	for prefix in prefixes:
		proteinfile.append(open(prefix + '_protein.txt', 'w', buffering = 500))
		dnafile.append(open(prefix + '_dna.txt', 'w', buffering = 500))
		dnanohififile.append(open(prefix + '_dnanohifi.txt', 'w', buffering = 500))
	
	return(proteinfile, dnafile, dnanohififile)
	
def close_outfiles(open_files):
	for idx in range(0, len(open_files)):
		open_files[0].close()
		
#Function to write to the correct files
#It requires the pdb string, the list of prefixes, the sequence, and the list of outfiles.
def write_sequences(pdbstring, prefix, seq, outfiles):
	#Figure out the index corresponding to the pdbstring.
	idx = prefix.index(pdbstring[2])
	#Write the sequence to that file.
	outfiles[idx].write(seq)
	
#Function to use codon_tools to reverse translate the protein sequence.
def codopt(protein_string):
	#Import modules
	from Bio.Seq import Seq
	from Bio.Alphabet import generic_protein
	from codon_tools import FopScorer, CodonOptimizer
	
	#Define a custom codon table.
	usable_codons = {  'A':['GCA', 'GCC', 'GCG'], 'R':['CGT', 'CGC'], 'N':['AAC', 'AAT'], 'D':['GAC', 'GAT'], 'C':['TGC', 'TGT' ], 'Q':['CAA', 'CAG'], 'E':['GAA'], 'G':['GGC', 'GGT'], 'H':['CAC', 'CAT'], 'I':['ATC', 'ATT'], 'L':['CTC', 'CTG', 'CTT', 'TTA', 'TTG'], 'F':['TTT', 'TTC'], 'P':['CCA', 'CCG', 'CCT'], 'S':['AGC', 'TCC','TCT'], 'T':[ 'ACT', 'ACC', 'ACG' ], 'Y':['TAC', 'TAT'], 'V':['GTA', 'GTC', 'GTG', 'GTT'], 'W':['TGG'], 'M':['ATG'], 'K':['AAA', 'AAG'], '*':['TAA', 'TAG', 'TGA']}
	
	#Define a codon scorer that uses the codon frequency table above.
	scorer = FopScorer(opt_codons = usable_codons)
	#Define an optimizer based on scorer
	opt = CodonOptimizer(scorer)
	
	#Create a sequence object from the protein_string
	protein_seq = Seq(protein_string, generic_protein)
	#Randomly reverse translate the sequence.
	randdnaseq = opt.random_reverse_translate(protein_seq)
	#Now codon optimize.
	dnaseq = opt.hillclimb(randdnaseq)
	
	#Assert that the dnaseq matches the original protein sequence.
	assert str(Seq.translate(dnaseq[0])) == str(protein_seq)
	
	#Convert dnaseq to a string and return.
	return(str(dnaseq[0]))

#Initialize pyrosetta.
from pyrosetta import *
from pyrosetta.rosetta import *
init()

#Open the poselist as read-only.
poselistfile = open(infile, 'r')

#Open the output files for writing
protein_outfiles, dna_outfiles, dnanohifi_outfiles = init_outfiles(outputprefixes)

#Convert to a CSV list of lists
poselist = []
for line in csv.reader(poselistfile):
	poselist.append(line)
	
#Remove the header
poselist = poselist[1:]

#For each item in the poselist...
for pdbname in poselist:
	#Open the pose based on name and path
	posepath = pdbname[1] + pdbname[0] + '.pdb'
	pose = pose_from_pdb(posepath)
	
	#Grab the domain extra sequence
	sequence = pose.chain_sequence(2)
	#Add the newline to protein sequence
	protseq = sequence + '\n'
	
	#Reverse translate sequence
	revtrans = codopt(sequence)
	
	#Some code here to convert to DNA and add stuff to the sequence
	dnaseq = dna5prime + revtrans + dna3prime + '\n'
	
	#Write the protein sequence to the appropriate files.
	write_sequences(pdbname, outputprefixes, protseq, protein_outfiles)
	write_sequences(pdbname, outputprefixes, dnaseq, dna_outfiles)
	write_sequences(pdbname, outputprefixes, revtrans+'\n', dnanohifi_outfiles)
	
#Close the output files
close_outfiles(protein_outfiles)
close_outfiles(dna_outfiles)
close_outfiles(dnanohifi_outfiles)
