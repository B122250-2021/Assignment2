import os,subprocess,shutil,sys
import numpy as np
import pandas as pd

#input test
os.system("mkdir protein_sequences && mkdir protein_sequences/individual_sequences && mkdir prosite_motifs")
tax_grp = input("Please enter the taxonomic group you want to investigate:\n\t")
prot_fam = input("Please enter the desired protein family:\n\t")
#Doing an esearch, and piping it to efetch

#need a prompt asking if you want partial or not
print('Fetching Accession values for {} proteins in {}...'.format(prot_fam,tax_grp))
print('Please quit now if these is not the desired query')
os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein] NOT PARTIAL"| efetch -format acc > protein_sequences/acc_numbers.txt')
#read the accession numbers to get the sequences, and print a warning if longer than 1000
#Ask the user if they want to seta limit
#efetch limit
acc_list=[]
for line in open("protein_sequences/acc_numbers.txt", "r"):
	acc_list.append(line.strip())
print('There are {} protein sequences that match this query.'.format(len(acc_list)))
if len(acc_list) > 1000:
	if input('There are more than 1000 sequences in this analysis, press y if you want to continue ') != 'y':
		print('Exiting analysis')
		sys.exit()
print('Fetching Protein Sequences...')
os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein] NOT PARTIAL" | efetch -format fasta > protein_sequences/input.fa')
print('Done')
#Need to figure out how many species I have
os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein] NOT PARTIAL" | efetch -format gbc | xtract -pattern INSDSeq -group INSDFeature -block INSDQualifier -element INSDQualifier_name INSDQualifier_value | grep {} | cut -f2 >  all_species.txt'.format('organism')) #used the .format to avoid sysntax errors with apostrophes
species=[]
for line in open("all_species.txt", "r"):
	species.append(line.strip())
non_redundant_species=set(species)
if input('There are {} species in your analysis, you can see them in all_species.txt. Press y if you wish to continue.\n\t'.format(len(non_redundant_species))) != 'y':
	print('Exiting analysis')
	sys.exit()
#I think I need to use BLAST and Clustalo as te tests before continuing to the next section
#use clustalo to align
print('Aligning sequences...')
os.system('clustalo -i protein_sequences/input.fa -o protein_sequences/aligned_sequences.fa')
print('Done. Producing conservation plot. Graph will appear on screen, but will also be saved as plotcon.svg')
os.system('plotcon -sequences protein_sequences/aligned_sequences.fa -winsize 4 -graph x11 -graph svg && plotcon -sequences protein_sequences/aligned_sequences.fa -winsize 4 -graph x11')
print('Done')
print('Searching for known protein motifs in our query sequences...')
#to scan for motifs, file needs to be split
foo = open("protein_sequences/input.fa", "r")
content = foo.read()
allsequences = content.split('>')
foo.close()
for sequence in allsequences[1:]:
	acc = sequence.split()[0] #retrieve the accession value
	current_file = open('protein_sequences/individual_sequences/' + acc + '.fa','w')
	current_file.write('>' + sequence)
#Scan prosite
#loop through files
directory = 'protein_sequences/individual_sequences/'
for filename in os.listdir(directory):
	file = os.path.join(directory, filename)
	os.system('patmatmotifs -sequence {} -outfile prosite_motifs/{}.patmatmotifs'.format(file,file[39:])) #file[39:] is the accession number
