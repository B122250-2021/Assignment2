import os,subprocess,shutil,sys
import numpy as np
import pandas as pd

#input test
os.system("mkdir protein_sequences")
tax_grp = input("Please enter the taxonomic group you want to investigate:\n\t")
prot_fam = input("Please enter the desired protein family:\n\t")
#Doing an esearch, and piping it to efetch
print('Fetching Accession values for {} proteins in {}...'.format(prot_fam,tax_grp))
print('Please quit now if these is not the desired query')
os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[protein] NOT PARTIAL"| efetch -format acc > protein_sequences/acc_numbers.txt')
#read the accession numbers to get the sequences, and print a warning if longer than 1000
#probs should add a warning to see if you want to continue
acc_list=[]
for line in open("protein_sequences/acc_numbers.txt", "r"):
	acc_list.append(line.strip())
print('There are {} protein sequences that match this query.'.format(len(acc_list)))
if len(acc_list) > 1000:
	if input('There are more than 1000 sequences in this analysis, press y if you want to continue ') != 'y':
		print('Exiting analysis')
		sys.exit()
print('Fetching Protein Sequences...')
for acc in acc_list:
	os.system('wget -q -O protein_sequences/' + acc + '.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="' + acc + '"&strand=1&rettype=fasta&retmode=text"')
print('Done')
#Need to figure out how many species I have
species_file = open("species.txt", "w")
all_species = []
for acc in acc_list:
	with open('protein_sequences/' + acc + '.fasta') as file:
		first_line = file.readline()
		first_line = first_line.split('[')
		species = first_line[1]
		species = species[0:-2]
		all_species.append(species)
		species_file.write(species+'\n')
non_redundant_species = set(all_species)
if input('There are {} species in your analysis, you can see them in species.txt. Press y if you wish to continue.\n\t'.format(len(non_redundant_species))) != 'y':
	print('Exiting analysis')
	sys.exit()
#I think I need to use BLAST and Clustalo as te tests before continuing to the next section

