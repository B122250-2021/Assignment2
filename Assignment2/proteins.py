import os,subprocess,shutil,sys,re,glob,ete3
import numpy as np
import pandas as pd
from ete3 import Tree
import argparse

#input test
os.system("mkdir protein_sequences && mkdir protein_sequences/individual_sequences && mkdir prosite_motifs")

parser=argparse.ArgumentParser(
    description='''proteins.py performs an analysis based on user input.\nThe user can provide a Taxonomic Group and a Protein Family of their choice.''',
    epilog="""Please message <my.email> to report any bugs.""")
parser.add_argument('--taxonomic_group','-t', type=str, help='Taxonomic group of your choice ')
parser.add_argument('--protein_family', '-p', type=str, help='Protein family of your choice')

args=parser.parse_args()
prot_fam = args.protein_family 
tax_grp = args.taxonomic_group

#If no arguments are provided
if tax_grp is None:
	tax_grp = input("No taxonomic group was selected. Please enter the taxonomic group you want to investigate:\n\t")
if prot_fam is None:
	prot_fam = input("No protein family was selected. Please enter the desired protein family:\n\t")

#Doing an esearch, and piping it to efetch

#need a prompt asking if you want partial or not
print('Fetching Accession values for {} proteins in {}...'.format(prot_fam,tax_grp))
print('Please quit now if these is not the desired query')
os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein] NOT PARTIAL"| efetch -format acc > protein_sequences/acc_numbers.txt')
#read the accession numbers to get the sequences, and print a warning if longer than 1000
#Ask the user if they want to set a limit
#efetch limit
acc_list=[]
for line in open("protein_sequences/acc_numbers.txt", "r"):
	acc_list.append(line.strip())
print('There are {} protein sequences that match this query.'.format(len(acc_list)))
if len(acc_list) > 1000:
	if input('There are more than 1000 sequences in this analysis, press y if you want to continue ') != 'y':
		print('Exiting analysis')
		sys.exit()

print('Done, Counting Species...')
#Need to figure out how many species I have
os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein] NOT PARTIAL" | efetch -format gbc | xtract -pattern INSDSeq -group INSDFeature -block INSDQualifier -element INSDQualifier_name INSDQualifier_value | grep {} | cut -f2 >  all_species.txt'.format('organism')) #used the .format to avoid sysntax errors with apostrophes
species=[]
for line in open("all_species.txt", "r"):
	species.append(line.strip())
non_redundant_species=set(species)
if input('There are {} species in your analysis, you can see them in all_species.txt. Press y if you wish to continue.\n\t'.format(len(non_redundant_species))) != 'y':
	print('Exiting analysis')
	sys.exit()

#create a dict with acc vs species, useful later on when building tree
#since the values are written to the file in the order that they are retrieved, acc and species name will be in order in their respective files
spp_list=[]
acc_spp = {}
for line in open("all_species.txt", "r"):
	spp_list.append(line.strip())
for a,s in zip(acc_list,spp_list):
	acc_spp[a] = s


print('Fetching Protein Sequences...')
os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein] NOT PARTIAL" | efetch -format fasta > protein_sequences/input.fa')

#I think I need to use BLAST and Clustalo as te tests before continuing to the next section
#use clustalo to align
print('Aligning sequences + building phylogenetic tree...')
os.system('clustalo --threads=64 -i protein_sequences/input.fa -o protein_sequences/aligned_sequences.fa --guidetree-out=clustalo_tree')


#Making tree file
#Now I need to change all acc values in the tree file for the format acc_spp, as a tree that only has the accession value is not very informative
newtree = open("newick_tree.txt", "a")
for line in open("clustalo_tree", "r"):
	if line.startswith('(') or line.startswith(')') or line.startswith(',') or line.startswith(';'):
		newtree.write(line)
	elif line.startswith('sp|'):
		acc = line.split('|')[1] #some sequences are found in this format
		spp = acc_spp.get(acc)
		spp = spp.replace(' ','_')
		newtree.write(acc + '_' + spp +':' + line.split(':')[1])
	else:
		acc = line.split(':')[0]
		spp = acc_spp.get(acc)
		spp = spp.replace(' ','_')
		newtree.write(acc + '_' + spp +':' + line.split(':')[1]) 
newtree.close()
sys.stdout = open("tree.txt", "w")
t = Tree("newick_tree.txt", format=1)
print(t)
sys.stdout.close()
#reopening stdout
sys.stdout = open("/dev/stdout", "w")
print('Done. File tree.nw contains Guide tree written in Newick format, tree.txt contains a visual representation of the infered tree.\n Producing conservation plot. Graph will appear on screen, but will also be saved as plotcon.svg')

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
	if acc.startswith('sp'):
		acc = acc.split('|')[2]
	current_file = open('protein_sequences/individual_sequences/' + acc + '.fa','w')
	current_file.write('>' + sequence)
	current_file.close()
#Scan prosite
#loop through files
directory = 'protein_sequences/individual_sequences/'
for filename in os.listdir(directory):
	file = os.path.join(directory, filename)
	os.system('patmatmotifs -sequence {} -outfile prosite_motifs/{}.patmatmotifs -auto'.format(file,file[39:])) #file[39:] is the accessioan number

#Extract motifs, and write them down to a file
files = glob.glob('prosite_motifs/*')
for file in files:
	motifs = open("motifs.txt", "a")
	elements = []
	with open(file,'r') as f:
		for line in f:
			pattern1 = 'Sequence'
			if re.search( pattern1, line):
				elements.append(line.split()[2]) #the accession value
			pattern2 = 'Motif'
			if re.search(pattern2,line):
				elements.append(line.split()[2]) # the motif
	for element in elements:
		motifs.write(element + '\t')
	motifs.write('\n')
print('Done. The motifs for each sequence can be found in motifs.txt')

#Option to perform tmap analysis
if input("Do you want to predict and plot possible transmembrane segments in the sequences (tmap)? If yes, press y:\n\t") =='y':
	os.system('tmap -sequences protein_sequences/aligned_sequences.fa -graph svg -outfile report.tmap && tmap -sequences protein_sequences/aligned_sequences.fa -graph x11 -outfile report.tmap')
