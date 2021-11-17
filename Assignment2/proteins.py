#!/usr/bin/python3

#Ask the user to select from the options:
#-use all species -use sequences most common species(if so, how many) - use sequences from least common -cycle one by one and decide if you use them
import os,subprocess,shutil,sys,re,glob,ete3,collections,argparse,fileinput
import numpy as np
import pandas as pd
from ete3 import Tree
from collections import Counter

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
if input('Do you want to include partial sequences? Press y if yes, any other key will discard partial sequences:\n\t') == 'y':
	partial = ''
else:
	partial = ' NOT PARTIAL'
#need a prompt asking if you want partial or not
print('Fetching Accession values for {} proteins in {}...'.format(prot_fam,tax_grp))
print('Please quit now if these is not the desired query')
os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein]' + partial +'"| efetch -format acc > protein_sequences/acc_numbers.txt')
#read the accession numbers to get the sequences, and print a warning if longer than 1000
#Ask the user if they want to set a limit
#efetch limit
acc_list=[]
for line in open("protein_sequences/acc_numbers.txt", "r"):
	acc_list.append(line.strip())
print('There are {} protein sequences that match this query.'.format(len(acc_list)))
if len(acc_list) == 0:
	print('There are no protein sequences that match your query. Please restart with different input arguments.')
	print('Exiting analysis')
	sys.exit()
elif len(acc_list) ==1:
	print('There is only one protein that matches your query, this programme can\'t use a single sequence as input. Please restart with different input arguments.')
	print('Exiting analysis')
	sys.exit()
if len(acc_list) > 1000:
	if input('There are more than 1000 sequences in this analysis, press y if you want to continue ') != 'y':
		print('Exiting analysis')
		sys.exit()

print('Done, Counting Species...')
#Need to figure out how many species I have
os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein]' + partial + '" | efetch -format gbc | xtract -pattern INSDSeq -group INSDFeature -block INSDQualifier -element INSDQualifier_name INSDQualifier_value | grep {} | cut -f2 >  all_species.txt'.format('organism')) #used the .format to avoid sysntax errors with apostrophes

species=[]
for line in open("all_species.txt", "r"):
	species.append(line.strip())
non_redundant_species=set(species)
#create a dict with acc vs species, useful later on when building tree
#since the values are written to the file in the order that they are retrieved, acc and species name will be in order in their respective files
spp_list=[]
acc_spp = {}
for line in open("all_species.txt", "r"):
	spp_list.append(line.strip())
for a,s in zip(acc_list,spp_list):
	acc_spp[a] = s
#Maybe add option for partial, and argument
print('Fetching Protein Sequences...')
os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein]' + partial + '" | efetch -format fasta > protein_sequences/all_sequences.fa')

#changing the format of the files that start with sp|, as pullseq will just skip them
def replaceAll(file,search,replace):
    for line in fileinput.input(file, inplace=1):
        if search in line:
            line = line.replace(search,replace)
        sys.stdout.write(line)

replaceAll('protein_sequences/all_sequences.fa','sp|','')
replaceAll('protein_sequences/all_sequences.fa','|',' ')

choice = input('\nThere are {} species in your analysis, you can see them in all_species.txt at the end of the analysis.\nFor now, please choose one of the following options, or press any other key to exit:\n\n\t-Use all sequences for the analysis (Will take a long time for large datasets) (Press 1)\n\n\t-Use sequences from the most common species (Press 2)\n\n\t-Use sequences from the least common species (Press 3)\n\n\t-Cycle one by one and chose which species you want (Not recommended for large datasets) (Press 4)\n\n\t'.format(len(non_redundant_species)))
if choice == '1' :
	print('All sequences will be used for the analysis. \nThe alignment may take a while, depending on the ammount of sequences')
	os.system('cp protein_sequences/all_sequences.fa protein_sequences/input.fa ')

elif choice == '2' :
	spp_no = int(input('These are the 10 most common species, please select how many species you want to use:\n\t{}\n'.format(Counter(species).most_common(10))))
	if spp_no > 10:
		print('The number you selected is greater than the displayed number, Exiting...')
		sys.exit()
	
	selected_species = []
	numbers = np.linspace(1,spp_no,spp_no)
	for x in numbers:
		selected_species.append(Counter(species).most_common(10)[int(x)-1][0])
	
	filtered_acc = []
	for species in selected_species:
		keys = [key for key,value in acc_spp.items() if value == species]
		for key in keys:
			filtered_acc.append(key)
	
	with open("protein_sequences/filtered_acc_numbers.txt", "w") as new_acc_values:
		for acc in filtered_acc:
			new_acc_values.write('{}\n'.format(acc))
	new_acc_values.close()
	os.system('/localdisk/data/BPSM/Assignment2/pullseq -i protein_sequences/all_sequences.fa -n protein_sequences/filtered_acc_numbers.txt > protein_sequences/input.fa')

elif choice == '3':
	spp_no = int(input('These are the 10 least common species, please select how many species you want to use:\n\t{}\n'.format(Counter(species).most_common()[-11:-1])))
	if spp_no > 10:
		print('The number you selected is greater than the displayed number, Exiting...')
	
	selected_species = []
	numbers = np.linspace(1,spp_no,spp_no)
	for x in numbers:
		selected_species.append(Counter(species).most_common()[-int(x)][0])
	filtered_acc = []
	for species in selected_species:
		keys = [key for key,value in acc_spp.items() if value == species]
		for key in keys:
			filtered_acc.append(key)
	
	with open("protein_sequences/filtered_acc_numbers.txt", "a") as new_acc_values:
		for acc in filtered_acc:
			new_acc_values.write('{}\n'.format(acc))
	new_acc_values.close()
	os.system('/localdisk/data/BPSM/Assignment2/pullseq -i protein_sequences/all_sequences.fa -n protein_sequences/filtered_acc_numbers.txt > protein_sequences/input.fa')

elif choice == '4':
	for current_species in non_redundant_species:
		if input('There are {} occurrances of {} in this dataset. Do you want to include the sequences in the analysis?\n Press y if yes, any other key will discard the sequence(s)'.format(species.count(current_species),current_species)) =='y':
			new_acc_values = open("protein_sequences/filtered_acc_numbers.txt", "a")
			keys = [key for key,value in acc_spp.items() if value == current_species]
			for key in keys:
				new_acc_values.write('{}\n'.format(key))
			new_acc_values.close()
	os.system('/localdisk/data/BPSM/Assignment2/pullseq -i protein_sequences/all_sequences.fa -n protein_sequences/filtered_acc_numbers.txt > protein_sequences/input.fa')

else:
	print('That wasn\'t one of the options...')
	print('Exiting analysis')
	sys.exit()

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
print('Done. File tree.nw contains Guide tree written in Newick format, tree.txt contains a visual representation of the infered tree.') 
winsize = input('Please specify a window size for the conservation plot (4 is recommended):\n\t')
if winsize.isnumeric() == False:
	print('That was not a number...\nExiting analysis')
	sys.exit()
if int(winsize) > 30:
	print('This window size is large, conservation plot might not be accurate')
if int(winsize) < 4:
	print('This window size is low, conservation plot might be noisy')

print('Producing conservation plot. Graph will appear on screen, but will also be saved as plotcon.svg')
os.system('plotcon -sequences protein_sequences/aligned_sequences.fa -winsize 4 -graph x11 -graph svg && plotcon -sequences protein_sequences/aligned_sequences.fa -winsize 4 -graph x11')
print('Done')
print('Searching for known protein motifs in the query sequences...')
#to scan for motifs, file needs to be split

inputseq = open("protein_sequences/input.fa", "r")
content = inputseq.read()
allsequences = content.split('>')
inputseq.close()
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
	os.system('patmatmotifs -sequence {} -outfile prosite_motifs/{}.patmatmotifs -auto'.format(file,file[39:])) #file[39:] is the accession number

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
	os.system('tmap -sequences protein_sequences/aligned_sequences.fa -graph svg -outfile report.tmap -auto && tmap -sequences protein_sequences/aligned_sequences.fa -graph x11 -outfile report.tmap -auto')
#iep analysis
if input('Would you like to calculate the isoelectric point of the chosen proteins? If yes, press y:\n\t') == 'y':
	os.system('iep -sequence protein_sequences/aligned_sequences.fa -outfile report.iep')
#hmoment analysis
#if input('Would you like to calculate and plot the hydrophobic moment for the chosen proteins? If yes, press y\n\t') == 'y':
#	os.system('hmoment -seqall protein_sequences/aligned_sequences.fa -graph svg -outfile report.hmoment  && hmoment -seqall protein_sequences/aligned_sequences.fa -graph x11 -outfile report.hmoment ')
#pepwindowall analysis
if input('Would you like to plot a supeimposed hydropathy plot for the chosen sequences (not recommended for large datasets)? If yes, press y\n\t') == 'y':
	os.system('pepwindowall -sequences protein_sequences/aligned_sequences.fa -window 4 -gxtitle="Residue Position" -gytitle="Hydropathy" -graph svg && pepwindowall -sequences protein_sequences/aligned_sequences.fa -window 4 -gxtitle="Residue Position" -gytitle="Hydropathy" -graph x11 ')
#interesting emboss apps: helixturnhelix, charge, pepstats?
