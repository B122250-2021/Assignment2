#!/usr/bin/python3

#Using pullseq
#Ask the user to select from the options:
#-use all species -use sequences most common species(if so, how many) - use sequences from least common -cycle one by one and decide if you use them

import os,subprocess,shutil,sys,re,glob,ete3,collections,argparse,fileinput
import numpy as np
import pandas as pd
from ete3 import Tree
from collections import Counter

parser=argparse.ArgumentParser(
    description='''proteins.py performs an analysis based on user input.\nThe user can provide a Taxonomic Group and a Protein Family of their choice.''',
    epilog="""Please message <my.email> to report any bugs.""")
parser.add_argument('--taxonomic_group','-t', type=str, help='Taxonomic group of your choice ')
parser.add_argument('--protein_family', '-p', type=str, help='Protein family of your choice')
parser.add_argument('--dir', '-d', type=str, help='Directory where program will be executed')
parser.add_argument('--winsize', '-w', type=str, help='Window size of conservation plot and hydropathy plot')
parser.add_argument('--partial', '-partial', type=str, help='y if you want to include partial sequences, anything else will exclude them')

args=parser.parse_args()
prot_fam = args.protein_family 
tax_grp = args.taxonomic_group
wd = args.dir
winsize = args.winsize
partial = args.partial

#If no arguments are provided
if wd is None:
	wd = input('Your current working directory is {}, please type the path where you want this analysis made:\n(the default is the current working directory, you can just press enter to stay in the current directory):\n\t'.format(os.getcwd()))
if tax_grp is None:
	tax_grp = input("No taxonomic group was selected. Please enter the taxonomic group you want to investigate:\n\t")
if prot_fam is None:
	prot_fam = input("No protein family was selected. Please enter the desired protein family:\n\t")
if partial is None:
	partial = input('Do you want to include partial sequences? Press y if yes, any other key will discard partial sequences:\n\t')
if partial == 'y':
	partial = ''
else:
	partial = ' NOT PARTIAL'

try:
	os.chdir(wd)
except FileNotFoundError:
	print("Directory: {} does not exist, defaulting to current directory".format(wd))
except NotADirectoryError:
	print("{} is not a directory, defaulting to current directory".format(wd))
os.system("mkdir protein_sequences && mkdir protein_sequences/individual_sequences && mkdir prosite_motifs reports")

#need a prompt asking if you want partial or not
print('Fetching Accession values for {} proteins in {}...'.format(prot_fam,tax_grp))
print('Please quit now if this is not the desired query')
#read the accession numbers to get the sequences, and print a warning if longer than 1000
#Ask the user if they want to set a limit
#efetch limit
os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein]' + partial + '" | efetch -format gbc | xtract -pattern INSDSeq -element INSDSeq_accession-version -element INSDSeq_organism | cut -f1,2 >  acc_spp.txt')
acc_list=[]
species=[]
acc_spp = {}
for line in open("acc_spp.txt", "r"):
	acc_list.append(line.split('\t')[0])
	species.append(line.split('\t')[1].strip('\n'))
	acc_spp[line.split('\t')[0]] = line.split('\t')[1].strip('\n')
non_redundant_species = set(species)
print('There are {} protein sequences, and {} species that match this query.'.format(len(acc_list),len(non_redundant_species)))


if len(acc_list) == 0:
	print('There are no protein sequences that match your query. Please restart with different input arguments.')
	print('Exiting analysis')
	sys.exit()
elif len(acc_list) ==1:
	print('There is only one protein that matches your query, this program can\'t use a single sequence as input. Please restart with different input arguments.')
	print('Exiting analysis')
	sys.exit()
if len(acc_list) > 1000:
	if input('There are more than 1000 sequences in this analysis, press y if you want to continue ') != 'y':
		print('Exiting analysis')
		sys.exit()

print('Fetching Protein Sequences...')
os.system('esearch -db protein -query "' + tax_grp + '[Organism] AND ' + prot_fam + '[Protein]' + partial + '" | efetch -format fasta > protein_sequences/all_sequences.fa')

#changing the format of the files that start with sp|, or tr|, as pullseq will just skip them, I encountered these two, but I am aware that there might be more
#which is why I added the next section which makes sure about the format and prints a warning if unrecognised

#Check the formats of the fasta file, to make sure that they will be analysed by pullseq. At this point, any sequence that doesn't start with the accession value will not be recognised by pullseq

for line in open("protein_sequences/all_sequences.fa", "r"):
	if line.startswith('>'):
		if any(acc in line for acc in acc_list) != True:
			print('WARNING! Unrecognised accession value\nThe following sequence:\n{}\nis either not in fasta format, or is in a format that is not supported in this program\n'.format(line))

choice = input('\nDo you want to filter your sequences? Please select one of the following options, or press any other key to exit:\n\n\t-Use all sequences for the analysis (Will take a long time for large datasets) (Press 1)\n\n\t-Use sequences from the most common species (Press 2)\n\n\t-Use sequences from the least common species (Press 3)\n\n\t-Cycle one by one and chose which species you want (Not recommended for large datasets) (Press 4)\n\n\t')
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
	
	
	for acc in filtered_acc:
		os.system('/localdisk/data/BPSM/Assignment2/pullseq -i protein_sequences/all_sequences.fa -g {} >> protein_sequences/input.fa'.format(acc))

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
	
	
	for acc in filtered_acc:
		os.system('/localdisk/data/BPSM/Assignment2/pullseq -i protein_sequences/all_sequences.fa -g {} >> protein_sequences/input.fa'.format(acc))
elif choice == '4':
	for current_species in non_redundant_species:
		if input('There are {} occurrances of {} in this dataset. Do you want to include the sequences in the analysis?\n Press y if yes, any other key will discard the sequence(s)'.format(species.count(current_species),current_species)) =='y':
			keys = [key for key,value in acc_spp.items() if value == current_species]
			for key in keys:
				os.system('/localdisk/data/BPSM/Assignment2/pullseq -i protein_sequences/all_sequences.fa -g {} >> protein_sequences/input.fa'.format(key))

else:
	print('That wasn\'t one of the options...')
	print('Exiting analysis')
	sys.exit()

#use clustalo to align
print('Aligning sequences + building phylogenetic tree...')
os.system('clustalo --threads=10 -i protein_sequences/input.fa -o protein_sequences/aligned_sequences.fa --guidetree-out=clustalo_tree')


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
if winsize is None:
	winsize = input('Please specify a window size for the conservation plot (4 is recommended). \nThe same window size will be used if you decide to perform more analysis:\n\t')
if winsize.isnumeric() == False:
	print('That was not a number...\nExiting analysis')
	sys.exit()
if int(winsize) > 30:
	print('WARNING:\nThis window size is large, conservation plot might not be accurate if sequences are short')
if int(winsize) < 4:
	print('WARNING:\nThis window size is low, conservation plot might be noisy')

print('Producing conservation plot. Graph will appear on screen, but will also be saved as plotcon.svg')
os.system('plotcon -sequences protein_sequences/aligned_sequences.fa -winsize {} -graph x11 -graph svg && plotcon -sequences protein_sequences/aligned_sequences.fa -winsize {} -graph x11'.format(winsize,winsize))
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
		acc = acc.split('|')[1]
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
	motifs = open("reports/motifs.txt", "a")
	elements = []
	with open(file,'r') as f:
		for line in f:
			pattern1 = 'Sequence'
			if re.search( pattern1, line):
				elements.append(file.split('/')[1][:-16]) #the accession value
			pattern2 = 'Motif'
			if re.search(pattern2,line):
				elements.append(line.split()[2]) # the motif
	for element in elements:
		motifs.write(element + '\t')
	motifs.write('\n')
print('Done. The motifs for each sequence can be found in reports/motifs.txt')

#Option to perform tmap analysis
if input("Do you want to predict and plot possible transmembrane segments in the sequences (tmap, small datasets might result in no hits)? If yes, press y:\n\t") =='y':
	os.system('tmap -sequences protein_sequences/aligned_sequences.fa -graph svg -outfile reports/report.tmap -auto && tmap -sequences protein_sequences/aligned_sequences.fa -graph x11 -outfile reports/report.tmap -auto')
#iep analysis
if input('Would you like to calculate the isoelectric point of the chosen proteins? If yes, press y:\n\t') == 'y':
	os.system('iep -sequence protein_sequences/aligned_sequences.fa -outfile reports/report.iep')
#hmoment analysis
if input('Would you like to calculate the hydrophobic moment for the chosen proteins? If yes, press y\n\t') == 'y':
	os.system('hmoment -seqall protein_sequences/aligned_sequences.fa  -outfile reports/report.hmoment ')
#pepwindowall analysis
if input('Would you like to plot a superimposed hydropathy plot for the chosen sequences (not recommended for large datasets)? If yes, press y\n\t') == 'y':
	os.system('pepwindowall -sequences protein_sequences/aligned_sequences.fa -window {} -gxtitle="Residue Position" -gytitle="Hydropathy" -graph svg && pepwindowall -sequences protein_sequences/aligned_sequences.fa -window {} -gxtitle="Residue Position" -gytitle="Hydropathy" -graph x11 '.format(winsize,winsize))

#Summarising all the results
acc_motif = {}
for line in open("reports/motifs.txt", "r"):
	line = line.strip() # remove newline
	motifs = ','.join(line.split('\t')[1:]) #Joining all the motifd that were found into comma-separated string
	acc_motif[line.split('\t')[0]] = motifs


df = pd.read_csv("acc_spp.txt", sep='\t', names=['Accession_number','Species','Motifs','Isoelectric_point','Max_Hmoment'])
for index, row in df.iterrows():
	df.loc[index,'Motifs'] = acc_motif[df.loc[index,'Accession_number']]


try:
	acc_iep = {}
	for line in open("reports/report.iep","r"):
		if line.startswith('IEP'):
			acc = line.split(' ')[2]
		elif line.startswith('Isoelectric'):
			iep = line.strip().split(' ')[3]
			acc_iep[acc] = iep
	for index, row in df.iterrows():
		try:
			df.loc[index,'Isoelectric_point'] = acc_iep[df.loc[index,'Accession_number']]
		except KeyError:
			df.loc[index,'Isoelectric_point'] = 'NaN'
except FileNotFoundError:
		print('Isoelectric point was not calculated, skipping from final report')

try:
	acc_hmoment = {}
	for line in open('reports/report.hmoment','r'):
		if line.startswith('HMOMENT'):
			acc = line.split(' ')[2]
		elif line.startswith('Window'):
			hmoment = line.strip().split(' ')[-1]
			acc_hmoment[acc] = hmoment
	for index, row in df.iterrows():
		try:
			df.loc[index,'Max_Hmoment'] = acc_hmoment[df.loc[index,'Accession_number']]
		except KeyError:
			df.loc[index,'Max_Hmoment'] = 'NaN'

except FileNotFoundError:
		print('Hydrophobic moment was not calculated, skipping from final report')

df.to_csv('output.tsv',sep='\t')
