import os,subprocess,shutil
#input test
os.system("mkdir protein_sequences")
tax_grp = input("Please enter the taxonomig group you want to investigate:\n\t")
prot_fam = input("Please enter the desired protein family:\n\t")
#Doing an esearch, and piping it to efetch
print('Fetching Accession values for {} proteins in {}'.format(prot_fam,tax_grp))
print('Please quit now if these is not the desired query')
os.system('esearch -db protein -query "' + prot_fam + ' AND ' + tax_grp + '"| efetch -format acc > protein_sequences/acc_numbers.txt')
#read the accession numbers to get the sequences, and print a warning if longer than 1000
#probs should add a warning to see if you want to continue
acc_list=[]
for line in open("protein_sequences/acc_numbers.txt", "r"):
	acc_list.append(line.strip())
if (len(acc_list) > 1000):
    if input('There are more than 1000 sequences in this analysis, press y if you want to continue ') != 'y':
        sys.exit()
print('Fetching Protein Sequences...')
for acc in acc_list:
	os.system('wget -q -O protein_sequences/' + acc + '.fasta "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id="' + acc + '"&strand=1&rettype=fasta&retmode=text" > protein_sequences/' + acc + '.txt')

