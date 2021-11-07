import os,subprocess,shutil
#input test
os.system("mkdir protein_sequences")
tax_grp = input("Please enter the taxonomig group you want to investigate:\n\t")
prot_fam = input("Please enter the desired protein family:\n\t")
#Doing an esearch, and piping it to efetch
print('Fetching Accession values for {} proteins in {}'.format(prot_fam,tax_grp))
print('Please quit now if these is not the desired query')
os.system('esearch -db protein -query "' + prot_fam + ' AND ' + tax_grp + '"| efetch -format acc > protein_sequences/acc_numbers.txt')

