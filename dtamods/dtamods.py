import argparse
import csv
import re

parser = argparse.ArgumentParser()
parser.add_argument("dta_file",help="The input dta-select")
parser.add_argument("fasta_file",help="The input protein fasta file")
parser.add_argument("out_file",help="The csv output file")
args = parser.parse_args()

f=open(args.dta_file)
f2=open(args.fasta_file)

peptides = dict()
name =''
modifications = {'(114.042927)':'Ub','(79.9663)':'Phos','(15.9949)':'Ox','(42.0106)':'Ac'}
for line in f2:
    if line[0] ==">":
        name = line.split('|')[0].split(':')[1].split('_')[0]  #strip out the extra crap and just return the transcript name
    else:
        peptides[name]=line
with open(args.out_file,'wb') as fout:
    writer = csv.writer(fout)
    for line in f:
        if line[:2] =="UC":
            name = line.split('|')[0].split(':')[1].split('_')[0] #strip out the extra crap and just return the transcript name
            continue
        array = line.rstrip('\n').split("\t")
        if len(array)<14:  #easy way to bypass the header and footer lines
            continue
        peptide = array[14]  #this will just be the peptide sequence
        if peptide.find('(')<0:
            continue
        #peptide = re.sub("\([^\)]+\)","",peptide)
        #print peptide
        m=re.search("\([^\)]+\)",peptide)
        mods = list()
        while m is not None:
            outpep = list()
            outpep.append(peptide[:m.start(0)])
            outpep.append(peptide[m.end(0):])
            peptide = ''.join(outpep)
            mods.append((m.start(0)-2,m.group(0)))
            m=re.search("\([^\)]+\)",peptide)
        #for m in re.finditer("\([^\)]+\)",peptide):
        #    outpep.append(peptide[:m.start(0)])
        #    outpep.append(peptide[m.end(0):])
            #print '%02d-%02d: %s' % (m.start(), m.end(), m.group(0)) 
        #m=re.search("\([^\)]+\)",peptide)
        peptide = re.sub("\.","",peptide)
        try:
            pos=peptides[name].find(peptide)
        except KeyError:
            continue
        for position,mod in mods:
            writer.writerow([name,peptides[name][pos+position],pos+position,modifications[mod]])
    
    