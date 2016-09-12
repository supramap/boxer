#!/bin/python -tt

import os
import re
import collections
import sys
import argparse
import shutil

'''
Takes output of TrimAll program and finds the best alignment. Alignment data is written to an output file.
'''



parser = argparse.ArgumentParser(description='This method takes in two directories, '
						'one with original files in phylip format and another separate '
						'directory that contains the files after the program (I don\'t '
						'know the name of it atm) has been run. The result is a file '
						'that shows data such as number of gaps, residues, taxa, and some '
						'other important data. '
						'(default is reAlignData.txt with information comparing the before '
						'and after runs. ')
parser.add_argument('-orig','--originalDirectory', metavar='originalDirectory',
                    action='store', 
                    help='Directory before the program ran. The files must be in PHY format, '
					'with taxa and chars the first line separated by a space, otherwise '
					'a NullPointerException will be thrown.',
					required=True)
parser.add_argument('-reAl','--realignedDirectory', metavar='realignedDirectory',
                    action='store',
                    help='The reAligned directory',
					required=True)	
parser.add_argument('--output', metavar='outputFile',
                    type=str, default='reAlignData.txt',
                    help='Output file that is written and holds following information: '
					'Number of gaps before and after, Number of chars before and after, '
					'and file names, all on the same line in a tab separated file.')
parser.add_argument('--taxacutoff',
                    help='Filters alignments below a minimum UNIQUE taxa threshold in the '
					'realigned files. NOTE: If the chosen taxacutoff is too high, the program '
					'will exit out AFTER the output file is generated, which can take a long time '
					'depending on number of files in the directories!!! Select a good taxacutoff FIRST '
					'in order to avoid wasting time!!!'
					,required=True)
parser.add_argument('-f', '--folder', metavar='outputFolder',
                    type=str, default='best_alignment',
					help='output folder, where the best alignment files are copied to. Files will '
					'not be copied here unless the "taxacutoff" flag is initialized.'
					'(default = best_alignment)')
parser.add_argument('--gappyness',
                    help='Filters alignments below a minimum gap percent threshold in the '
					'realigned files. NOTE: If the chosen gap percent is too high, the program '
					'will exit out AFTER the output file is generated. '
					,required=True)
args = parser.parse_args()


dirpath1=args.originalDirectory
dirpath2=args.realignedDirectory
output=args.output
outputdir=args.folder
taxaCO=args.taxacutoff
gapCO=args.gappyness


def getMismatch(dirpath1, dirpath2, output):
	'''
	The main engine of the program, that calls the getGapNum() and getCharNum() functions.
	The is also the function that writes the file. It goes through both the original and the realigned file directories
	and outputs data from each of files, to an output file.
	NOTE: If you want to change the output, THIS is the function you would change.
	'''
	
	firsttime=True
	outFile = open(output, "w")
	outFile.write("#uniqueTaxa"+"\t"+"origLength"+"\t"+"numGapsO"+"\t"+"taxa*char"+"\t"+"original"+"\t"+"uniqueTaxa"+"\t"+"reAlLength"+"\t"+"gapNumsR"+"\t"+"taxa*char"+"\t"+"gapPercent"+"\t"+"realign"+"\n")
	for file1 in os.listdir(dirpath1):
		tokstring1=str(file1).split(".")
		for file2 in os.listdir(dirpath2):
			tokstring2=str(file2).split(".")
			if tokstring2[0]==tokstring1[0]:
				phyfile1=dirpath1+"/"+file1
				phyfile2=dirpath2+"/"+file2
				with open(phyfile1, 'r') as f1, open(phyfile2, 'r') as f2:
					info2=f2.readline().split(" ")
					info1=f1.readline().split(" ")
					gapOrig=getGapNum(file1,dirpath1)
					uniqueTaxaO=countUniqueTaxa(file1,dirpath1)
					outFile.write(str(uniqueTaxaO)+"\t"+info1[2].strip("\n")+"\t"+str(gapOrig)+"\t"+str(int(info1[1].strip())*int(info1[2].strip()))+"\t"+file1)
					outFile.write("\t")
					gapReal=getGapNum(file2,dirpath2)
					gapPerc=((float(gapReal)/float(int(info2[1])*int(info2[2])))*100)
					uniqueTaxaR=countUniqueTaxa(file2,dirpath2)
					if info1[1]==info2[1]:
						outFile.write(str(uniqueTaxaR)+"\t"+info2[1].strip()+"\t"+str(gapReal)+"\t"+str(int(info2[1].strip())*int(info2[2].strip()))+"\t"+str(gapPerc)+"\t"+file2)
						outFile.write("\n")
					elif info1[2]==info2[2]:
						outFile.write(str(uniqueTaxaR)+"\t"+info2[2].strip()+"\t"+str(gapReal)+"\t"+str(int(info2[1].strip())*int(info2[2].strip()))+"\t"+str(gapPerc)+"\t"+file2)
						outFile.write("\n")
					else:
						outFile.write(str(uniqueTaxaR)+"\t"+info2[2].strip()+"\t"+str(gapReal)+"\t"+str(int(info2[1].strip())*int(info2[2].strip()))+"\t"+str(gapPerc)+"\t"+file2)
						outFile.write("\n")
					if firsttime==True:
						firsttime=False
		outFile.write("####"+"\n")
	outFile.close()

	
def getGapNum(file,dirpath):
	'''
	Counts the total number of gaps in the file
	'''
	
	dashcount=0
	firstline=True
	phyfile=dirpath+"/"+file
	for line in open(phyfile,"r"):
		if firstline==True:
			firstline=False
		else:
			for x in line:
				if x=="-":
					dashcount=dashcount+1
	return dashcount


def countUniqueTaxa(file,dirpath):
	'''
	Counts the number of unique taxa in each realign file.
	NOTE: bj22|555 is equal to bj22|666. This function considers uniqueness BEFORE the pipe!!!
	'''
	
	tok=''
	dupFreeList=set()
	firstline=True
	phyfile=dirpath+"/"+file
	for line in open(phyfile,"r"):
		if firstline==True:
			firstline=False
		else:
			if re.match("^[a-z]", line):
				tok=line.split("|")
				dupFreeList.add(tok[0])
	return len(dupFreeList)
	

def cutOff(dirpath2,taxaCO,output,outputdir,gapCO):
	'''
	Filters by minimum number of taxa and gap percent, then moves the file to the best alignment folder (or whatever is set in the -f flag)
	'''	
	
	taxaCutoff=int(taxaCO)
	fortaxa,storedname,firstime,tuplist,tok,tok7,tok8 =True,'',True,[],'','',''
	for line in open(output):
		if not line.startswith("#"):
				tok=line.split("\t")
				taxaNum=tok[5]
				gapPerc=tok[9]
				taxa=taxaNum
				gap=gapPerc
				if int(taxa)>int(taxaCutoff) and float(gap)<float(gapCO):	
					tup=()
					if tok[4] != storedname and firstime==False: 
						if firstime==False:
							tok=line.split("\t")
							reAlname=tok[10]
							tok7=tok[7]
							tok8=tok[8]
							mathForRealign=(float(tok7)/int(tok8))
							tup=(mathForRealign,reAlname)
							stored = min(tuplist, key=lambda x: x[0])
							nametomove=stored[1]
							if (firstsamp==True):
								tuplist[:]=[]
								phyfile=dirpath2+"/"+nametomove.strip("\n")
								shutil.copy(phyfile,outputdir)
								storedname=tok[4]
								firstime=True
								firstsamp=False
								tuplist.append(tup)
					else:
						tok=line.split("\t")
						reAlname=tok[10]
						tok7=tok[7]
						tok8=tok[8]
						mathForRealign=(float(tok[7])/int(tok[8]))
						tup=(mathForRealign,reAlname)
						tuplist.append(tup)
						firstime=False
		if line.startswith("#"):
			firstsamp=True
		if firstime==False and int(taxa)>int(taxaCutoff):				
			storedname=tok[4]
	if not tok7:
		print "#####ERROR####"
		print "Chosen taxa cut-off is higher than the maximum unique taxa number. Choose a lower cut-off!!!"
		print "Exiting....."
		exit()
	mathForRealign=(float(tok7)/int(tok8))
	tup=(mathForRealign,reAlname)
	tuplist.append(tup)
	stored = min(tuplist, key=lambda x: x[0])
	nametomove=stored[1]
	phyfile=dirpath2+"/"+nametomove.strip("\n")
	shutil.copy(phyfile,outputdir)

	
'''
Executing the program...
'''

if not os.path.exists(args.folder):
	os.makedirs(args.folder)
		

getMismatch(dirpath1,dirpath2,output)

if args.taxacutoff:
	cutOff(dirpath2,taxaCO,output,outputdir,gapCO)
