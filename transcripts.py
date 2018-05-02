
from Bio import SeqIO
import statistics as stat
import pandas as pd
import argparse

def get_args():
	parser = argparse.ArgumentParser(description="Python seminar assignment. Will parse transcriptome data.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--ifasta', type=str, help='Name of the peptide file', required=True)
	parser.add_argument('-b', '--blastx', type=str, help='Name of the blastx file', required=True)
	parser.add_argument('-p', '--prefix', type=str, help='A prefix to label the output', required=True)

	args = parser.parse_args()
	INFASTA = args.ifasta
	BLASTX = args.blastx
	PREFIX = args.prefix

	return INFASTA, BLASTX, PREFIX
INFASTA, BLASTX, PREFIX = get_args()

print("Input fasta = " + INFASTA + ".")
print("Input blastx = " + BLASTX + ".")
print("Prefix = " + PREFIX + ".")

####PHASE 1
ID_LIST = []
LENLIST = []
IDLIST = []

with open(PREFIX + '.stats.txt', 'w') as STATS:


	for RECORD in SeqIO.parse(INFASTA, 'fasta'):
		LENLIST.append(len(RECORD.seq))
			
	print('Total bases = ' + str(sum(LENLIST)) + '.' )
	STATS.write('Total bases = ' + str(sum(LENLIST)) + '.\n')
	print('Total transcripts = ' + str(len(LENLIST)) + '.')
	STATS.write('Total transcripts = ' + str(len(LENLIST)) + '.\n')
	print('Mean length of transcripts = ' + str(stat.mean(LENLIST)) + '.')
	STATS.write('Mean length of transcripts = ' + str(stat.mean(LENLIST)) + '.\n')
	print('Median length of transcripts = ' + str(stat.median(LENLIST)) + '.')
	STATS.write('Median length of transcripts = ' + str(stat.median(LENLIST)) + '.\n')
	
	LENLIST = sorted(LENLIST)
	HALF = sum(LENLIST)/2
	TOTALENTRIES = len(LENLIST)
	for ENTRY in range(1, TOTALENTRIES):
		RUNNINGTOTAL = sum(LENLIST[0:ENTRY])
		if RUNNINGTOTAL >= HALF:
			print('N50 = ' + str(LENLIST[ENTRY-1]) + '.')
			STATS.write('N50 = ' + str(LENLIST[ENTRY-1]) + '.\n')
			break

	for RECORD in SeqIO.parse(INFASTA, 'fasta'):
		IDLIST.append(RECORD.id.rsplit('_',1))

	IDS = pd.DataFrame(IDLIST, columns = ['uni', 'iso'])
	COUNTUNI = len(IDS['uni'].unique())
	print('There are ' + str(COUNTUNI) + ' unigenes in the file.')
	STATS.write('There are ' + str(COUNTUNI) + ' unigenes in the file.\n')

####PHASE 2
#Task 1 - Number of transcripts with BLAST hits in a reference proteome
	#Read in blastx file as dataframe
	BLASTHITS = pd.read_table(BLASTX, sep='\t', usecols=[1,2], names=['ID', 'db'])
	NUMHITS = len(BLASTHITS.ID.unique())
	print('Transcripts with hits = ' + str(NUMHITS) + '.')

#Task 2 - Collapse factor (average number of transcripts that match the same protein sequence)	
	#Convert to a dataframe with the database hits and how many times they're hits
	DBCOUNTS = BLASTHITS.db.value_counts().reset_index().rename(columns={'index': 'BLASTHITS', 0: 'count'})
	#keep only the counts greater than 1
	DBCOUNTSKEEP = DBCOUNTS[DBCOUNTS['db'] > 1]
	
	#find the mean of the counts columns
	COLLAPSE = DBCOUNTSKEEP['db'].mean()
	print('Collapse = ' + str(COLLAPSE) + '.')
	
	
	


