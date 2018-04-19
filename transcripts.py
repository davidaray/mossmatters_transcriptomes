
from Bio import SeqIO
import statistics as stat
import pandas as pd
import argparse

def get_args():
	parser = argparse.ArgumentParser(description="Python seminar assignment. Will parse transcriptome data.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--ifasta', type=str, help='Name of the peptide file', required=True)
	parser.add_argument('-b', '--blastx', type=str, help='Name of the blastx file', required=True)

	args = parser.parse_args()
	INFASTA = args.ifasta
	BLASTX = args.blastx
	
	return INFASTA, BLASTX
#	return INFASTA, PREFIX
INFASTA, BLASTX = get_args()

print("Input fasta = " + INFASTA + ".")
print("Input blastx = " + BLASTX + ".")

ID_LIST = []
LENLIST = []
IDLIST = []

with open('stats.txt', 'w') as STATS:


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

	for RECORD in SeqIO.parse(INFASTA, 'fasta'):
		IDLIST.append(RECORD.id.rsplit('_',1))

IDS = pd.DataFrame(IDLIST, columns = ['uni', 'iso'])
COUNTUNI = len(IDS['uni'].unique())
print('There are ' + str(COUNTUNI) + ' unigenes in the file.')





