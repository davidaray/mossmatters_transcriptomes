
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

STATS = open('stats.txt', 'w')

#TD = SeqIO.to_dict(SeqIO.parse(INFASTA, "fasta"))

LENLIST = []
IDLIST = []
ID_LIST = []
for RECORD in SeqIO.parse(INFASTA, "fasta"):
	LENLIST.append(len(RECORD.seq))
	IDLIST.append(RECORD.id)
	for LINE in IDLIST:
		LINE = LINE.rsplit("_",1)
		ID_LIST.append(LINE)
			
#print(LIST)
print('Total bases = ' + str(sum(LENLIST)))
STATS.write('Total bases = ' + str(sum(LENLIST)))
print('Total transcripts = ' + str(len(LENLIST)))
STATS.write('Total transcripts = ' + str(len(LENLIST)))
print('Mean length of transcripts = ' + str(stat.mean(LENLIST)))
STATS.write('Mean length of transcripts = ' + str(stat.mean(LENLIST)))
print('Median length of transcripts = ' + str(stat.median(LENLIST)))
STATS.write('Median length of transcripts = ' + str(stat.median(LENLIST)))


IDS = pd.DataFrame(ID_LIST, columns = ['uni', 'iso'])
#print(IDS)
COUNTUNI = len(IDS['uni'].unique())
print('There are ' + str(COUNTUNI) + 'unigenes in the file.')





