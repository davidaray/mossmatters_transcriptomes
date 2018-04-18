
from Bio import SeqIO
import statistics as stat
import argparse

def get_args():
	parser = argparse.ArgumentParser(description="Python seminar assignment. Will parse transcriptome data.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i', '--ifasta', type=str, help='Name of the peptide file', required=True)
	parser.add_argument('-b', '--blastx', type=str, help='Name of the blastx file', required=True)

	args = parser.parse_args()
	INFASTA = args.ifasta
	BLASTP = args.blastp
	
	return INFASTA, BLASTP
#	return INFASTA, PREFIX
INFASTA, BLASTP = get_args()

print("Input fasta = " + INFASTA + ".")
print("Input blastp = " + BLASTP + ".")

STATS = open('stats.txt', 'w')

#TD = SeqIO.to_dict(SeqIO.parse(INFASTA, "fasta"))

for RECORD in SeqIO.to_dict(SeqIO.parse(INFASTA, "fasta")):
   LIST.append(len(RECORD.seq))
   
print('Total bases = ' + str(sum(LIST)))
STATS.write('Total bases = ' + str(sum(LIST)))
print('Total transcripts = ' + str(len(LIST)))
STATS.write('Total transcripts = ' + str(len(LIST)))
print('Mean length of transcripts = ' + str(stat.mean(LIST)))
STATS.write('Mean length of transcripts = ' + str(stat.mean(LIST)))
print('Median length of transcripts = ' + str(stat.median(LIST)))
STATS.write('Median length of transcripts = ' + str(stat.median(LIST)))



