from Bio import SeqIO
import sys
import pandas as pd

def get_length(path):
    chrom = []
    chrom_length = []

    for seq_record in SeqIO.parse(path,"fasta"):
      chrom.append(str(seq_record.id))
      chrom_length.append(str(len(seq_record)))
    chr_length = pd.DataFrame({'chr': chrom,'start': 0,'end': chrom_length})
    return(chr_length)
