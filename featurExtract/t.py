import gffutils
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO


def genome_dict(genome_fasta_path):
    genome = dict()
    for record in SeqIO.parse(genome_fasta_path, 'fasta'):
        genome[record.id] = record.seq
    return genome

def CDS(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    genome_path = args
    genome = genome_dict(genome_path)
    cds_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','CDS'])
    db = gffutils.FeatureDB('gff.db', keep_order=True)
    index = 0
    for t in db.features_of_type('mRNA', order_by='start'):
        seq = ''
        for c in db.children(t, featuretype='three_prime_UTR', order_by='start'):
            print(c)
            s = c.sequence(args, use_strand=False)
            seq += s
        seq = Seq(seq)
        if t.strand == '-':
            seq = seq.reverse_complement()
        
        print(seq)
        cds_seq.loc[index] = [t.id,t.chrom,t.start,t.end,t.strand,seq]
        index += 1
        if index == 2:
            break
    cds_seq.to_csv('cds.csv',sep=',',index=False)


def uORF(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    uORF_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','Type','uORF'])
    db = gffutils.FeatureDB('gff.db', keep_order=True)
    index = 0
    for t in db.features_of_type('mRNA', order_by='start'):
        print(t)
        t_seq = t.sequence(args, use_strand=True)
        cdna = ''
        for c in db.children(t, featuretype='CDS', order_by='start'):
            print(c)
            s = c.sequence(args, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
            cdna += s
        cdna = Seq(cdna)
        if t.strand == '-':
            cdna = cdna.reverse_complement()
        print(t_seq)
        print(cdna)
        index += 1
        if index == 2:
            break
        
uORF('ath.fa')
