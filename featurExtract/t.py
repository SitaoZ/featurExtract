import gffutils
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
from extract_uORF import uorf 

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
    for t in db.features_of_type('intron', order_by='start'):
        # primary transcript (pt) 是基因组上的转录本的序列，
        # 有的会包括intron，所以提取的序列和matural transcript 不一致
        #print(t)
        #print(t.id)
        pt = t.sequence(args, use_strand=True)
        # matural transcript (mt)
        # exon 提取的是转录本内成熟的MRNA的序列,即外显子收尾相连的matural transcript
        mt = ''
        for e in db.children(t, featuretype='intron', order_by='start'):
            s = e.sequence(args, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
            mt += s
        mt = Seq(mt)
        if t.strand == '-':
            mt = mt.reverse_complement()
        # CDS 提取的是编码区 ATG ... TAG
        cdna = ''
        for c in db.children(t, featuretype='intron', order_by='start'):
            #print(c)
            s = c.sequence(args, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
            cdna += s
        cdna = Seq(cdna)
        if t.strand == '-':
            cdna = cdna.reverse_complement()
        #print(pt)
        #print(mt)
        #print(cdna)
        #uORF_dict = uorf(t.id, mt, cdna)
        uORF_dict = uorf(t.id, t.chrom, t.strand, mt, cdna)
        #print(t.id)
        print(pt)
        #for key in uORF_dict:
        #    print(key,len(uORF_dict[key]))
        #    for it in uORF_dict[key]:
        #        print(it)
        #print(uORF_dict)
        index += 1
        if index == 8:
        #if index == 2:
            break
        
        
uORF('ath.fa')
