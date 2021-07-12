# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

def utr(db, genome, transcript_id ,output):
    '''
    parameters:
     db : database create by gffutils 
     genome: genome fasta file path 
     transcript_id: transcript id
     output: output file path  
    '''
    # header
    utr_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','UTR5','UTR3'])
    if not transcript_id:
        # all UTR in genome 
        index = 0
        for t in db.features_of_type('mRNA', order_by='start'):
            seq3, seq5 = '', ''
            # utr3
            for c in db.children(t, featuretype='three_prime_UTR', order_by='start'):
                s = c.sequence(genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                seq3 += s
            # utr5
            for c in db.children(t, featuretype='five_prime_UTR', order_by='start'):
                s = c.sequence(genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                seq5 += s
            seq3 = Seq(seq3)
            seq5 = Seq(seq5)
            if t.strand == '-':
                seq3 = seq3.reverse_complement()
                seq5 = seq5.reverse_complement()

            utr_seq.loc[index] = [t.id,t.chrom,t.start,t.end,t.strand,seq5,seq3]
            index += 1
        utr_seq.to_csv(output, sep=',', index=False)
    else:
        # return a specific transcript
        out = [] 
        for t in db.features_of_type('mRNA', order_by='start'):
            print(t.id)
            if transcript_id in t.id:
                seq3, seq5 = '', ''
                # utr3
                for c in db.children(t, featuretype='three_prime_UTR', order_by='start'):
                    s = c.sequence(genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                    seq3 += s
                # utr5
                for c in db.children(t, featuretype='five_prime_UTR', order_by='start'):
                    s = c.sequence(genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                    seq5 += s
                seq3 = Seq(seq3)
                seq5 = Seq(seq5)
                if t.strand == '-':
                    seq3 = seq3.reverse_complement()
                    seq5 = seq5.reverse_complement()
                seq3Record = SeqRecord(seq3,id=transcript_id, description='utr3 length=%d'%(len(seq3)))
                seq5Record = SeqRecord(seq5,id=transcript_id, description='utr5 length=%d'%(len(seq5)))
                out.append(seq3Record)
                out.append(seq5Record)
                break 
        SeqIO.write(out, output, "fasta")
