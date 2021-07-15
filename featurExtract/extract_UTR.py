# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict



def utr(args):
    '''
    parameters:
        args: parse from argparse
    return:
        file write to file or stdout
    '''
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    # header
    utr_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','UTR5','UTR3'])
    if not args.transcript:
        # all UTR in genome 
        index = 0
        for t in db.features_of_type('mRNA', order_by='start'):
            seq3, seq5 = '', ''
            # utr3
            for c in db.children(t, featuretype='three_prime_UTR', order_by='start'):
                s = c.sequence(args.genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                seq3 += s
            # utr5
            for c in db.children(t, featuretype='five_prime_UTR', order_by='start'):
                s = c.sequence(args.genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                seq5 += s
            seq3 = Seq(seq3)
            seq5 = Seq(seq5)
            if t.strand == '-':
                seq3 = seq3.reverse_complement()
                seq5 = seq5.reverse_complement()

            utr_seq.loc[index] = [t.id,t.chrom,t.start,t.end,t.strand,seq5,seq3]
            index += 1
        utr_seq.to_csv(args.output, sep=',', index=False)
    else:
        # return a specific transcript
        out = [] 
        for t in db.features_of_type('mRNA', order_by='start'):
            if args.transcript in t.id:
                seq3, seq5 = '', ''
                # utr3
                for c in db.children(t, featuretype='three_prime_UTR', order_by='start'):
                    s = c.sequence(args.genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                    seq3 += s
                # utr5
                for c in db.children(t, featuretype='five_prime_UTR', order_by='start'):
                    s = c.sequence(args.genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                    seq5 += s
                seq3 = Seq(seq3)
                seq5 = Seq(seq5)
                if t.strand == '-':
                    seq3 = seq3.reverse_complement()
                    seq5 = seq5.reverse_complement()
                seq3Record = SeqRecord(seq3,id=args.transcript, description='strand %s utr3 length=%d'%(t.strand, len(seq3)))
                seq5Record = SeqRecord(seq5,id=args.transcript, description='strand %s utr5 length=%d'%(t.strand, len(seq5)))
                out.append(seq3Record)
                out.append(seq5Record)
                if args.print:
                    SeqIO.write(out, sys.stdout, "fasta")
                else:
                    SeqIO.write(out, args.output, "fasta")
                break 



def utr_genbank(args):
    '''
    parameters:
        parse from argparse
    return:
        elements write to a file or stdout
    '''
    genbank = args.genbank
