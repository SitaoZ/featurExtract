# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
from io import StringIO 
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, deque
from featurExtract.utils.util import utr3_type, utr5_type, mRNA_type, parse_output



_CSV_HEADER = ['TranscriptID','Chrom', \
               'Start','End','Strand', \
               'UTR5','UTR5_Region','UTR5_Count', \
               'UTR3','UTR3_Region','UTR3_Count']

def sub_utr(params):
    t, database, genome, style, output_format, utr3_t, utr5_t = params
    db = gffutils.FeatureDB(database, keep_order=True)
    utr_seq_list = deque()
    seq3, seq5 = '', ''
    seq3_count, seq5_count = 0, 0
    seq3_region, seq5_region = [], []
    seq3_gff_lines = deque()
    seq5_gff_lines = deque()
    # utr3
    for c in db.children(t, featuretype=utr3_t, order_by='start'):
        seq3_gff_lines.append(c)
        s = c.sequence(genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
        seq3 += s
        seq3_count += 1
        seq3_region.append('-'.join(map(str,[c.start,c.end])))
    # utr5
    for c in db.children(t, featuretype=utr5_t, order_by='start'):
        seq5_gff_lines.append(c)
        s = c.sequence(genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
        seq5 += s
        seq5_count += 1
        seq5_region.append('-'.join(map(str,[c.start,c.end])))
    seq3 = Seq(seq3)
    seq5 = Seq(seq5)
    if t.strand == '-':
        seq3 = seq3.reverse_complement()
        seq5 = seq5.reverse_complement()
    if output_format == 'gff':
        for line in seq5_gff_lines:
            utr_seq_list.append(line)
        for line in seq3_gff_lines:
            utr_seq_list.append(line)
    elif output_format == 'fasta':
        if len(seq5) != 0:
            utr5_id = ".".join(map(str,[t.id, '|'.join(seq5_region), '5UTR']))
            desc='strand:%s UTR=5UTR length=%d'%(t.strand,len(seq5))
            utr5Record = SeqRecord(seq5, id=utr5_id.replace('transcript:',''), description=desc)
            utr_seq_list.append(utr5Record)
        if len(seq3) != 0:
            utr3_id = ".".join(map(str,[t.id, '|'.join(seq3_region), '3UTR']))
            desc='strand:%s UTR=3UTR length=%d'%(t.strand,len(seq3))
            utr3Record = SeqRecord(seq3, id=utr3_id.replace('transcript:',''), description=desc)
            utr_seq_list.append(utr3Record)
    else:
        it = [t.id.replace('transcript:',''), t.chrom, t.start, t.end, t.strand,\
              seq5, '|'.join(seq5_region), seq5_count,\
              seq3, '|'.join(seq3_region), seq3_count]
        utr_seq_list.append(dict(zip(_CSV_HEADER, it)))
    return utr_seq_list
    

def utr(args):
    '''
    parameters:
        args: parse from argparse
    return:
        file write to file or stdout
    '''
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    feature_types = db.featuretypes()
    # header
    # utr_seq = pd.DataFrame(columns=_CSV_HEADER)
    # dist fastest way 
    if args.rna_feature == 'mRNA':
        mRNA_str = mRNA_type(args.style)
    else:
        mRNA_str = tuple(x for x in list(feature_types) if 'RNA' in x or 'transcript' in x)
    utr3_t = utr3_type(args.style)
    utr5_t = utr5_type(args.style)
    if not args.transcript:
        # all UTR in genome 
        param_list = [(t, args.database, args.genome, args.style, args.output_format, utr3_t, utr5_t) for t in db.features_of_type(mRNA_str, order_by='start')]
        with Pool(processes=args.process) as p:
            utr_seq_list = list(tqdm(p.imap(sub_utr, param_list), total=len(param_list), ncols = 80, desc='UTR Processing:'))
        utr_seq_list = [d for de in utr_seq_list if de != None for d in de]
        # output
        parse_output(args, utr_seq_list)
    else:
        # return a specific transcript
        for t in db.features_of_type(mRNA_str, order_by='start'):
            if args.transcript in t.id:
                params = (t, args.database, args.genome, args.style, args.output_format, utr3_t, utr5_t)
                utr_seq_list = sub_utr(params)
                # output 
                parse_output(args, utr_seq_list)
                break 



def utr_genbank(args):
    '''
    parameters:
        parse from argparse
    return:
        elements write to a file or stdout
    '''
    genbank = args.genbank
