# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, deque
from featurExtract.utils.util import mRNA_type, parse_output
from featurExtract.database.database import create_db

_CSV_HEADER = ['TranscriptID','Chrom','Start','End','Strand','Exon', 'Exon_length', 'Exon_index']

def sub_exon(params):
    """
    parameters:
        params:
    return:
        exon_seq_list, a deque class 
    """
    t, database, genome, output_format = params
    db = gffutils.FeatureDB(database, keep_order=True)
    exon_seq_list = deque()
    exon_index = 0
    for e in db.children(t, featuretype='exon', order_by='start'):
        if output_format == 'gff':
            exon_seq_list.append(e)
            continue
        exon = e.sequence(genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
        exon = Seq(exon)
        if t.strand == '-':
            exon = exon.reverse_complement()
        exon_index += 1
        # fasta
        if output_format == 'fasta':
            exonRecord = SeqRecord(exon,id=t.id,
                                description='strand %s exon %d start %d end %d length=%d'%(t.strand,
                                exon_index, e.start, e.end, len(exon)))
            exon_seq_list.append(exonRecord)
    
        else:
            it = [t.id.replace('transcript:',''), t.chrom, e.start, e.end, t.strand,\
                  exon, len(exon), exon_index]
            exon_seq_list.append(dict(zip(_CSV_HEADER, it)))
    return exon_seq_list    



def get_exon(args):
    '''
    parameters:
     
    '''
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    feature_types = db.featuretypes()
    if args.rna_feature == 'mRNA':
        mRNA_str = mRNA_type(args.style)
    else:
        mRNA_str = tuple(x for x in list(feature_types) if 'RNA' in x or 'transcript' in x)
    if args.transcript:
        for t in db.features_of_type(mRNA_str, order_by='start'):
            if args.transcript in t.id:
                param = (t, args.database, args.genome, args.output_format)
                exon_seq_list = sub_exon(param)
                parse_output(args, exon_seq_list)
                break 
    else:
        param_list = [(t, args.database, args.genome, args.output_format) for t in db.features_of_type(mRNA_str, order_by='start')]
        # param_list = param_list[:10]
        with Pool(processes=args.process) as p:
            exon_seq_list = list(tqdm(p.imap(sub_exon, param_list), total=len(param_list), ncols = 80, desc='Exon Processing:'))
        exon_seq_list = [d for de in exon_seq_list if de != None for d in de]
        # output
        parse_output(args, exon_seq_list)

def get_exon_gb(args):
    exons = []
    for record in create_db(args.genbank):
        index = 1
        print(type(record))
        for feature in record.features:
            if feature.type == 'exon': # CDS promoter UTR 
                exon_seq = feature.extract(record.seq)
                # feature.strand 会将FeatureLocation -1的反向互补
                # 判断提取的和已知的是否一致
                exon_seq_record = SeqRecord(exon_seq,id='%s exon %d'%(record.id, index),
                                 description='strand %s length %d'%(feature.strand, len(exon_seq)))
                exons.append(exon_seq_record)
                index += 1 
    if args.print:
        SeqIO.write(exons, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(exons, args.output, 'fasta')
 
