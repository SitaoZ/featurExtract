# -*- coding: utf-8 -*-
import re
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from tqdm import tqdm
from Bio.Seq import Seq
from collections import deque
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from featurExtract.utils.util import parse_output
from featurExtract.database.database import genome_dict

# https://pythonhosted.org/gffutils/autodocs/gffutils.Feature.html
# 1-based coordinates

_CSV_HEADER = ['GeneID','Chrom','Start','End','Strand','Promoter']

def sub_promoter(params):
    """
    parameters:
        params: args for multiple processing 
    return:
        deque list 
    """
    promoter_seq_list = deque()
    g, database, genome, promoter_length, utr5_upper_length, output_format = params
    genomeDict = genome_dict(genome) # load fasta 
    db = gffutils.FeatureDB(database, keep_order=True)
    if g.strand == "+":
        gene_seq = g.sequence(genome, use_strand=True)
        p_start = g.start - promoter_length - 1 # 往前数,和terminator不同需要-1, transform 1-based to 0-based
        if p_start < 0 :
            p_start = g.start - 1
        p_end = g.start + utr5_upper_length - 1      # include first base of TSS, transform 1-based to 0-based
        promoter = genomeDict[g.chrom][p_start:p_end].seq # pyfaidx
    else:
        # g.strand == "-"
        gene_seq = g.sequence(genome, use_strand=True)
        p_start = g.end - utr5_upper_length - 1
        p_end = g.end + promoter_length - 1
        promoter = genomeDict[g.chrom][p_start:p_end].reverse.complement.seq
    p_start_in_genome = p_start
    p_end_in_genome = p_end
    if output_format == 'fasta':
        desc='strand:%s start:%d end:%d length=%d'%(g.strand, p_start_in_genome,
                                                 p_end_in_genome, len(promoter))
        promoter = Seq(promoter) # str to Seq
        promoterRecord = SeqRecord(promoter, id=g.id.replace('gene:',''), description=desc)
        promoter_seq_list.append(promoterRecord)
    else:
        # default csv
        it = [g.id.replace('gene:',''), g.chrom, p_start_in_genome, \
              p_end_in_genome, g.strand, promoter]
        promoter_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
    return promoter_seq_list

    
def get_promoter(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    feature_types = db.featuretypes()
    promoter_seq_list = deque()
    promoter_record = []
    index = 0 
    if not args.gene:
        # loop all gene
        param_list = [(g, args.database, args.genome, args.promoter_length, args.utr5_upper_length, args.output_format) for g in db.all_features(featuretype='gene', order_by="seqid")]
        with Pool(processes=args.process) as p:
            promoter_seq_list = list(tqdm(p.imap(sub_promoter, param_list), total=len(param_list), ncols = 80, desc='CDS Processing:'))
        promoter_seq_list = [d for de in promoter_seq_list if de != None for d in de]
        parse_output(args, promoter_seq_list)
    else:
        for g in db.all_features(featuretype='gene', order_by="seqid"):
            if args.gene in g.id:
                param = (g, args.database, args.genome, args.promoter_length, args.utr5_upper_length, args.output_format)
                promoter_seq_list = sub_promoter(param)
                parse_output(args, promoter_seq_list)
                break
