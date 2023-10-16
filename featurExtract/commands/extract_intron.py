# -*- coding: utf-8 -*-
import sys
import gffutils
import itertools
import pandas as pd 
from tqdm import tqdm
from Bio import SeqIO
from Bio.Seq import Seq
from multiprocessing import Pool
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, deque
from featurExtract.database.database import genome_dict
from featurExtract.utils.util import utr3_type, utr5_type, mRNA_type, parse_output

_CSV_HEADER = ['TranscriptID', 'Intron', 'Chrom','Start_genome', \
               'End_genome', 'Strand', 'Sequence']

def sub_intron(params):
    '''
    parameters:
        params: argparse
    return:
        intron_seq_list, a deque
    '''
    t, database, genome, style, output_format = params
    genomeDict = genome_dict(genome)
    db = gffutils.FeatureDB(database, keep_order=True) # load database
    intron_seq_list = deque()
    exons = list(itertools.chain(*[[e.start,e.end] for e in db.children(t, featuretype='exon', order_by='start')])) # exon position include start and end
    if len(exons)/2 == 1:
        # only one exon
        return None
    else:
        left = exons[1:-1]
        intron_index = 0
        for i in range(0, len(left), 2):
            s = left[i]   # intron start
            e = left[i+1] # intron end, shoukd be minus 1 
            if s+1 != e:
                # exon 与exon 不相连，就是intron的位置
                # 注意索引的点，
                intron_index += 1
                intron_seq= genomeDict[t.chrom][s:e-1]
                if t.strand == '-':
                    intron_seq = intron_seq.reverse
                    intron_seq = intron_seq.complement
                    intron_rank = int(len(exons)/2 - intron_index)
                else:
                    intron_rank = intron_index
                if output_format == 'gff':
                    # 1       araport11       gene    3631    5899    .       +       .       ID=gene:AT1G01010;
                    intron_gff = [t.chrom, 'featurExtract', 'intron', s+1, e-1, '.', t.strand, '.', f'ID={t.id};Name=intron{intron_rank};']
                    intron_gff_line = list(map(str, intron_gff))
                    intron_seq_list.append("\t".join(intron_gff_line))
                    continue
                elif output_format == 'fasta':
                    intron_seq= Seq(str(intron_seq))
                    # 1 based position
                    intronRecord = SeqRecord(intron_seq, id=t.id, 
                                             description='strand %s intron%d start %d end %d length=%d exons=%d'%(
                                             t.strand, intron_rank, s+1, e-1, len(intron_seq), len(exons)/2)
                                            )
                    intron_seq_list.append(intronRecord)
                else:
                    intron_seq= str(intron_seq)
                    it = [t.id.replace('transcript:',''), f'intron{intron_rank}',t.chrom, s+1, e-1, t.strand, intron_seq]
                    intron_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
        return intron_seq_list

def get_intron(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database 
    genomeDict = genome_dict(args.genome) # load fasta
    # assert GTF or GFF
    if args.rna_feature == 'mRNA':
        mRNA_str = mRNA_type(args.style)
    else:
        mRNA_str = tuple(x for x in list(feature_types) if 'RNA' in x or 'transcript' in x)

    intron_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','Exon']) # header 
    if args.transcript:
        # return a specific transcript
        for t in tqdm(db.features_of_type(mRNA_str, order_by='start'),\
                      total = len(list(db.features_of_type(mRNA_str, order_by='start'))), \
                      ncols = 80, desc = "Intron Processing"):
            if args.transcript in t.id:
                param = (t, args.database, args.genome, args.style, args.output_format)
                intron_seq_list = sub_intron(param)
                parse_output(args, intron_seq_list)
                break
    else:
        param_list = [(t, args.database, args.genome, args.style, args.output_format) for t in db.features_of_type(mRNA_str, order_by='start')]
        with Pool(processes=args.process) as p:
            intron_seq_list = list(tqdm(p.imap(sub_intron, param_list), total=len(param_list), ncols = 80, desc='Intron Processing:'))
        intron_seq_list = [d for de in intron_seq_list if de != None for d in de]
        parse_output(args, intron_seq_list)
