# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from featurExtract.database.database import create_db


def get_gene(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    gene_seq = pd.DataFrame(columns=['GeneID','Chrom','Start','End','Strand','Seq'])
    index = 0
    if not args.gene:
        for g in db.features_of_type('gene', order_by='start'):
            seq = g.sequence(args.genome, use_strand=False)
            seq = Seq(seq)
            if g.strand == '-':
                seq = seq.reverse_complement()
            gene_seq.loc[index] = [g.id,g.chrom,g.start,g.end,g.strand,seq]
            index += 1
        gene_seq.to_csv(args.output, sep=',', index=False)
    else:
        for g in db.features_of_type('gene', order_by='start'):
            if args.gene in g.id:
                seq = g.sequence(args.genome, use_strand=False)
                seq = Seq(seq)
                if g.strand == '-':
                    seq = seq.reverse_complement()
                geneRecord = SeqRecord(seq, id=g.id, 
                             description='strand %s start %d end %d length=%d'%(
                                          g.strand, g.start, g.end, len(seq))
                                      )
                if args.print:
                    SeqIO.write([geneRecord], sys.stdout, "fasta") 
                else:
                    SeqIO.write([geneRecord], args.output, "fasta") 
                break 

def get_gene_gb(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    gene = []
    for record in create_db(args.genbank):
        for feature in record.features:
            if feature.type == 'gene': # CDS promoter UTR 
                gene_seq = ''
                for part in feature.location.parts:
                    gene_seq += part.extract(record.seq)
                # part.strand 会将FeatureLocation -1的反向互补
                # 判断提取的和已知的是否一致
                gene_id = feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else 'Null'
                gene_seq_record = SeqRecord(gene_seq, id='gene:%s'%(gene_id), 
                                            description='strand %s length %d'%(
                                            feature.strand, len(gene_seq))
                                           )
                gene.append(gene_seq_record)
    if args.print and args.format == 'dna':
        SeqIO.write(gene, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(gene, args.output, 'fasta')
