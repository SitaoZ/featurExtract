# -*- coding: utf-8 -*-
import os
import argparse
import gffutils
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def create(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    fn = args.gff
    db = gffutils.create_db(fn, dbfn='gff.db', force=True, keep_order=True,\
        disable_infer_genes=True, disable_infer_transcripts=True,\
        merge_strategy='merge', sort_attribute_values=True)
    return db

def genome_dict(genome_fasta_path):
    '''
    parameters:
     genome_fasta_path: genome reference of organism
    return:
     genome fasta dict
    '''
    genome = dict()
    for record in SeqIO.parse(genome_fasta_path, 'fasta'):
        genome[record.id] = record.seq
    return genome


def UTR(args):
    '''
    parameters:
     args: arguments from argparse 
    '''
    utr_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','UTR5','UTR3'])
    db = gffutils.FeatureDB('gff.db', keep_order=True)
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



def uORF(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    uORF_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','Type','uORF'])
    db = gffutils.FeatureDB('gff.db', keep_order=True)
    index = 0
    for t in db.features_of_type('mRNA', order_by='start'):
        t_seq = t.sequence(args.genome, use_strand=True)
        cds = ''
        for c in db.children(t, featuretype='CDS', order_by='start'):
            s = c.sequence(args.genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
            cds += s
        cds = Seq(cds)
        if t.strand == '-':
            cds = cds.reverse_complement()
        # uORF parse
        
def CDS(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    cds_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','CDS'])
    db = gffutils.FeatureDB('gff.db', keep_order=True)
    index = 0
    for t in db.features_of_type('mRNA', order_by='start'):
        seq = ''
        for c in db.children(t, featuretype='CDS', order_by='start'):
            s = c.sequence(args.genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
            seq += s
        seq = Seq(seq)
        if t.strand == '-':
            seq= seq.reverse_complement()
        cds_seq.loc[index] = [t.id,t.chrom,t.start,t.end,t.strand,seq]
        index += 1
    cds_seq.to_csv(args.output, sep=',', index=False)

def promoter(args):
    promoter_length = args.length
    utr_head_length = args.utr_head
    genome_path = args.genome
    gff_path = args.gff
    genome = genome_dict(genome_path)
    #db = create_db(gff_path)
    db = gffutils.FeatureDB('gff.db', keep_order=True)
    index = 0
    promoter_seq = pd.DataFrame(columns=['GeneID','Chrom','Start','End','Strand','Promoter'])
    for f in db.all_features(featuretype='gene', order_by="seqid"):
        chrom = f.chrom
        geneid = f.id
        strand = f.strand
        if f.strand == "+":
            gene_seq = f.sequence(genome_path, use_strand=True)
            p_start = f.start - int(promoter_length)
            if p_start < 0 :
                continue
            p_end = f.start + utr_head_length
            promoter = genome[f.chrom][p_start:p_end]
            # print(promoter)
        elif f.strand == "-":
            gene_seq = f.sequence(genome_path, use_strand=True)
            p_start = f.end - utr_head_length
            p_end = f.end + int(promoter_length)
            promoter = genome[f.chrom][p_start:p_end].reverse_complement()
            # print(promoter)
        p_start_in_genome = p_start
        p_end_in_genome = p_end
        promoter_seq.loc[index] = [geneid,chrom,p_start_in_genome,p_end_in_genome,strand,promoter]
        index += 1
    promoter_seq.to_csv(agrs.output, sep=',', index=False)

parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(help='sub-command help')
# create subcommand 
parser_create = subparsers.add_parser('create', help='create annotation database')
parser_create.add_argument('-g', '--gff', type=str, help='genome annotation file')
parser_create.set_defaults(func=create)
# promoter subcommand
parser_promoter = subparsers.add_parser('promoter', help='extract promoter in genome or gene')
parser_promoter.add_argument('-l', '--length', type=int, help='promoter length before TSS')
parser_promoter.add_argument('-u', '--utr_head', type=int, help='utr5 length after TSS')
parser_promoter.add_argument('-f', '--genome', type=str, help='genome fasta')
parser_promoter.add_argument('-o', '--output', type=str, help = 'output csv file path')
parser_promoter.add_argument('-v', '--version', help = 'promoterExtract version', action = "store_true")
parser_promoter.set_defaults(func=promoter)

# UTR subcommand
parser_utr = subparsers.add_parser('UTR', help='extract untranslated region sequence in genome or gene')
parser_utr.add_argument('-o', '--output', type=str, help='output file path')
parser_utr.set_defaults(func=UTR)

# UTR subcommand
parser_uORF = subparsers.add_parser('uORF', help='extract upper stream open reading sequence in genome or gene')
parser_uORF.add_argument('-o', '--output', type=str, help='output file path')
parser_uORF.set_defaults(func=uORF)

# CDS subcommand
parser_cds = subparsers.add_parser('CDS', help='extract coding sequence in genome or gene')
parser_cds.add_argument('-o', '--output', type=str, help='output file path')
parser_utr.set_defaults(func=CDS)

args = parser.parse_args()
print('runing ...')
args.func(args)
print('finished ...')



'''
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--length', type=int, help='promoter length before TSS')
    parser.add_argument('-u', '--utr_head', type=int, help='length after TSS')
    parser.add_argument('-f', '--genome', type=str, help='genome fasta')
    parser.add_argument('-g', '--gff', type=str, help='genome annotation file')
    parser.add_argument('-o', '--outdir', type=str, help='output directory')
    args = parser.parse_args()
    get_promoter(args.length, args.utr_head, args.genome, args.gff, args.outdir)
'''                                             
