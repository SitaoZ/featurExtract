# -*- coding: utf-8 -*-
import os
import time
import argparse
import gffutils
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from featurExtract.extract_UTR import utr
from featurExtract.extract_uORF import get_uorf
from featurExtract.extract_dORF import get_dorf
from featurExtract.extract_CDS import get_cds
from featurExtract.extract_promoter import get_promoter
from featurExtract.extract_exon import get_exon
from featurExtract.extract_intron import get_intron
from featurExtract.extract_gene import get_gene
from featurExtract.extract_IGR import get_IGR
from featurExtract.extract_cdna import get_cdna
from featurExtract.extract_mRNA import get_mRNA


def create(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    fn = args.genomefeature
    database_id = args.output_prefix +'.'+ args.file_type
    db = gffutils.create_db(fn, dbfn=database_id, force=True, keep_order=True,\
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
    utr(args)


def uORF(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_uorf(args)
    
 
def CDS(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_cds(args)


def dORF(args):
    '''
    parameters:
     args: arugmensts from argparse
    '''
    get_dorf(args)


def exon(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_exon(args)

def intron(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_intron(args)


def promoter(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_promoter(args)

def gene(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_gene(args)

def mRNA(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_mRNA(args)

def cdna(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_cdna(args)

def IGR(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_IGR(args)

parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(help='sub-command help')
# create subcommand 
parser_create = subparsers.add_parser('create', help='create annotation database')
parser_create.add_argument('-f', '--file_type', choices=['GFF','GTF'],
                           help='genome annotation file')
parser_create.add_argument('-g', '--genomefeature', type=str, required=True, 
                           help='genome annotation file')
parser_create.add_argument('-o', '--output_prefix', type=str, required=True, 
                           help='database absolute path')
parser_create.set_defaults(func=create)

# promoter subcommand
parser_promoter = subparsers.add_parser('promoter', help='extract promoter in genome or gene')
parser_promoter.add_argument('-d', '--database', type=str, required=True, 
                             help='database generated by subcommand create')
parser_promoter.add_argument('-f', '--genome', type=str, required=True,
                             help='genome fasta path')
parser_promoter.add_argument('-g', '--gene', type=str, 
                             help='specific gene; if not given, return whole genes')
parser_promoter.add_argument('-l', '--promoter_length', type=int, default=0,
                             help='promoter length before TSS')
parser_promoter.add_argument('-u', '--utr5_upper_length', type=int, default=0,
                             help='utr5 length after TSS')
parser_promoter.add_argument('-o', '--output', type=str, 
                             help = 'output file path')
parser_promoter.add_argument('--output_format', type=str, choices=['csv','fasta'], 
                             help = 'output format')
parser_promoter.add_argument('-p', '--print', action="store_true", 
                             help = 'output to stdout')
parser_promoter.set_defaults(func=promoter)

# gene subcommand 
parser_gene = subparsers.add_parser('gene', help='extract gene in genome or gene')
parser_gene.add_argument('-d', '--database', type=str, required=True,
                         help='database generated by subcommand create')
parser_gene.add_argument('-f', '--genome', type=str, required=True,
                         help='genome fasta')
parser_gene.add_argument('-g', '--gene', type=str, 
                         help='specific gene; if not given, return whole genes')
parser_gene.add_argument('-o', '--output', type=str, 
                         help = 'output file path')
parser_gene.add_argument('-p', '--print', action="store_true", 
                         help='output to stdout')
parser_gene.set_defaults(func=gene)

# mRNA subcommand 
parser_mRNA = subparsers.add_parser('mRNA', help='extract mature messager RNA in genome or gene')
parser_mRNA.add_argument('-d', '--database', type=str, required=True, 
                         help='database generated by subcommand create')
parser_mRNA.add_argument('-f', '--genome', type=str, required=True, 
                         help='genome fasta')
parser_mRNA.add_argument('-t', '--transcript', type=str, 
                         help='specific transcript; if not given, return whole transcripts')
parser_mRNA.add_argument('-o', '--output', type=str, 
                         help = 'output file path')
parser_mRNA.add_argument('--output_format', type=str, choices=['csv','fasta'], 
                         help = 'output format')
parser_mRNA.add_argument('-p', '--print', action="store_true", 
                         help='output to stdout')
parser_mRNA.add_argument('-u', '--upper', action="store_true", 
                         help='upper CDS and lower utr')
parser_mRNA.add_argument('-s', '--style', choices=['GFF','GTF'], 
                         help = 'GTF database or GFF database')
parser_mRNA.set_defaults(func=mRNA)


# cdna subcommand 
parser_cdna = subparsers.add_parser('cdna', help='extract cdna (or refMrna) in genome or gene')
parser_cdna.add_argument('-d', '--database', type=str, required=True, 
                         help='database generated by subcommand create')
parser_cdna.add_argument('-f', '--genome', type=str, required=True, 
                         help='genome fasta')
parser_cdna.add_argument('-t', '--transcript', type=str, 
                         help='specific transcript; if not given, return whole transcripts')
parser_cdna.add_argument('-o', '--output', type=str, 
                         help = 'output file path')
parser_cdna.add_argument('--output_format', type=str, choices=['csv','fasta'], 
                         help = 'output format')
parser_cdna.add_argument('-p', '--print', action="store_true", 
                         help='output to stdout')
parser_cdna.add_argument('-u', '--upper', action="store_true", 
                         help='upper CDS and lower utr')
parser_cdna.add_argument('-s', '--style', choices=['GFF','GTF'], 
                         help = 'GTF database or GFF database')
parser_cdna.set_defaults(func=cdna)


# IGR subcommand 
parser_IGR = subparsers.add_parser('IGR', help='extract IGR in genome or gene')
parser_IGR.add_argument('-d', '--database', type=str, required=True,
                        help='database generated by subcommand create')
parser_IGR.add_argument('-f', '--genome', type=str, required=True, 
                        help='genome fasta')
parser_IGR.add_argument('-l', '--IGR_length', type=int, default=100,
                        help='IGR length threshold')
parser_IGR.add_argument('-o', '--output', type=str, 
                        help = 'output fasta file path')
parser_IGR.add_argument('-p', '--print', action="store_true", 
                        help='output to stdout')
parser_IGR.add_argument('-s', '--style', choices=['GFF','GTF'], 
                        help = 'GTF database only contain \
                       protein genes, while GFF database contain protein genes and nocoding genes')
parser_IGR.set_defaults(func=IGR)

# UTR subcommand
parser_utr = subparsers.add_parser('UTR', help='extract untranslated region sequence in genome or gene')
parser_utr.add_argument('-d', '--database', type=str, required=True, 
                        help='database generated by subcommand create')
parser_utr.add_argument('-f', '--genome', type=str, required=True,
                        help='genome fasta file')
parser_utr.add_argument('-t', '--transcript', type=str, 
                        help='specific transcript id; if not given, \
                        whole transcript will return')
parser_utr.add_argument('-o', '--output', type=str, 
                        help='output file path')
parser_utr.add_argument('-p', '--print', action="store_true", 
                        help='output to stdout')
parser_utr.add_argument('-s', '--style', choices=['GFF','GTF'], 
                        help = 'GTF database or GFF database')
parser_utr.set_defaults(func=UTR)

# uORF subcommand
parser_uORF = subparsers.add_parser('uORF', help='extract upper stream open reading sequence in genome or gene')
parser_uORF.add_argument('-d', '--database', type=str, required=True, 
                         help='database generated by subcommand create')
parser_uORF.add_argument('-f', '--genome', type=str, required=True, 
                         help='genome fasta')
parser_uORF.add_argument('-t', '--transcript', type=str, 
                         help='specific transcript id; if not given, \
                               whole transcript will return')
parser_uORF.add_argument('-o', '--output', type=str, 
                         help='output file path')
parser_uORF.add_argument('-s', '--style', choices=['GFF','GTF'], 
                         help = 'GTF database or GFF database')
parser_uORF.set_defaults(func=uORF)

# CDS subcommand
parser_cds = subparsers.add_parser('CDS', help='extract coding sequence in genome or gene')
parser_cds.add_argument('-d', '--database', type=str, required=True, 
                        help='database generated by subcommand create')
parser_cds.add_argument('-f', '--genome', type=str, required=True,
                        help='genome fasta')
parser_cds.add_argument('-t', '--transcript', type=str, 
                        help='specific transcript id; if not given, \
                        whole transcript will return')
parser_cds.add_argument('-o', '--output', type=str, 
                        help='output file path')
parser_cds.add_argument('-p', '--print', action="store_true", 
                        help='output to stdout')
parser_cds.add_argument('-s', '--style', choices=['GFF','GTF'], 
                        help = 'GTF database or GFF database')
parser_cds.set_defaults(func=CDS)

# dORF subcommand 
parser_dORF = subparsers.add_parser('dORF', help='extract down stream open reading frame sequence in a genome or gene')
parser_dORF.add_argument('-d', '--database', type=str, required=True,
                         help='database generated by subcommand create')
parser_dORF.add_argument('-f', '--genome', type=str, required=True, 
                         help='genome fasta')
parser_dORF.add_argument('-t', '--transcript', type=str, 
                         help='specific transcript id; if not given, \
                               whole transcript will return')
parser_dORF.add_argument('-o', '--output', type=str, 
                         help='output file path')
parser_dORF.add_argument('-p', '--print', action="store_true", 
                         help='output to stdout')
parser_dORF.add_argument('-s', '--style', choices=['GFF','GTF'], 
                         help = 'GTF database or GFF database')
parser_dORF.set_defaults(func=dORF)


# exon 
parser_exon = subparsers.add_parser('exon', help='extract exon sequence for a given transcript')
parser_exon.add_argument('-d', '--database', type=str, required=True, 
                         help='database generated by subcommand create')
parser_exon.add_argument('-f', '--genome', type=str, required=True, 
                         help='genome fasta')
parser_exon.add_argument('-t', '--transcript', type=str, required=True, 
                         help='specific transcript id; needed')
parser_exon.add_argument('-o', '--output', type=str, 
                         help='output file path')
parser_exon.add_argument('-p', '--print', action="store_true", 
                         help='output to stdout')
parser_exon.add_argument('-s', '--style', choices=['GFF','GTF'], 
                         help = 'GTF database or GFF database')
parser_exon.set_defaults(func=exon)

# intron 
parser_intron = subparsers.add_parser('intron', help='extract exon sequence for a given transcript')
parser_intron.add_argument('-d', '--database', type=str, required=True, 
                           help='database generated by subcommand create')
parser_intron.add_argument('-f', '--genome', type=str, required=True, 
                           help='genome fasta')
parser_intron.add_argument('-t', '--transcript', type=str, 
                           help='specific transcript id; needed')
parser_intron.add_argument('-o', '--output', type=str, 
                           help='output file path')
parser_intron.add_argument('-p', '--print', action="store_true", 
                           help='output to stdout')
parser_intron.add_argument('-s', '--style', choices=['GFF','GTF'], 
                           help = 'GTF database or GFF database')
parser_intron.set_defaults(func=intron)


args = parser.parse_args()
print('[%s runing ...]'%(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime())))
args.func(args)
print('[%s finished ...]'%(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime())))

