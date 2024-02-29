# -*- coding: utf-8 -*-
import os
import gc
import sys
import time
import argparse
import gffutils
import pandas as pd
import multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq
import _pickle as cPickle
from Bio.SeqRecord import SeqRecord
from featurExtract.utils.util import record
from featurExtract.utils.util import gff_feature_dict
from featurExtract.utils.util import gtf_feature_dict
from featurExtract.commands.extract_utr import get_utr
from featurExtract.commands.extract_uorf import get_uorf
from featurExtract.commands.extract_dorf import get_dorf
from featurExtract.commands.extract_cds import get_cds
from featurExtract.commands.extract_promoter import get_promoter
from featurExtract.commands.extract_terminator import get_terminator
from featurExtract.commands.extract_exon import get_exon
from featurExtract.commands.extract_intron import get_intron
from featurExtract.commands.extract_gene import get_gene
from featurExtract.commands.extract_igr import get_igr
from featurExtract.commands.extract_transcript import get_transcript
from featurExtract.commands.extract_mrna import get_mrna
from featurExtract.commands.feature_stat import get_stat

gc.disable()

def create2(args):
    '''
    parameters:
     args: arguments from argparse
     This is a time- and memory-intense procedure, 
     but it needs to be done only once for a given genome.
    return:
     a db object
    '''
    fn = args.genomefeature
    database_id = args.output_prefix +'.'+ args.file_type
    db = gffutils.create_db(fn, dbfn=database_id, force=True, keep_order=True,
        disable_infer_genes=True, disable_infer_transcripts=True,
        merge_strategy='merge', sort_attribute_values=True)
    return db

def create(args):
    '''
    parameters:
     args: arguments from argparse
     This is a time- and memory-intense procedure, 
     but it needs to be done only once for a given genome.
    return:
     a dict object
    '''
    if args.style == 'gff':
        db, t2g = gff_feature_dict(args.genomefeature, args.style)
    elif args.style == 'gtf':
        db, t2g = gtf_feature_dict(args.genomefeature, args.style)
    else:
        sys.exit(1)
        print('The file format should be gff ot gtf.') 
    database_id = os.path.join(args.output, args.prefix +'.'+ args.style)
    with open(database_id, 'wb') as f:
        # fastest way 
        cPickle.dump((db, t2g), f, protocol=-1)

def stat(args):
   '''
   parameters:
     args: arguments from argparse
    return:
     a genome statistics to stdout
   '''
   get_stat(args)
   
def utr(args):
    '''
    parameters:
     args: arguments from argparse 
    '''
    get_utr(args)


def uorf(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_uorf(args)
    
 
def cds(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_cds(args)


def dorf(args):
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

def terminator(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_terminator(args)

def gene(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_gene(args)

def mrna(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_mrna(args)

def transcript(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_transcript(args)

def igr(args):
    '''
    parameters:
     args: arguments from argparse
    '''
    get_igr(args)

def motif(args):
    '''
    parameters:
     args: arguments from argparse
    '''

parser = argparse.ArgumentParser()
subparsers = parser.add_subparsers(help='sub-command help')
# create subcommand 
parser_create = subparsers.add_parser('create', help='create annotation database')
parser_create.add_argument('-g', '--genomefeature', type=str, required=True, 
                           help='genome annotation file, gff or gtf')
parser_create.add_argument('-o', '--output', type=str, required=True, 
                           help='database output dir path')
parser_create.add_argument('-p', '--prefix', type=str, required=True,
                           help='database prefix')
parser_create.add_argument('-s', '--style', choices=['gff','gtf'],
                           help='genome annotation file format')
parser_create.set_defaults(func=create)

# stat subcommand
parser_stat = subparsers.add_parser('stat', help='database stat')
parser_stat.add_argument('-d', '--database', type=str, required=True,
                           help='database created from creat command')
parser_stat.add_argument('-g', '--genome', type=str, required=True,
                             help='genome fasta path')
parser_stat.add_argument('-o', '--output', type=str, required=True,
                           help='stat output')
parser_stat.add_argument('-s', '--style', choices=['gff','gtf'],
                           help='genome annotation file format')
parser_stat.set_defaults(func=stat)


# promoter subcommand
parser_promoter = subparsers.add_parser('promoter', help='extract promoter sequence')
parser_promoter.add_argument('-d', '--database', type=str, required=True, 
                             help='database generated by subcommand create')
parser_promoter.add_argument('-f', '--output_format', type=str, 
                             choices=['csv','fasta'], default='csv',
                             help = 'output format')
parser_promoter.add_argument('-g', '--genome', type=str, required=True,
                             help='genome fasta path')
parser_promoter.add_argument('-i', '--gene', type=str, 
                             help='specific gene (optional); if not given, return whole genes')
parser_promoter.add_argument('-l', '--promoter_length', type=int, default=100,
                             help='promoter length before TSS (default: 100)')
parser_promoter.add_argument('-o', '--output', type=str, 
                             help = 'output file path')
parser_promoter.add_argument('-p', '--process', type=int, default=4,
                             help='number of promoter extract process, (default: 4)')
parser_promoter.add_argument('-u', '--utr5_upper_length', type=int, default=10,
                             help='5\' utr length after TSS (default: 10)')
parser_promoter.add_argument('-v', '--print', action="store_true",
                             help = 'output to stdout')
parser_promoter.set_defaults(func=promoter)


# terminator 
parser_terminator = subparsers.add_parser('terminator', help='extract terminator sequence')
parser_terminator.add_argument('-d', '--database', type=str, required=True,
                             help='database generated by subcommand create')
parser_terminator.add_argument('-f', '--output_format', type=str, choices=['csv','fasta'],
                             help = 'output format')
parser_terminator.add_argument('-g', '--genome', type=str, required=True,
                             help='genome fasta path')
parser_terminator.add_argument('-i', '--gene', type=str,
                             help='specific gene (optional); if not given, return whole genes')
parser_terminator.add_argument('-l', '--terminator_length', type=int, default=100,
                             help='terminator length (default: 100)')
parser_terminator.add_argument('-o', '--output', type=str,
                             help = 'output file path')
parser_terminator.add_argument('-u', '--utr3_lower_length', type=int, default=10,
                             help='3\' length (default: 10)')
parser_terminator.add_argument('-v', '--print', action="store_true", 
                             help = 'output to stdout')
parser_terminator.set_defaults(func=terminator)

# gene subcommand 
parser_gene = subparsers.add_parser('gene', help='extract gene sequence')
parser_gene.add_argument('-d', '--database', type=str, required=True,
                         help='database generated by subcommand create')
parser_gene.add_argument('-f', '--output_format', type=str, choices=['csv','fasta', 'gff', 'gtf'],
                             help = 'output format')
parser_gene.add_argument('-g', '--genome', type=str, required=True,
                         help='genome fasta')
parser_gene.add_argument('-i', '--gene', type=str, 
                         help='specific gene (optional); if not given, return whole genes')
parser_gene.add_argument('-o', '--output', type=str, 
                         help = 'output file path')
parser_gene.add_argument('-p', '--print', action="store_true", 
                         help='output to stdout')
parser_gene.add_argument('-s', '--style', choices=['gff','gtf'],
                         help = 'gtf database or gff database')
parser_gene.set_defaults(func=gene)

# mrna subcommand 
parser_mrna = subparsers.add_parser('mrna', help='extract messager RNA')
parser_mrna.add_argument('-d', '--database', type=str, required=True, 
                         help='database generated by subcommand create')
parser_mrna.add_argument('-f', '--output_format', type=str, choices=['csv','fasta'], default='csv',
                         help = 'output format')
parser_mrna.add_argument('-g', '--genome', type=str, required=True, 
                         help='genome fasta')
parser_mrna.add_argument('-i', '--transcript', type=str, 
                         help='specific transcript (optional); if not given, return whole transcripts')
parser_mrna.add_argument('-o', '--output', type=str, 
                         help = 'output file path')
parser_mrna.add_argument('-p', '--print', action="store_true", 
                         help='output to stdout')
parser_mrna.add_argument('-s', '--style', choices=['gff','gtf'],
                         help = 'gtf database or gff database')
parser_mrna.add_argument('-u', '--upper', action="store_true", 
                         help='upper cds and lower utr')
parser_mrna.set_defaults(func=mrna)


# transcript subcommand 
parser_transcript = subparsers.add_parser('transcript', help='extract transcript')
parser_transcript.add_argument('-d', '--database', type=str, required=True, 
                         help='database generated by subcommand create')
parser_transcript.add_argument('-f', '--output_format', type=str, choices=['csv','fasta'], default='csv',
                         help = 'output format')
parser_transcript.add_argument('-g', '--genome', type=str, required=True, 
                         help='genome fasta')
parser_transcript.add_argument('-i', '--transcript', type=str, 
                         help='specific transcript (optional); if not given, return whole transcripts')
parser_transcript.add_argument('-o', '--output', type=str, 
                         help = 'output file path')
parser_transcript.add_argument('-p', '--process', type=int, default=4,
                         help='number of cDNA extract process, (default: 4)')
parser_transcript.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
                         help='The type of RNA for extract transcript, (default: mrna)')
parser_transcript.add_argument('-s', '--style', choices=['gff','gtf'],
                         help = 'gtf database or gff database')
parser_transcript.add_argument('-u', '--upper', action="store_true", 
                         help='upper cds and lower utr')
parser_transcript.add_argument('-v', '--print', action="store_true",
                         help='output to stdout')
parser_transcript.set_defaults(func=transcript)


# igr subcommand 
parser_igr = subparsers.add_parser('igr', help='extract intergenic region')
parser_igr.add_argument('-d', '--database', type=str, required=True,
                        help='database generated by subcommand create')
parser_igr.add_argument('-f', '--output_format', type=str, choices=['csv','fasta', 'gff'], default='csv',
                         help = 'output format')
parser_igr.add_argument('-g', '--genome', type=str, required=True, 
                        help='genome fasta')
parser_igr.add_argument('-l', '--igr_length', type=int, default=100,
                        help='igr length threshold')
parser_igr.add_argument('-o', '--output', type=str, 
                        help = 'output fasta file path')
parser_igr.add_argument('-p', '--process', type=int, default=4,
                         help='number of igr extract process, (default: 4)')
parser_igr.add_argument('-s', '--style', choices=['gff','gtf'], 
                        help = 'gtf database only contain \
                       protein genes, while gff database contain protein genes and nocoding genes')
parser_igr.add_argument('-v', '--print', action="store_true",
                         help='output to stdout')
parser_igr.set_defaults(func=igr)

# utr subcommand
parser_utr = subparsers.add_parser('utr', help='extract untranslated region sequence')
parser_utr.add_argument('-d', '--database', type=str, required=True, 
                        help='database generated by subcommand create')
parser_utr.add_argument('-f', '--output_format', choices=['csv','fasta','gff'], default='csv',
                        help='output format (default: csv)')
parser_utr.add_argument('-g', '--genome', type=str, required=True,
                        help='genome fasta file')
parser_utr.add_argument('-i', '--transcript', type=str, 
                        help='specific transcript (optional); if not given, \
                        return whole transcripts')
parser_utr.add_argument('-o', '--output', type=str, 
                        help='output file path')
parser_utr.add_argument('-p', '--process', type=int, default=4,
                         help='number of utr extract process, (default: 4)')
parser_utr.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
                         help='The type of RNA for extract utr, (default: mrna)')
parser_utr.add_argument('-s', '--style', choices=['gff','gtf'], 
                        help = 'gtf database or gff database')
parser_utr.add_argument('-v', '--print', action="store_true", 
                        help='output to stdout. -v and -o option are mutually exclusive')
parser_utr.set_defaults(func=utr)

# uorf subcommand
parser_uorf = subparsers.add_parser('uorf', help='extract upperstream open reading sequence')
parser_uorf.add_argument('-d', '--database', type=str, required=True, 
                         help='database generated by subcommand create')
parser_uorf.add_argument('-f', '--output_format', choices=['csv','fasta','gff'], default='csv',
                        help='output format (default: csv)')
parser_uorf.add_argument('-g', '--genome', type=str, required=True,
                         help='genome fasta')
parser_uorf.add_argument('-i', '--transcript', type=str, 
                         help='specific transcript (optional); if not given, \
                               return whole transcripts')
parser_uorf.add_argument('-l', '--length', type=int, default=6,
                         help='uorf length, (default: 6)')
parser_uorf.add_argument('-m', '--schematic_without_intron', action='store_true',
                         help='schematic figure file for uorf, cds and transcript without intron')
parser_uorf.add_argument('-n', '--schematic_with_intron', action='store_true',
                         help='schematic figure file for uorf, cds and transcript with intron')
parser_uorf.add_argument('-o', '--output', type=str, 
                         help='output file path')
parser_uorf.add_argument('-p', '--process', type=int, default=4,
                         help='number of uorf extract process, (default: 4)')
parser_uorf.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
                         help='The type of RNA for uorf extraction (default: mrna)')
parser_uorf.add_argument('-s', '--style', choices=['gff','gtf'], 
                         help = 'gtf database or gff database')
parser_uorf.add_argument('-v', '--print', action="store_true",
                        help='output to stdout. -v and -o option are mutually exclusive')
parser_uorf.set_defaults(func=uorf)

# cds subcommand
parser_cds = subparsers.add_parser('cds', help='extract coding sequence')
parser_cds.add_argument('-d', '--database', type=str, required=True, 
                        help='database generated by subcommand create')
parser_cds.add_argument('-f', '--output_format', type=str, choices=['csv','fasta','gff'], default='csv',
                         help = 'output format')
parser_cds.add_argument('-g', '--genome', type=str, required=True,
                        help='genome fasta')
parser_cds.add_argument('-i', '--transcript', type=str, 
                        help='specific transcript (optional); if not given, \
                        return whole transcripts')
parser_cds.add_argument('-o', '--output', type=str, 
                        help='output file path')
parser_cds.add_argument('-p', '--process', type=int, default=4,
                        help='number of cds extract process, (default: 4)')
parser_cds.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
                         help='The type of RNA for extract cds, (default: mrna)')
parser_cds.add_argument('-s', '--style', choices=['gff','gtf'], 
                        help = 'gtf database or gff database')
parser_cds.add_argument('-v', '--print', action="store_true",
                        help='output to stdout. -v and -o option are mutually exclusive')
parser_cds.set_defaults(func=cds)

# dorf subcommand 
parser_dorf = subparsers.add_parser('dorf', help='extract downstream open reading frame sequence')
parser_dorf.add_argument('-d', '--database', type=str, required=True,
                         help='database generated by subcommand create')
parser_dorf.add_argument('-f', '--output_format', choices=['csv','fasta','gff'], default='csv',
                        help='output format')
parser_dorf.add_argument('-g', '--genome', type=str, required=True, 
                         help='genome fasta')
parser_dorf.add_argument('-i', '--transcript', type=str, 
                         help='specific transcript (optional); if not given, \
                               return whole transcripts')
parser_dorf.add_argument('-l', '--length', type=int, default=6,
                         help='dorf length, (default: 6)')
parser_dorf.add_argument('-m', '--schematic_without_intron', action='store_true',
                         help='schematic figure file for dorf, cds and transcript without intron')
parser_dorf.add_argument('-n', '--schematic_with_intron', action='store_true',
                         help='schematic figure file for dorf, cds and transcript with intron')
parser_dorf.add_argument('-o', '--output', type=str, 
                         help='output file path')
parser_dorf.add_argument('-p', '--process', type=int, default=4,
                         help='number of dorf extract process, (default: 4)')
parser_dorf.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
                         help='The type of RNA for dorf extraction (default: mrna)')
parser_dorf.add_argument('-s', '--style', choices=['gff','gtf'], 
                         help = 'gtf database or gff database')
parser_dorf.add_argument('-v', '--print', action="store_true",
                        help='output to stdout. -v and -o option are mutually exclusive')
parser_dorf.set_defaults(func=dorf)


# exon 
parser_exon = subparsers.add_parser('exon', help='extract exon sequence')
parser_exon.add_argument('-d', '--database', type=str, required=True, 
                         help='database generated by subcommand create')
parser_exon.add_argument('-f', '--output_format', choices=['csv','fasta','gff'], default='csv',
                        help='output format')
parser_exon.add_argument('-g', '--genome', type=str, required=True, 
                         help='genome fasta')
parser_exon.add_argument('-i', '--transcript', type=str, 
                         help='specific transcript (optional);  if not given, \
                               return whole transcripts')
parser_exon.add_argument('-o', '--output', type=str, 
                         help='output file path')
parser_exon.add_argument('-p', '--process', type=int, default=4,
                         help='number of exon extract process, (default: 4)')
parser_exon.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
                         help='The type of RNA for exon extraction (default: mrna)')
parser_exon.add_argument('-s', '--style', choices=['gff','gtf'], 
                         help = 'gtf database or gff database')
parser_exon.add_argument('-v', '--print', action="store_true", 
                         help='output to stdout')
parser_exon.set_defaults(func=exon)

# intron 
parser_intron = subparsers.add_parser('intron', help='extract intron sequence')
parser_intron.add_argument('-d', '--database', type=str, required=True, 
                           help='database generated by subcommand create')
parser_intron.add_argument('-f', '--output_format', choices=['csv','fasta','gff'], default='csv',
                        help='output format')
parser_intron.add_argument('-g', '--genome', type=str, required=True, 
                           help='genome fasta')
parser_intron.add_argument('-i', '--transcript', type=str, 
                           help='specific transcript (optional);  if not given, \
                               return whole transcripts')
parser_intron.add_argument('-o', '--output', type=str, 
                           help='output file path')
parser_intron.add_argument('-p', '--process', type=int, default=4, 
                           help='number of exon extract process, (default: 4)')
parser_intron.add_argument('-r', '--rna_feature', choices=['mrna', 'all'], default='mrna',
                         help='The type of RNA for intron extraction (default: mrna)')
parser_intron.add_argument('-s', '--style', choices=['gff','gtf'], 
                           help = 'gtf database or gff database')
parser_intron.add_argument('-v', '--print', action="store_true",
                         help='output to stdout')
parser_intron.set_defaults(func=intron)

# motif 
parser_motif = subparsers.add_parser('motif', help='extract sequence motf from input fasta')
parser_motif.add_argument('-i', '--input', type=str,
                           help='input sequence in fasta format')
parser_motif.add_argument('-o', '--output', type=str,
                           help='output file path')
parser_motif.add_argument('-p', '--process', type=int, default=4,
                           help='number of exon extract process, (default: 4)')
parser_motif.set_defaults(func=motif)

args = parser.parse_args()
#print('[%s runing ...]'%(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime())))
args.func(args)
#print('[%s finished ...]'%(time.strftime("%a %b %d %H:%M:%S %Y", time.localtime())))

gc.enable()
