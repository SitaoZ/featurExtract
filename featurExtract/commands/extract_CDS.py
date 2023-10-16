# -*- coding: utf-8 -*-
import os
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
from featurExtract.database.database import create_db
from featurExtract.utils.util import add_stop_codon, mRNA_type, parse_output

'''
https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
'''
def stop_codon(db, transcript, genome):
    s = ''
    for codon in db.children(transcript, featuretype='stop_codon', order_by='start'):
        s = codon.sequence(genome, use_strand=False) # 不反向互补,序列全部连接后再互补
    return s


_CSV_HEADER = ['TranscriptID','Chrom','Start','End','Strand','CDS']

def sub_cds(params):
    t, database, genome, style, output_format = params
    db = gffutils.FeatureDB(database, keep_order=True)
    cds_seq_list = deque()
    seq = ''
    cds_start_transcript = 0
    cds_end_transcript = 0
    cds_gff_lines = deque()
    for c in db.children(t, featuretype='CDS', order_by='start'):
        cds_gff_lines.append(c)
        # 不反向互补，对于负链要得到全部的cds后再一次性反向>互补
        s = c.sequence(genome, use_strand=False)
        seq += s
        if not cds_start_transcript :
            cds_start_transcript = c.start - t.start
        cds_end_transcript += len(s)
    cds_end_transcript = cds_end_transcript + cds_start_transcript
    
    if style == 'GTF':
        stop_codon_seq = stop_codon(db, t, genome)
        seq = add_stop_codon(seq, t.strand, stop_codon_seq)
    seq = Seq(seq)
    if t.strand == '-':
        seq = seq.reverse_complement()
    # csv output
    if output_format == 'gff':
        for line in cds_gff_lines:
            cds_seq_list.append(line)
    # defalut fasta
    elif output_format == 'fasta':
        desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(t.strand,
                                                              t.start,
                                                              t.end,
                                                              len(seq),
                                                              cds_start_transcript,
                                                              cds_end_transcript)
        cdsRecord = SeqRecord(seq, id=t.id.replace('transcript:',''), description=desc)
        cds_seq_list.append(cdsRecord)
    else:
        it = [t.id,t.chrom,t.start,t.end,t.strand,seq]
        cds_seq_list.append(dict((_CSV_HEADER[i], it[i]) for i in range(len(_CSV_HEADER))))
    return cds_seq_list

def get_cds(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    # print(args)
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    feature_types = db.featuretypes()
    # assert GTF or GFF
    if args.rna_feature == 'mRNA':
        mRNA_str = mRNA_type(args.style)
    else:
        mRNA_str = tuple(x for x in list(feature_types) if 'RNA' in x or 'transcript' in x)
    # loop all transcript in genome
    if not args.transcript:
        param_list = [(t, args.database, args.genome, args.style, args.output_format) for t in db.features_of_type(mRNA_str, order_by='start')]
        with Pool(processes=args.process) as p:
            cds_seq_list = list(tqdm(p.imap(sub_cds, param_list), total=len(param_list), ncols = 80, desc='CDS Processing:'))
        cds_seq_list = [d for de in cds_seq_list if de != None for d in de]
        parse_output(args, cds_seq_list)
    else:
        # only one transcript
        for t in db.features_of_type(mRNA_str, order_by='start'):
            if args.transcript in t.id:
                param = (t, args.database, args.genome, args.style, args.output_format)
                cds_seq_list = sub_cds(param) 
                parse_output(args, cds_seq_list)
                break 


# GenBank 基因组文件，是每条染色体一个LOCUS

def get_cds_gb(args):
    '''
    '''
    cds = []
    proteins = []
    for record in create_db(args.genbank):
        for feature in record.features:
            if feature.type == 'CDS': # CDS promoter UTR 
                cds_seq = ''
                #for part in feature.location.parts:
                #    cds_seq += part.extract(record.seq)
                cds_seq = feature.location.extract(record).seq
                if len(cds_seq)%3 != 0:
                    continue # reject not triple 
                # part.strand 会将FeatureLocation -1的反向互补
                # geneid or locus_tag 
                if 'locus_tag' in feature.qualifiers:
                    gene_id = feature.qualifiers['locus_tag'][0]
                elif 'gene' in feature.qualifiers:
                    gene_id = feature.qualifiers['gene'][0]
                else:
                    gene_id = "Null"
                if 'transl_table' in feature.qualifiers:
                    table_translate = feature.qualifiers['transl_table'][0] 
                else:
                    table_translate = 1
                if str(cds_seq)[:3] in ['AAT','ATA','GTG','TTG','ATT','ACG','TCA','AGG']:
                    # ATA, GTG and TTG (Yokobori et al. 1999). 
                    # ATT is the start codon for the CytB gene
                    # in Halocynthia roretzi (Gissi and Pesole, 2003).
                    pep = Seq(feature.qualifiers['translation'][0])
                else:
                    pep = cds_seq.translate(table=table_translate, cds=True)
                #pep = cds_seq.translate(to_stop=True)
                # 判断提取的和已知的是否一致
                if 'translation' in feature.qualifiers:
                    assert(str(pep) == feature.qualifiers['translation'][0])
                if 'protein_id' in feature.qualifiers:
                    protein_id = feature.qualifiers['protein_id'][0] 
                else:
                    protein_id = 'Null'
                cds_seq_record = SeqRecord(cds_seq,id='gene:%s protein:%s'%(gene_id, protein_id), 
                                 description='strand %s length %d'%(feature.strand, len(cds_seq)))
                pep_record = SeqRecord(pep,id='gene:%s protein:%s'%(gene_id, protein_id),
                             description='strand %s length %d'%(feature.strand, len(pep)))
                cds.append(cds_seq_record)
                proteins.append(pep_record)
                break 
        break     
    if args.print and args.format == 'dna':
        SeqIO.write(cds, sys.stdout, 'fasta')
    elif args.print and args.format == 'protein':
        SeqIO.write(proteins, sys.stdout, 'fasta')
    elif args.output and args.format == 'dna':
        SeqIO.write(cds, args.output, 'fasta')
    elif args.output and args.format == 'protein':
        SeqIO.write(proteins, args.output, 'fasta')
