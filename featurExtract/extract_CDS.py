# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from featurExtract.database import create_db


'''
https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
'''
def stop_codon(db, transcript, genome):
    s = ''
    for codon in db.children(transcript, featuretype='stop_codon', order_by='start'):
        s = codon.sequence(genome, use_strand=False) # 不反向互补,序列全部连接后再互补
    return s

def add_stop_codon(seq, strand, stop_codon_seq):
    if strand == '+':
        seq += stop_codon_seq
    elif strand == '-':
        seq = stop_codon_seq + seq # 负链 stop codon 位置靠前
    else:
        seq += stop_codon_seq
    return seq 

def mRNA_type(file_format):
    if file_format == 'GFF':
        return 'mRNA'

    elif file_format == 'GTF':
        return 'transcript'
    else:
        sys.stderr.write("parameter -s should be assign \n")
        sys.exit(1)

def get_cds(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    cds_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','CDS'])
    index = 0
    # assert GTF or GFF
    mRNA_str = mRNA_type(args.style)
    if not args.transcript:
        for t in db.features_of_type(mRNA_str, order_by='start'):
            seq = ''
            for c in db.children(t, featuretype='CDS', order_by='start'):
                # 不反向互补，对于负链要得到全部的cds后再一次性反向>互补
                s = c.sequence(args.genome, use_strand=False)
                seq += s
            if args.style == 'GTF':
                stop_codon_seq = stop_codon(db, t, args.genome)
                seq = add_stop_codon(seq, t.strand, stop_codon_seq)
            seq = Seq(seq)
            if t.strand == '-':
                seq = seq.reverse_complement()
            cds_seq.loc[index] = [t.id,t.chrom,t.start,t.end,t.strand,seq]
            index += 1
            if index == 100:
                break 
        cds_seq.to_csv(args.output, sep=',', index=False)
    else:
        for t in db.features_of_type(mRNA_str, order_by='start'):
            if args.transcript in t.id:
                seq = ''
                for c in db.children(t, featuretype='CDS', order_by='start'):
                    # 不反向互补，对于负链要得到全部的cds后再一次性>反向互补
                    s = c.sequence(args.genome, use_strand=False)
                    seq += s
                if args.style == 'GTF':
                    stop_codon_seq = stop_codon(db, t, args.genome)
                    seq = add_stop_codon(seq, t.strand, stop_codon_seq)
                seq = Seq(seq)
                if t.strand == '-':
                    seq= seq.reverse_complement()
                cds_seq.loc[index] = [t.id,t.chrom,t.start,t.end,t.strand,seq]
                index += 1
                desc='strand %s start %d end %d length=%d'%(t.strand,t.start,t.end,len(seq))
                cdsRecord = SeqRecord(seq, id=t.id, description=desc)
                if args.print:
                    SeqIO.write([cdsRecord], sys.stdout, "fasta") 
                else:
                    SeqIO.write([cdsRecord], args.output, "fasta") 
                break 



def get_cds_gb(args):
    cds = []
    proteins = []
    for record in create_db(args.genbank):
        for feature in record.features:
            if feature.type == 'CDS': # CDS promoter UTR 
                cds_seq = ''
                for part in feature.location.parts:
                    cds_seq += part.extract(record.seq)
                # part.strand 会将FeatureLocation -1的反向互补
                if 'transl_table' in feature.qualifiers:
                    table_translate = feature.qualifiers['transl_table'][0] 
                else:
                    table_translate = 1
                if str(cds_seq)[:3] in ['ATA','GTG','TTG','ATT','ACG']:
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
                if 'gene' in feature.qualifiers:
                    gene_id = feature.qualifiers['gene'][0] 
                else:
                    gene_id = feature.qualifiers['locus_tag'][0]
                cds_seq_record = SeqRecord(cds_seq,id='gene:%s protein:%s'%(gene_id, protein_id), 
                                 description='strand %s length %d'%(feature.strand, len(cds_seq)))
                pep_record = SeqRecord(pep,id='gene:%s protein:%s'%(gene_id, protein_id),
                             description='strand %s length %d'%(feature.strand, len(pep)))
                cds.append(cds_seq_record)
                proteins.append(pep_record)
                #break 
    if args.print and args.format == 'dna':
        SeqIO.write(cds, sys.stdout, 'fasta')
    elif args.print and args.format == 'protein':
        SeqIO.write(proteins, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(proteins, args.output, 'fasta')
