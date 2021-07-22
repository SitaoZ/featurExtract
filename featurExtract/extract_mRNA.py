# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from featurExtract.database import create_db
from featurExtract.util import utr3_type, utr5_type, mRNA_type

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


def seq_upper_lower(seq,start,end):
    '''
    parameter:
     seq: sequence 
     start: start codon position 
     end: stop codon position 
    '''
    utr5 = seq[:start-1].lower()
    coding = seq[start-1:end].upper()
    utr3 = seq[end:].lower()
    return utr5+coding+utr3
    
def anchor_CDS(db, genome, transcript, style):
    '''
    parameters: 
     db : database create by gffutils
     transcript: gffutils transcript feature type
     style: database file type
    return:
     start_codon position 
     stop_codon positiopn 
     cds length 
    '''
    if style == 'GTF':
        cds = ''
        for c in db.children(transcript, featuretype='CDS', order_by='start'):
            cds += c.sequence(genome, use_strand=False)
        cds = Seq(cds)
        if transcript.strand == '-':
            cds = cds.reverse_complement()
        cds_len = len(cds) + 3 # add stop codon length
        if transcript.strand == '-':
            for i in db.children(transcript, featuretype='start_codon', order_by='start'):
                start_codon_s = i.end
                start_codon_e = i.start
            for i in db.children(transcript, featuretype='stop_codon', order_by='start'):
                stop_codon_s = i.end
                stop_codon_e = i.start    
        else:
            # contain + .
            for i in db.children(transcript, featuretype='start_codon', order_by='start'):
                start_codon_s = i.start
                start_codon_e = i.end
            for i in db.children(transcript, featuretype='stop_codon', order_by='start'):
                stop_codon_s = i.start
                stop_codon_e = i.end
        
        return start_codon_s, stop_codon_e, cds_len
    elif style == 'GFF':
        cds = ''
        for c in db.children(transcript, featuretype='CDS', order_by='start'):
            cds += c.sequence(genome, use_strand=False)
        cds = Seq(cds)
        if transcript.strand == '-':
            cds = cds.reverse_complement()
        cds_len = len(cds) # cds full length
        if transcript.strand == '-':
            for i in db.children(transcript, featuretype='five_prime_UTR', order_by='start'):
                # occasionally a transcript has more then one three_prime_UTR
                # the first five_prime_UTR position should be saved in minus strand
                start_codon_s = i.start - 1
                start_codon_e = i.start - 3
                break
            for i in db.children(transcript, featuretype='three_prime_UTR', order_by='start'):
                # occasionally a transcript has more then one three_prime_UTR
                # the last save (position large)
                stop_codon_s = i.end - 1
                stop_codon_e = i.end - 3
        else:
            # contain + .
            for i in db.children(transcript, featuretype='five_prime_UTR', order_by='start'):
                # occasionally a transcript has more then one three_prime_UTR
                # the last five_prime_UTR position should be saved in plus strand
                start_codon_s = i.end + 1
                start_codon_e = i.end - 2
            for i in db.children(transcript, featuretype='three_prime_UTR', order_by='start'):
                # occasionally a transcript has more then one three_prime_UTR
                # the first three_prime_UTR position should be saved in plus strand
                stop_codon_s = i.start - 2
                stop_codon_e = i.start + 1
        return start_codon_s, stop_codon_e, cds_len

def plus_strand():
    pass

def minus_strand():
    pass


def get_mRNA(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    mrna_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','CDS'])
    index = 0
    # assert GTF or GFF
    mRNA_str = mRNA_type(args.style)
    utr3_t = utr3_type(args.style)
    utr5_t = utr5_type(args.style)
    
    if not args.transcript:
        for t in db.features_of_type(mRNA_str, order_by='start'):
            #start_codon, stop_codon = anchor_CDS(db, args.genome, t, args.style)
            utr5, cds , utr3 = '', '', ''
            for u in db.children(t, featuretype=utr5_t, order_by='start'):
                utr5 += u.sequence(args.genome, use_strand=False)
            for u in db.children(t, featuretype=utr3_t, order_by='start'):
                utr3 += u.sequence(args.genome, use_strand=False)
            for c in db.children(t, featuretype='CDS', order_by='start'):
                # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                cds = c.sequence(args.genome, use_strand=False)
            seq = utr5 + cds + utr3
            seq = Seq(seq)
            if t.strand == '-':
                seq = seq.reverse_complement()
            mrna_seq.loc[index] = [t.id.replace('transcript:',''),t.chrom,t.start,t.end,t.strand,seq]
            index += 1
        mrna_seq.to_csv(args.output, sep=',', index=False)
    else:
        for t in db.features_of_type(mRNA_str, order_by='start'):
            if args.transcript in t.id:
                #start_codon_s, stop_codon_e, cds_len = anchor_CDS(db, args.genome, t, args.style)
                utr5, cds , utr3 = '', '', ''
                for u in db.children(t, featuretype=utr5_t, order_by='start'):
                    utr5 += u.sequence(args.genome, use_strand=False)
                for u in db.children(t, featuretype=utr3_t, order_by='start'):
                    utr3 += u.sequence(args.genome, use_strand=False)
                for c in db.children(t, featuretype='CDS', order_by='start'):
                    # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                    cds = c.sequence(args.genome, use_strand=False)
                seq = utr5 + cds + utr3
                seq = Seq(seq)
                if t.strand == '-':
                    seq= seq.reverse_complement()
                index += 1
                if t.strand == "-":
                    desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(t.strand,t.start,t.end,len(seq),
                                                                      len(utr5) ,
                                                                      len(utr5) + cds_len - 1)
                else:
                    desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(t.strand,t.start,t.end,len(seq),
                                                                      len(utr5) ,
                                                                      len(utr5) + cds_len - 1)
                if args.upper:
                    seq = utr5.lower() + cds.upper() + utr3.lower()
                mrnaRecord = SeqRecord(seq, id=t.id.replace('transcript:',''), description=desc)
                if args.print:
                    SeqIO.write([mrnaRecord], sys.stdout, "fasta") 
                else:
                    SeqIO.write([mrnaRecord], args.output, "fasta") 
                break 


# GenBank 基因组文件，是每条染色体一个LOCUS

def get_mRNA_gb(args):
    '''
    '''
    cdna = []
    for record in create_db(args.genbank):
        for feature in record.features:
            if feature.type == 'exon': # CDS promoter UTR 
                cdna_seq = ''
                #for part in feature.location.parts:
                #    cds_seq += part.extract(record.seq)
                cdna_seq = feature.location.extract(record).seq
                # part.strand 会将FeatureLocation -1的反向互补
                # geneid or locus_tag 
                if 'locus_tag' in feature.qualifiers:
                    gene_id = feature.qualifiers['locus_tag'][0]
                elif 'gene' in feature.qualifiers:
                    gene_id = feature.qualifiers['gene'][0]
                else:
                    gene_id = "Null"
                # protein id  
                if 'protein_id' in feature.qualifiers:
                    protein_id = feature.qualifiers['protein_id'][0]
                else:
                    protein_id = 'Null'

                cdna_seq_record = SeqRecord(cds_seq,id='gene:%s protein:%s'%(gene_id, protein_id),
                                 description='strand %s length %d'%(feature.strand, len(cds_seq)))
                cdna.append(cds_seq_record)
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
