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


def get_cdna(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    cdna_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','CDS'])
    index = 0
    # assert GTF or GFF
    mRNA_str = mRNA_type(args.style)
    if not args.transcript:
        for t in db.features_of_type(mRNA_str, order_by='start'):
            start_codon, stop_codon = anchor_CDS(db, args.genome, t, args.style)
            first_exon_position = ''
            seq = ''
            for e in db.children(t, featuretype='exon', order_by='start'):
                # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                s = e.sequence(args.genome, use_strand=False)
                seq += s
                if t.strand == '-':
                    if e.start >= start_codon_s :
                        atg2firstexon += len(s)
                    elif e.start < start_codon_s < e.end:
                        truncated = e.end - start_codon_s + 1 # 反向减法
                        atg2firstexon += truncated
                else:
                    # contain + .
                    if start_codon_s >= e.end:
                        atg2firstexon += len(s)
                    elif e.end > start_codon_s > e.start:
                        truncated = start_codon_s - e.start + 1
                        atg2firstexon += truncated
            seq = Seq(seq)
            if t.strand == '-':
                seq = seq.reverse_complement()
            cdna_seq.loc[index] = [t.id.replace('transcript:',''),t.chrom,t.start,t.end,t.strand,seq]
            index += 1
        cdna_seq.to_csv(args.output, sep=',', index=False)
    else:
        for t in db.features_of_type(mRNA_str, order_by='start'):
            if args.transcript in t.id:
                start_codon_s, stop_codon_e, cds_len = anchor_CDS(db, args.genome, t, args.style)
                atg2firstexon = 0
                seq = ''
                for e in db.children(t, featuretype='exon', order_by='start'):
                    # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                    s = e.sequence(args.genome, use_strand=False)
                    seq += s
                    
                    if t.strand == '-':
                        if e.start >= start_codon_s :
                            atg2firstexon += len(s) 
                        elif e.start < start_codon_s < e.end:
                            truncated = e.end - start_codon_s + 1 # 反向减法
                            atg2firstexon += truncated
                    else:
                        # contain + .
                        if start_codon_s >= e.end:
                            atg2firstexon += len(s)
                        elif e.end > start_codon_s > e.start:
                            truncated = start_codon_s - e.start + 1
                            atg2firstexon += truncated
                seq = Seq(seq)
                if t.strand == '-':
                    seq= seq.reverse_complement()
                index += 1
                if t.strand == "-":
                    desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(t.strand,t.start,t.end,len(seq),
                                                                      atg2firstexon ,
                                                                      atg2firstexon + cds_len - 1)
                else:
                    desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(t.strand,t.start,t.end,len(seq),
                                                                      atg2firstexon ,
                                                                      atg2firstexon + cds_len - 1)
                if args.upper:
                    seq = seq_upper_lower(seq,atg2firstexon,atg2firstexon + cds_len - 1)
                cdnaRecord = SeqRecord(seq, id=t.id.replace('transcript:',''), description=desc)
                if args.print:
                    SeqIO.write([cdnaRecord], sys.stdout, "fasta") 
                else:
                    SeqIO.write([cdnaRecord], args.output, "fasta") 
                break 


# GenBank 基因组文件，是每条染色体一个LOCUS

def get_cdna_gb(args):
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
