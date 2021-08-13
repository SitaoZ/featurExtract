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
from featurExtract.util import stop_codon_seq, add_stop_codon

'''
https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
'''

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
            utr5, cds , utr3 = '', '', ''
            for u in db.children(t, featuretype=utr5_t, order_by='start'):
                utr5 += u.sequence(args.genome, use_strand=False)
            for u in db.children(t, featuretype=utr3_t, order_by='start'):
                utr3 += u.sequence(args.genome, use_strand=False)
            for c in db.children(t, featuretype='CDS', order_by='start'):
                # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                cds += c.sequence(args.genome, use_strand=False)
            if args.style == 'GTF':
                stop_codon_s = stop_codon_seq(db, t, args.genome)
                cds = add_stop_codon(cds, t.strand, stop_codon_s)
            if args.upper:
                if t.strand == '-':
                    seq = utr3.lower() + cds.upper() + utr5.lower()
                else:
                    seq = utr5.lower() + cds.upper() + utr3.lower()
            else:
                if t.strand == '-':
                    seq = utr3 + cds + utr5
                else:
                    seq = utr5 + cds + utr3
            
            seq = Seq(seq)
            if t.strand == '-':
                seq = seq.reverse_complement()
            mrna_seq.loc[index] = [t.id.replace('transcript:',''),t.chrom,t.start,t.end,t.strand,seq]
            index += 1
            if index == 200:
                break 
        mrna_seq.to_csv(args.output, sep=',', index=False)
    else:
        for t in db.features_of_type(mRNA_str, order_by='start'):
            if args.transcript in t.id:
                utr5, cds , utr3 = '', '', ''
                for u in db.children(t, featuretype=utr5_t, order_by='start'):
                    utr5 += u.sequence(args.genome, use_strand=False)
                for u in db.children(t, featuretype=utr3_t, order_by='start'):
                    utr3 += u.sequence(args.genome, use_strand=False)
                for c in db.children(t, featuretype='CDS', order_by='start'):
                    # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                    cds += c.sequence(args.genome, use_strand=False)
                if args.style == 'GTF':
                    stop_codon_s = stop_codon_seq(db, t, args.genome)
                    cds = add_stop_codon(cds, t.strand, stop_codon_s)

                if args.upper:
                    if t.strand == '-':
                        seq = utr3.lower() + cds.upper() + utr5.lower()
                    else:
                        seq = utr5.lower() + cds.upper() + utr3.lower()
                else:
                    if t.strand == '-':
                        seq = utr3 + cds + utr5
                    else:
                        seq = utr5 + cds + utr3
                seq = Seq(seq)
                if t.strand == '-':
                    seq= seq.reverse_complement()
                index += 1
                if t.strand == "-":
                    desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(t.strand,t.start,t.end,len(seq),
                                                                      len(utr5)+1 ,
                                                                      len(utr5) + len(cds))
                else:
                    desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(t.strand,t.start,t.end,len(seq),
                                                                      len(utr5)+1 ,
                                                                      len(utr5) + len(cds))
                mrnaRecord = SeqRecord(seq, id=t.id.replace('transcript:', ''), description=desc)
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
