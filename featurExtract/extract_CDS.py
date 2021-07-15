# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from featurExtract.database import create_db

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
    if not args.transcript:
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
    else:
        for t in db.features_of_type('mRNA', order_by='start'):
            if args.transcript in t.id:
                seq = ''
                for c in db.children(t, featuretype='CDS', order_by='start'):
                    s = c.sequence(args.genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                    seq += s
                seq = Seq(seq)
                if t.strand == '-':
                    seq= seq.reverse_complement()
                cds_seq.loc[index] = [t.id,t.chrom,t.start,t.end,t.strand,seq]
                index += 1
                cdsRecord = SeqRecord(seq, id=t.id, description='strand %s start %d end %d length=%d'%(t.strand, t.start, t.end, len(seq)))
                if args.print:
                    SeqIO.write([cdsRecord], sys.stdout, "fasta") 
                else:
                    #cds_seq.to_csv(args.output, sep=',', index=False)
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
                table_translate = feature.qualifiers['transl_table'][0] if 'transl_table' in feature.qualifiers else 1
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
                protein_id = feature.qualifiers['protein_id'][0] if 'protein_id' in feature.qualifiers else 'Null'
                gene_id = feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else 'Null'
                cds_seq_record = SeqRecord(cds_seq, id='gene:%s protein_id:%s'%(gene_id, protein_id), description='strand %s length %d'%(feature.strand, len(cds_seq)))
                pep_record = SeqRecord(pep, id='gene:%s protein_id:%s'%(gene_id, protein_id), description='strand %s length %d'%(feature.strand, len(pep)))
                cds.append(cds_seq_record)
                proteins.append(pep_record)
                #break 
    if args.print and args.format == 'dna':
        SeqIO.write(cds, sys.stdout, 'fasta')
    elif args.print and args.format == 'protein':
        SeqIO.write(proteins, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(proteins, args.output, 'fasta')
