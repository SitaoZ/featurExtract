# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from tqdm import tqdm 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from featurExtract.utils.util import mRNA_type
from featurExtract.database.database import create_db

_CSV_HEADER = ['TranscriptID','Chrom','Start','End','Strand','Exon']

def get_exon(args):
    '''
    parameters:
     
    '''
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    # exon_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','Exon']) # header
    exon_seq_list = []
    mRNA_str = mRNA_type(args.style)
    if args.transcript:
        # return a specific transcript
        out = [] 
        for t in tqdm(db.features_of_type(mRNA_str, order_by='start'), \
                      total = len(list(db.features_of_type(mRNA_str, order_by='start'))), \
                      ncols = 80, desc = "Exon Processing:"):
            if args.transcript in t.id:
                # exon
                exon_index = 1
                for e in db.children(t, featuretype='exon', order_by='start'):
                    exon = e.sequence(args.genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                    exon = Seq(exon)
                    if t.strand == '-':
                        exon = exon.reverse_complement()
                    exonRecord = SeqRecord(exon,id=args.transcript, 
                                 description='strand %s exon %d start %d end %d length=%d'%(t.strand, 
                                             exon_index, e.start, e.end, len(exon)))
                    out.append(exonRecord)
                    exon_index += 1
                break 
        if not args.print:
            SeqIO.write(out, args.output, "fasta")
        elif args.output:
            SeqIO.write(out, args.output, "fasta")
        else:
            SeqIO.write(out, sys.stdout, "fasta")
    else:
        whole_exons = []
        for t in tqdm(db.features_of_type(mRNA_str, order_by='start'), \
                      total = len(list(db.features_of_type(mRNA_str, order_by='start'))), \
                      ncols = 80, desc = "Exon Processing:"):
            exon_index = 1
            for e in db.children(t, featuretype='exon', order_by='start'):
                exon = e.sequence(args.genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                exon = Seq(exon)
                if t.strand == '-':
                    exon = exon.reverse_complement()
                exonRecord = SeqRecord(exon,id=t.id,
                             description='strand %s exon %d start %d end %d length=%d'%(t.strand,
                                        exon_index, e.start, e.end, len(exon)))
                whole_exons.append(exonRecord)
                exon_index += 1
        if not args.print:
            SeqIO.write(whole_exons, args.output, "fasta")
        elif args.output:
            SeqIO.write(whole_exons, args.output, "fasta")
        else:
            SeqIO.write(whole_exons, sys.stdout, "fasta")
        

def get_exon_gb(args):
    exons = []
    for record in create_db(args.genbank):
        index = 1
        print(type(record))
        for feature in record.features:
            if feature.type == 'exon': # CDS promoter UTR 
                exon_seq = feature.extract(record.seq)
                # feature.strand 会将FeatureLocation -1的反向互补
                # 判断提取的和已知的是否一致
                exon_seq_record = SeqRecord(exon_seq,id='%s exon %d'%(record.id, index),
                                 description='strand %s length %d'%(feature.strand, len(exon_seq)))
                exons.append(exon_seq_record)
                index += 1 
    if args.print:
        SeqIO.write(exons, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(exons, args.output, 'fasta')
 
