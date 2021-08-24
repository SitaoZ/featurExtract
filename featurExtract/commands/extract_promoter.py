# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from featurExtract.database.database import genome_dict

# https://pythonhosted.org/gffutils/autodocs/gffutils.Feature.html
# 1-based coordinates

def get_promoter(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    genomeDict = genome_dict(args.genome) # load fasta 
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    promoter_seq = pd.DataFrame(columns=['GeneID','Chrom','Start','End','Strand','Promoter'])
    promoter_record = []
    index = 0 
    if not args.gene:
        for g in db.all_features(featuretype='gene', order_by="seqid"):
            if g.strand == "+":
                gene_seq = g.sequence(args.genome, use_strand=True)
                p_start = g.start - int(args.promoter_length) - 1 # 往前数,和terminator不同需要-1, transform 1-based to 0-based
                if p_start < 0 :
                    continue
                p_end = g.start + args.utr5_upper_length - 1      # include first base of TSS, transform 1-based to 0-based
                promoter = genomeDict[g.chrom][p_start:p_end]
            elif g.strand == "-":
                gene_seq = g.sequence(args.genome, use_strand=True)
                p_start = g.end - args.utr5_upper_length - 1 
                p_end = g.end + int(args.promoter_length) - 1
                promoter = genomeDict[g.chrom][p_start:p_end].reverse_complement()
            p_start_in_genome = p_start
            p_end_in_genome = p_end
            if args.output_format == 'csv':
                promoter_seq.loc[index] = [g.id.replace('transcript:',''), g.chrom, p_start_in_genome,
                                                                  p_end_in_genome, g.strand, promoter]
                index += 1
            # default fasta
            else:
                desc='strand:%s start:%d end:%d length=%d'%(g.strand, p_start_in_genome,
                                                         p_end_in_genome, len(promoter))
                promoterRecord = SeqRecord(promoter, id=g.id.replace('transcript:',''), description=desc)
                promoter_record.append(promoterRecord)
        # csv
        if args.output_format == 'csv':
            promoter_seq.to_csv(args.output, sep=',', index=False)
        # default fasta
        else:
            SeqIO.write(promoter_record, args.output, "fasta")
    else:
        for g in db.all_features(featuretype='gene', order_by="seqid"):
            if args.gene in g.id:
                if g.strand == "+":
                    gene_seq = g.sequence(args.genome, use_strand=True)
                    p_start = g.start - int(args.promoter_length) - 1
                    if p_start < 0 :
                        continue
                    p_end = g.start + args.utr5_upper_length - 1
                    promoter = genomeDict[g.chrom][p_start:p_end]
                elif g.strand == "-":
                    gene_seq = g.sequence(args.genome, use_strand=True)
                    p_start = g.end - args.utr5_upper_length - 1
                    p_end = g.end + int(args.promoter_length) - 1
                    promoter = genomeDict[g.chrom][p_start:p_end].reverse_complement()
                p_start_in_genome = p_start
                p_end_in_genome = p_end
                promoter_seq.loc[index] = [args.gene, g.chrom, p_start_in_genome,
                                          p_end_in_genome, g.strand, promoter]
                index += 1
                promoterSeq = SeqRecord(promoter,id=args.gene, 
                              description='chrom:%s strand:%s gene:%d-%d promoter:%d-%d length=%d'%(
                              g.chrom, g.strand, g.start, g.end ,
                              p_start_in_genome, p_end_in_genome, len(promoter)))
                if args.print:
                    if args.output_format == 'csv':
                        promoter_seq.to_csv(sys.stdout, sep=',', index=False)
                    else:
                        SeqIO.write(promoterSeq, sys.stdout, "fasta")
                else:
                    if args.output_format == 'csv':
                        promoter_seq.to_csv(args.output, sep=',', index=False)
                    else:
                        SeqIO.write(promoterSeq, args.output, "fasta")
                break
