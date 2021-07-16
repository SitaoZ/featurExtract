# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from featurExtract.database import genome_dict


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
    index = 0 
    if not args.gene:
        for g in db.all_features(featuretype='gene', order_by="seqid"):
            if g.strand == "+":
                gene_seq = g.sequence(args.genome, use_strand=True)
                p_start = g.start - int(args.promoter_length)
                if p_start < 0 :
                    continue
                p_end = g.start + args.utr5_upper_length
                promoter = genomeDict[g.chrom][p_start:p_end]
                # print(promoter)
            elif g.strand == "-":
                gene_seq = g.sequence(args.genome, use_strand=True)
                p_start = g.end - args.utr5_upper_length
                p_end = g.end + int(args.promoter_length)
                promoter = genomeDict[g.chrom][p_start:p_end].reverse_complement()
                # print(promoter)
            p_start_in_genome = p_start
            p_end_in_genome = p_end
            promoter_seq.loc[index] = [g.id,g.chrom,p_start_in_genome,p_end_in_genome,g.strand,promoter]
            index += 1
        promoter_seq.to_csv(args.output, sep=',', index=False)
    else:
        for g in db.all_features(featuretype='gene', order_by="seqid"):
            if args.gene in g.id:
                if g.strand == "+":
                    gene_seq = g.sequence(args.genome, use_strand=True)
                    p_start = g.start - int(args.promoter_length)
                    if p_start < 0 :
                        continue
                    p_end = g.start + args.utr5_upper_length
                    promoter = genomeDict[g.chrom][p_start:p_end]
                    # print(promoter)
                elif g.strand == "-":
                    gene_seq = g.sequence(args.genome, use_strand=True)
                    p_start = g.end - args.utr5_upper_length
                    p_end = g.end + int(args.promoter_length)
                    promoter = genomeDict[g.chrom][p_start:p_end].reverse_complement()
                p_start_in_genome = p_start
                p_end_in_genome = p_end
                promoter_seq.loc[index] = [g.id,g.chrom,p_start_in_genome,p_end_in_genome,g.strand,promoter]
                index += 1
                promoterSeq = SeqRecord(promoter,id=args.gene, 
                              description='chrom %s strand %s promoter start %d end %d length=%d'%(
                              g.chrom, g.strand, p_start_in_genome, p_end_in_genome, len(promoter)))
                if args.print:
                    SeqIO.write(promoterSeq, sys.stdout, "fasta")
                else:
                    SeqIO.write(promoterSeq, args.output, "fasta")

                '''
                if args.print:
                    print(">{} {} {} {} {}".format(g.id,g.chrom,p_start_in_genome,p_end_in_genome,g.strand))
                    print(promoter)
                else:
                    promoter_seq.to_csv(args.output, sep=',', index=False)
                '''
                break
