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

def get_terminator(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    genomeDict = genome_dict(args.genome) # load fasta 
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    terminator_seq = pd.DataFrame(columns=['GeneID','Chrom','Start','End','Strand','Promoter'])
    terminator_record = []
    index = 0 
    if not args.gene:
        for g in db.all_features(featuretype='gene', order_by="seqid"):
            if g.strand == "+":
                gene_seq = g.sequence(args.genome, use_strand=True)
                t_start = g.end - args.utr3_upper_length     # transform 1-based to 0-based
                t_end = g.end + int(args.terminator_length)  # transform 1-based to 0-based
                if t_start < 0 :
                    continue
                terminator = genomeDict[g.chrom][t_start:t_end]
            elif g.strand == "-":
                gene_seq = g.sequence(args.genome, use_strand=True)
                t_start = g.start - int(args.terminator_length)
                t_end = g.start + args.utr3_upper_length 
                terminator = genomeDict[g.chrom][t_start:t_end].reverse_complement()
            t_start_in_genome = t_start
            t_end_in_genome = t_end
            if args.output_format == 'csv':
                terminator_seq.loc[index] = [g.id.replace('transcript:',''), g.chrom, t_start_in_genome,
                                                                  t_end_in_genome, g.strand, terminator]
                index += 1
            # default fasta
            else:
                desc='chrom:%s strand:%s gene:%d-%d terminator:%d-%d length=%d'%(g.chrom, g.strand, 
                                                         g.start, g.end,
                                                         p_start_in_genome,
                                                         p_end_in_genome, 
                                                         len(terminator))
                terminatorRecord = SeqRecord(terminator, id=g.id.replace('transcript:',''), description=desc)
                terminator_record.append(terminatorRecord)
        # csv
        if args.output_format == 'csv':
            terminator_seq.to_csv(args.output, sep=',', index=False)
        # default fasta
        else:
            SeqIO.write(terminator_record, args.output, "fasta")
    else:
        for g in db.all_features(featuretype='gene', order_by="seqid"):
            if args.gene in g.id:
                if g.strand == "+":
                    gene_seq = g.sequence(args.genome, use_strand=True)
                    t_start = g.end - args.utr3_upper_length # 往后数，和promoter不一样，不需要-1
                    if t_start < 0 :
                        continue
                    t_end = g.end + int(args.terminator_length)
                    terminator = genomeDict[g.chrom][t_start:t_end]
                elif g.strand == "-":
                    gene_seq = g.sequence(args.genome, use_strand=True)
                    t_start = g.start - int(args.terminator_length)
                    t_end = g.start + args.utr3_upper_length 
                    terminator = genomeDict[g.chrom][t_start:t_end].reverse_complement() # index 0-based
                t_start_in_genome = t_start
                t_end_in_genome = t_end
                terminator_seq.loc[index] = [args.gene, g.chrom, t_start_in_genome,
                                          t_end_in_genome, g.strand, terminator]
                index += 1
                terminatorSeq = SeqRecord(terminator,id=args.gene, 
                              description='chrom:%s strand:%s gene:%d-%d terminator:%d-%d length=%d'%(
                                  g.chrom, g.strand, 
                                  g.start, g.end,
                                  t_start_in_genome + 1 , 
                                  t_end_in_genome + 1, 
                                  len(terminator)))
                if args.print:
                    if args.output_format == 'csv':
                        terminator_seq.to_csv(sys.stdout, sep=',', index=False)
                    else:
                        SeqIO.write(terminatorSeq, sys.stdout, "fasta")
                else:
                    if args.output_format == 'csv':
                        terminator_seq.to_csv(args.output, sep=',', index=False)
                    else:
                        SeqIO.write(terminatorSeq, args.output, "fasta")
                break
