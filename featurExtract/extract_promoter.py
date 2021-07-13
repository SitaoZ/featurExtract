# -*- coding: utf-8 -*-
import sys
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict



def get_promoter(args, db, genomeDict, genome, gene_id, promoter_length, utr_head_length, output):
    '''
    parameters:
     db: database generated by gffutils 
     genome: genome fasta 
     gene_id: gene id 
     output: output file
    '''
    promoter_seq = pd.DataFrame(columns=['GeneID','Chrom','Start','End','Strand','Promoter'])
    index = 0 
    if not gene_id:
        for g in db.all_features(featuretype='gene', order_by="seqid"):
            if g.strand == "+":
                gene_seq = g.sequence(genome, use_strand=True)
                p_start = g.start - int(promoter_length)
                if p_start < 0 :
                    continue
                p_end = g.start + utr_head_length
                promoter = genomeDict[g.chrom][p_start:p_end]
                # print(promoter)
            elif g.strand == "-":
                gene_seq = g.sequence(genome, use_strand=True)
                p_start = g.end - utr_head_length
                p_end = g.end + int(promoter_length)
                promoter = genomeDict[g.chrom][p_start:p_end].reverse_complement()
                # print(promoter)
            p_start_in_genome = p_start
            p_end_in_genome = p_end
            promoter_seq.loc[index] = [g.id,g.chrom,p_start_in_genome,p_end_in_genome,g.strand,promoter]
            index += 1
        promoter_seq.to_csv(output, sep=',', index=False)
    else:
        for g in db.all_features(featuretype='gene', order_by="seqid"):
            if gene_id in g.id:
                if g.strand == "+":
                    gene_seq = g.sequence(genome, use_strand=True)
                    p_start = g.start - int(promoter_length)
                    if p_start < 0 :
                        continue
                    p_end = g.start + utr_head_length
                    promoter = genomeDict[g.chrom][p_start:p_end]
                    # print(promoter)
                elif g.strand == "-":
                    gene_seq = g.sequence(genome, use_strand=True)
                    p_start = g.end - utr_head_length
                    p_end = g.end + int(promoter_length)
                    promoter = genomeDict[g.chrom][p_start:p_end].reverse_complement()
                p_start_in_genome = p_start
                p_end_in_genome = p_end
                promoter_seq.loc[index] = [g.id,g.chrom,p_start_in_genome,p_end_in_genome,g.strand,promoter]
                index += 1
                promoterSeq = SeqRecord(promoter,id=gene_id, description='chrom %s strand %s promoter start %d end %d length=%d'%(g.chrom, g.strand, p_start_in_genome, p_end_in_genome, len(promoter)))
                if args.print:
                    SeqIO.write(promoterSeq, sys.stdout, "fasta")
                else:
                    SeqIO.write(promoterSeq, output, "fasta")

                '''
                if args.print:
                    print(">{} {} {} {} {}".format(g.id,g.chrom,p_start_in_genome,p_end_in_genome,g.strand))
                    print(promoter)
                else:
                    promoter_seq.to_csv(output, sep=',', index=False)
                '''
                break
