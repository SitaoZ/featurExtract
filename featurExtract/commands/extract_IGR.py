# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from featurExtract.database.database import create_db, genome_dict

def get_IGR(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    # GTF noncoding are annotated transcript,while GFF noncoding ared annotated gene
    if not args.style:
        sys.stderr.write("parameter -s should be assign \n")
        sys.exit()
    genomeDict = genome_dict(args.genome) # load fasta 
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    gene_seq = pd.DataFrame(columns=['GeneID','Chrom','Start','End','Strand','Seq'])
    index = 0
    seq_out = []
    gene_plus = []
    gene_minus = []
    for g in db.features_of_type('gene', order_by='start'):
        #seq = g.sequence(args.genome, use_strand=False)
        #seq = Seq(seq)
        start = g.start
        end = g.end
        if g.strand == '+':
            gene_plus.append((start, end, g.chrom, g.id, '+'))
        elif g.strand == '-':
            gene_minus.append((start, end, g.chrom, g.id, '-'))
            #seq = seq.reverse_complement()
            #gene_seq.loc[index] = [g.id,g.chrom,g.start,g.end,g.strand,seq]
            #index += 1
        else:
            sys.stderr.write("No strand indicated %d-%d. Assuming +\n" % (start, end))
            gene_plus.append((start, end, 1))

    # 注意：这里循环时，需要捕获起一个基因的终止和下一个基因的起始，
    # 所以gene_plus[1:]是从第二位置开始循环
    for i, pospair in enumerate(gene_plus[1:]):
        last_end = gene_plus[i][1]
        last_id = gene_plus[i][3]
        this_start = pospair[0]
        chrom = pospair[2]
        this_id = pospair[3]
        strand = pospair[4]
        # 处理正链 threshold length
        if this_start - last_end >= args.IGR_length:
            intergenic_seq = genomeDict[chrom][last_end:this_start]
            intergenic_record = SeqRecord(intergenic_seq, 
                                id="%s--%s"%(last_id, this_id), 
                                description='strand:%s start:%d end:%d length=%d'%(
                                strand, last_end+1, this_start, len(intergenic_seq)))
            seq_out.append(intergenic_record)
    for i,pospair in enumerate(gene_minus[1:]):
        last_end = gene_minus[i][1]
        last_id = gene_minus[i][3]
        this_start = pospair[0]
        chrom = pospair[2]
        this_id = pospair[3]
        strand = pospair[4]
        # 处理负链
        if this_start - last_end >= args.IGR_length:
            intergenic_seq = genomeDict[chrom][last_end:this_start]
            if strand == '-':
                intergenic_seq.reverse_complement()
            intergenic_record = SeqRecord(intergenic_seq, 
                                id="%s--%s"%(last_id, this_id), 
                                description='strand:%s start:%d end:%d length=%d'%(
                                strand, last_end+1, this_start, len(intergenic_seq)))
            seq_out.append(intergenic_record)
    
    if args.print:
        SeqIO.write(seq_out, sys.stdout, "fasta") 
    else:
        SeqIO.write(seq_out, args.output, "fasta") 

def get_IGR_gb(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    gene = []
    for record in create_db(args.genbank):
        for feature in record.features:
            if feature.type == 'gene': # CDS promoter UTR 
                gene_seq = ''
                for part in feature.location.parts:
                    gene_seq += part.extract(record.seq)
                # part.strand 会将FeatureLocation -1的反向互补
                # 判断提取的和已知的是否一致
                gene_id = feature.qualifiers['gene'][0] if 'gene' in feature.qualifiers else 'Null'
                gene_seq_record = SeqRecord(gene_seq, id='gene:%s'%(gene_id), 
                                  description='strand %s length %d'%(feature.strand, len(gene_seq)))
                gene.append(gene_seq_record)
    if args.print and args.format == 'dna':
        SeqIO.write(gene, sys.stdout, 'fasta')
    elif args.output:
        SeqIO.write(gene, args.output, 'fasta')
