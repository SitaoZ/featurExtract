# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from featurExtract.database.database import genome_dict
from featurExtract.utils.util import utr3_type, utr5_type, mRNA_type

def intron(genomeDict, chrom, transcript_id, strand, exons):
    '''
    '''
    introns = [] 
    left = exons[1:-1]
    intron_index = 1
    for i in range(0,len(left), 2):
        s = left[i]   # intron start
        e = left[i+1] # intron end, shoukd be minus 1 
        if s+1 != e:
            # exon 与exon 不相连，就是intron的位置
            # 注意索引的点，
            intron_seq= genomeDict[chrom][s:e-1]
            #intron_seq= Seq(intron_seq)
            # 1 based position
            intronRecord = SeqRecord(intron_seq, id=transcript_id, 
                                     description='strand %s intron %d start %d end %d length=%d'%(
                                     strand, intron_index, s+1, e-1, len(intron_seq)
                                        )
                                    )
            introns.append(intronRecord)
            intron_index += 1
    return introns 
        

def get_intron(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database 
    genomeDict = genome_dict(args.genome) # load fasta
    mRNA_str = mRNA_type(args.style)
    intron_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','Exon']) # header 
    if args.transcript:
        # return a specific transcript
        for t in db.features_of_type(mRNA_str, order_by='start'):
            introns = ''
            if args.transcript in t.id:
                exons = [] # exon position include start and end
                for e in db.children(t, featuretype='exon', order_by='start'):
                    exons.append(e.start)
                    exons.append(e.end)
                introns = intron(genomeDict, t.chrom, args.transcript, t.strand, exons)
                break 
        if not args.print:
            SeqIO.write(introns, args.output, "fasta")
        else:
            SeqIO.write(introns, sys.stdout, "fasta")
    else:
        whole_introns = []
        for t in db.features_of_type(mRNA_str, order_by='start'):
            if t.id:
                exons = [] # exon position include start and end
                for e in db.children(t, featuretype='exon', order_by='start'):
                    exons.append(e.start)
                    exons.append(e.end)
                introns = intron(genomeDict, t.chrom, t.id, t.strand, exons)
                whole_introns.extend(introns)
        if not args.print:
            SeqIO.write(whole_introns, args.output, "fasta")
        else:
            SeqIO.write(whole_introns, sys.stdout, "fasta") 
        
