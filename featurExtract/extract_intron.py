# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

def intron(genomeDict, chrom, transcript_id, strand, exons):
    '''
    '''
    introns = [] 
    left = exons[1:-1]
    intron_index = 1
    for i in range(0,len(left), 2):
        s = left[i]
        e = left[i+1]
        if s+1 != e:
            # exon 与exon 不相连，就是intron的位置
            intron_seq= genomeDict[chrom][s-1:e]
            #intron_seq= Seq(intron_seq)
            intronRecord = SeqRecord(intron_seq, id=transcript_id, description='strand %s intron %d start %d end %d length=%d'%(strand, intron_index, s, e, len(intron_seq)))
            introns.append(intronRecord)
            intron_index += 1
    return introns 
        

def get_intron(args, db, genome, genomeDict, transcript_id ,output):
    '''
    parameters:
     
    '''
    # header
    intron_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','Exon'])
    if transcript_id:
        # return a specific transcript
        for t in db.features_of_type('mRNA', order_by='start'):
            introns = ''
            if transcript_id in t.id:
                exons = [] # exon position include start and end
                for e in db.children(t, featuretype='exon', order_by='start'):
                    exons.append(e.start)
                    exons.append(e.end)
                introns = intron(genomeDict, t.chrom, transcript_id, t.strand, exons)
                break 
        if not args.print:
            SeqIO.write(introns, output, "fasta")
        else:
            SeqIO.write(introns, sys.stdout, "fasta")
    else:
        print('transcript id needed!')
        sys.exit(1)
        
