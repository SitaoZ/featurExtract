# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict

def get_exon(args, db, genome, transcript_id ,output):
    '''
    parameters:
     
    '''
    # header
    exon_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','Exon'])
    if transcript_id:
        # return a specific transcript
        out = [] 
        for t in db.features_of_type('mRNA', order_by='start'):
            if transcript_id in t.id:
                # exon
                exon_index = 1
                for e in db.children(t, featuretype='exon', order_by='start'):
                    exon = e.sequence(genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                    exon = Seq(exon)
                    if t.strand == '-':
                        exon = exon.reverse_complement()
                    exonRecord = SeqRecord(exon,id=transcript_id, description='strand %s exon %d start %d end %d length=%d'%(t.strand, exon_index, e.start, e.end, len(exon)))
                    out.append(exonRecord)
                    exon_index += 1
                break 
        if not args.print:
            SeqIO.write(out, output, "fasta")
        else:
            SeqIO.write(out, sys.stdout, "fasta")
    else:
        print('transcript id needed!')
        sys.exit(1)
        
