# -*- coding: utf-8 -*-
import sys
import pandas as pd 
from Bio.Seq import Seq
from collections import defaultdict

# table header 
['TranscriptID','Chrom', 'Strand','CDS Interval','uORF Start','uORF End','uORF Type','uORF Length','uORF']

def dorf(transcript_id, chrom, strand, matural_transcript, coding_sequence):
    '''
    parameters:
     transcript_id:      transcript id
     matural_transcript: a Seq type (Biopython) from mature transcript without intron
     coding_sequence:    a Seq type (Biopython) from coding sequence start with ATG , 
                         end with TAA, TGA, TAG
    return:
     down stream open reading frame in utr3
    '''
    dORF_dict = defaultdict(list)
    stop_codon_list = ['TAA','TAG','TGA']
    # start_codon means the first base position in matural_transcript
    start_codon = matural_transcript.index(coding_sequence)
    # stop_codon means the last base position in matural transcript
    cds_len = len(coding_sequence)
    stop_codon = start_codon + cds_len
    cds_intervel = str(start_codon) + '-' + str(stop_codon)
    mt_len = len(matural_transcript)
    utr5 = matural_transcript[:start_codon]
    utr5_len = len(utr5)
    utr3 = matural_transcript[stop_codon:]
    utr3_len = len(utr3)
    for i in range(0, utr3_len-3):
        # start codon find 
        if utr3[i:i+3] == "ATG":
            for j in range(i+3, utr3_len, 3):
            # stop codon find
                if utr3[j:j+3] in stop_codon_list:
                    # type1 uORF  upstream; not unique 
                    dORF = utr3[i:j+3]
                    out1 = [transcript_id, chrom, strand, cds_intervel, stop_codon+i+1, stop_codon+j+3, 'dORF', len(dORF), dORF]
                    if not dORF_dict.get('dORF'):
                        dORF_dict['dORF'] = [out1]
                    else:
                        dORF_dict['dORF'].append(out1)
                    break 
    return dORF_dict


def get_dorf(db, genome, transcript_id, output):
    '''
    parameters:
    
    '''
    dORF_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Strand','CDS Interval','dORF Start','dORF End','Type','dORF Length','dORF'])
    index = 0
    if not transcript_id:
        for t in db.features_of_type('mRNA', order_by='start'):
            # primary transcript (pt) 是基因组上的转录本的序列，
            # 有的会包括intron，所以提取的序列和matural transcript 不一致
            # print(t)
            pt = t.sequence(genome, use_strand=True)
            # matural transcript (mt)
            # exon 提取的是转录本内成熟的MRNA的序列,即外显子收尾相连的matural transcript
            mt = ''
            for e in db.children(t, featuretype='exon', order_by='start'):
                s = e.sequence(genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                mt += s
            mt = Seq(mt)
            if t.strand == '-':
                mt = mt.reverse_complement()
            # CDS 提取的是编码区 ATG ... TAG
            cds = ''
            for c in db.children(t, featuretype='CDS', order_by='start'):
                #print(c)
                s = c.sequence(genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                cds += s
            cds = Seq(cds)
            if t.strand == '-':
                cds = cds.reverse_complement()
            dORF_dict = dorf(t.id, t.chrom, t.strand, mt, cds)
            for key in sorted(dORF_dict.keys()):
                #print(key,len(uORF_dict[key]))
                for it in dORF_dict[key]:
                    #print(it)
                    dORF_seq.loc[index] = it
                    index += 1
        dORF_seq.to_csv(output, sep=',', index=False)
    else:
        for t in db.features_of_type('mRNA', order_by='start'):
            # primary transcript (pt) 是基因组上的转录本的序列，
            # 有的会包括intron，所以提取的序列和matural transcript 不一致
            # print(t)
            if transcript_id in t.id:
                pt = t.sequence(genome, use_strand=True)
                # matural transcript (mt)
                # exon 提取的是转录本内成熟的MRNA的序列,即外显子收尾相连的matural transcript
                mt = ''
                for e in db.children(t, featuretype='exon', order_by='start'):
                    s = e.sequence(genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                    mt += s
                mt = Seq(mt)
                if t.strand == '-':
                    mt = mt.reverse_complement()
                # CDS 提取的是编码区 ATG ... TAG
                cds = ''
                for c in db.children(t, featuretype='CDS', order_by='start'):
                    #print(c)
                    s = c.sequence(genome, use_strand=False) # 不反向互补，对于负链要得到全部的cds后再一次性反向互补
                    cds += s
                cds = Seq(cds)
                if t.strand == '-':
                    cds = cds.reverse_complement()
                dORF_dict = dorf(t.id, t.chrom, t.strand, mt, cds)
                for key in sorted(dORF_dict.keys()):
                    #print(key,len(uORF_dict[key]))
                    for it in dORF_dict[key]:
                        #print(it)
                        dORF_seq.loc[index] = it
                        index += 1
                dORF_seq.to_csv(output, sep=',', index=False)
                break 
