# -*- coding: utf-8 -*-
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from featurExtract.utils.util import utr3_type, utr5_type, mRNA_type
from featurExtract.utils.util import stop_codon_seq, add_stop_codon
from featurExtract.utils.visualize import visual_transcript, visual_transcript_genome

# table header 
['TranscriptID','Chrom', 'Strand','CDS Interval','uORF Start','uORF End','uORF Type','uORF Length','uORF']

class GFF(object):
    '''
    a parse gff from line of uorf
    # uorf_out format "TranscriptID,Chrom,Strand,CDS Interval,uORF Start,uORF End,uORF Type,uORF Length,uORF"
    '''
    def __init__(self, exons, uorf_out):
        self.exons = exons
        self.uorf_out = uorf_out
        self.t_id = self.uorf_out[0].replace('transcript:','uORF.')
        self.chrom = self.uorf_out[1] 
        self.strand = self.uorf_out[2]
        self.start = self.uorf_out[4]
        self.end = self.uorf_out[5]
        self.uorf_type = self.uorf_out[6]
        self.uorf_id = ".".join(map(str,[self.t_id, self.start, self.end, self.uorf_type]))
    def parse(self):
        cds_interval = self.uorf_out[3]
        uORF_start, uORF_end, uORF_type = self.uorf_out[4:7]
        # build a class for uROF genome location 
        ugl = uORF_genome_location(self.exons, uORF_start, uORF_end, self.strand, cds_interval)
        genome_start, genome_end = ugl.parse()
        features = ugl.split_feature() # uORF exons blocks, genome location 
        features_gff_lines = self.feature_gff_format(features)
        #if strand == '-':
            # exchange location, because start less then end in gff 
        #    genome_start, genome_end = genome_end, genome_start
        # gff annotation line https://genome.ucsc.edu/FAQ/FAQformat#format3
        # col9 : id start end type 
        uorf_gff_line = "\t".join([self.chrom, "featurExtract", "uORF", \
                   str(genome_start), str(genome_end), ".", self.strand, "0", "ID="+self.uorf_id])
        
        return uorf_gff_line,features_gff_lines
        
    def uorf_exons_location(self):
        '''
        extract uorf exons blocks in genome, no reverse 
        '''
        cds_interval = self.uorf_out[3]
        uORF_start, uORF_end, uORF_type = self.uorf_out[4:7]
        # build a class for uROF genome location 
        ugl = uORF_genome_location(self.exons, uORF_start, uORF_end, self.strand, cds_interval)
        genome_start, genome_end = ugl.parse()
        features = ugl.split_feature()
        return features
    
    def feature_gff_format(self, features):
        '''
        formate the split feature into gff line
        '''
        gff_lines = []
        if self.strand == "+":
            for i,exon in enumerate(features):
                gff_9_col = "ID="+self.uorf_id+":cds:"+str(i)+";Parent="+self.uorf_id
                gff_line = "\t".join([self.chrom, "featurExtract", "CDS", \
                   str(exon[0]), str(exon[1]), ".", self.strand, "0", gff_9_col])
                gff_lines.append(gff_line)
            return gff_lines
        else:
            # gff annotation line https://genome.ucsc.edu/FAQ/FAQformat#format3
            features = features[::-1] # reverse, ??????????????????
            for i,exon in enumerate(features):
                cds_index = len(features) - i 
                gff_9_col = "ID="+self.uorf_id+":cds:"+str(cds_index)+";Parent="+self.uorf_id
                gff_line = "\t".join([self.chrom, "featurExtract", "CDS", \
                   str(exon[0]), str(exon[1]), ".", self.strand, "0", gff_9_col])
                gff_lines.append(gff_line)
            return gff_lines
        

class uORF_genome_location(object):
    '''
    a parse for uORF location from transcript map to genome 
    uORF????????????????????????????????????????????????,??????????????????????????????????????????????????????,
    ??????uORF???????????????????????????????????????????????????+???????????????exon???????????????,
    uORF???????????????????????????????????????????????????+???????????????exon???????????????.
    '''
    def __init__(self, exons, uorf_start, uorf_end, strand, cds_interval):
        self.exons = exons
        self.uorf_start = uorf_start
        self.uorf_end = uorf_end
        self.strand = strand 
        # strand affect the utr5 location, reverse -
        if self.strand == '-':
            self.exons = self.exons[::-1]
        self.cds_interval = cds_interval # for 5utr and 3utr
        if self.cds_interval:
            cdss = self.cds_interval.split('-')
            self.cds_start = cdss[0]
            self.cds_end = cdss[1]
    def parse(self):
        '''
        uORF transcript location to genome location 
        '''
        genome_start = self.location(self.uorf_start)
        #genome_start = self.exon_index(self.uorf_start)
        genome_end = self.location(self.uorf_end)
        #genome_end = self.exon_index(self.uorf_end)
        
        if self.strand == '-':
            # exchange location, because start less then end in gff 
            genome_start, genome_end = genome_end, genome_start 
        return genome_start, genome_end
        
    def location(self, position):
        '''
        used for uORF start and end mapping to genome location  
        Params:
         - position: the uORF start and end position in cDNA(mRNA)
        Return:
         - the uORF start or uORF end posiiton in genome 
        ''' 
        exon_len_sum = 0
        genome_location = 0
        for i, exon in enumerate(self.exons):
            exon_len = exon[1] - exon[0] + 1 # 1 based; length should add 1
            prior_exon_len = exon_len_sum    # pay attention to the coordinate, assignment should be prior to sum
            exon_len_sum += exon_len 
            if position <= prior_exon_len:
                # uORF???????????????????????????
                if self.strand == '+':
                    genome_location = exon[0] + position - 1 # 1 based
                else:
                    genome_location = exon[1] - position + 1
                return genome_location
            # ????????????????????????????????????????????????
            elif prior_exon_len < position <= exon_len_sum:
                # UTR ???????????????, uORF ??????????????????????????????
                # length ?????????????????????????????????????????????????????????????????????
                length = position - prior_exon_len
                if self.strand == '+':
                    genome_location = exon[0] + length - 1 # 1 based 
                else:
                    genome_location = exon[1] - length + 1
                return genome_location
    def exon_index(self, position):
        '''
        Position locate the exon in mRNA
        Params:
         - position: the index in cDNA(mRNA)
        Return:
         - the exon index containing position
        '''
        exon_len_sum = 0
        for i, exon in enumerate(self.exons):
            exon_len = exon[1] - exon[0] + 1 # 1 based; length should add 1
            prior_exon_len = exon_len_sum    # pay attention to the coordinate, assignment should be prior to sum
            exon_len_sum += exon_len
            if prior_exon_len < position <= exon_len_sum:
                return i
    def split_feature(self):
        '''
        Split the uORF in feature including exon and intron 
        '''
        uorf_g_start, uorf_g_end = self.parse()
        # true position in genome
        if self.strand == '-':
            uorf_g_start, uorf_g_end = uorf_g_end, uorf_g_start
        u_start_exon_index, u_end_exon_index = self.exon_index(self.uorf_start), self.exon_index(self.uorf_end)
        # cds_g_start, cds_g_end = self.location((self.cds_start), self.location(self.cds_end)
        if u_start_exon_index == u_end_exon_index:
            # uORF start and end in same exon, no intron
            # type1, type2, type3 uORF
            if self.strand == '+':
                return [[uorf_g_start,uorf_g_end]]
            else:
                return [[uorf_g_end,uorf_g_start]]
        elif u_start_exon_index < u_end_exon_index:
            # have at least once intron, CDS feature collect
            # type1, type2, type3 uORF
            if self.strand == '+':
                s = [[uorf_g_start, self.exons[u_start_exon_index][1]]]
                e = [[self.exons[u_end_exon_index][0], uorf_g_end]]
            else:
                s = [[self.exons[u_start_exon_index][0], uorf_g_start]]
                e = [[uorf_g_end, self.exons[u_end_exon_index][1]]]
            
            inters = []
            for i,inter in enumerate(self.exons[u_start_exon_index+1:u_end_exon_index]):
                if len(inter) != 0:
                    inters.append(inter)
            features = s + inters + e # ????????????
            return features
        else:
            sys.exit('exon index error')

class uORF(object):
    '''
    uORF finder using object-oriented
    '''
    START_CODON = 'ATG'
    STOP_CODON_LIST = ['TAA','TAG','TGA']
    
    def __init__(self, transcript_id, chrom, strand, matural_transcript, coding_sequence):
        self.transcript_id = transcript_id # transcript id, str
        self.chrom = chrom                 # chromosome id ,str
        self.strand = strand               # transcript strand, str
        self.mt = matural_transcript       # a Seq type (Biopython) from mature transcript without intron
        self.cds = coding_sequence         # a Seq type (Biopython) from coding sequence start with ATG ; end with TAA, TGA, TAG
    def transcript_location(self):
        '''
        transcript start and end 
        1d list regradless of strand  
        '''    
        return [0,len(self.mt)]
    def cds_location_transcript(self):
        '''
        extract the cds start and end in transcript
        1d list regardless of strand 
        '''
        start_codon_position = self.mt.index(self.cds)
        cds_start = start_codon_position + 1
        cds_end = start_codon_position + len(self.cds)
        return [cds_start, cds_end]
        
    def uorf_location_transcript(self):
        '''
        extract uorfs start and end in transcript
        2d list regardless strand 
        ''' 
        d2 = []
        uORF_dict = self.uorf_parse()
        for uorf_type in uORF_dict.keys():
            for it in uORF_dict[uorf_type]:
                uorf_start, uorf_end = it[4], it[5]
                d2.append([uorf_start, uorf_end])
        return d2 
        
    def uorf_parse(self):
        '''
        uorf finder main program 
        return a dict
        '''
        uORF_dict = defaultdict(list)
        # start_codon_position means the first base position in matural_transcript
        start_codon_position = self.mt.index(self.cds)
        # stop_codon_position means the last base position in matural transcript
        stop_codon_position = start_codon_position + len(self.cds)
        # CDS on mt coordinate
        cds_interval = str(start_codon_position + 1) + '-' + str(stop_codon_position) # transcript coordinate 1 based 
        mt_len = len(self.mt)
        cds_len = len(self.cds)
        utr5 = self.mt[:start_codon_position]
        utr5_len = len(utr5)
        utr3 = self.mt[stop_codon_position:]
        utr3_len = len(utr3)
        
        # UTR5 range for ATG searching
        for i in range(0, utr5_len-3+1):
            # start codon find 
            if self.mt[i:i+3] == uORF.START_CODON:
                for j in range(i+3, stop_codon_position, 3):
                # stop codon find
                    if self.mt[j:j+3] in uORF.STOP_CODON_LIST and j < utr5_len:
                        # type1 uORF  upstream; not unique 
                        type1_uORF = self.mt[i:j+3]
                        out1 = [self.transcript_id, self.chrom, self.strand, \
                                cds_interval, i+1, j+3, 'type1', len(type1_uORF), type1_uORF]
                        if not uORF_dict.get('type1_uORF'):
                            uORF_dict['type1_uORF'] = [out1]
                        else:
                            uORF_dict['type1_uORF'].append(out1)
                        # only use the first stop codon, so break must
                        break
                    if self.mt[j:j+3] in uORF.STOP_CODON_LIST and j + 3 > utr5_len and j + 3 < stop_codon_position:
                        # type2 uORF across; the overlap region is triple or not; not unique
                        type2_uORF = self.mt[i:j+3]
                        out2 = [self.transcript_id, self.chrom, self.strand, \
                                cds_interval, i+1, j+3, 'type2', len(type2_uORF), type2_uORF]
                        if not uORF_dict.get('type2_uORF'):
                            uORF_dict['type2_uORF'] = [out2]
                        else:
                            uORF_dict['type2_uORF'].append(out2)
                        break
                    if self.mt[j:j+3] in uORF.STOP_CODON_LIST and j + 3 > utr5_len and j+3 == stop_codon_position and (utr5_len - i)%3 == 0:
                        # type3 uORF; N extention ??????; not unique 
                        type3_uORF = self.mt[i:utr5_len+cds_len]
                        out3 = [self.transcript_id, self.chrom, self.strand, \
                                cds_interval, i+1, j+3, 'type3', len(type3_uORF), type3_uORF]
                        if not uORF_dict.get('type3_uORF'):
                            uORF_dict['type3_uORF'] = [out3]
                        else:
                            uORF_dict['type3_uORF'].append(out3)
                        break
        return uORF_dict
            


def uorf_procedure_oriented(transcript_id, chrom, strand, matural_transcript, coding_sequence):
    '''
    uorf finder using procedure-oriented
    parameters:
     transcript_id:      transcript id
     matural_transcript: a Seq type (Biopython) from mature transcript without intron
     coding_sequence:    a Seq type (Biopython) from coding sequence start with ATG , 
                         end with TAA, TGA, TAG
    return:
     upper stream open reading frame
    '''
    uORF_dict = defaultdict(list)
    stop_codon_list = ['TAA','TAG','TGA']
    # start_codon means the first base position in matural_transcript
    start_codon = matural_transcript.index(coding_sequence)
    # print('xxxx',start_codon)
    # stop_codon means the last base position in matural transcript
    cds_len = len(coding_sequence)
    stop_codon = start_codon + cds_len
    # print('yyyy',cds_len)
    cds_interval = str(start_codon+1) + '-' + str(stop_codon) # transcript coordinate 1 based 
    mt_len = len(matural_transcript)
    utr5 = matural_transcript[:start_codon]
    utr5_len = len(utr5)
    utr3 = matural_transcript[stop_codon:]
    # UTR5 range 
    for i in range(0, utr5_len-3+1):
        # start codon find 
        if matural_transcript[i:i+3] == "ATG":
            for j in range(i+3,stop_codon,3):
            # stop codon find
                if matural_transcript[j:j+3] in stop_codon_list and j < utr5_len:
                    # type1 uORF  upstream; not unique 
                    type1_uORF = matural_transcript[i:j+3]
                    out1 = [transcript_id, chrom, strand, cds_interval, i+1, j+3, 'type1', len(type1_uORF), type1_uORF]
                    if not uORF_dict.get('type1_uORF'):
                        uORF_dict['type1_uORF'] = [out1]
                    else:
                        uORF_dict['type1_uORF'].append(out1)
                    break 
                if matural_transcript[j:j+3] in stop_codon_list and j + 3 > utr5_len and j + 3 < stop_codon:
                    # type2 uORF across; the overlap region is triple or not; not unique
                    type2_uORF = matural_transcript[i:j+3]
                    out2 = [transcript_id, chrom, strand, cds_interval, i+1, j+3, 'type2', len(type2_uORF), type2_uORF]
                    if not uORF_dict.get('type2_uORF'):
                        uORF_dict['type2_uORF'] = [out2]
                    else:
                        uORF_dict['type2_uORF'].append(out2)
                    break 
                if matural_transcript[j:j+3] in stop_codon_list and j + 3 > utr5_len and j+3 == stop_codon and (utr5_len - i)%3 == 0:
                    # N extention ??????; not unique 
                    type3_uORF = matural_transcript[i:utr5_len+cds_len]
                    out3 = [transcript_id, chrom, strand, cds_interval, i+1, j+3, 'type3', len(type3_uORF), type3_uORF]
                    if not uORF_dict.get('type3_uORF'):
                        uORF_dict['type3_uORF'] = [out3]
                    else:
                        uORF_dict['type3_uORF'].append(out3)
                    break
    return uORF_dict



def get_uorf(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    db = gffutils.FeatureDB(args.database, keep_order=True)
    # csv
    uORF_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Strand','CDS Interval','uORF Start','uORF End','uORF Type','uORF Length','uORF'])
    # gff
    uORF_gff = open(args.output,'w')
    # fasta
    uorf_record = []
    index = 0
    mRNA_str = mRNA_type(args.style) 
    # not specify the transcript id
    if not args.transcript:
        for t in db.features_of_type(mRNA_str, order_by='start'):
            # primary transcript (pt) ???????????????????????????????????????
            # ???????????????intron???????????????????????????matural transcript ?????????
            # print(t)
            pt = t.sequence(args.genome, use_strand=True)
            # matural transcript (mt)
            # exon ?????????????????????????????????MRNA?????????,???????????????????????????matural transcript
            mt = ''
            exons_list = []
            for e in db.children(t, featuretype='exon', order_by='start'):
                exons_list.append([e.start, e.end])
                s = e.sequence(args.genome, use_strand=False) # ????????????????????????????????????????????????cds???????????????????????????
                mt += s
            mt = Seq(mt)
            if t.strand == '-':
                mt = mt.reverse_complement()
            # CDS ????????????????????? ATG ... TAG
            cds = ''
            for c in db.children(t, featuretype='CDS', order_by='start'):
                s = c.sequence(args.genome, use_strand=False) # ????????????????????????????????????????????????cds???????????????????????????
                cds += s
            if args.style == 'GTF':
                sc_seq = stop_codon_seq(db, t, args.genome)
                cds = add_stop_codon(cds, t.strand, sc_seq)
            cds = Seq(cds)
            if t.strand == '-':
                cds = cds.reverse_complement()
            uORF_dict = uorf_procedure_oriented(t.id, t.chrom, t.strand, mt, cds)
            # loop for type1, type2 and type3
            for key in sorted(uORF_dict.keys()):
                # print(key,len(uORF_dict[key]))
                for it in uORF_dict[key]:
                    # csv output format 
                    if args.output_format == 'csv':
                        uORF_seq.loc[index] = it
                        index += 1
                    elif args.output_format == 'fasta':
                        uorf_id = ".".join(map(str,[it[0], it[4], it[5], it[6]]))
                        seq = it[8] # it[8] is s Seq type
                        desc='strand:%s uORF=%d-%d %s length=%d CDS=%s-%s'%(it[2],
                                                                      it[4],
                                                                      it[5],
                                                                      it[6],
                                                                      len(seq),
                                                                      it[3].split('-')[0],
                                                                      it[3].split('-')[1])
                        uorfRecord = SeqRecord(seq, id=uorf_id.replace('transcript:',''), description=desc)
                        uorf_record.append(uorfRecord)
                        
                    elif args.output_format == 'gff':
                        gff = GFF(exons_list, it)
                        uorf_gff_line,features_gff_lines = gff.parse() 
                        # print(uorf_gff_line)
                        uORF_gff.write(uorf_gff_line+"\n") 
                        for feature_line in features_gff_lines:
                            # print(feature_line)
                            uORF_gff.write(feature_line+"\n")
        # output file 
        if args.output_format == 'csv':
            uORF_seq.to_csv(args.output, sep=',', index=False)
        elif args.output_format == 'gff':
            uORF_gff.close()
        elif args.output_format == 'fasta':
            with open(args.output,'w') as handle:
                SeqIO.write(uorf_record, handle, "fasta")
    else:
        for t in db.features_of_type(mRNA_str, order_by='start'):
            # primary transcript (pt) ???????????????????????????????????????
            # ???????????????intron???????????????????????????matural transcript ?????????
            # id specify
            if args.transcript in t.id:
                pt = t.sequence(args.genome, use_strand=True)
                # matural transcript (mt)
                # exon ?????????????????????????????????MRNA?????????,???????????????????????????matural transcript
                mt = ''
                exons_list = []
                for e in db.children(t, featuretype='exon', order_by='start'):
                    exons_list.append([e.start, e.end]) # 1 based location 
                    s = e.sequence(args.genome, use_strand=False) # ????????????????????????????????????????????????cds???????????????????????????
                    mt += s
                mt = Seq(mt)
                mt_len = len(mt)
                if t.strand == '-':
                    mt = mt.reverse_complement()
                # intron for uORF locate position 
                for i in db.children(t, featuretype='intron', order_by='start'):
                    pass
                # CDS ????????????????????? ATG ... TAG
                cds = ''
                cds_list = []
                for c in db.children(t, featuretype='CDS', order_by='start'):
                    cds_list.append([c.start, c.end])
                    s = c.sequence(args.genome, use_strand=False) # ????????????????????????????????????????????????cds???????????????????????????
                    # print('cds len i:', len(c))
                    cds += s
                if args.style == 'GTF':
                    sc_seq = stop_codon_seq(db, t, args.genome)
                    cds = add_stop_codon(cds, t.strand, sc_seq)
                cds = Seq(cds)
                cds_len = len(cds)
                if t.strand == '-':
                    cds = cds.reverse_complement()
                # procedure-oriented 
                #uORF_dict = uorf_procedure_oriented(t.id, t.chrom, t.strand, mt, cds)
                # objected-oriented
                uORF_ = uORF(t.id, t.chrom, t.strand, mt, cds)
                uORF_dict = uORF_.uorf_parse()
                uorf_location_genome = [] # for schematic on genome # 3D list 
                for key in sorted(uORF_dict.keys()):
                    # print(key,len(uORF_dict[key]))
                    for it in uORF_dict[key]:
                        # csv output format 
                        if args.output_format == 'csv':
                            uORF_seq.loc[index] = it
                            index += 1
                        elif args.output_format == 'fasta':
                            uorf_id = ".".join(map(str,[it[0], it[4], it[5], it[6]]))
                            seq = it[8] # it[8] is a Seq type
                            desc='strand:%s uORF=%d-%d %s length=%d CDS=%s-%s'%(it[2],
                                                                        it[4],
                                                                        it[5],
                                                                        it[6],
                                                                        len(seq),
                                                                        it[3].split('-')[0],
                                                                        it[3].split('-')[1])
                            uorfRecord = SeqRecord(seq, id=uorf_id.replace('transcript:',''), description=desc)
                            uorf_record.append(uorfRecord)
                        elif args.output_format == 'gff':
                            gff = GFF(exons_list, it)
                            uorf_gff_line,features_gff_lines = gff.parse()
                            ex_locations = gff.uorf_exons_location()
                            uorf_location_genome.append(ex_locations) # for schematic on genome
                            uORF_gff.write(uorf_gff_line+"\n")
                            for feature_line in features_gff_lines:
                                uORF_gff.write(feature_line+"\n")
                # figure the schametic 
                # without intron 
                if args.schematic_without_intron:
                    figure = visual_transcript(args.schematic_without_intron, 
                                               args.transcript, 
                                               uORF_.transcript_location(),
                                               uORF_.cds_location_transcript(),
                                               uORF_.uorf_location_transcript())
                    figure.draw()
                # with intron
                if args.schematic_with_intron:
                    figure = visual_transcript_genome(t.strand, 
                                                      args.schematic_with_intron, 
                                                      args.transcript, 
                                                      exons_list, 
                                                      cds_list,
                                                      uorf_location_genome,# 3D list 
                                                      )
                    figure.draw() 
                # file out
                if args.output_format == 'csv':
                    uORF_seq.to_csv(args.output, sep=',', index=False)
                elif args.output_format == 'gff':
                    uORF_gff.close()
                elif args.output_format == 'fasta':
                    with open(args.output,'w') as handle:
                        SeqIO.write(uorf_record, handle, "fasta")
                break 
