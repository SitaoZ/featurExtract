# -*- coding: utf-8 -*-
import os
import sys
import gffutils
import pandas as pd 
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict
from featurExtract.database.database import create_db
from featurExtract.utils.util import add_stop_codon,mRNA_type

'''
https://www.ncbi.nlm.nih.gov/genbank/genomes_gff/
'''
def stop_codon(db, transcript, genome):
    s = ''
    for codon in db.children(transcript, featuretype='stop_codon', order_by='start'):
        s = codon.sequence(genome, use_strand=False) # 不反向互补,序列全部连接后再互补
    return s


def get_cds(args):
    '''
    parameters:
        args: parse from argparse
    return:
        elements write to a file or stdout
    '''
    print(args)
    db = gffutils.FeatureDB(args.database, keep_order=True) # load database
    cds_seq = pd.DataFrame(columns=['TranscriptID','Chrom','Start','End','Strand','CDS'])
    cds_record = []
    index = 0
    # assert GTF or GFF
    mRNA_str = mRNA_type(args.style)
    if not args.transcript:
        for t in db.features_of_type(mRNA_str, order_by='start'):
            seq = ''
            cds_start_transcript = 0
            cds_end_transcript = 0
            for c in db.children(t, featuretype='CDS', order_by='start'):
                # 不反向互补，对于负链要得到全部的cds后再一次性反向>互补
                s = c.sequence(args.genome, use_strand=False)
                seq += s
                if not cds_start_transcript :
                    cds_start_transcript = c.start - t.start
                cds_end_transcript += len(s)
            cds_end_transcript = cds_end_transcript + cds_start_transcript
            
            if args.style == 'GTF':
                stop_codon_seq = stop_codon(db, t, args.genome)
                seq = add_stop_codon(seq, t.strand, stop_codon_seq)
            seq = Seq(seq)
            if t.strand == '-':
                seq = seq.reverse_complement()
            # csv output
            if args.output_format == 'csv':
                cds_seq.loc[index] = [t.id,t.chrom,t.start,t.end,t.strand,seq]
                index += 1
            # defalut fasta
            else:
                desc='strand:%s start:%d end:%d length=%d CDS=%d-%d'%(t.strand,
                                                                      t.start,
                                                                      t.end,
                                                                      len(seq),
                                                                      cds_start_transcript,
                                                                      cds_end_transcript)
                cdsRecord = SeqRecord(seq, id=t.id.replace('transcript:',''), description=desc)
                cds_record.append(cdsRecord)
            # break 
        # csv
        if args.output_format == 'csv':
            cds_seq.to_csv(args.output, sep=',', index=False)
        # default fasta
        elif args.output_format == 'fasta': 
            with open(args.output,'w') as handle:
                SeqIO.write(cds_record, handle, "fasta") 
        else:
            sys.exit('output format not be specified')
    else:
        for t in db.features_of_type(mRNA_str, order_by='start'):
            if args.transcript in t.id:
                seq = ''
                for c in db.children(t, featuretype='CDS', order_by='start'):
                    # 不反向互补，对于负链要得到全部的cds后再一次性>反向互补
                    s = c.sequence(args.genome, use_strand=False)
                    seq += s
                if args.style == 'GTF':
                    stop_codon_seq = stop_codon(db, t, args.genome)
                    seq = add_stop_codon(seq, t.strand, stop_codon_seq)
                seq = Seq(seq)
                if t.strand == '-':
                    seq= seq.reverse_complement()
                cds_seq.loc[index] = [t.id,t.chrom,t.start,t.end,t.strand,seq]
                index += 1
                desc='strand %s start %d end %d length=%d'%(t.strand,t.start,t.end,len(seq))
                cdsRecord = SeqRecord(seq, id=t.id, description=desc)
                if args.print:
                    SeqIO.write([cdsRecord], sys.stdout, "fasta") 
                else:
                    SeqIO.write([cdsRecord], args.output, "fasta") 
                break 


# GenBank 基因组文件，是每条染色体一个LOCUS

def get_cds_gb(args):
    '''
    '''
    cds = []
    proteins = []
    for record in create_db(args.genbank):
        for feature in record.features:
            if feature.type == 'CDS': # CDS promoter UTR 
                cds_seq = ''
                #for part in feature.location.parts:
                #    cds_seq += part.extract(record.seq)
                cds_seq = feature.location.extract(record).seq
                if len(cds_seq)%3 != 0:
                    continue # reject not triple 
                # part.strand 会将FeatureLocation -1的反向互补
                # geneid or locus_tag 
                if 'locus_tag' in feature.qualifiers:
                    gene_id = feature.qualifiers['locus_tag'][0]
                elif 'gene' in feature.qualifiers:
                    gene_id = feature.qualifiers['gene'][0]
                else:
                    gene_id = "Null"
                if 'transl_table' in feature.qualifiers:
                    table_translate = feature.qualifiers['transl_table'][0] 
                else:
                    table_translate = 1
                if str(cds_seq)[:3] in ['AAT','ATA','GTG','TTG','ATT','ACG','TCA','AGG']:
                    # ATA, GTG and TTG (Yokobori et al. 1999). 
                    # ATT is the start codon for the CytB gene
                    # in Halocynthia roretzi (Gissi and Pesole, 2003).
                    pep = Seq(feature.qualifiers['translation'][0])
                else:
                    pep = cds_seq.translate(table=table_translate, cds=True)
                #pep = cds_seq.translate(to_stop=True)
                # 判断提取的和已知的是否一致
                if 'translation' in feature.qualifiers:
                    assert(str(pep) == feature.qualifiers['translation'][0])
                if 'protein_id' in feature.qualifiers:
                    protein_id = feature.qualifiers['protein_id'][0] 
                else:
                    protein_id = 'Null'
                cds_seq_record = SeqRecord(cds_seq,id='gene:%s protein:%s'%(gene_id, protein_id), 
                                 description='strand %s length %d'%(feature.strand, len(cds_seq)))
                pep_record = SeqRecord(pep,id='gene:%s protein:%s'%(gene_id, protein_id),
                             description='strand %s length %d'%(feature.strand, len(pep)))
                cds.append(cds_seq_record)
                proteins.append(pep_record)
                break 
        break     
    if args.print and args.format == 'dna':
        SeqIO.write(cds, sys.stdout, 'fasta')
    elif args.print and args.format == 'protein':
        SeqIO.write(proteins, sys.stdout, 'fasta')
    elif args.output and args.format == 'dna':
        SeqIO.write(cds, args.output, 'fasta')
    elif args.output and args.format == 'protein':
        SeqIO.write(proteins, args.output, 'fasta')
