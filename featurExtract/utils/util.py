import sys

def stop_codon_seq(db, transcript, genome):
    s = ''
    for codon in db.children(transcript, featuretype='stop_codon', order_by='start'):
        s = codon.sequence(genome, use_strand=False) # 不反向互补,序列全部连接后再互补
    return s

def add_stop_codon(seq, strand, stop_codon_seq):
    '''
    use for GTF seq add stop codon in CDS extracing  
    parameter:
     seq: sequence string 
     strand : ranscript strand 
     stop_codon_seq : stop codong sequence (3 bases)
    return:
     whole CDS
    '''
    if strand == '+':
        seq += stop_codon_seq
    elif strand == '-':
        seq = stop_codon_seq + seq # 负链 stop codon 位置靠前
    else:
        # no clear strand info
        seq += stop_codon_seq
    return seq

def mRNA_type(file_format):
    '''
    parameter:
     file_format: the file format of genome annotation file 
                  for creating database
    return: 
     mRNA type string
    '''
    if file_format == 'GFF':
        return 'mRNA'
        
    elif file_format == 'GTF':
        return 'transcript'
    else:
        sys.stderr.write("parameter -s should be assign \n")
        sys.exit(1)

def utr3_type(file_format):
    '''
    parameter:
     file_format: the file format of genome annotation file 
                  for creating database
    return: 
     utr3 type string
    '''
    if file_format == 'GFF':
        return 'three_prime_UTR'
    elif file_format == 'GTF':
        return 'three_prime_utr'
    else:
        sys.stderr.write("parameter -s should be assign \n")
        sys.exit(1)

def utr5_type(file_format):
    '''
    parameter:
     file_format: the file format of genome annotation file 
                  for creating database
    return:
     utr5 type string
    '''
    if file_format == 'GFF':
        return 'five_prime_UTR'
    elif file_format == 'GTF':
        return 'five_prime_utr'
    else:
        sys.stderr.write("parameter -s should be assign \n")
        sys.exit(1)

def seq_upper_lower(seq,start,end):
    '''
    parameter:
     seq: sequence 
     start: start codon position 
     end: stop codon position 
    return:
     seq with upper and lower bases
    '''
    utr5 = seq[:start].lower()
    coding = seq[start:end].upper()
    utr3 = seq[end:].lower()
    return utr5+coding+utr3
