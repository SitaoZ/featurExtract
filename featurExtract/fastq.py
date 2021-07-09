#/usr/bin/python3

import os,sys
import re
import gzip
from collections import defaultdict

"""
Author: Zhu Sitao
Date : 2018-3-21
Dest: a fastq class; used in python 3.6
"""

RANGES = {
    'sanger': (33, 75),
    'illumina-1.8': (33, 79),
    'solexa': (59, 106),
    'phred64': (64, 106),
}

class fastq_record(object):
    """a class for fastq record
    """
    def __init__(self,name,seq,desc,qual):
        self.read_name = name
        self.seq = seq
        self.desc = desc
        self.qual = qual

    def read_id(self):
        return self.read_name
    def seq(self):
        return self.seq
    def desc(self):
        return self.desc
    def qual(self):
        return self.qual
    def fastq_parse(self):
        """ parse to a list """
        return '\n'.join([self.read_name,self.seq,self.desc,self.qual])

class Fastq(object):
    """a class deal fastq file
        new fastq name: @HISEQ:310:C5MH9ANXX:1:1101:3517:2043 2:N:0:TCGGTCAC
        old fastq name: @FCD1PB1ACXX:4:1101:1799:2201#GAAGCACG/2
    """
    def __init__(self,path):
        self.path = path
        self.patternNew = re.compile(r'(?P<id>\d):N:(\d):(?P<index>\S+)',re.VERBOSE)
        self.patternOld = re.compile(r'\@\S+\#(?P<index>\S+?)\/(?P<id>\d)',re.VERBOSE)

    def readFastq(self):
        """ Return fastq file line by line """
        if self.path.endswith("gz"):
            FH = gzip.open(self.path,'r')
        else:
            FH = open(self.path, 'r')
        for line in FH:
            yield line
        FH.close()
    def fastq_to_dict(self):
        """ return dict{fastq_id:}"""
        fqdict = {}
        flag = 0
        record = []
        for line in self.readFastq():
            line = line.strip()
            flag += 1
            if flag == 1:
                record.append(line)  # ID
            if flag == 2:
                record.append(line)  # sequence
            if flag == 3:
                record.append(line)  # description
            if flag == 4:
                record.append(line)  # quality
                fqdict[record[0]] = fastq_record(record[0],record[1],record[2],record[3])
                record.clear() # clear list
                flag = 0       # break
        return fqdict


    def qualitySystem(self):
        """ Return quality system """
        lineCount = 0
        checklines = ''
        for line in self.readFastq():
            lineCount += 1
            if lineCount > 100:
                break
            else:
                if lineCount % 4 == 0:
                    checklines += line.strip()
        c = [ord(x) for x in checklines]
        mi,ma = min(c),max(c)
        for format,v in RANGES.items():
            m1,m2 = v
            if mi >= m1 and ma <= m2:
                return format
        MAX = 0
        MIN = 100
        for one in map(ord,checklines):
            if one > MAX:
                MAX = one
            if one < MIN:
                MIN = one
        if MAX <= 74 and MIN < 59:
            # sanger and Ilumina 1.8+
            return "Phred+33"
        elif MAX > 73 and MIN >=64:
            # Illumina 1.3+ and Illumina 1.5+
            return "Phred+64"
        elif MAX >73 and MIN >= 59 and MIN < 64:
            # Solexa+64
            return "Solexa+64"
        else:
            return "Unknown score encoding"

    def qualsFormat(self):
        """ Return quality system """
        lineCount = 0
        checklines = ''
        for line in self.readFastq():
            lineCount += 1
            if lineCount > 100:
                break
            else:
                if lineCount % 4 == 0:
                    checklines += line.strip()
        c = [ord(x) for x in checklines]
        mi,ma = min(c),max(c)
        for format,v in RANGES.items():
            m1,m2 = v
            if mi >= m1 and ma <= m2:
                return format

    def indexSequence(self):
        """ Return index list """
        indexDict = defaultdict(int)
        indexList = []
        for line in self.readFastq():
            matchNew = self.patternNew.search(line)
            matchOld = self.patternOld.search(line)
            if matchNew:
                indexSequence = matchNew.group("index")
                indexDict[indexSequence] += 1
            if matchOld:
                indexSequence = matchOld.group("index")
                indexDict[indexSequence] += 1
        for key in indexDict.keys():
            indexList.append(key)
        return indexList

    def pairOrSingel(self):
        """ Return the fastq type: single end or paired ends """
        count, singel, pair = 0, 0, 0
        for line in self.readFastq():
            matchNew = self.patternNew.search(line)
            matchOld = self.patternOld.search(line)
            if matchNew:
                count += 1
                if matchNew.group("id") == "1":
                    singel += 1
                else:
                    pair += 1
            if matchOld:
                count += 1
                if matchOld.group("id") == '1':
                    singel += 1
                else:
                    pair += 1

        if count == singel:
            return "singel end"
        if count == pair :
            return "pair end"


    def to_fasta(self,output_fasta):
        """ Transformate fastq into fasta """
        OUT=open(output_fasta,'w')
        flag = 0
        for line in self.readFastq():
            flag += 1
            if flag == 1:
                OUT.writelines(">"+line) ## ID
            if flag == 2:
                OUT.writelines(line) ## sequence
            if flag == 4:
                flag = 0             ## break
        OUT.close()




