import re
from fasta import Fasta,fasta_record
"""
Author: Zhu Sitao
Date: 2018-3-21
Dest: The readGFF.py is use to read GFF3,Python version is 3.6
GFF is a standard file format for storing genomic features in a text file. 
GFF stands for Generic Feature Format. GFF files are plain text, 9 column, tab-delimited files. 
GFF databases also exist. They use a schema custom built to represent GFF data.
GFF is frequently used in GMOD for data exchange and representation of genomic data.
http://gmod.org/wiki/GFF3

Fields
Fields must be tab-separated. Also, all but the final field in each feature line must contain a value; "empty" columns should be denoted with a '.'

seqname - name of the chromosome or scaffold; chromosome names can be given with or without the 'chr' prefix. Important note: the seqname must be one used within Ensembl, i.e. a standard chromosome name or an Ensembl identifier such as a scaffold ID, without any additional content such as species or assembly. See the example GFF output below.
source - name of the program that generated this feature, or the data source (database or project name)
feature - feature type name, e.g. Gene, Variation, Similarity
start - Start position of the feature, with sequence numbering starting at 1.
end - End position of the feature, with sequence numbering starting at 1.
score - A floating point value.
strand - defined as + (forward) or - (reverse).
frame - One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on..
attribute - A semicolon-separated list of tag-value pairs, providing additional information about each feature.
"""





class gff_record(object):
	""" a parse for gff line """
	def __init__(self,line):
		self.line = line
		self._attrs = line.strip().split()
	def seqid(self):
		return self._attrs[0]
	def source(self):
		return self._attrs[1]
	def feature(self):
		return self._attrs[2]
	def start(self):
		return self._attrs[3]
	def end(self):
		return self._attrs[4]
	def score(self):
		return self._attrs[5]
	def strand(self):
		return self._attrs[6]
	def frame(self):
		return self._attrs[7]
	def attribute(self):
		return self._attrs[8]
	def region(self):
		return eval(self.end()) - eval(self.start())
	def out(self):
		return self.line





class GFF(object):
	""" Read GFF(general feature format) """
	def __init__(self,gffPath):
		self.gffPath = gffPath
		self.gff = self.all_gff()

		self._class = None
		self._mRNA = dict()
		self._CDS = dict()
		self.version3 = self.gff_version()
		self.geneid_pattern = re.compile(r'geneID=(?P<id>\S+);?')
		self.mRNA_pattern = re.compile(r'ID=(?P<id>\S+);')
		self.cds_pattern = re.compile(r'Parent=(?P<id>\S+);')


	def __iter__(self):
		""" Supports traversal with a for loop"""
		return iter(self._gff)
	def _handle(self):
		""" return a gff handle """
		if 'gz' in self.gffPath:
			fh= gzip.open(self.gffPath, 'rb')
		else:
			fh = open(self.gffPath, 'r')
		return fh

	def readGFF(self):
		""" Return gff file line by line through generator """
		for line in self._handle():
			yield line
	def all_gff(self):
		"""return a dict by rownum"""
		allgff=dict()
		for num,line in enumerate(self.readGFF()):
			if not line.startswith("#"):
				allgff[num]=gff_record(line)
		return allgff

	def gff_version(self):
		""" return gff version """
		#version = int(self._handle().readline().strip().split()[1])
		version = 3
		return True if version == 3 else False

	def annotation_chr_list(self):
		""" Return a chromosome list """
		chr_list = []
		for line in self.readGFF():
			if line.startswith("#"):
				continue
			else:
				chr_list.append(line.strip().split()[0])
		chromosomes = list(set(chr_list))
		chromosomes.sort(key=chr_list.index)
		return chromosomes

	def annotation_source(self):
		""" Return annotation source class list """
		source_list = []
		for line in self.readGFF():
			if line.startswith("#"):
				continue
			else:
				source_list.append(line.strip().split()[1])
		source = list(set(source_list))
		source.sort(key=source_list.index)
		return source

	def annotation_type(self):
		""" Return annotation type list """
		type_list = []
		for line in self.readGFF():
			if line.startswith("#"):
				continue
			else:
				type_list.append(line.strip().split()[2])
		types = list(set(type_list))
		types.sort(key=type_list.index)
		return types


	def gene_id_list(self):
		""" Return gene id list """
		geneid = []
		for line in self.readGFF():
			line = line.strip()
			match = self.geneid_pattern.search(line)
			if match:
				id = match.group("id")
				geneid.append(id)
		geneid_uniq = list(set(geneid))
		geneid_uniq_sort = sorted(geneid_uniq)
		return geneid_uniq_sort
	def geneCount(self):
		""" Return gene total number """
		geneid = self.gene_id_list()
		return len(geneid)

	def mrna_id_list(self):
		""" Return transcript id list """
		mrnaid = []
		for line in self.readGFF():
			line = line.strip()
			match = self.mRNA_pattern.search(line)
			if match:
				id = match.group("id")
				mrnaid.append(id)
		mrnaid_uniq = list(set(mrnaid))
		mrnaid_uniq_sort = sorted(mrnaid_uniq)
		return mrnaid_uniq_sort

	def mrnaCount(self):
		mrnaid = self.mrna_id_list()
		return len(mrnaid)

	def exonCount(self):
		exon = 0
		for line in self.readGFF():
			if line.startswith("#"):
				continue
			else:
				if line.strip().split()[2].lower() == "exon":
					exon += 1
		return exon


	def mRNA_fasta(self,genomefasta_path,mRNAfasta_path):
		genome = Fasta(genomefasta_path)
		genome.readFasta()
		genomeDict = genome._fasta
		forward = {}
		for record in self.all_gff().values():
			start,end = int(record.start()),int(record.end())
			if record.feature() == "mRNA":
				key = self.mRNA_pattern.search(record.attribute()).group("id")
				seq = genomeDict[record.seqid()].get_seq()[start-1:end]
				if record.strand() == "+":
					forward[key] = seq
				else:
					seq = seq.replace('A','{A}').replace('T','{T}').replace('C','{C}').replace('G','{G}')
					seq = seq.format(A='T', T='A', C='G', G='C')[::-1]
					forward[key] = seq

		mRNAfasta = open(mRNAfasta_path,'w')
		for key in forward:
			fa = fasta_record(key,forward[key])
			mRNAfasta.writelines(fa.fasta_parse())
		mRNAfasta.close()


	def cds_fasta(self,genomefasta_path,cds_outpath):
		genome = Fasta(genomefasta_path)
		genome.readFasta()
		genomeDict = genome._fasta
		forward = {}
		reverse = {}
		for record in self.all_gff().values():
			start,end = map(int,[record.start(),record.end()])
			if record.feature() == "mRNA":
				seq2 = ''
				seq3 = ''
			elif record.feature() == "CDS":
				key = self.cds_pattern.search(record.attribute()).group("id")
				out = genomeDict[record.seqid()].get_seq()[start-1:end]
				seq2 += out
				seq3 += out
				if record.strand() == "+":
					forward[key] = seq2
				else:
					reverse[key] = seq3
		for key in reverse.keys():
			seq = reverse[key]
			seq = seq.replace('A','{A}').replace('T','{T}').replace('C','{C}').replace('G','{G}')
			seq = seq.format(A='T', T='A', C='G', G='C')[::-1]
			forward[key] = seq
		CDS = open(cds_outpath,'w')
		for key in sorted(forward.keys()):
			fa = fasta_record(key,forward[key])
			CDS.writelines(fa.fasta_parse())
		CDS.close()

	def sortGFF(self):
		""" a """
		pass

	def gff2gtf(self,gtf_path):
		""" transform gff to gtf """
		GTF=open(gtf_path,'w')
		for line in self._handle():
			# skip comment lines that start with '#'
			if line[0] != "#":
				data = line.strip().split('\t') # split line into columns by tab
				ID = ''
				if data[2] == 'gene':
					ID = data[-1].split('ID=')[-1].split(';')[0]
				else:
					ID = data[-1].split('Parent=')[-1].split(';')[0]
				#modify last column
				data[-1] = 'gene_id "'+ID+'"; transcript_id "' + ID
				GTF.writelines('\t'.join(data))
		GTF.close()


