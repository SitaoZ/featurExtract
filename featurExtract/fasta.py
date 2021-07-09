import re,os
from random import sample
from itertools import groupby

"""
Author: Zhu Sitao
Date : 2018-3-21
Dest: a DNA fasta class; used in python 3.6
"""



BASES = ['T', 'C', 'A', 'G']
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
AMINO_ACIDS = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
CODON_TABLE = dict(zip(CODONS, AMINO_ACIDS))


class fasta_record(object):
	""" a atom of fasta file """
	def __init__(self,id,seq,form="DNA"):
		self.id = id
		self.seq = seq
		self.form = form
	def get_id(self):
		return self.id
	def get_seq(self):
		return self.seq
	def reverse(self,seq):
		"""reverse fasta seq"""
		return seq[::-1]
	def complement(self,seq):
		"""reverse fasta seq"""
		COMPLEMENT_TRANS = str.maketrans('TAGCtagc', 'ATCGATCG')
		return seq.translate(COMPLEMENT_TRANS)
	def reverse_complement(self,seq):
		""" return reverse and complement seq"""
		return self.complement(self.reverse(seq))
	def translate(self,seq):
		seq = seq.lower().replace('\n', '').replace(' ', '')
		peptide = ''
		for i in range(0, len(seq), 3):
			codon = seq[i: i + 3]
			amino_acid = CODON_TABLE.get(codon, '!')
			if amino_acid != '!':  # end of seq
				peptide += amino_acid
		return peptide
	def fasta_parse(self,number=60):
		""" a parse for a fasta """
		seqOut = ''
		for i in range(1, len(self.seq) // number + 1):
			start = (i - 1) * number
			end = i * number
			seqOut += self.seq[start:end] + "\n"
		left = len(self.seq) % number
		remainder = self.seq[-left:] if left != 0 else ""
		totalOut = ">" + self.id + "\n" + seqOut + remainder + "\n"
		return totalOut





class Fasta(object):
	""" Fasta class """

	def __init__(self, filePath):
		self._fasta = dict()
		self.path = filePath
		#self.readFasta()

	def readFasta(self):
		""" Read Fasta file and load in a dict ,normal method
			使用groupby 将文本文件做成一个生成器，生成器没有把所有值存在内存中，而是在运行时生成的值，可以快速访问大文件。
			生成器你只能对其迭代一次。
		"""
		fh = open(self.path)
		faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
		for header in faiter:
			header = header.__next__()[1:].strip()  # [1:] 为了去除 > 符号 header is a class 'itertools._grouper'
			header = header.split()[0]  # header.split()[0] 对名称进行简化，当前的做法是保存全部名称
			seq = "".join(s.strip() for s in faiter.__next__()) # faiter is a class 'generator'
			self._fasta[header] = fasta_record(header, seq)  # self._fasta[header] = seq.upper()
		return self._fasta


	def __iter__(self):
		""" Supports traversal with a for loop ,for ID loop """
		for record in self._fasta.values():
			yield record

	def __getitem__(self, index):
		""" Index one fasta ID """
		outPut = self._fasta[index]
		return outPut


	def __len__(self):
		""" Total sequence length of the fasta file """
		totalLength = 0
		for record in self._fasta.values():
			totalLength += len(record.seq)
		return totalLength

	def fasta_key(self):
		""" return a list of fasta name"""
		return self._fasta.keys()


	def fasta_sequence(self):
		""" return a list of sequence"""
		return self._fasta.values()


	def item_count(self):
		""" return numbers fasta names"""
		item_total_number = len(self.fasta_key())
		return item_total_number


	def std_out(self, number=60):
		""" printf sequence to stdout """
		for record in self.fasta_sequence():
			print(record.fasta_parse(number=number))


	def gc_rate(self, output_txt):
		""" Show each fasta GC and GC rate"""
		Stat = open(output_txt, 'w')
		Stat.writelines("ID\tGC\tGCrate(%)\tN\tNrate(%)\n")
		totalGC, totalN = 0, 0
		for record in self.fasta_sequence():
			seqIn = record.seq
			GC = seqIn.count("G") + seqIn.count("C")
			N = seqIn.count("N")
			totalGC += GC
			totalN += N
			GCrate = GC / len(seqIn)
			Nrate = N / len(seqIn)
			Stat.writelines("%s\t%d\t%.4f\t%d\t%.4f\n" % (ID, GC, GCrate, N, Nrate))
		totalGCrate = totalGC / len(self)
		totalNrate = totalN / len(self)
		Stat.writelines("total\t%d\t%.4f\t%d\t%.4f" % (totalGC, totalGCrate, totalN, totalNrate))
		Stat.close()


	def stat_length(self):
		""" Statistic of each id length """
		# LenDict = defaultdict(int)
		LenDict = dict()
		for ID,record in self._fasta.items():
			LenDict[ID] = len(record.seq)
		return LenDict


	@classmethod
	def basename(cls):
		""" class method """
		return os.path.basename(cls)


	@staticmethod
	def print_working_directory():
		""" static method ,no default parameter """
		return os.getcwd()


	def _max_min(self, type):
		sort_dict = sorted(self.stat_length().items(), key=lambda d: d[1], reverse=True)
		if type == 'max':
			max_item = sort_dict[0]
			return max_item
		elif type == 'min':
			min_item = sort_dict[self.item_count() - 1]
			return min_item
		else:
			return "%s not max or min" % type


	def extract_item(self, type, outfile):
		""" extract max or min length fasta item"""
		if type.lower() not in ['max', 'min']:
			raise KeyError(str(type) + " not contain,should be min or max")
		with open(outfile, 'w') as OUT:
			key, value = self._max_min(type)
			record = self._fasta[key]
			if type == 'max':
				OUT.writelines(record.fasta_parse())
			elif type == 'min':
				OUT.writelines(record.fasta_parse())
			else:
				print("%s not correct,should be max or min" % type)


	def random_sample(self, number, outpath):
		""" ramdom sample from the all fasta file"""
		OUT = open(outpath, 'w')
		key_list = self.fasta_key()
		sample_list = sample(key_list, number)
		for key in sample_list:
			record = self._fasta[key]
			OUT.writelines(record.fasta_parse())
		OUT.close()




