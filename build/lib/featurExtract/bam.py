import re,os

"""
Author :Zhu Sitao
Date: 2018-3-23
Dest: read bam file by samtools ;samtolls should be installed in your ~/.bashrc
      python version is 3.6
"""


class Bam(object):
	""" a readBam class """
	# class global variable , each funcation could use the variable
	patternAnnotation = re.compile(r'^@')

	def __init__(self,filePath):
		self.path = filePath
		self.sam = self.readBam()

	def readBam(self):
		""" Read bam from samtools view """
		samtools = os.popen('which samtools').read().strip()
		sam = os.popen('%s view -h %s' % (samtools,self.path)).read()
		iterm = sam.split("\n")
		itermFilter = [x for x in iterm if x !='']
		return itermFilter
		#for i in iterm:
		#	yield i

	def bamHeader(self):
		""" Return bam header to stdout """
		for line in self.sam:
			if self.patternAnnotation.match(line) and line != '':
				print (line)
	def mappingRate(self):
		""" Return bam mapping rate """
		clean_reads, clean_bases = 0, 0
		mapped_reads, mapped_bases = 0, 0
		mapping_rate, uniq_rate = 0, 0
		dup_reads, dup_rate = 0, 0
		mismatch_bases, mismatch_rate = 0, 0
		uniq_reads, uniq_bases = 0, 0
		mappedPattern = re.compile(r'XC:i:(\d+)')
		uniqPattern = re.compile(r'X0:i:(\d+)')
		mismatchPattern1 = re.compile(r'XM:i:(\d+)')
		mismatchPattern2 = re.compile(r'MD:Z:(\S+)')

		bamfh = os.popen("%s view %s" % (self.samtools, bam))
		for line in bamfh.readlines():
			tmp = line.strip().split()
			flag = int(tmp[1])
			if (flag & 0x100):
				continue
			else:
				clean_reads += 1
				clean_bases += len(tmp[9])
				if not (flag & 0x4):
					mapped_reads += 1
					if mappedPattern.search(line):
						mapped_bases += int(mappedPattern.search(line).group(1))
					else:
						mapped_bases += len(tmp[9])
					if flag & 0x400 == 0:
						if uniqPattern.search(line):
							if int(uniqPattern.search(line).group(1)) == 1:
								uniq_reads += 1
							uniq_bases += len(tmp[9])
						elif int(tmp[4]) > 0:
							uniq_reads += 1
							uniq_bases += len(tmp[9])
						if mismatchPattern1.search(line):
							mismatch_bases += int(mismatchPattern1.search(line).group(1))
						elif mismatchPattern2.search(line):
							mismatch_bases += self.CoutMismathNo(line)
					else:
						dup_reads += 1
		mismatch_rate = mismatch_bases / mapped_bases
		mapping_rate = mapped_reads / clean_reads
		dup_rate = dup_reads / mapped_reads
		uniq_rate = uniq_reads / mapped_reads
		name = os.path.basename(bam)
		sample = name.split('.')[0]
		print("Sample\t%s\n" % (sample))
		print("Clean reads\t%s\n" % (clean_reads))
		print("Clean bases (Mb)\t%.2f\n" % (clean_bases / 1000000))
		print("Mapping rate (%%)\t%.2f\n" % (100 * mapping_rate))
		print("Unique rate (%%)\t%.2f\n" % (100 * uniq_rate))
		print("Duplicate rate (%%)\t%.2f\n" % (100 * dup_rate))
		print("Mismatch rate (%%)\t%.2f\n" % (100 * mismatch_rate))

	def readLength(self):
		""" Read length """
		for line in self.sam:
			if not self.patternAnnotation.match(line) and line != '':
				readLength = len(line.split('\t')[9])
				return readLength

	def cleanReads(self):
		""" Clean reads """
		cleanReads = 0
		for line in self.sam:
			if line !='' and not self.patternAnnotation.match(line):
				cleanReads += 1
		return cleanReads

	def cleanBases(self):
		""" Clean bases """
		readLength = self.readLength()
		cleanBases = 0
		for line in self.sam:
			if line != '' and not self.patternAnnotation.match(line):
				cleanBases += readLength
		return cleanBases
				
	def depth(self):
		""" Bam deepth """
		pass

	def get_insersize(self):
		pass


	def extract(self,chrName):
		""" extract given ID name bam file """
		for line in self.sam:
			if line != '' and not self.patternAnnotation.match(line):
				chrID = line.split("\t")[2]
				if chrName == chrID:
					print (line)
