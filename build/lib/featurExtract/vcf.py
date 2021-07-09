import re,os
import sys
import gzip
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import pandas as pd

"""
Author: Zhu Sitao
Date : 2018-3-21
Dest: a DNA fasta class; used in python 3.6
VCF is a text file format (most likely stored in a compressed manner). 
It contains meta-information lines, a header line, and then data lines 
each containing information about a position in the genome.
http://www.internationalgenome.org/wiki/Analysis/vcf4.0
http://samtools.github.io/hts-specs/VCFv4.2.pdf
"""

class Header(object):
	""" a class for vcf header """
	def __init__(self,vcfpath):
		self.path = vcfpath
		self.info_pattern = re.compile(r'''\#\#INFO=<
		            		ID=(?P<id>[^,]+),\s*
		            		Number=(?P<number>-?\d+|\.|[AGR])?,\s*
		            		Type=(?P<type>Integer|Float|Flag|Character|String),\s*
		            		Description="(?P<desc>[^"]*)"
		            		(?:,\s*Source="(?P<source>[^"]*)")?
		            		(?:,\s*Version="?(?P<version>[^"]*)"?)?
		            		>''', re.VERBOSE)
		self.filter_pattern = re.compile(r'''\#\#FILTER=<
		            		ID=(?P<id>[^,]+),\s*
		            		Description="(?P<desc>[^"]*)"
		            		>''', re.VERBOSE)
		self.alt_pattern = re.compile(r'''\#\#ALT=<
		            		ID=(?P<id>[^,]+),\s*
		            		Description="(?P<desc>[^"]*)"
		            		>''', re.VERBOSE)
		self.format_pattern = re.compile(r'''\#\#FORMAT=<
		  			ID=(?P<id>.+),\s*
		            		Number=(?P<number>-?\d+|\.|[AGR]),\s*
		            		Type=(?P<type>.+),\s*
		            		Description="(?P<desc>.*)"
		            		>''', re.VERBOSE)
		self.contig_pattern = re.compile(r'''\#\#contig=<
		            		ID=(?P<id>[^>,]+)
		            		(,.*length=(?P<length>-?\d+))?
		            		.*
		            		>''', re.VERBOSE)
		self.head_pattern = re.compile(r'\#CHROM\s*POS', re.VERBOSE)
		self.meta_pattern = re.compile(r'''##(?P<key>.+?)=(?P<val>.+)''')
		self.referPatter = re.compile(r'##reference=file:\/\/(\S+)')
		self.samples = self.getSampleID()
	def readHeader(self):
		""" read header of vcf"""
		with open(self.path,'r') as F:
			for line in F:
				if line.startswith('#') and len(line) !=0:
					yield line
				else:
					break
	def getChromLength(self):
		""" Return a chromosome length dict """
		chromLenDict = dict()
		for line in self.readHeader():
			match = self.contig_pattern.match(line)
			if match:
				chromLenDict[match.group('id')] = match.group('length')
		return chromLenDict

	def getReferPath(self):
		""" Return reference genome path """
		for line in self.readHeader():
			match = self.referPatter.search(line)
			if match:
				return match.group(1)


	def getReferLength(self):
		""" Calculate reference genome length """
		total_length = 0
		chromLenDict = self.getChromLength()
		for key in chromLenDict.keys():
			total_length += int(chromLenDict[key])
		return total_length

	def getSampleID(self):
		""" Return population vcf file sample ids """
		sampleList = []
		for line in self.readHeader():
			match = self.head_pattern.match(line)
			if match:
				temp = line.strip().split()
				for sample in range(9,len(temp)):
					sampleList.append(temp[sample])
		return sampleList


	def header_to_file(self,out_path):
		""" write a header to a file """
		F = open(out_path,'w')
		for line in self.readHeader():
			F.write(line)
		F.close()

class Recorder(object):
	""" a recorder for vcf """
	def __init__(self, CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO,FORMAT,CALLS):
		#: A ``str`` with the chromosome name
		self.CHROM = CHROM
		#: An ``int`` with a 1-based begin position
		self.POS = POS
		#: An ``int`` with a 0-based begin position
		self.begin = eval(POS) - 1
		#: An ``int`` with a 0-based end position
		self.end = None  # XXX
		#: A list of the semicolon-separated values of the ID column
		#self.ID = list(ID)
		self.ID = ID
		#: A ``str`` with the REF value
		self.REF = REF
		#: A list of alternative allele records of type :py:class:`AltRecord`
		#self.ALT = list(ALT)
		self.ALT = ALT
		#: The quality value, can be ``None``
		self.QUAL = QUAL
		#: A list of strings for the FILTER column
		self.FILTER = FILTER
		#: An OrderedDict giving the values of the INFO column, flags are
		#: mapped to ``True``
		self.INFO = INFO
		#: A list of strings for the FORMAT column.  Optional, must be given if
		#: and only if ``calls`` is also given.
		self.FORMAT = FORMAT or ()
		#: A list of genotype :py:class:`Call` objects.  Optional, must be given if
		#: and only if ``FORMAT`` is also given.
		self.CALLS = CALLS or ()


	def is_snv(self):
		""" Return ture is a snv"""
		return len(self.REF) == 1 and all(a.type == "SNV" for a in self.ALT)
	def record_to_line(self):
		line = [self.CHROM,self.POS,self.ID,self.REF,self.ALT,self.QUAL,self.FILTER,self.INFO,self.FORMAT]
		for call in self.CALLS:
			line.append(call)
		return "\t".join(line)

	def filter_by_quality(self, cutoff=30):
		""" filter vcf by information """
		if self.QUAL > cutoff:
			return self.record_to_line()

	def filter_by_depth(self,down,up):
		""" filter vcf by depth <4 or >200"""
		info = Info(self.INFO)
		if info.DP <= up or info.DP >= down:
			return self.record_to_line()


class AltRecord(object):
	""" A class  for alternative allele record """
	def __init__(self,allele=None):
		self.ale = allele
		self.ale_list= self.allele_list()
	def allele_list(self):
		""" return all allele item include the ref,
		which corresponding to the genotype of each sample call"""
		if ',' in self.ale:
			return self.ale.split(',')
		else:
			return self.ale



class Call(object):
	""" a class for single call """
	def __init__(self, sample, data):
		#: the name of the sample for which the call was made
		self.sample = sample
		#: an OrderedDict with the key/value pair information from the
		#: call's data
		self.data = data
		#: the allele numbers (0, 1, ...) in this calls or None for no-call
		self.gt_alleles = self.parser_call()
	def parser_call(self):
		genotype_index = self.data.split(":")[0].split("/")
		return genotype_index

class Info(object):
	""" a class for vcf INFO column """
	def __init__(self,data):
		# An integer,Allele count in genotypes, for each ALT allele, in the same order as listed
		self.data = data
		self.DP =self.depth()
		self.AC = self.allele_count()
		self.AF = self.allele_frequency()
	def depth(self):
		""" return a Combined depth across samples"""
		DP_pattern = re.compile(r'DP=(?P<value>\d+);')
		DP_value = DP_pattern.search(self.data).group('value')
		return DP_value
	def allele_count(self):
		""" return a allele count """
		AC_pattern = re.compile(r'AC=(?P<value>\d+);')
		AC_value = AC_pattern.search(self.data).group('value')
		return AC_value
	def allele_frequency(self):
		""" return an allele frequency """
		AF_pattern = re.compile(r'AF=(?P<value>-?\d+\.?\d*e?-?\d*)?;')
		AF_value = AF_pattern.search(self.data).group('value')
		return AF_value








class VCF(object):
	""" VCF class """
	def __init__(self,vcfpath):
		self.path = vcfpath
		self.basename = os.path.basename(vcfpath)

	def readVCF(self):
		""" Return a vcf file line by line """
		if not os.path.exists(self.path):
			sys.exit("%s file not exist\n"%(self.path))
		if self.path.endswith(".gz"):
			F = gzip.open(self.path)
		else:
			F = open(self.path,'r')
		for line in F:
			if not line.startswith('#'):
				array = line.strip().split()
				CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT= array[0:9]
				CALLS = array[9:len(array)] # a list
				recorder = Recorder(CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,CALLS)
				yield recorder
		F.close()

	def __len__(self):
		""" Return numbers of items in self """
		totalLength = 0
		for i in self.readVCF():
			totalLength += 1
		return totalLength

	def Homo_Hete_Sites_Stat(self):
		""" Stat homozygosis and heterozygosis variation sites"""
		All_samples = dict()
		sampleids = Header(self.path).samples
		for sample in sampleids:
			All_samples[sample] = [0,0,0,0] # HomozygousSitesStat, heterozygousSitesStat, missingSitesStat, nonSitesStat
		for recorder in self.readVCF():
			alleles = AltRecord(allele=recorder.ALT).ale_list
			for i in range(len(sampleids)):
				GT = Call(sampleids[i],recorder.CALLS[i]).gt_alleles
				if GT[0] == "." and GT[1] == ".":
					All_samples[sampleids[i]][3] += 1
					#nonSitesStat += 1
				elif "*" in GT[0] or "*" in GT[1]:
					All_samples[sampleids[i]][2] += 1
					#missingSitesStat += 1
				else:
					if GT[0] == GT[1]:
						All_samples[sampleids[i]][0] += 1
						#HomozygousSitesStat += 1
					else:
						All_samples[sampleids[i]][1] += 1
						#heterozygousSitesStat += 1
		#Table header
		print ("Type\t"+"\t".join(sampleids))
		Type_sites = ["Homozygous   Sites", "Heterozygous Sites", "MissingSites(contain *)", "NonSites(contain .)"]
		for type in Type_sites:
			output = [type]
			for sample in sampleids:
				output.append(str(All_samples[sample][Type_sites.index(type)]))
			print ("\t".join(output))
			output = []

	def getHomoSites(self,HomoSites_outpath):
		""" Return homozygosis sites , including header line and variation line"""
		HomoSites =open(HomoSites_outpath,'w')
		# header
		header = Header()
		for line in header.readHeader(self.path):
			HomoSites.write(line)
		# record
		for recorder in self.readVCF():
			genotype = recorder.CALLS[0]
			GT = genotype.split(":")[0].split("/")
			## homo
			if GT[0] == "." and GT[1] == ".":
				continue
			elif "*" in GT[0] or "*" in GT[1]:
				continue
			else:
				if GT[0] == GT[1]:
					line = recorder.record_to_line()
					HomoSites.write(line+"\n")
		HomoSites.close()


	def getHeteSites(self,HeteSites_outpath):
		"""Return heterozygosis sites , including header and variation line """
		HeteSites = open(HeteSites_outpath, 'w')
		header = Header()
		for line in header.readHeader(self.path):
			HeteSites.write(line)
		for recorder in self.readVCF():
			genotype = recorder.CALLS[0]
			GT = genotype.split(":")[0].split("/")
			## homo
			if GT[0] == "." and GT[1] == ".":
				continue
			elif "*" in GT[0] or "*" in GT[1]:
				continue
			else:
				if GT[0] != GT[1]:
					line = recorder.record_to_line()
					HeteSites.write(line+"\n")
		HeteSites.close()

	def reduceNonSites(self):
		pass

	def extract_snp(self,snp_path):
		""" extract snp into a new vcf """
		F = open(snp_path,'w')
		for recorder in self.readVCF():
			if len(recorder.REF) == len(recorder.ALT):
				F.write(recorder.record_to_line())
		F.close()



	def extract_indel(self,indel_path):
		""" extract indel into a new vcf """
		F =open(indel_path,'w')
		for recorder in self.readVCF():
			if ',' in recorder.ALT:
				if len(recorder.REF) != len(recorder.ALT.split(',')[0]) :
					F.write(recorder.record_to_line())
			else:
				if len(recorder.REF) != len(recorder.ALT):
					F.write(recorder.record_to_line())
		F.close()

	def stat_indel_length(self):
		""" return a indel length dict for draw indel length histogram """
		indel_length = {"delete":{},"insert":{}}

		for recorder in self.readVCF():
			if len(recorder.REF) == 1 and len(recorder.ALT) == 1:
				continue
			else:
				if len(recorder.ALT) - len(recorder.REF) > 0: # insert
					delete_lengthinsert_length = len(recorder.ALT) - len(recorder.REF)
					if insert_length in insert_length["insert"].keys():
						indel_length["insert"][length] += 1
					else:
						indel_length["insert"][length] = 1
				else:
					delete_length = len(recorder.REF) - len(recorder.ALT)
					if delete_length in delete_length["delete"].keys():
						indel_length["delete"][length] += 1
					else:
						indel_length["delete"][length] = 1
		return indel_length
	def indel_length_plot(self):
		insert_df = pd.DataFrame(self.indel_length_plot()['insert'].items())
		insert_df.columns = ['length','count']
		delete_df = pd.DataFrame(self.indel_length_plot()['delete'].items())
		base_plot = (ggplot(insert_df,aes(x='length',y='count')))+geom_bar(stat='identity')
		




	def vcf2genotype(self,genotype_outpath):
		""" Transformate population vcf to genotype """
		IUPAC = {'-' : '-', 'AA' : 'A', 'TT' : 'T', 'CC' : 'C', 'GG' : 'G', 'AG' : 'R',
				'GA' : 'R', 'CT' : 'Y', 'TC' : 'Y', 'AC' : 'M', 'CA' : 'M', 'GT' : 'K',
				'TG' : 'K', 'GC' : 'S', 'CG' : 'S', 'AT' : 'W', 'TA' : 'W',
				'ATC' : 'H', 'ACT' : 'H', 'TAC' : 'H', 'TCA' : 'H', 'CAT' : 'H', 'CTA' : 'H',
				'GTC' : 'B', 'GCT' : 'B', 'CGT' : 'B', 'CTG' : 'B', 'TCG' : 'B', 'TGC' : 'B',
				'GAC' : 'V', 'GCA' : 'V', 'CGA' : 'V', 'CAG' : 'V', 'AGC' : 'V', 'ACG' : 'V',
				'GAT' : 'D', 'GTA' : 'D', 'AGT' : 'D', 'ATG' : 'D', 'TGA' : 'D', 'TAG' : 'D', 
				'ATCG' : 'N', 'ATGC' : 'N', 'ACTG' : 'N', 'ACGT' : 'N', 'AGTC' : 'N', 'AGCT' : 'N', 'TACG' : 'N',
				'TAGC' : 'N', 'TCAG' : 'N', 'TCGA' : 'N', 'TGAC' : 'N', 'TGCA' : 'N', 'CATG' : 'N', 'CAGT' : 'N',
				'CTAG' : 'N', 'CTGA' : 'N', 'CGAT' : 'N', 'CGTA' : 'N', 'GATC' : 'N', 'GACT' : 'N', 'GTAC' : 'N',
				'GTCA' : 'N', 'GCAT' : 'N', 'GCTA' : 'N'}
				#this hash was use for converting the genotype to a single letter according to th IUPAC.
		OUT = open(genotype_outpath,"w")
		header = "Chr\tPosi\tRef"
		for sample in self.getSampleID():
			sampleTab = "\t"+sample
			header += sampleTab
		OUT.write(header+"\n")

		for line in self.readVCF():
			temp,alts = [],[]
			if len(line) != 0 and not line.startswith("#"):
				temp = line.strip().split("\t")
				OUT.write(temp[0]+"\t"+temp[1]+"\t"+temp[3]+"\t") #output the chromosome number, positions and the allele of reference.
				if "," in temp[4]:
					alts = temp[4].split(',')
				else:
					alts.append(temp[4])
				alts.insert(0,temp[3])
				for num in range(9,len(temp)):
					geno = ''
					if re.compile(r'\.\/\.').match(temp[num]):
						geno = "-"
					else:
						aa = re.compile(r'(\d)\/(\d).*').match(temp[num])
						allele1 = alts[int(aa.group(1))]
						allele2 = alts[int(aa.group(2))]
						if "*" in allele1 or "*" in allele2:
							geno = "-"
						else:
							geno = allele1+allele2
					OUT.write(IUPAC[geno]+"\t")
				OUT.write("\n")
		OUT.close()


	def ts_tv(self):
		""" transitions and transversions for snp """
		ts,tv = 0,0
		ts_flag = ["AG","GA","CT","TC"]
		tv_flag = ["AC","CA","GT","TG","GC","CG","AT","TA"]
		for recorder in self.readVCF():
			ref = recorder.REF
			ale = recorder.ALT
			if ref.upper()+ale.split(',')[0].upper() in ts_flag:
				ts += 1
			else:
				tv += 1
		return ts,tv,'%.6f'%(ts/float(tv))

	def compare_vcf(self,other):
		list1 = []
		list2 = []
		for record in self.readVCF():
			chrom,pos,ref,ale = record.CHROM,record.POS,record.REF,record.ALT
			list1.append(chrom+pos+ref+ale)
		for record in other.readVCF():
			chrom, pos, ref, ale = record.CHROM, record.POS, record.REF, record.ALT
			list2.append(chrom + pos + ref + ale)
		set1 = set(list1)
		set2 = set(list2)
		plt.figure(figsize=(4,4))
		venn2(subsets = [set1,set2],set_labels=(self.basename,other.basename),set_colors=('r','g'))
		plt.savefig("{t1}_vs_{t2}.png".format(t1=self.basename,t2=other.basename))





