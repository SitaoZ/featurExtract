# Overview

The featurExtract is python package for bioinformatics.  
The packages contains two executable command programs.  
The first executable program is featurExtract including  
ten subroutines termed create, gene, promoter, UTR, uORF,  
CDS, dORF, exon, intron, intergenic. The create subroutine is   
used for creating database. The promoter subroutine is used  
for extracting promoter sequence. uORF subroutine is used  
for extracting upstream open reading frames sequence. UTR  
subroutine is used for extracting untranslated region sequence.  
CDS subroutine is used for extracting coding sequence.intergenic  
subroutine is used for extracting intergenic sequence between two  
genes. The second executable program is genBankExtract including  
four subroutines termed gene, CDS, rRNA, tRNA.


## Brief introduction of featurExtract package

### Install
Two way offer to install featurExtract module.

#### install command line

```bash
pip install featurExtract
# other
git clone https://github.com/SitaoZ/featurExtract.git
cd featurExtract
python setup.py install
```

#### Requirements

python >= 3.7.6 [python](https://www.python.org/)  
pandas >= 1.2.4 [pandas](https://pandas.pydata.org/docs/)  
gffutils >= 0.10.1 [gffutils](https://pythonhosted.org/gffutils/)  
setuptools >= 49.2.0 [setuptools](https://pypi.org/project/setuptools/)  
biopython >= 1.78 [biopython](https://biopython.org/wiki/Documentation/)  

### Usage
featurExtract is designed for GFF and GTF file  
and GenBankExtract is suited for GenBank file. 

#### featurExtract

```bash
# gff or gtf database 
which featurExtract
featurExtract -h 
featurExtract create -h 
featurExtract promoter -h 
featurExtract UTR -h 
featurExtract uORF -h 
featurExtract CDS -h 
featurExtract dORF -h
featurExtract exon -h
featurExtract intron -h
featurExtract intergenic -h
```

#### genBankExtract

```bash 
# GenBank database
which genBankExtract
genBankExtract -h
genBankExtract gene -h
genBankExtract CDS  -h
genBankExtract rRNA -h
genBankExtract tRNA -h
```
### Examples

#### featurExtract

```bash
# step 1 create database
featurExtract create -f GFF -g ath.gff3 -o ath
# step 2 command
# promoter whole genome
featurExtract promoter -d ath.GFF -f ath.fa -l 200 -u 100 -o promoter.csv --output_format fasta
# promoter one gene to stdout 
featurExtract promoter -d ath.GFF -f ath.fa -l 200 -u 100 -g AT1G01010 -p --output_format fasta
featurExtract UTR -d ath.GFF -f ath.fa -o UTR.csv -s GFF
featurExtract uORF -d ath.GFF -f ath.fa -o uORF.csv -s GFF
featurExtract CDS -d ath.GFF -f ath.fa -o CDS.csv -s GFF
featurExtract mRNA -d ath.GFF -f ath.fa -o mRNA.fasta -s GFF --output_format fasta
featurExtract exon -d ath.GFF -f ath.fa -t AT1G01010.1 -p -s GFF
featurExtract intron -d ath.GFF -f ath.fa -t AT1G01010.1 -p -s GFF
```
    
#### genBankExtract

```bash 
# GenBank step 3
genBankExtract gene -g NC_000932.gb -f dna -p  
genBankExtract CDS  -g NC_000932.gb -f dna -p 
genBankExtract rRNA -g NC_000932.gb -f dna -p
genBankExtract tRNA -g NC_000932.gb -f dna -p
```
    
