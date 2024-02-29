# Overview

featurExtract is a Python package for bioinformatics, containing two command programs. The first, featurExtract, includes ten subroutines: create, gene, promoter, UTR, uORF, CDS, dORF, exon, intron and intergenic. The create subroutine is used to create a database, while the promoter subroutine is used to extract promoter sequences. The uORF subroutine extracts upstream open reading frames sequences, and the UTR subroutine extracts untranslated region sequences. The CDS subroutine extracts coding sequences and the intergenic subroutine extracts intergenic sequences between two genes. The second command program, genBankExtract, includes four subroutines: gene, CDS, rRNA and tRNA.


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

featurExtract -h

Program:  featurExtract (tools for genomic feature extract)
Version:  0.2.6.0
Contact:  Sitao Zhu <zhusitao1990@163.com>
Usage  :  featurExtract <command> [parameters] 
Command: 
          create        create GFF/GTF database
          stat          database statistics
          cds           extract CDS sequence
          dorf          extract dORF sequence
          exon          extract exon sequence
          gene          extract gene sequence
          intron        extract intron sequence
          igr           extract intergenic region
          mrna          extract mRNA sequence
          promoter      extract promoter sequence
          terminator    extract terminator sequence
          transcript    extract transcript sequence
          uorf          extract uORF sequence
          utr           extract 5/3UTR sequence
```

-create

```bash
featurExtract create -h

usage: featurExtract create [-h] -g GENOMEFEATURE -o OUTPUT -p PREFIX
                            [-s {gff,gtf}]

optional arguments:
  -h, --help            show this help message and exit
  -g GENOMEFEATURE, --genomefeature GENOMEFEATURE
                        genome annotation file, gff or gtf
  -o OUTPUT, --output OUTPUT
                        database output dir path
  -p PREFIX, --prefix PREFIX
                        database prefix
  -s {gff,gtf}, --style {gff,gtf}
                        genome annotation file format
```

-stat
```bash
featurExtract stat -h

usage: featurExtract stat [-h] -d DATABASE -g GENOME -o OUTPUT [-s {gff,gtf}]

optional arguments:
  -h, --help            show this help message and exit
  -d DATABASE, --database DATABASE
                        database created from creat command
  -g GENOME, --genome GENOME
                        genome fasta path
  -o OUTPUT, --output OUTPUT
                        stat output
  -s {gff,gtf}, --style {gff,gtf}
                        genome annotation file format
```


```
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
    
