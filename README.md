# Overview

The featurExtract is python package for bioinformatics. 
The packages contains two executable command programs.
The first executable program is featurExtract including 
nine subroutines termed create, gene, promoter, UTR, uORF,
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
    **install command line** <br>
    ```bash
    pip install featurExtract
    # other
    git clone https://github.com/SitaoZ/featurExtract.git
    cd featurExtract; python setup.py install
    ```
python>=3.7.6',
                        'argparse>=1.1', 
                        'pandas>=1.2.4', 
                        'gffutils>=0.10.1',
                        'setuptools>=49.2.0',
                        'biopython>=1.78']

    **Requirements** <br>
    ```
    python >= 3.7.6 [python](https://www.python.org/)
    pandas >= 1.2.4 [pandas](https://pandas.pydata.org/docs/)
    gffutils >= 0.10.1 [gffutils](https://pythonhosted.org/gffutils/)
    setuptools >= 49.2.0 [setuptools](https://pypi.org/project/setuptools/)
    biopython >= 1.78 [biopython](https://biopython.org/wiki/Documentation)
    ```
### Usage
    featurExtract is designed for GFF and GTF file and 
    GenBankExtract is suited for GenBank file. 
    **featurExtract** <br> 
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
    **genBankExtract** <br>
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

    **featurExtract** <br>
    ```bash
    # step 1 
    featurExtract create -g ath.gff3 
    # step 2 command
    featurExtract promoter -l 200 -u 100 -f ath.fa -o promoter.csv
    featurExtract UTR  -o UTR.csv
    featurExtract uORF -o uORF.csv
    featurExtract CDS  -o CDS.csv
    featurExtract exon -f ath.fa -t AT1G01010.1 -p 
    featurExtract intron -f ath.fa -t AT1G01010.1 -p  
    ```
    
    **genBankExtract** <br>
    ```bash 
    # GenBank step 3
    genBankExtract gene -g NC_000932.gb -f dna -p  
    genBankExtract CDS  -g NC_000932.gb -f dna -p 
    genBankExtract rRNA -g NC_000932.gb -f dna -p
    genBankExtract tRNA -g NC_000932.gb -f dna -p
    ```
    
