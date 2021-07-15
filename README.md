# Overview

The featurExtract is python package for bioinformatics. 
The packages contains five subcommands.
The create subcommand is used for creating database.
promoter subcomand is used for extracting promoter sequence.
uORF subcomand is used for extracting upstream open reading frames sequence.
UTR subcomand is used for extracting untranslated region sequence.
CDS subcomand is used for extracting coding sequence.


## Brief introduction of format package

1. **Install** <br>
    ```bash
    pip install featurExtract
    # other
    git clone https://github.com/SitaoZ/featurExtract.git
    cd featurExtract; python setup.py install
    ```

2. **Usage** <br>
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
    
    # GenBank database
    which genBankExtract
    genBankExtract -h
    genBankExtract gene -h
    genBankExtract CDS  -h
    genBankExtract rRNA -h
    genBankExtract tRNA -h
    ```

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
    # GenBank step 3
    genBankExtract gene -g NC_000932.gb -f dna -p  
    genBankExtract CDS  -g NC_000932.gb -f dna -p 
    genBankExtract rRNA -g NC_000932.gb -f dna -p
    genBankExtract tRNA -g NC_000932.gb -f dna -p
    ```
    
