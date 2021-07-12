import argparse

class Parameter():
    def __init__(self):
        self._parser = argparse.ArgumentParser()
    def parse(self):
        parser = argparse.ArgumentParser()
        subparsers = parser.add_subparsers(help='sub-command help')
        # create subcommand 
        parser_create = subparsers.add_parser('create', help='create annotation database')
        parser_create.add_argument('-g', '--gff', type=str, help='genome annotation file')
        parser_create.set_defaults(func=create)
        # promoter subcommand
        parser_promoter = subparsers.add_parser('promoter', help='extract promoter in genome or gene')
        parser_promoter.add_argument('-g', '--gene', type=str, help='specific gene; if not given, return whole genes')
        parser_promoter.add_argument('-l', '--length', type=int, help='promoter length before TSS')
        parser_promoter.add_argument('-u', '--utr_head', type=int, help='utr5 length after TSS')
        parser_promoter.add_argument('-f', '--genome', type=str, help='genome fasta')
        parser_promoter.add_argument('-o', '--output', type=str, help = 'output csv file path')
        parser_promoter.add_argument('-p', '--print', action="store_true", help='boolean type, stdin')
        parser_promoter.add_argument('-v', '--version', help = 'promoterExtract version', action = "store_true")
        parser_promoter.set_defaults(func=promoter)
        
        # UTR subcommand
        parser_utr = subparsers.add_parser('UTR', help='extract untranslated region sequence in genome or gene')
        parser_utr.add_argument('-f', '--genome', type=str, help='genome fasta file')
        parser_utr.add_argument('-t', '--transcript', type=str, help='specific transcript id; if not given, \
                               whole transcript will return')
        parser_utr.add_argument('-o', '--output', type=str, help='output file path')
        parser_utr.set_defaults(func=UTR)
        
        # uORF subcommand
        parser_uORF = subparsers.add_parser('uORF', help='extract upper stream open reading sequence in genome or gene')
        parser_uORF.add_argument('-f', '--genome', type=str, help='genome fasta')
        parser_uORF.add_argument('-t', '--transcript', type=str, help='specific transcript id; if not given, \
                               whole transcript will return')
        parser_uORF.add_argument('-o', '--output', type=str, help='output file path')
        parser_uORF.set_defaults(func=uORF)
        
        # CDS subcommand
        parser_cds = subparsers.add_parser('CDS', help='extract coding sequence in genome or gene')
        parser_cds.add_argument('-f', '--genome', type=str, help='genome fasta')
        parser_cds.add_argument('-t', '--transcript', type=str, help='specific transcript id; if not given, \
                               whole transcript will return')
        parser_cds.add_argument('-o', '--output', type=str, help='output file path')
        parser_cds.add_argument('-p', '--print', action="store_true", help='boolean type; stdin')
        parser_cds.set_defaults(func=CDS)
        
        self._parser = parser
        return self._parser

