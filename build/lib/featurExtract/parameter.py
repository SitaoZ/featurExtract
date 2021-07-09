import argparse

class Parameter():
    def __init__(self):
        self._parser = argparse.ArgumentParser()
    
    def parse(self):
        subparsers = self._parser.add_subparsers(help='sub-command help')
        # create subcommand 
        parser_create = subparsers.add_parser('create', help='create help')
        parser_create.add_argument('-g', '--gff', type=str, help='genome annotation file')
        parser_create.set_defaults(func='create')
        # extract subcommand
        parser_extract = subparsers.add_parser('extract', help='extract help')
        parser_extract.add_argument('-l', '--length', type=int, help='promoter length before TSS')
        parser_extract.add_argument('-u', '--utr_head', type=int, help='length after TSS')
        parser_extract.add_argument('-f', '--genome', type=str, help='genome fasta')
        parser_extract.add_argument('-g', '--gff', type=str, help='genome annotation file')
        parser_extract.add_argument('-o', '--output', type=str, help = 'output csv file path')
        parser_extract.add_argument('-v', '--version', help = 'promoterExtract version', action = "store_true")
        parser_extract.set_defaults(func='extract')
        #
        args = self._parser.parse_args()
        #args.func(args)

        #if args.version:
        #    print("promorerExtract version 0.9.3")
        #    exit(1)

        #self.length = args.length
        #self.utr_head = args.utr_head
        #self.genome = args.genome
        self.gff = args.gff
        #self.output = args.output
        
        #if self.length is None:
        #    raise  Exception("Error: promoter length is required")
        #elif self.utr_head is None:
        #    raise Exception("Error: utr length is required")
        #elif self.genome is None:
        #    raise Exception("Error: genome fasta is required")
        #elif self.gff is None:
        #    raise Exception("Error: gff is required")
        #elif self.output is None:
        #    raise Exception("Error: output is required")

        return self
