# -*- coding: utf-8 -*-
import sys

def sub_usage(args):
    if len(args) == 1:
        print("          \033[1;33;40m\nUsage  :\033[1m\033[1;35;40m  %s\033[1m" %s (args[0]) )
    elif len(args) == 2:
        print("          \033[1;35;40m%s\033[0m    \033[1;32;40m%s\033[0m" % (args[0], args[1]) )

def main_usage():
    print("\n\033[1;33;40mProgram: \033[0m\033[1;35;40m featurExtract \033[1;31;40m(pipeline for genome feature extract)\033[0m")
    print("\033[1;33;40mVersion: \033[0m\033[1;32;40m 0.1.5\033[0m")
    print("\033[1;33;40mContact: \033[0m\033[1;32;40m Sitao Zhu <zhusitao1990@163.com>\033[0m")
    print("\033[1;33;40mUsage  : \033[0m\033[1;35;40m featurExtract\033[0m \033[1;31;40m<command>\033[0m")
    print("\033[1;33;40mCommand: \033[0m")
    sub_usage(["create  ", "create database for GTF or GFF"])
    sub_usage(["promoter", "extract promoter for genome or gene"])
    sub_usage(["UTR     ", "extract UTR for genome or gene"])
    sub_usage(["uORF    ", "extract uORF for genome or gene"])
    sub_usage(["CDS     ", "extract CDS for genome or gene"])
    sub_usage(["dORF    ", "extract dORF for genome or gene"])
    sub_usage(["exon    ", "extract exon for transcript"])
    sub_usage(["intron  ", "extract intron for transcript"])
    
    sys.exit(1)


def main():
    if len(sys.argv) == 1:
        main_usage()
    elif len(sys.argv) >= 2:
        if sys.argv[1] in ['create','promoter','UTR','uORF','CDS','dORF','exon','intron']:
            # import 就执行feature_extract()
            # 安装后，系统存在featurExtract包
            from featurExtract import feature_extract
            #import feature_extract
        else:
            main_usage()

if __name__ == "__main__":
    main()
