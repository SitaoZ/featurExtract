#coding:utf-8
import os
import sys
import argparse

"""
https://genome.ucsc.edu/FAQ/FAQformat.html#format5
The multiple alignment format stores a series of multiple alignments in a format that is easy to parse and relatively easy to read
MAF格式通常用于记录体细胞突变。

MAF本来并不是一个常见的文本文件格式，只是因为癌症研究实在是太热门了，对它的理解也变得需求旺盛起来了。

一、MAF的说明

这些文件应该使用下面描述的突变注释格式（MAF）进行格式化。另外下文中有文件命名规范。

以下几种类型的体细胞突变会在MAF文件中出现：

错义突变及无义突变
剪接位点，其定义为剪接位点2 bp以内的SNP
沉默突变
与基因的编码区、剪接位点或遗传元件目标区域重叠的引物。
移码突变
*调控区突变
大部分MAF提交提交的是原始数据。这些原始数据中在体细胞中标记的位点与已知的变异类型相重合的。为避免有可能出现的细胞系污染，MAF规定了一定的下细胞过滤标准。根据现行政策，可开放获取MAF资料应满足：

包括所有已验证的体细胞突变名称
包括与编码区域或剪接位点重叠的所有未验证的体细胞突变名称
排除所有其他类型的突变（即非体细胞突变、不在编码区域或剪接位点的未验证体的细胞突变以及未在dbSNP、COSMIC或OMIM中注释为体细胞的dbSNP位点）
我们提交给DCC MAF存档的数据包括两种：Somatic MAF（named .somatic.maf）的开放访问数据以及不经过筛选的包含原始数据的Protected MAF（named.protected.maf）。所有数据将使用MAF标准进行格式化。
MAF文件可能有两种格式 ，可能是47列，或者120列，第一行一般都是 头文件，注释着每一列的信息，的确，信息量有点略大。如下：
     1  Hugo_Symbol
     2  Entrez_Gene_Id
     3  Center
     4  NCBI_Build
     5  Chromosome
     6  Start_Position
     7  End_Position
     8  Strand
     9  Consequence
    10  Variant_Classification
    11  Variant_Type
    12  Reference_Allele
    13  Tumor_Seq_Allele1
    14  Tumor_Seq_Allele2
    15  dbSNP_RS
    16  dbSNP_Val_Status
    17  Tumor_Sample_Barcode
    18  Matched_Norm_Sample_Barcode
    19  Match_Norm_Seq_Allele1
    20  Match_Norm_Seq_Allele2
    21  Tumor_Validation_Allele1
    22  Tumor_Validation_Allele2
    23  Match_Norm_Validation_Allele1
    24  Match_Norm_Validation_Allele2
    25  Verification_Status
    26  Validation_Status
    27  Mutation_Status
    28  Sequencing_Phase
    29  Sequence_Source
    30  Validation_Method
    31  Score
    32  BAM_File
    33  Sequencer
    34  t_ref_count
    35  t_alt_count
    36  n_ref_count
    37  n_alt_count
    38  HGVSc
    39  HGVSp
    40  HGVSp_Short
    41  Transcript_ID
    42  RefSeq
    43  Protein_position
    44  Codons
    45  Hotspot
    46  cDNA_change
    47  Amino_Acid_Change
     1  Hugo_Symbol
     2  Entrez_Gene_Id
     3  Center
     4  NCBI_Build
     5  Chromosome
     6  Start_Position
     7  End_Position
     8  Strand
     9  Variant_Classification
    10  Variant_Type
    11  Reference_Allele
    12  Tumor_Seq_Allele1
    13  Tumor_Seq_Allele2
    14  dbSNP_RS
    15  dbSNP_Val_Status
    16  Tumor_Sample_Barcode
    17  Matched_Norm_Sample_Barcode
    18  Match_Norm_Seq_Allele1
    19  Match_Norm_Seq_Allele2
    20  Tumor_Validation_Allele1
    21  Tumor_Validation_Allele2
    22  Match_Norm_Validation_Allele1
    23  Match_Norm_Validation_Allele2
    24  Verification_Status
    25  Validation_Status
    26  Mutation_Status
    27  Sequencing_Phase
    28  Sequence_Source
    29  Validation_Method
    30  Score
    31  BAM_File
    32  Sequencer
    33  Tumor_Sample_UUID
    34  Matched_Norm_Sample_UUID
    35  HGVSc
    36  HGVSp
    37  HGVSp_Short
    38  Transcript_ID
    39  Exon_Number
    40  t_depth
    41  t_ref_count
    42  t_alt_count
    43  n_depth
    44  n_ref_count
    45  n_alt_count
    46  all_effects
    47  Allele
    48  Gene
    49  Feature
    50  Feature_type
    51  One_Consequence
    52  Consequence
    53  cDNA_position
    54  CDS_position
    55  Protein_position
    56  Amino_acids
    57  Codons
    58  Existing_variation
    59  ALLELE_NUM
    60  DISTANCE
    61  TRANSCRIPT_STRAND
    62  SYMBOL
    63  SYMBOL_SOURCE
    64  HGNC_ID
    65  BIOTYPE
    66  CANONICAL
    67  CCDS
    68  ENSP
    69  SWISSPROT
    70  TREMBL
    71  UNIPARC
    72  RefSeq
    73  SIFT
    74  PolyPhen
    75  EXON
    76  INTRON
    77  DOMAINS
    78  GMAF
    79  AFR_MAF
    80  AMR_MAF
    81  ASN_MAF
    82  EAS_MAF
    83  EUR_MAF
    84  SAS_MAF
    85  AA_MAF
    86  EA_MAF
    87  CLIN_SIG
    88  SOMATIC
    89  PUBMED
    90  MOTIF_NAME
    91  MOTIF_POS
    92  HIGH_INF_POS
    93  MOTIF_SCORE_CHANGE
    94  IMPACT
    95  PICK
    96  VARIANT_CLASS
    97  TSL
    98  HGVS_OFFSET
    99  PHENO
   100  MINIMISED
   101  ExAC_AF
   102  ExAC_AF_Adj
   103  ExAC_AF_AFR
   104  ExAC_AF_AMR
   105  ExAC_AF_EAS
   106  ExAC_AF_FIN
   107  ExAC_AF_NFE
   108  ExAC_AF_OTH
   109  ExAC_AF_SAS
   110  GENE_PHENO
   111  FILTER
   112  CONTEXT
   113  src_vcf_id
   114  tumor_bam_uuid
   115  normal_bam_uuid
   116  case_id
   117  GDC_FILTER
   118  COSMIC
   119  MC3_Overlap
   120  GDC_Validation_Status
"""
class maf(object):
    """ a class for maf """
    def __init__(self):
        pass
