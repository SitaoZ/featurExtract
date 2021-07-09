#coding:utf=8

import os,re

"""
Author: Zhu Sitao
Date: 2019-4-3
Description: The BED format consists of one line per feature, each containing 3-12 columns of data,plus optional track definition lines.
 BED 文件(Browser Extensible Data)格式是ucsc 的genome browser的一个格式 ,提供了一种灵活的方式来定义的数据行，以用来描述注释信息.
 BED文件结构：
-------------------------------------------------------------必须有以下3列------------------------------------------------------------------------
    chrom :即染色体号
    chromStart :即feature在染色体上起始位置 。在染色体上最左端坐标是0
    chromEnd :即feature在染色体上的终止位置。例如一个染色体前100个碱基定义为chromStart=0,  chromEnd=100, 跨度为 0-99.
----------------------------------------------------------------可选9列-------------------------------------------------------------------------------
    name ：feature的名字 ，在基因组浏览器左边显示；
    score ：在基因组浏览器中显示的灰度设定，值介于0-1000；
    strand ：定义链的方向，‘‘+” 或者”-”
    thickStart ：起始位置(例如，基因起始编码位置）
    thickEnd ：终止位置（例如：基因终止编码位置）　
    itemRGB ：是一个RGB值的形式, R, G, B (eg. 255, 0,0), 如果itemRgb设置为‘On”, 这个RBG值将决定数据的显示的颜色。
    blockCount ：BED行中的block数目，也就是外显子数目
    blockSize：用逗号分割的外显子的大小, 这个item的数目对应于BlockCount的数目
    blockStarts ：用逗号分割的列表, 所有外显子的起始位置，数目也与blockCount数目对应
"""

class bedParse(object):
    """ a class for one line of a bed file """
    def __init__(self,chrom,start,end,name=None,score=None,strand=None,thickStart=None,thickEnd=None,itemRgb=None,blockCount=None,blockSize=None,blockStarts=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand
        self.thickStart = thickStart
        self.thickEnd = thickEnd
        self.itemRgb = itemRgb
        self.blockCount = blockCount
        self.blockSize = blockSize
        self.blockStarts = blockStarts

    def length(self):
        """get length of a line """
        return self.end - self.start


class bed(object):
    """ a bed class """
    def __init__(self,bed_path):
        self.bedPath = bed_path
        self.basename = os.path.basename(bed_path)
        self.dirname = os.path.dirname(bed_path)
    def readBed(self):
        """ a method for read bed file """
        for line in open(self.bedPath,'r').readline():
            yield line.strip()



            


