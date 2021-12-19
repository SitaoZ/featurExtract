import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle, FancyArrowPatch
import matplotlib as mpl
mpl.rcParams['axes.formatter.useoffset'] = False

class visual_transcript(object):
    '''
    Visualize the uORF in transcript location (relative position) 
    Parameters:
     - transcript_location: a 1D list postion of transcript
     - cds_location: a 1D list CDS postion of transcript
     - uorf_location: a 2D list position of uORFs
    Return:
     - a figure object
    '''
    def __init__(self, output, transcript_id, transcript_location, cds_location, uorf_location):
        self.output = output
        self.t_id = transcript_id
        self.transcript = transcript_location
        self.CDS = cds_location
        self.uORFs = uorf_location
        self.xt = self.transcript[0] # real position 
        self.yt = 10                 # relative position
        self.xc = self.CDS[0]
        self.yc = 10
        self.height = 1
        self.scale = 0.5 
        self.tip = 1
        
    def draw(self):
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        # transcript 
        ax.add_patch(Rectangle((self.xt, self.yt + 0.5*self.scale * self.height), \
                    self.transcript[1] - self.xt, \
                    self.scale * self.height, color='#DA70D6')) # xy, width, height
        plt.text(self.transcript[1]*1.01, \
                 self.yt, self.t_id, \
                 fontsize='small')
        # CDS, CDS overlap with transcript 
        ax.add_patch(Rectangle((self.xc, self.yt), \
                    self.CDS[1] - self.xc, \
                    self.height, color='#FF8C00')) # lightcoral
        plt.text((self.CDS[1] + self.CDS[0])/2, \
                  self.yc + 2, \
                  "CDS", fontsize='small')
        # UTR5
        if self.CDS[0] != self.transcript[0]:
            # 5 primer utr exist 
            plt.text((self.xt + self.xc)/2, \
                      self.yc + 2, \
                      "5' UTR", fontsize='small')
        # UTR3
        if self.CDS[1] != self.transcript[1]:
            # 3 primer utr exist 
            plt.text((self.transcript[1] + self.CDS[1])/2, \
                      self.yc + 2, \
                      "3' UTR", fontsize='small')
        # uORF
        for i,uorf in enumerate(self.uORFs):
            i = i + 1
            xu = uorf[0]
            yu = self.yc - 2 * i * self.height
            ax.add_patch(Rectangle((xu, yu), \
                         uorf[1] - xu, 0.8 * self.height, color='#43CD80'))
            plt.text(uorf[1]*1.01, yu, 'uORF%d'%(i), fontsize='small')
        plt.xlim(self.transcript[0] - 10, self.transcript[1] + 10)
        plt.ylim(-2 * len(self.uORFs),10 + 2 * len(self.uORFs))
        plt.axis('off')
        plt.savefig(f'{self.output}.pdf')
        plt.show()


# a = visual_transcript([[1,100],\
#                        [10,80],\
#                        [[3,9],[6,44]]])
# a.draw()


class visual_transcript_genome(object):
    '''
    visualize a transcript on genome, absolute position
    Parameters:
     - output: output file path, str
     - transcript_id: transcript id, str
     - exon_location: exon location blocks in genome, 2D list 
     - cds_location: CDS location blocks in genome, 2D list 
     - uorf_location: uORF location blocks in genome, 2D list 
    '''
    def __init__(self, strand, output, transcript_id, exon_location, cds_location, uorf_location):
        self.strand = strand
        self.output = output 
        self.transcript_id = transcript_id
        self.exon_location = exon_location # 2D list 
        self.cds_location = cds_location   # 2D list 
        self.uorf_location = uorf_location # 3D list  
        self.height = 1
        self.yt = 10     # relative position
        self.yc = 10
        self.scale = 0.5
        self.tip = 1
        self.transcript_intron_location = self.transcript_intron_locate()
    
    def transcript_intron_locate(self):
        '''
        intron of transcript in genome coordinate, same to genome location 
        '''
        intron_list = []
        intron_start = 0
        intron_end = 0
        # no intron 
        if len(self.exon_location) == 1:
            return None
        # intron exist 
        else:
            for i,exon in enumerate(self.exon_location):
                i = i + 1
                if i == 1 :
                    # first exon 
                    intron_start = exon[1] # genome location should add 1 exon[1]+1; visualization ignore this
                    intron_list.append(intron_start)
                    continue 
                elif i == len(self.exon_location):
                    # last exon 
                    intron_end = exon[0]
                    intron_list.append(intron_end)
                else:
                    # inter exon
                    intron_list.append(exon[0])
                    intron_list.append(exon[1])
            array = np.array(intron_list)
            intron_location_2d = array.reshape(-1, 2)
            #if self.strand == '+':
            #    genome_intron_location_2d = intron_location_2d
            #else:
            #    genome_intron_location_2d = intron_location_2d[::-1]
            return intron_location_2d

    def cds_intron_location(self):
        '''
        intron of cds in genome coordinate, same to genome location 
        '''
        intron_list = []
        intron_start = 0
        intron_end = 0
        # no intron 
        if len(self.cds_location) == 1:
            return None
        # intron exist 
        else:
            for i,exon in enumerate(self.cds_location):
                i = i + 1
                if i == 1 :
                    # first exon 
                    intron_start = exon[1] # genome location should add 1 exon[1]+1; visualization ignore this
                    intron_list.append(intron_start)
                    continue
                elif i == len(self.cds_location):
                    # last exon 
                    intron_end = exon[0]
                    intron_list.append(intron_end)
                else:
                    # inter exon
                    intron_list.append(exon[0])
                    intron_list.append(exon[1])
            array = np.array(intron_list)
            intron_location_2d = array.reshape(-1, 2)
            #if self.strand == '+':
            #    genome_intron_location_2d = intron_location_2d
            #else:
            #    genome_intron_location_2d = intron_location_2d[::-1]
            return intron_location_2d 
    
    def uorf_intron_location(self):
        '''
         uorf intron location  
        '''
        intron_list_all = []
        intron_start = 0
        intron_end = 0
        # 2D list 
        for uorf in self.uorf_location:
            # wotyout intron 
            if len(uorf) == 1:
                intron_list_all.append([])
            else:
                # with intron 
                if self.strand == '-':
                    uorf = uorf[::-1]
                intron_list = []
                for i,exon in enumerate(uorf):
                    i = i + 1
                    if i == 1 :
                        # first exon 
                        intron_start = exon[1] # genome location should add 1 exon[1]+1; visualization ignore this
                        intron_list.append(intron_start)
                        continue
                    elif i == len(uorf):
                        # last exon 
                        intron_end = exon[0]
                        intron_list.append(intron_end)
                    else:
                        # inter exon
                        intron_list.append(exon[0])
                        intron_list.append(exon[1])
                array = np.array(intron_list)
                intron_location_2d = array.reshape(-1, 2)
                intron_list_all.append(intron_location_2d)    
                #if self.strand == '+':
                #    genome_intron_location_2d = intron_location_2d
                #else:
                #    genome_intron_location_2d = intron_location_2d[::-1]
        return intron_list_all
        
    def get_arrow_ab(self, exon):
        '''
        get arrow posA and posB 
        '''
        colors = {'transcript': '#FF8C00',
                  'cds':        '#FF8C00',
                  'uorf':       '#DA70D6'}
        if self.strand == '+':
            pos_a = exon[0]
            pos_b = exon[1]
        else:
            pos_a = exon[1]
            pos_b = exon[0]
        return pos_a, pos_b
         
    def draw (self):
        #fig = plt.figure(figsize=(15, 6), dpi=300)
        fig = plt.figure()
        ax = fig.add_axes([0.1, 0.1, 0.6, 0.75])
        # transcript intron line for transcript
        for intron in self.transcript_intron_location:
            if len(intron) != 0 :
                ax.plot([intron[0], intron[1]], \
                        [self.yt + self.scale * self.height, \
                        self.yt + self.scale * self.height], \
                        solid_capstyle='butt', dash_capstyle = 'butt', linewidth=0.7,color='grey')
                # plt.axvline(x=intron[1])
        # transcript exon rectangle 
        for exon in self.exon_location:
            # Rectangle
            #ax.add_patch(Rectangle((exon[0], self.yt), \
            #            exon[1]-exon[0]-1, self.height, color='#FF8C00')) # xy, width(exon length), height
            # FancyArrowPatch
            pos_a, pos_b = self.get_arrow_ab(exon)
            arrow = FancyArrowPatch(
               posA=(pos_a, self.yt + self.scale * self.height), \
               posB=(pos_b, self.yt + self.scale * self.height), fc='#FF8C00', ec='#FF8C00',
               arrowstyle='simple, head_width=5, head_length=5, tail_width=5', capstyle='butt')
            ax.add_artist(arrow)
            # plt.text((exon[0]+exon[1])/2, self.yt, 'exon', fontsize='small')
        # transcript text 
        plt.text(self.exon_location[-1][1] * 1.01, \
                 self.yt + self.scale * self.height, \
                 self.transcript_id, fontsize=5)
        
        # cds intron line, different with transcript
        for intron in self.cds_intron_location():
            if len(intron) != 0 :
                ax.plot([intron[0], intron[1]], \
                        [self.yt - 2 + self.scale * self.height, \
                         self.yt - 2 + self.scale * self.height], \
                         solid_capstyle='butt', linewidth=0.7, color='grey')
        # cds exon rectangle
        for cds in self.cds_location:
            # Rectangle
            #ax.add_patch(Rectangle((cds[0], self.yt - 2), \
            #            cds[1] - cds[0], self.scale * self.height, color='#DA70D6')) # xy, width(exon length), height  
            # FancyArrowPatch
            pos_a, pos_b = self.get_arrow_ab(cds)
            arrow = FancyArrowPatch(
                 posA=(pos_a, self.yt - 2 + self.scale * self.height), \
                 posB=(pos_b, self.yt - 2 + self.scale * self.height), fc='#DA70D6', ec='#DA70D6',
                 arrowstyle='simple, head_width=5, head_length=5, tail_width=5', capstyle='butt')
            ax.add_artist(arrow)
        # cds text 
        plt.text(self.cds_location[-1][1] * 1.01, \
                 self.yt - 2 + self.scale * self.height, \
                 self.transcript_id + '.CDS', fontsize=5)
        
        # uORF
        for i,uorf in enumerate(self.uorf_location):
            i = i + 1
            # uorf intron, different with transcript
            if len(self.uorf_intron_location()[i-1]) != 0:
                for intron in self.uorf_intron_location()[i-1]:
                    # intron 1D list 
                    if len(intron) != 0 :
                        ax.plot([intron[0], intron[1]], \
                            [self.yt - 2 - i + self.scale**2 * self.height , \
                            self.yt - 2 - i  + self.scale**2 * self.height], \
                            solid_capstyle='butt', linewidth=0.7,color='grey')
                    else:
                        continue
            # urof exon rectangle
            for j,u in enumerate(uorf):
                # Rectangle
                ax.add_patch(Rectangle((u[0], self.yt  - 2 - i), \
                           u[1]-u[0], self.scale * self.height, capstyle='butt', color='#43CD80'))
                #pos_a, pos_b = self.get_arrow_ab(u)
                #arrow = FancyArrowPatch(
                #     posA=(pos_a, self.yt - 2 - i + self.scale * self.height), \
                #     posB=(pos_b, self.yt - 2 - i + self.scale * self.height), fc='#43CD80', ec='#43CD80',
                #     arrowstyle='simple, head_width=1, head_length=1, tail_width=5', capstyle='butt')
                #ax.add_artist(arrow)
            # uorf text 
            if self.strand == '+':
                uorf_text_position = uorf[-1][1] * 1.01
            else:
                uorf_text_position = uorf[0][1] * 1.01
            plt.text(uorf_text_position, self.yt - 2 - i, \
                     self.transcript_id + '.uORF', fontsize=5)
        # arrow show the strand
        #if self.strand == '+':
        #    arrow = FancyArrowPatch(
        #        posA=(self.exon_location[0][0], 15), posB=(self.exon_location[-1][1], 15), fc=None, ec='#1f77b4',
        #        arrowstyle='simple, head_width=0.5, head_length=0.5, tail_width=0.2')
        #else:
        #    arrow = FancyArrowPatch(
        #        posA=(self.exon_location[-1][1], 15), posB=(self.exon_location[0][0], 15), fc=None, ec='#1f77b4',
        #        arrowstyle='simple, head_width=10, head_length=10, tail_width=5')

        ax.add_artist(arrow)
        # plt.axvline(x=self.exon_location[-1][1])
        plt.xlim(self.exon_location[0][0] - 100, self.exon_location[-1][1] + 500)
        plt.ylim(-2 * len(self.uorf_location), 10 + 2 * len(self.uorf_location))
        plt.axis('off')
        plt.savefig(f'{self.output}.pdf')
        plt.show()

# a = visual_transcript_genome([[1,10],[20,30],[34,40],[44,55]], \
#                              [[5,10],[20,30],[34,40],[44,50]], \
#                              [[[2,10],[20,24]],[[2,10],[20,30],[34,36]],[[5,10],[20,30]]])
# a.draw()
