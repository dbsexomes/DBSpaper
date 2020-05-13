import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pysam
import numpy
import math
import bedcheck 

def read_tx(txfile):
    txlist = []
    with open(txfile, 'rt') as stream:
        tx = None
        for line in stream:
            if line[:1] == '>':
                if tx:
                    txlist.append(tx)
                    
                title = line[1:].rstrip()
                tx = {'title': title,
                      'intervals': []}
            else:
                region, start, end = line.split('\t')
                tx['intervals'].append((region, int(start), int(end)))
        
    if tx:
        txlist.append(tx)
    return txlist

class CoveragePlot(object):

    def __init__(self, filename, bamfiles, txfile, fontsize=6,
                 numrow=10, numcol=2, bedfile='', genelist=()):
        self.output = PdfPages(filename)
        self._bams = [pysam.Samfile(name, 'rb') for name in bamfiles]
        self.txfile = txfile
        self.fontsize = fontsize
        self.numrow = numrow
        self.numcol = numcol
        self.genelist = genelist
        self.rnum = 0
        self.cnum = 0
        self.fig = None
        self.title = ''
        if bedfile:
            self.bedregion = bedcheck.BedFileChecker(bedfile)
        #for idx, bam in enumerate(bamfiles):
        #    print idx, bam

    def newfig(self):
        if self.fig:
            plt.gca()
            plt.clf()
        else:
            #self.fig = plt.figure(figsize=(8.27, 11.69), dpi=600)
            self.fig = plt.figure(figsize=(11, 8.5), dpi=600)
            #self.fig = plt.figure(figsize=(8.5, 6.5), dpi=600)
            plt.subplots_adjust(left=0.033, right=0.99, bottom=0.01, top=0.985)
        self.rnum = self.cnum = 0
        plt.rc('font', size=self.fontsize)
        return self.fig

    def initplot(self, title, start, end, glen):
        plt.subplot2grid((self.numrow, self.numcol), (self.rnum, self.cnum))
        
        gene, chrom, strand, coords = title.split(':')
        #plt.title(gene, fontsize=self.fontsize)
        plt.text((end - start)/2, 3.15, gene, horizontalalignment='center',fontsize=self.fontsize)
        if self.cnum == 0:
            plt.ylabel('log10(Coverage)', fontsize=self.fontsize)
        #plt.label('CDS + UTR bases = %d' % glen, fontsize=self.fontsize)
        #print "start, end", start, end
        lspace = int((end - start) * 0.03)
        plt.text(lspace, 2.75, '%s: %s' % (chrom, coords), fontsize=2.2)
        plt.text(lspace, 2.60, 'bases=%d: strand=%s' % (glen,strand), fontsize=2.2)
        plt.ylim(-0.5, 3.0)
        plt.xlim(0, end - start + 1)
        plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')

    def draw_legend1(self):
        plt.subplot2grid((self.numrow, self.numcol), (self.rnum, self.cnum))
        
        #plt.title("Legend1", fontsize=self.fontsize)
        plt.text(500, 3.13, 'Legend1', horizontalalignment='center',fontsize=self.fontsize)
        plt.ylim(-0.5, 3.0)
        plt.xlim(0, 1000)
        plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
        plt.tick_params(axis='y', which='both', right='off', left='off', labelleft='off')
        #p1 = plt.Rectangle((0, 0), 1, 1, fc="0.8", alpha=0.6)
        p1 = plt.Rectangle((0, 0), 1, 1, fc="0.3", alpha=0.7)
        p2 = plt.Rectangle((0, 0), 1, 1, fc="0.75", alpha=0.7)
        #p3 = plt.Rectangle((0, 0), 1, 1, fc="0.4", alpha=0.6)
        p3 = plt.Rectangle((0, 0), 1, 1, fc="0.9", alpha=0.6)
        p4 = plt.Rectangle((0, 0), 1, 1, fc="blue", alpha=0.8)
        p5 = plt.Rectangle((0, 0), 1, 1, fc='white')
        legtxt = ["median coverage", "20th percentile coverage", "minimum coverage", "Nimblegen capture", "Outside capture"]
        plt.legend([p1, p2, p3, p4, p5], legtxt, loc=10)
        leg = plt.gca().get_legend()
        lt = leg.get_texts()
        plt.setp(lt, fontsize=3)

    def draw_legend2(self):
        plt.subplot2grid((self.numrow, self.numcol), (self.rnum, self.cnum))
        
        #plt.title("Legend2", fontsize=self.fontsize)
        plt.text(500, 3.13, 'Legend2', horizontalalignment='center',fontsize=self.fontsize)
        plt.ylim(-0.5, 3.0)
        plt.xlim(0, 1000)
        plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
        plt.tick_params(axis='y', which='both', right='off', left='off', labelleft='off')
        extra = plt.Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
        p4 = plt.Rectangle((0, 0), 1, 1, fc="red", alpha=0.7)
        p5 = plt.Rectangle((0, 0), 1, 1, fc="#FF6400", alpha=0.7)
        p6 = plt.Rectangle((0, 0), 1, 1, fc="#DAB419", alpha=0.7)
        p7 = plt.Rectangle((0, 0), 1, 1, fc="#9ACD32", alpha=0.7)
        p8 = plt.Rectangle((0, 0), 1, 1, fc="green", alpha=0.7)
        #legtxt = ["median coverage", "20th percentile coverage", "minimum coverage"]
        #legtxt = ["> 20% exon < 5x", "> 20% exon b/w 5 and 15", "> 20% exon b/w 15 and 25", "> 20% exon b/w 25 and 35"]
        legtxt = ["> 15% less than 10x "] # red
        legtxt.append("5 - 15% below 10x") # orange
        legtxt.append("95% at 10x or less than 5% below 10x") # brown
        legtxt.append("95% at 15x and less than 2% below 10x") # light green
        legtxt.extend(["95% at 20x"]) # dark green
        plt.legend([p4, p5, p6, p7, p8], legtxt, loc=10)
        leg = plt.gca().get_legend()
        lt = leg.get_texts()
        plt.setp(lt, fontsize=2.5)


        
    def txplot(self, chrom, start, end, strand):
        med, mins, p5 = self.getcov(chrom, start, end)
        if strand == '-':
            med.reverse()
            mins.reverse()
            p5.reverse()
        for idx, v in enumerate(p5):
            if v <= 5.0:
                plt.plot(idx, med[idx], 'r.')
        plt.plot(range(1, len(med) + 1), med, 'b')
        plt.plot(range(1, len(med) + 1), mins, 'k')

    def exondata(self, chrom, start, end, strand, emed, emins, ep5, erawmed, ebed, glen):
        rawmed, med, mins, p5, bed = self.getcov(chrom, start, end)
        #if glen:
        glen.extend([0] * 50)
        
        if strand == '-':
            med.reverse()
            mins.reverse()
            p5.reverse()
            rawmed.reverse()
            bed.reverse()
        glen.extend(med)
        #print "extending ", chrom, start, end
        emed.append(med)
        emins.append(mins)
        ep5.append(p5)
        erawmed.append(rawmed)
        ebed.append(bed)

    def plotgene(self, emed, emins, ep5, erawmed, ebed, glen):
        #for idx, v in enumerate(p5):
        #    if v <= 5.0:
        #        plt.plot(idx, med[idx], 'r.')
        #plt.plot(range(1, len(med) + 1), med, 'b')
        #plt.plot(range(1, len(med) + 1), mins, 'k')
        start = 50
        covcol = 'red'
        for idx, med in enumerate(emed):
            if not med:
                continue
            if idx == 0 or idx == len(emed) - 1:
                cfdepth = -0.25
            else:
                cfdepth = -0.5
            fbk = {'lw':0.0, 'edgecolor':None}
            v15 = numpy.percentile(erawmed[idx], 15)
            v5 = numpy.percentile(erawmed[idx], 5)
            exlen = len(erawmed[idx])
            p35 = []
            p20 = []
            p15 = []
            p10 = []
            pl10 = []
            p5 = []
            for i in erawmed[idx]:
                if i >= 35: p35.append(i)
                if i >= 20: p20.append(i)
                if i >= 15: p15.append(i)
                if i >= 10: p10.append(i)
                if i < 10: pl10.append(i)
                if i >= 5: p5.append(i)                                        
            p35 = len(p35) * 100.0/exlen
            p20 = len(p20) * 100.0/exlen
            p15 = len(p15) * 100.0/exlen
            p10 = len(p10) * 100.0/exlen
            pl10 = len(pl10) * 100.0/exlen
            p5 = len(p5) * 100.0/exlen
            #print "p5, p10, p15, p35", p5, p10, p15, p35
            if p20 >= 95:
                covcol = 'green'
            elif (p15 >= 95 and pl10 <=2) or (p20 >= 80 and pl10 <= 2):
                covcol = '#9ACD32'
            elif p10 >= 95 or (pl10 <= 5 and pl10 > 0):
                covcol = '#DAB419'
            elif (pl10 < 15 and pl10 > 5):
                covcol = '#FF6400'
            else:
                covcol = 'red'

            cov = len(med) * [cfdepth]
            #plt.fill_between(range(start, len(med) + start), med, 0, facecolor='#6495ED', alpha=0.5, **fbk)
            #plt.fill_between(range(start, len(med) + start), med, 0, facecolor='0.8', alpha=0.6, **fbk)
            plt.fill_between(range(start, len(med) + start), med, 0, facecolor='0.30', alpha=0.7, **fbk)
            #plt.fill_between(range(start, len(ep5[idx]) + start), ep5[idx], 0, facecolor='#4169E1', alpha=0.2, **fbk)
            plt.fill_between(range(start, len(ep5[idx]) + start), ep5[idx], 0, facecolor='0.75', alpha=0.3, **fbk)
            #plt.fill_between(range(start, len(emins[idx]) + start), emins[idx], 0, facecolor='0.4', alpha=0.4, **fbk)
            plt.fill_between(range(start, len(emins[idx]) + start), emins[idx], 0, facecolor='0.9', alpha=0.3, **fbk)
            #print "covcol", covcol
            plt.fill_between(range(start, len(emins[idx]) + start), cov, 0, facecolor=covcol, alpha=0.6, **fbk)
            capx, capy, noncapx, noncapy = self.capreg(ebed[idx], start - 1)
            for capidx, capxel in enumerate(capx):
                plt.fill_between(capxel, capy[capidx], facecolor='blue', alpha=0.8, **fbk)
            for noncapidx, noncapxel in enumerate(noncapx):
                plt.fill_between(noncapxel, noncapy[noncapidx], facecolor='white', **fbk)
            if idx == 0 or idx == (len(emed) - 2) :
                start += len(med)
            else:
                start += len(med) + 50
        plt.plot(range(1, len(glen) + 1), len(glen) * [0], color='black')
        # scale bar
        sbpos =  int(len(glen) * 0.97)
        plt.plot(range(sbpos - 500, sbpos), [2.75] * 500, color='black')
        plt.plot([sbpos - 500, sbpos - 500], [2.7, 2.75], color='black')
        plt.plot([sbpos, sbpos], [2.7, 2.75], color='black')
        plt.text(sbpos - 250, 2.80, '500', horizontalalignment='center',fontsize=3)   
        

    def capreg(self, bed, globalpos):
        capx = []
        capy = []
        tmpcx = []
        tmpcy = []
        noncapx = []
        noncapy = []
        tmpncx = []
        tmpncy = []
        for idx, el in enumerate(bed):
            globalpos += 1
            
            if el:
                tmpcx.append(globalpos)
                tmpcy.append(-0.1)
                if tmpncx:
                    noncapx.append(tmpncx)
                    noncapy.append(tmpncy)
                    tmpncx = []
                    tmpncy = []
            else:
                if tmpcx:
                    capx.append(tmpcx)
                    capy.append(tmpcy)
                    tmpcx = []
                    tmpcy = []
                tmpncx.append(globalpos)
                tmpncy.append(-0.1)
        capx.append(tmpcx)
        capy.append(tmpcy)
        noncapx.append(tmpncx)
        noncapy.append(tmpncy)        
        
        return capx, capy, noncapx, noncapy   

    def plot(self):
        ntx = 0
        #print "numrow, numcol", self.numrow, self.numcol
        fig = self.newfig()
        self.draw_legend1()
        self.cnum += 1
        self.draw_legend2()
        self.cnum += 1
        for g in self.genelist:
            for tx in read_tx(self.txfile):
                #gene, _, chrom, strand = tx['title'].split(':')
                gene, chrom, strand, _ = tx['title'].split(':')
                if g != gene:
                    continue
                self.title = tx['title']
                #if self.genelist and gene not in self.genelist:
                #    continue
                chrom = chrom.replace('chr', '')
                if '_' in chrom: continue
                #fig = self.newfig()
                emed = []
                emins = []
                ep5 = []
                erawmed = []
                ebed = []
                glen = []
                for reg, start, end in tx['intervals']:
                    #print "reg start end", reg, start, end
                    self.exondata(chrom, start, end, strand, emed, emins, ep5, erawmed, ebed, glen)
                glen = glen[:(len(glen) - 50)]
                totglen = 0
                for m in emed: totglen += len(m)
                

                self.initplot(tx['title'], 0, len(glen), totglen)
                self.plotgene(emed, emins, ep5, erawmed, ebed, glen)
                #self.output.savefig(fig)
                self.cnum += 1
                if self.cnum == self.numcol:
                    self.cnum = 0
                    self.rnum += 1
                    if self.rnum == self.numrow:
                        #plt.tight_layout()
                        self.output.savefig(fig)
                        fig = self.newfig()            
                
                #self.output.savefig(fig)
                ntx += 1
                print 'Processed %d transcripts' % ntx
        self.output.savefig(fig)
        self.output.close()
        return

    def getcov(self, chrom, start, end):
        dall = [[] for i in xrange(start, end)]
        for bam in self._bams:
            d = [0] * (end - start)
            for read in bam.fetch(chrom, start, end):
                if read.is_duplicate or \
                      read.is_qcfail or read.is_secondary or \
                      read.is_unmapped or read.mapq <= 0 or \
                      len(read.seq) != len(read.qual) or read.mapq == 255:
                    continue
                for pos in read.positions:
                    if pos >= start and pos < end:
                        d[pos - start] += 1
            for idx, v in enumerate(d):
                dall[idx].append(v)
        rawmed = [numpy.median(el) for el in dall]
        medians = self.tolog(rawmed)
        mins = [min(el) for el in dall]
        minidx = [el.index(min(el)) for el in dall] 
        mins = self.tolog(mins)
        percentile = [numpy.percentile(el, 20) for el in dall]
        percentile = self.tolog(percentile)
        bed = []
        for pos in range(start, end + 1):
            if self.bedregion.inbed(chrom, pos):
                bed.append(0.1)
            else:
                bed.append(0)
        return rawmed, medians, mins, percentile, bed

    def tolog(self, v):
        return [math.log10(el) if el >= 1 else 0 for el in v]


if __name__ == "__main__":
    import sys
    import argparse

    parser = argparse.ArgumentParser(description='Coverage Plots for Genes')
    pa = parser.add_argument
    pa('-o', '--output', dest='output', type=str, required=True,
       help='Output file name for plots')
    pa('-s', '--fontsize', dest='fontsize', type=float, default=3.7,
       help='font size to use in plots')
    pa('-r', '--numrow', dest='numrow', type=int, default=8,
       help='Number of rows in each page')
    pa('-c', '--numcol', dest='numcol', type=int, default=5,
       help='Number of columns in each page')
    pa('-t', '--txfile', dest='txfile', type=str, required=True,
       help='Name of file containing transcript definitions')
    pa('-g', '--genes', dest='genes', type=str,
       help='Name of file containing gene names')
    pa('-b', '--bedfile', dest='bedfile', type=str, default= '',
       help='Capture bed file')
    pa('-f', '--filelist', dest='filelist', type=str, default= '',
       help='bam file list')
    #pa('bamfiles', type=str, nargs='+',
    #   help='BAM files to be used for coverage')

    args = parser.parse_args()
    if args.genes:
        genelist = set((line.strip() for line in open(args.genes, 'rt')))
        genelist = list(genelist)
        genelist.sort()
    else:
        genelist = ()
    bamfiles = []
    if args.filelist:
       for el in open(args.filelist):
           bamfiles.append(el.strip())
    c = CoveragePlot(args.output, bamfiles, args.txfile, args.fontsize,
                     args.numrow, args.numcol, args.bedfile, genelist)
    c.plot()
















