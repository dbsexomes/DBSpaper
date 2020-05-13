from __future__ import division
from collections import defaultdict, Counter
from ostats import OStats
import pysam
import os
import sys
from multiprocessing import Pool
from bed import readfile


def bediter(fname, np):
    with open(fname, 'r') as stream:
        for line in stream:
            if line[:1] == '#':
                continue
            args = line.rstrip().split('\t')
            chr = args[0]
            args[1] = int(args[1])
            args[2] = int(args[2])
            if np:
                if chr[:3] != 'chr':
                    args[0] = 'chr' + chr
            else:
                if chr[:3] == 'chr':
                    args[0] = chr[3:]
            if len(args) == 3:
                args.append('')
            yield args[:4]


class BedCov(object):

    def __init__(self, bedfile, bamfile, globalonly=False, hqual=False, uniq=False):
        self.bamfile = bamfile
        self.bedfile = bedfile
        self.gstats = OStats()
        self.globalonly = globalonly
        self.hqual = hqual
        self.uniq = uniq
        self._results = []

    def run(self):
        self._bam = pysam.Samfile(self.bamfile, 'rb')
        np = False
        for seq in self._bam.references:
            if seq[:3] == 'chr':
                np = True
                break
        self.biter = bediter(self.bedfile, np)
        ogene = None
        if self.globalonly:
            self.rung()
            return self.reportg(self.gstats)

        for item in self.biter:
            gene = item[-1].split(':')[0]
            item[-1] = gene
            if gene != ogene:
                if ogene:
                    gs = self._doivs(intervals)
                    self.reportg(gs, intervals[-1])
                    self.gstats.update(gs)
                intervals = [item]
                ogene = gene
            else:
                intervals.append(item)
        if intervals:
            gs = self._doivs(intervals)
            self.reportg(gs, intervals[-1])
            self.gstats.update(gs)
        self.reportg(self.gstats)

    def rung_old(self):
        nrec = 0
        bam = self._bam
        for chr, start, end, _ in self.biter:
            clen = end - start
            counts = [0] * clen
            for pc in bam.pileup(chr, start, end, stepper='samtools'):
                pos = pc.pos - start
                if pos >= 0 and pos < clen:
                    counts[pos] += pc.n
            self.gstats.update(counts)

    def rung(self):
        intervals = 0
        for rec in readfile(self.bedfile):
            intervals += 1
            if intervals % 10000 == 0: print "processed.. ", intervals
            chrom = rec.chrom
            #if self.need_prefix and chrom[:3] != 'chr':
            #    chrom = 'chr' + chrom

            start = rec.start
            end = rec.end
            p = self._bam.fetch(bytes(chrom), start, end)
            blen = end - start
            counts = [0] * blen
            readcounts = defaultdict(set)
            preads = 0
            #print "chrom, start, end, blen", chrom, start, end, blen

            for read in p:
                if self.hqual:
                    if read.is_duplicate or \
                          read.is_qcfail or read.is_secondary or \
                          read.is_unmapped or read.mapq < 20 or \
                          (not read.is_proper_pair) or \
                          len(read.seq) != len(read.qual) or read.mapq == 255:
                        continue
                #preads += 1
                #self.nreads += 1
                for pos in read.positions:
                    if pos >= start and pos < end:
                        #print 'fetch pos', pos
                        if self.uniq:
                            if read.qname not in readcounts[pos - start]:
                                counts[pos - start] += 1
                                readcounts[pos - start].add(read.qname)
                        else:
                            counts[pos - start] += 1
            
            #print "counts"
            #print counts
            self.gstats.update(counts)
            
    def _doivs(self, ivs):
        gstats = OStats()
        for iv in ivs:
            cs = self.dobedcov(iv)
            self.report(cs, iv)
            gstats.update(cs)
        return gstats

    def dobedcov(self, interval):
        chr, start, end, rest = interval
        counts = [0] * (end - start)
        clen = len(counts)
        pc = self._bam.pileup(chr, start, end, stepper='samtools')
        for rec in pc:
            pos = rec.pos - start
            if pos >= clen:
                break
            if pos >= 0:
                counts[pos] = rec.n

        return OStats(counts)

    def valid_alignment(self, rec):
        if not rec.is_proper_pair:
            return False
        if rec.is_duplicate:
            return False
        if rec.is_secondary:
            return False
        if rec.mapq <= 0 or rec.mapq >= 255:
            return False
        if rec.is_unmapped:
            return False
        if rec.rlen != len(rec.qual):
            return False
        return True

    def report(self, cs, interval):
        chr, start, end, rest = interval
        bincounts = cs.binned_counts((1, 5, 10, 15, 20, 25, 30, 45), True)
        bc = '\t'.join(['%s' % (round(k, 5)) \
                            for i, (j, k) in bincounts])
        if not rest:
            rest = str(end - start)
        mode = ','.join(str(round(el, 0)) for el in cs.mode)
        row = '%s:%s-%s' % (chr, start, end)
        res = '\t'.join((row, rest, str(cs.nobs),
                         str(round(cs.mean, 1)), str(round(cs.sdev, 1)),
                         str(round(cs.median, 1)), mode, bc))
        self._results.append(res)

    def reportg(self, cs, iv=None):
        if iv:
            r = [iv[0], iv[-1]]
        else:
            r = ['', 'All']
        r.append(str(cs.nobs))
        r.append(str(round(cs.mean, 1)))
        r.append(str(round(cs.sdev, 1)))
        r.append(str(round(cs.median, 1)))
        mode = ','.join(map(str, [round(el, 1) for el in cs.mode]))
        r.append(mode)
        bincounts = cs.binned_counts((1, 5, 10, 15, 20, 25, 30, 45), True)
        bc = '\t'.join(['%s' % (round(k, 5)) \
                            for i, (j, k) in bincounts])
        r.append(bc)
        self._results.append('\t'.join(r))


def run(item):
    bc = BedCov(item[0], item[1], item[2], item[3], item[4])
    bc.run()
    return (item[1], bc._results)

def main(bamfiles, bedfile, nprocs=1, hqual=False, uniq=False, prefix=''):
    gonly = 1
    args = [(bedfile, bam, gonly, hqual, uniq) for bam in bamfiles]
    if nprocs == 1:
        r = [run(args[0])]
    else:
        p = Pool(nprocs)
        r = p.map(run, args)
    r.sort()
    if hqual:
        hqstr = '_hq'
    else:
        hqstr = ''
    if prefix:
        outfile = prefix + '_cov' + hqstr + '.tsv'
    else:
        outfile = 'exomes_' + str(int(len(bamfiles))) + '_cov' + hqstr + '.tsv'
    
    out = open(outfile, 'w')    
    header1 = '#sample\tmean_cov\tsdev_cov\tmed_cov\tmode_cov\t1x\t5x\t10x\t15x\t20x\t25x\t30x\t45x\n'
    out.write(header1)
    for samplefile, result in r:
        sample = os.path.basename(samplefile).split('.')[0]
        row = [sample] + result[0].split('\t')[3:]
        tmp = '\t'.join(row)
        out.write(tmp + '\n')
     

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Generate statistics for overall coverage in bed region')
    parser.add_argument('bamfiles', type=str, nargs='*', help='BAM files')
    parser.add_argument('-f', dest='bamfilelist', default='', help='bamfilelist')
    parser.add_argument('-b', dest='bedfile', help='Bed file')
    parser.add_argument('-n', dest='nprocs', default=1, help='number of processors')    
    parser.add_argument('-q', dest='hqual', action='store_true', help='number of processors')    
    parser.add_argument('-u', dest='uniq', action='store_true', help='count only unique fragments')    
    parser.add_argument('-o', dest='prefix', default='', help='outfile prefix')    
    
    options = parser.parse_args()
    if not options.bamfiles and not options.bamfilelist:
        print 'No bamfiles input'
        sys.exit(1)
    bamfiles = []
    if options.bamfilelist:
        for el in open(options.bamfilelist): bamfiles.append(el.strip())
    else:
        bamfiles = options.bamfiles
    main(bamfiles, options.bedfile, int(options.nprocs), options.hqual, options.uniq, options.prefix)


