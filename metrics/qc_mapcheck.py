"""
Various checks on mapped reads
"""
from __future__ import division
import pysam
from collections import defaultdict
import ostats
import os
import sys
from multiprocessing import Pool
from bedcheck import *

class Mapcheck(object):

    __slots__ = ['bamfile', 'capfile', 'gencode', 'mapq', 'nfreads', 'nbreads', 
                  'capreads', 'gcreads', 'cgreads', 'lqbases', 'npairs', 'ndups', 'ndc',
                 'nunmapped_r', 'nunmapped_p', 'nfr', 'nff', 'nbadfr',
                 'nlowmapq_r', 'nlowmapq_p',
                 'nunpaired', 'fragsize', 'hqpairs', 'hqpairs_pp', 'hqpairs_ot']

    def __init__(self, bamfile, capfile, gencode, mapq=20):
        self.bamfile = bamfile
        self.capfile = capfile
        self.gencode = gencode
        self.mapq = mapq

        # initializations
        self.nfreads = 0
        self.nbreads = 0
        self.capreads = 0
        self.gcreads = 0
        self.cgreads = 0
        self.lqbases = 0
        self.npairs = 0
        self.ndups = 0
        self.ndc = 0
        self.nunmapped_r = 0
        self.nunmapped_p = 0
        self.nfr = 0
        self.nff = 0
        self.nbadfr = 0
        self.nlowmapq_r = 0
        self.nlowmapq_p = 0
        self.nunpaired = 0
        self.hqpairs = 0
        self.hqpairs_pp = 0
        self.hqpairs_ot = 0
        self.fragsize = defaultdict(int)

    def run(self):
        cache = {}
        smapq = self.mapq
        _isize = self._isize
        capregion = BedFileChecker(self.capfile)
        gcregion = BedFileChecker(self.gencode)
        sf = pysam.Samfile(self.bamfile, 'rb')

        with pysam.Samfile(self.bamfile, 'rb') as stream:
            for rec in stream:
                if not rec.is_secondary:
                    self.nfreads += 1
                self.nbreads += 1
                cap = False
                gc = False
                if not rec.is_secondary:
                    self.lqbases += self.numlqb(rec)
                    if (not rec.is_duplicate) and (rec.mapq >= 20 and rec.mapq < 255):
                        if self.inbed(capregion, rec, sf):
                            cap = True
                            self.capreads += 1
                        #if self.inbed(gcregion, rec, sf):
                        #    gc = True
                        #    self.gcreads += 1
                        #if cap and gc:
                        #   self.cgreads += 1
                
                if rec.qname not in cache:
                    cache[rec.qname] = rec
                    if len(cache) >= 10000000:
                        raise ValueError('May not be a paired-end file')
                    continue
                prec = cache.pop(rec.qname)
                self.npairs += 1

                do_isize = True
                if rec.is_duplicate:
                    self.ndups += 2
                    do_isize = False

                nun = rec.is_unmapped + prec.is_unmapped
                if nun:
                    self.nunmapped_r += nun
                    self.nunmapped_p += 1
                    do_isize = False

                if rec.tid != prec.tid:
                    self.ndc += 1
                    do_isize = False

                try:
                    nlow = (rec.mapq < smapq) + (prec.mapq < smapq)
                except:
                    pass
                else:
                    if nlow:
                        self.nlowmapq_r += nlow
                        self.nlowmapq_p += 1
                        do_isize = False

                if not do_isize: continue
                self.hqpairs += 1
                if rec.is_proper_pair:
                    self.hqpairs_pp += 1
                    if cap and self.inbed(capregion, prec, sf):
                        self.hqpairs_ot += 1
                A = rec.is_reverse
                B = prec.is_reverse
                if (A and B) or (not A and not B):
                    self.nff += 1
                else:
                    self.nfr += 1
                    if (A and not B):
                        _isize(prec, rec)
                    else:
                        _isize(rec, prec)

        self.nunpaired = len(cache)
    
    def inbed(self, bedregion, read, sf):
        for pos in read.positions:
            if bedregion.inbed(sf.getrname(read.tid), pos):
                return True
        
        return False
    
    def numlqb(self, read):
        lqb = 0
        for l in read.query_qualities:
            if l < 20: lqb += 1
        return lqb

    def _isize(self, r1, r2):

        hasindel = False
        for i, j in r1.cigar:
            if i not in (0, 7, 8):
                hasindel = True
                break
        if hasindel:
            return
        for i, j in r2.cigar:
            if i not in (0, 7, 8):
                hasindel = True
                break
        if not hasindel:

            # r1.pos       r1.aend
            # +------------>
            #      <---------------+
            #      r2.pos          r2.aend

            if r1.pos > r2.pos:
                self.nbadfr += 1
            else:
                tlen = r2.aend - r1.pos
                self.fragsize[tlen] += 1

    def stats(self):
        s = ostats.OStats()
        s.update(self.fragsize)
        return s


    def report(self):
        s = self.stats()
        r = ['%d' % self.nfreads,
             '%d' % self.nbreads,
             '%d(%0.1f%%)' % (self.capreads, self.capreads * 100.0 / self.nfreads,),
             '%d(%0.1f%%)' % (self.gcreads, self.gcreads * 100.0 / self.nfreads,),
             '%d(%0.1f%%)' % (self.cgreads, self.cgreads * 100.0 / self.nfreads,),
             '%d' % self.npairs,
             '%d(%0.1f%%)' % (self.lqbases, self.lqbases * 100.0/(self.nfreads * s.minobs),),
             '%d(%0.1f%%)' % (self.ndups, self.ndups * 100.0 / self.nfreads,),
             '%d(%0.1f%%)' % (self.nunmapped_r, self.nunmapped_r * 100.0/self.nfreads,),
             '%d(%0.2f%%)' % (self.ndc, self.ndc * 100.0/self.nfreads,),
             '%d(%0.1f%%)' % (self.nlowmapq_r, self.nlowmapq_r * 100.0/(self.nfreads - self.nunmapped_r),),
             '%d(%0.1f%%)' % (self.hqpairs_pp, self.hqpairs_pp * 100.0/(self.nfreads/2),),
             '%d(%0.1f%%)' % (self.hqpairs_ot, self.hqpairs_ot * 100.0/(self.nfreads/2),),
#             '%d(%0.1f%%)' % (self.nfr, self.nfr * 100.0/ (self.npairs - self.nunmapped_p - self.ndups//2),),
             '%d(%0.1f%%)' % (self.nfr, self.nfr * 100.0/self.hqpairs,),
             '%d' % self.nff,
             '%d' % s.minobs,
             '%d' % s.maxobs,
             '%d' % s.median,
             '%d' % s.mean,
             '%0.2f' % s.sdev,]

        return '\t'.join(r)


def mapcheck(infiles):
    fname, capfile, gencode = infiles
    m = Mapcheck(fname, capfile, gencode)
    m.run()
    return m

def format(rawdata):
    formreport = []
    formreports.append(rawdata[0])
    formreports.append(rawdata[1])

def updtglob(gvals, sample):
    gvals.nfreads += sample.nfreads
    gvals.nbreads += sample.nbreads
    gvals.capreads += sample.capreads
    gvals.gcreads += sample.gcreads
    gvals.cgreads += sample.cgreads
    gvals.npairs += sample.npairs
    gvals.ndups += sample.ndups
    gvals.ndc += sample.ndc
    gvals.nunmapped_r += sample.nunmapped_r
    gvals.nunmapped_p += sample.nunmapped_p
    gvals.nfr += sample.nfr
    gvals.nff += sample.nff
    gvals.nbadfr += sample.nbadfr
    gvals.nlowmapq_r += sample.nlowmapq_r
    gvals.nlowmapq_p += sample.nlowmapq_p
    gvals.nunpaired += sample.nunpaired
    gvals.hqpairs += sample.hqpairs
    gvals.hqpairs_pp += sample.hqpairs_pp
    gvals.hqpairs_ot += sample.hqpairs_ot
    for el in sample.fragsize:
        gvals.fragsize[el] += sample.fragsize[el]

def hist(fname, fragsize):
    out = open(fname + '_hist', 'w')
    k = fragsize.keys()
    k.sort()
    for el in k:
        out.write(str(el) + '\t' + str(fragsize[el]) + '\n')

def main(samples, capfile, gencode, nprocs=1, prefix=''):
    nprocs = int(nprocs)
    p = Pool(nprocs)
    infiles = []
    for s in samples:
        infiles.append((s, capfile, gencode))
    results = p.map(mapcheck, infiles)
    gstats = ostats.OStats()
    gvals = Mapcheck('global.bam', capfile, gencode)
    if prefix:
        outfile = prefix + '_mapcheckdtld.tsv'
    else:
        outfile = 'exomes_' + str(len(samples)) + '_mapcheck.tsv'    
    out = open(outfile, 'w')    
    header1 = ['#sample', 'reads_uniq', 'reads_bam', 'reads_cap', 'reads_gc', 
                'reads_cap_gc', 'paired_reads', 'bqlt20', 'duplicates', 'unmapped',
                'rpdc', 'reads_mqlt20', 'umhq', 'umhq_cap', 'ros', 'rss',
               'shortest_fs', 'longest_fs', 'median_fs', 'mean_fs', 'sd_fs']
    
    tmp = '\t'.join(header1)
    out.write(tmp + '\n')
    for r in results:
        fname = os.path.basename(r.bamfile)
        out.write('%s\t%s\n' % (fname.split('.')[0], r.report()))
        
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Generate statistics for mapped reads')
    parser.add_argument('bamfiles', type=str, nargs='+', help='BAM files')
    parser.add_argument('-n', dest='nprocs', default=1, help='number of processors')
    parser.add_argument('-b', dest='capfile', help='capture file')
    parser.add_argument('-g', dest='gencode', help='gencode file')
    parser.add_argument('-o', dest='prefix', default='', help='outfile prefix')    
    
    options = parser.parse_args()
    main(options.bamfiles, options.capfile, options.gencode, options.nprocs, options.prefix)








