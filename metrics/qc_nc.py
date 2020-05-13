from __future__ import division
import sys
import vcf
import normalize
import os
from bedcheck import *

CMAPS = []
NUCS = 'ACGT'
COMP = dict(zip(NUCS, NUCS[::-1]))
for i in range(len(NUCS)):
    b1 = NUCS[i]
    for j in range(len(NUCS)):
        b2 = NUCS[j]
        if b1 == b2: continue
        CMAPS.append('%s>%s' % (b1, b2))
CMAPS.sort()

def is_ma_snp(var):
    for el in var:
        if len(el) == 1:
            return True
    return False

def freq(kgaf):
    if not kgaf:
        return False
    for af in kgaf:
        if float(af) >= 0.001:
            return True
    return False

def is_mnp(ref, alt):
    for el in alt:
        if len(ref) == len(el):
            return True
    return False
            
def main(fname, bedfile, outfile='', gqcutoff=30.0, freqcutoff=0.001):

    v = vcf.VCFParser(fname)
    if bedfile:
        bedregion = BedFileChecker(bedfile)
    samples = v.samples
    if not outfile:
        outfile = os.path.basename(fname).split('.')[0] + '_nc.tsv'
    out = open(outfile, 'w')
    snp = {}
    for key in samples:
        snp[key] = {}
        snp[key]['hq'] = [0, 0]
        snp[key]['sites'] = [0, 0]
        for c in CMAPS:
            snp[key][c] = [0, 0]
        snp[key]['hom'] = [0, 0]
    for rec in v:
        if 'PASS' not in rec.filter: continue
        v.parseinfo(rec)
        if bedfile and not bedregion.inbed(rec.chrom, rec.pos): continue
        v.parsegenotypes(rec)
        for alidx, allele in enumerate(rec.alt):
            if allele == "*": continue
            normp, normref, normalt = normalize.normalize(int(rec.pos), rec.ref, allele)
            for p, r, a1 in normalize.primitives(int(normp), normref, normalt):
                if len(r) == 1 and len(a1) == 1:
                    field = snp
                else:
                    continue
                nc = r + '>' + a1
                for sid in samples:
                    s = getattr(rec, sid)
                    hom = False
                    if s.GT != './.' and s.GT != '0/0':
                        aa = s.GT.split('/')
                        if aa[0] == aa[1]:
                            hom = True
                        if 'PASS' in rec.filter and s.GQ >= gqcutoff:
                            if rec.info.KGAF and float(rec.info.KGAF[alidx]) >= freqcutoff:
                                field[sid]['sites'][0] += 1
                                if hom: field[sid]['hom'][0] += 1
                            else:
                                field[sid]['sites'][1] += 1
                                if hom: field[sid]['hom'][1] += 1
                            gts = s.GT.split('/')
                            #if rec.info.DB137 or rec.info.KGDB:
                            for gidx, g in enumerate(gts):
                                if alidx + 1 == int(g):
                                    if 'PASS' in rec.filter and s.GQ >= gqcutoff:
                                        #if rec.info.DB137 or rec.info.KGDB:
                                        if rec.info.KGAF and float(rec.info.KGAF[alidx]) >= freqcutoff:
                                            field[sid]['hq'][0] += 1
                                            field[sid][nc][0] += 1
                                        else:
                                            field[sid]['hq'][1] += 1
                                            field[sid][nc][1] += 1

    hdr = ['#sample', 'hq_snps']
    
    for c in CMAPS:
        hdr.append(c)
    hdr.append('sites')
    hdr.append('hom')
    hdr.append('het/hom')
    out.write('\t'.join(hdr) + '\n')
    
    for s in samples:
        tot = snp[s]['hq'][0] + snp[s]['hq'][1]
        l = [s, tot]
        for c in CMAPS:
            ctot = snp[s][c][0] + snp[s][c][1]
            l.append(round(ctot/tot, 4))
        sites = snp[s]['sites'][0] + snp[s]['sites'][1]
        hom = snp[s]['hom'][0] + snp[s]['hom'][1]
        l.append(sites)
        l.append(hom)
        try:
            l.append(round((sites - hom)/hom, 2))
        except ZeroDivisionError:
            l.append(0)
        l = [str(el) for el in l]
        out.write('\t'.join(l) + '\n')

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Variant statistics in the capture region')
    parser.add_argument('-i', dest='infile', help='input vcf file')
    parser.add_argument('-b', dest='bedfile', help='Bed file')
    parser.add_argument('-q', dest='mingq', default=30.0, help='min genotype quality cutoff (default 30.0)')
    parser.add_argument('-f', dest='freq', default=0.001, help='freq cutoff')
    parser.add_argument('-o', dest='outfile', default='', help='output file')
    
    options = parser.parse_args()
    main(options.infile,options.bedfile, options.outfile, float(options.mingq), float(options.freq))


