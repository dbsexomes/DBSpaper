from __future__ import division
import sys
import vcf
import argparse
import normalize
import os
from bedcheck import *

S, V = 0, 1
TsTvDict = {
    'AT': V, 'AG': S, 'AC': V,
    'GT': V, 'GA': S, 'GC': V,
    'CA': V, 'CG': V, 'CT': S,
    'TA': V, 'TG': V, 'TC': S,
}

def parsesamples(samplefile):
    slist = []
    f = open(samplefile, 'rU')
    for l in f: slist.append(l.strip())
    return slist

def indel(alt):
    for altallele in alt:
        if len(altallele) != 1:
            return True
    return False

def freq(kgaf, freqcutoff):
    if not kgaf:
        return False
    for af in kgaf:
        if float(af) >= float(freqcutoff):
            return True
    return False

def get_family(sid):
    if sid[-1].isdigit():
        return sid[:-2]
    else:
        return sid[:-1]
         

def shared(rec, samples, alidx):
    families = set()
    numsamples = []
    for sid in samples:
        s = getattr(rec, sid)
        if s.GT not in ['./.', '.|.', '0/0', '0|0', '.', '0']:
            if '/' in s.GT:
                sep = '/'
            else:
                sep = '|'

            gts = s.GT.split(sep)
            #if rec.info.DB137 or rec.info.KGDB:
            for gidx, g in enumerate(gts):
                if alidx + 1 == int(g) and s.GQ >= 30:
                    # allele matched
                    family = get_family(sid)
                    families.add(family)
                    numsamples.append(sid)
                    #print sid, s.GT, s.GQ
    #print "families", families
    if len(families) > 2:
        return True
    else:
        return False
        
def convert_gt(ref, alt, gt):
    if '/' in gt:
        int_alleles = gt.split('/')
    else:
        int_alleles = gt.split('|')
    all_allel = [ref] + alt
    sample_alleles = []
    for i in int_alleles:
        for idx, el in enumerate(all_allel):
            if idx == int(i):
                sample_alleles.append(el)
     
    # remove extra allele for homozygous
    if sample_alleles[0] == sample_alleles[1] and sample_alleles[0] != ref:
        sample_alleles = [sample_alleles[0]]
    return sample_alleles
        
def main(fname, bedfile, outfile='', gqcutoff=30.0, freqcutoff=0.001, detailed=False):
    
    v = vcf.VCFParser(fname)
    samples = v.samples

    tstvf = {}
    tstvr = {}
    tstvn1 = {}
    tstvn2 = {}
    tstvfhq = {}
    tstvrhq = {}
    tstvn1hq = {}
    tstvn2hq = {}
    for key in samples:
        tstvf[key] = [0, 0]
        tstvr[key] = [0, 0]
        tstvn1[key] = [0, 0]
        tstvn2[key] = [0, 0]
        tstvfhq[key] = [0, 0]
        tstvrhq[key] = [0, 0]
        tstvn1hq[key] = [0, 0]
        tstvn2hq[key] = [0, 0]
    if not outfile:
        outfile = os.path.basename(fname).split('.')[0] + '_tstv.tsv'
    out = open(outfile, 'w')
    if bedfile:
        bedregion = BedFileChecker(bedfile)
    for rec in v:
        #if indel(rec.alt) or  len(rec.ref) != 1:
        if bedfile and not bedregion.inbed(rec.chrom, rec.pos): continue
        if not 'PASS' in rec.filter:
            continue
        for alidx, allele in enumerate(rec.alt):
            if allele == "*": continue
            normp, normref, normalt = normalize.normalize(int(rec.pos), rec.ref, allele)
            #print "normp", normp, normref, normalt
            for p, r, a1 in normalize.primitives(int(normp), normref, normalt):
                #print "p r a1", p, r, a1
                if len(r) == 1 and len(a1) == 1:
                    # Its a SNP, so process
                    v.parseinfo(rec)
                    v.parsegenotypes(rec)
                    numnovel = 0

                    shared_variant = False
                    if rec.info.KGDB_AF:
                        pass
                    else:
                        if shared(rec, samples, alidx):
                            if 'PASS' in rec.filter:
                                shared_variant = True
                                #print rec.chrom, rec.pos, rec.ref, rec.alt, rec.info.AF
				
                    if rec.info.KGDB_AF:
                        if float(rec.info.KGDB_AF[alidx]) >= freqcutoff:
                            field = tstvf
                            if 'PASS' in rec.filter:
                                fieldhq = tstvfhq
                        else:
                            field = tstvr
                            if 'PASS' in rec.filter:
                                fieldhq = tstvrhq
                                
                    else:
                        # novel SNPs
                        if not shared_variant:
                            field = tstvn1
                            if 'PASS' in rec.filter:
                                fieldhq = tstvn1hq
                        else:
                            field = tstvn2
                            if 'PASS' in rec.filter:
                                fieldhq = tstvn2hq

                    #v.parsegenotypes(rec) # call this if you need to query Genotype data
                    sv = TsTvDict[r + a1]
                    for sid in samples:
                        s = getattr(rec, sid)
                        if s.GT not in ['./.', '.|.', '0/0', '0|0', '.', '0']:
                            if '/' in s.GT:
                                sep = '/'
                            else:
                                sep = '|'
                            gts = s.GT.split(sep)

                            for gidx, g in enumerate(gts):
                                if alidx + 1 == int(g):
                                    field[sid][sv] += 1
                                    if 'PASS' in rec.filter and s.GQ >= gqcutoff:
                                        fieldhq[sid][sv] += 1

 
    if detailed:
        out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('#sample', 'tstv_all', 'tstv_fr', 'tstv_rare', \
            'tstv_rarenovel', 'tstv_novel1', 'tstv_novel2', 'ts_freq', 'tv_freq', 'ts_rare', 'tv_rare', \
            'ts_novel1', 'tv_novel1', 'tv_novel2', 'tv_novel2'))
    else:
        out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('#sample', 'tstv_all', 'tstv_freq', 'tstv_rare', 'tstv_rarenovel', 'tstv_novel1', 'tstv_novel2'))
    
    #print "sample id, tsk, tvk, tsn, tvn"
    for s in samples:
        tsf, tvf = tstvfhq[s]
        tsr, tvr = tstvrhq[s]
        tsn1, tvn1 = tstvn1hq[s]
        tsn2, tvn2 = tstvn2hq[s]
        
        #print "sample id, tsk, tvk, tsn, tvn"
        #print '%s\t%d\t%d\t%d\t%d' % (s, tsk, tvk, tsn, tvn)
        try:
            tstvall = (tsf + tsr + tsn1 + tsn2)/(tvf + tvr + tvn1 + tvn2)
        except ZeroDivisionError:
            tstvall = 0
        try:
            tstvf = tsf/tvf
        except ZeroDivisionError:
            tstvf = 0
        try:
            tstvr = tsr/tvr
        except ZeroDivisionError:
            tstvr = 0
        try:
            tstvrn = (tsr + tsn1)/(tvr + tvn1)
        except ZeroDivisionError:
            tstvrn = 0

        try:
            tstvn1 = tsn1/tvn1
        except ZeroDivisionError:
            tstvn1 = 0        
        try:
            tstvn2 = tsn2/tvn2
        except ZeroDivisionError:
            tstvn2 = 0              
        try:
            if detailed:
                out.write('%s\t%g\t%g\t%g\t%g\t%g\t%g\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n' % (s, round(tstvall, 3), round(tstvf, 3), round(tstvr, 3), round(tstvrn, 3), round(tstvn1, 3), round(tstvn2, 3), 
                            tsf, tvf, tsr, tvr, tsn1, tvn1, tsn2, tvn2))
            else:
                out.write('%s\t%g\t%g\t%g\t%g\t%g\t%g\n' % (s, round(tstvall, 3), round(tstvf, 3), round(tstvr, 3), round(tstvrn, 3), round(tstvn1, 3), round(tstvn2, 3)))
        except ZeroDivisionError:
                if detailed:
                    out.write('%s\t%g\t%g\t%g\n' % (s, 0, 0, 0))
                else:
                    out.write('%s\t%g\t%g\t%g\n' % (s, 0, 0, 0))

    #out.write('\n#High quality SNPs in Capture region marked PASS and GQ >= %d for %s\n' % (gqcutoff, os.path.abspath(fname)))
    #if bedfile:
    #    out.write('#restricted to regions in: %s\n' % os.path.abspath(bedfile)) 
    #out.write('\n\n#Transition to transversion ratio for 1000 genomes (P3) across Nimblegen capture region = 2.51\n')
    #out.write('#ts/tv - Number of Transitions/Number of Transversions\n')
    #out.write('#frequent - SNPs with freq greater than %g\n' % freqcutoff)
    #out.write('#rare - SNPs with freq less than %g\n' % freqcutoff)
    
if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Calculate Ti/Tv ratios for SNPs in the capture region')
    parser.add_argument('-i', dest='infile', help='input vcf file')
    parser.add_argument('-b', dest='bedfile', help='Bed file')
    parser.add_argument('-q', dest='mingq', default=30.0, help='min genotype quality cutoff (default 30.0)')
    parser.add_argument('-f', dest='freq', default=0.001, help='freq cutoff')
    parser.add_argument('-d', dest='dtld', action='store_true', help='detailed report')
    parser.add_argument('-o', dest='outfile', default='', help='output file')
    
    options = parser.parse_args()
    main(options.infile,options.bedfile, options.outfile, float(options.mingq), float(options.freq), options.dtld)








