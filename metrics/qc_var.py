from __future__ import division
import sys
import vcf
import normalize
import os
from bedcheck import *

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
                    
def main(fname, bedfile, outfile='', gqcutoff=30.0, freqcutoff=0.001):

    v = vcf.VCFParser(fname)
    if bedfile:
        bedregion = BedFileChecker(bedfile)
    samples = v.samples
    if not outfile:
        outfile = os.path.basename(fname).split('.')[0] + '_var.tsv'
    out = open(outfile, 'w')
    snp = {}
    indel = {}
    mnp = {}
    for key in samples:
        snp[key] = {}
        indel[key] = {}
        mnp[key] = {}
        snp[key]['tot'] = [0, 0, 0, 0]
        snp[key]['hq'] = [0, 0, 0, 0]
        indel[key]['tot'] = [0, 0, 0, 0]
        indel[key]['hq'] = [0, 0, 0, 0]
        mnp[key]['tot'] = [0, 0, 0, 0]
        mnp[key]['hq'] = [0, 0, 0, 0]

    for rec in v:
        v.parseinfo(rec)
        if bedfile and not bedregion.inbed(rec.chrom, rec.pos): continue
        v.parsegenotypes(rec)
        #print rec.chrom, rec.pos, rec.ref, rec.alt
        for alidx, allele in enumerate(rec.alt):
            if allele == "*": continue
            normp, normref, normalt = normalize.normalize(int(rec.pos), rec.ref, allele)
            #print "normp", normp, normref, normalt
            for p, r, a1 in normalize.primitives(int(normp), normref, normalt):
                #print "p r a1", p, r, a1
                if len(r) == 1 and len(a1) == 1:
                    field = snp
                else:
                    field = indel
                # check if its a novel variant and shared by many and hence fp
                shared_variant = False
                if rec.info.KGDB_AF:
                    pass
                else:
                    if shared(rec, samples, alidx):
                        if 'PASS' in rec.filter:
                            #if rec.info.KGDB_AF:
                                #print 'kgaf', rec.info.KGDB_AF[alidx]
                            #print "shared false positive", rec.chrom, rec.pos, rec.ref, rec.alt, alidx, rec.filter
                            shared_variant = True
                
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
                            if alidx + 1 == int(g):
                                if rec.info.KGDB_AF:
                                    if float(rec.info.KGDB_AF[alidx]) >= freqcutoff:
                                        field[sid]['tot'][0] += 1
                                    else:
                                        field[sid]['tot'][1] += 1
                                else:
                                    # novel 
                                    if not shared_variant:
                                        field[sid]['tot'][2] += 1
                                    else:
                                        field[sid]['tot'][3] += 1
                                
                                
                                if 'PASS' in rec.filter and s.GQ >= gqcutoff:
                                    #if rec.info.DB137 or rec.info.KGDB:
                                    if rec.info.KGDB_AF:
                                        if float(rec.info.KGDB_AF[alidx]) >= freqcutoff:
                                            field[sid]['hq'][0] += 1
                                        else:
                                            field[sid]['hq'][1] += 1
                                    else:
                                        #novel
                                        if not shared_variant:
                                            r = [rec.chrom, rec.pos, sid, s.GT]
                                            r = [str(rr) for rr in r]
                                            #out.write('\t'.join(r) + '\n')
                                            field[sid]['hq'][2] += 1
                                        else:
                                            field[sid]['hq'][3] += 1
    

    out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('#sample', 'snps_all_freq', 'snps_all_rare', 'snps_all_novel1', 'snps_all_novel2', \
                            'indels_all_freq', 'indels_all_rare', 'indels_all_novel1', 'indels_all_novel2', \
                            'snps_hq_freq', 'snps_hq_rare', 'snps_hq_novel1', 'snps_hq_novel2', 'snps_hq_rarenovel',\
                            'indels_hq_freq', 'indels_hq_rare', 'indels_hq_novel1', 'indels_hq_novel2', 'indels_hq_rarenovel'))
    
    for s in samples:
        out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (s, snp[s]['tot'][0], snp[s]['tot'][1], snp[s]['tot'][2], snp[s]['tot'][3],\
                               indel[s]['tot'][0], indel[s]['tot'][1],indel[s]['tot'][2], indel[s]['tot'][3],  \
                                snp[s]['hq'][0], snp[s]['hq'][1], snp[s]['hq'][2], snp[s]['hq'][3], snp[s]['hq'][1] + snp[s]['hq'][2],\
                                indel[s]['hq'][0], indel[s]['hq'][1], indel[s]['hq'][2], indel[s]['hq'][3], indel[s]['hq'][1] + indel[s]['hq'][2]))
    
    
    #if bedfile:
    #    out.write('\n#Variant statistics from %s over bed file %s\n' % (os.path.abspath(fname), bedfile))
    #else:
    #    out.write('\n#Variant statistics from %s :\n' % os.path.abspath(fname))
    #out.write('\n#Filtered - Marked PASS and with GQ >= %d\n' % gqcutoff)
    #out.write('#frequent - SNPs with freq greater than %g\n' % freqcutoff)
    #out.write('#rare - SNPs with freq less than %g\n' % freqcutoff)

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





