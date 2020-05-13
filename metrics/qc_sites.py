import vcf
from bedcheck import *
import sys

def sites(infile, outfile, bedfile=''):
    
    numconf = 0
    allsites = 0
    nocall = 0
    out = open(outfile, 'w')
    v = vcf.VCFParser(infile)
    samples = v.samples
    samplestats = {}
    for s in samples:
        samplestats[s] = {'all':0, 'nc':0, 'lqref':0, 'lqsnp':0, 'confref':0, 'confsnp':0}
    if bedfile:
        bedregion = BedFileChecker(bedfile)
    for rec in v:
        v.parseinfo(rec)
        if rec.info.END:
            end = rec.info.END
        else:
            end = rec.pos
        v.parsegenotypes(rec)
        for p in range(rec.pos, end + 1):
            if bedfile and not bedregion.inbed(rec.chrom, p - 1):
                # subtract 1 when going from vcf to check in bed
                continue
            #print rec.chrom, rec.pos, rec.info.END, (rec.info.END - rec.pos + 1)
            for el in samples:
                s = getattr(rec, el)
                samplestats[el]['all'] += 1
                if s.GT == './.':
                    samplestats[el]['nc'] += 1
                    continue
                if s.GT == '0/0':
                    if s.DP == 0 or (s.GQ == 0 and s.DP <= 5):
                        samplestats[el]['nc'] += 1
                        continue
                    if s.GQ >= 30:
                        samplestats[el]['confref'] += 1
                    else:
                        samplestats[el]['lqref'] += 1
                else:
                    # variants
                    if s.DP == 0 or (s.GQ == 0 and s.DP <= 5):
                        samplestats[el]['nc'] += 1
                        continue
                    if s.GQ >= 30:
                        samplestats[el]['confsnp'] += 1
                    else:
                        samplestats[el]['lqsnp'] += 1

    hdr = '#sample all_sites conf_ref lq_ref conf_snp lq_snp no_call hq_fr lq_fr'.split()
    tmp = '\t'.join(hdr)
    out.write(tmp + '\n')
    for el in samples:
        ss = samplestats[el]
        tmp = [el, ss['all'], ss['confref'], ss['lqref'], ss['confsnp'], \
                        ss['lqsnp'], ss['nc'], (ss['confref'] + ss['confsnp'] * 1.0)/ss['all'],
                        (ss['lqref'] + ss['lqsnp'] + ss['nc']) * 1.0/ss['all']]
        tmp = [str(el) for el in tmp]
        tmp = '\t'.join(tmp)
        out.write(tmp + '\n')

if __name__=='__main__':
    
    infile = sys.argv[1]
    outfile = sys.argv[2]
    bedfile = sys.argv[3]
    sites(infile, outfile, bedfile)


