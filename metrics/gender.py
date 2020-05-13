import vcf
from collections import defaultdict
import sys

def in_PAR(chrom, pos):
    PAR = {'X': ((60000, 2699520), (154931044, 155260560)),
           'Y': ((10000, 2649520), (59034049, 59363566))}
    for a, b in PAR[chrom]:
        if a <= pos <= b:
            return True
    return False
    
def check_gender(infile):
    
    v = vcf.VCFParser(infile)
    samples = v.samples
    stats = defaultdict(list)
    for s in samples:
        stats[s] = [0, 0, 0, 0]
    for rec in v:
        if 'PASS' not in rec.filter:
            continue
        
        if not (rec.chrom == 'Y' or rec.chrom == 'X'):
            continue
        v.parseinfo(rec)
        if not rec.info.NimbleGenCaptureCore:
            continue
        v.parsegenotypes(rec)
        for sample in samples:
            
            s = getattr(rec, sample)
            if s.GT[0] != '.' and s.GT != '0/0' and s.GQ > 30:
                if in_PAR(rec.chrom, rec.pos):
                    continue
                if rec.chrom == 'Y':
                    stats[sample][1] += 1
                    if s.GT[0] != s.GT[2]:
                        stats[sample][3] += 1
                        
                if rec.chrom == 'X':
                    stats[sample][0] += 1
                    if s.GT[0] != s.GT[2]:
                        stats[sample][2] += 1
    print '\t'.join('x_snps y_snps xhet xhet_fr y_het yhet_fr'.split())
    for s in samples:
        try:
            tmp = [s, stats[s][0], stats[s][1], stats[s][2], stats[s][2] * 1.0/stats[s][0], 
                stats[s][3], stats[s][3] * 1.0/stats[s][1]]
        except ZeroDivisionError:
            xhet_fr = 0
            yhet_fr = 0
            if stats[s][0] == 0:
                xhet_fr = 0
            else:
                xhet_fr = stats[s][2] * 1.0/stats[s][0]
            if stats[s][1] == 0:
                yhet_fr = 0
            else:
                yhet_fr = stats[s][3] * 1.0/stats[s][1]
            tmp = [s, stats[s][0], stats[s][1], stats[s][2], xhet_fr, stats[s][3], yhet_fr]
        tmp = [str(el) for el in tmp]
        print '\t'.join(tmp)

check_gender(sys.argv[1])
    
        
        
    
    


