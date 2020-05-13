import glob
import sys

#import xhets
xhets = set()
for el in open(sys.argv[1]):
    if el.startswith('x_snps'): continue
    args = el.strip().split()
    xhets.add(tuple(args))

cov_gender = set()
for el in glob.glob('*sex_cov.tsv'):
    for l in open(el):
        if l.startswith('#'): continue
        args = l.strip().split()
        cov_gender.add(tuple(args))

merged = []
for el in cov_gender:
    for el2 in xhets:
        if el2[0] == el[0]:
            r1 = round(float(el[1])/float(el[2]), 6)
            r2 = round(float(el[2])/float(el[1]), 6)
           
            merged.append([el[0], el[1], el[2], el[3], el[4], el[5], el[6], el[7], el[8], el[9], el[10], el2[4], el[11]])
            break
hdr = '#sample	Chrom_4	Chrom_15	Chrom_X	Chrom_Y	Chrom_X/4	Chrom_X/15	Chrom_Y/4	Chrom_Y/15	Dist_Female	Dist_Male Xhet/Xtot Sex'.split()
print '\t'.join(hdr)

merged.sort()
for el in merged:
    print '\t'.join(el)
        

