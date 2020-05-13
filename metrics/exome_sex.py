from __future__ import division
from collections import Counter, defaultdict
import pysam
import os
import math


def get_intervals(bedfile):
    intervals = defaultdict(list)
    PAR = {'X': ((60000, 2699520), (154931044, 155260560)),
           'Y': ((10000, 2649520), (59034049, 59363566))}
    with open(bedfile, 'rt') as stream:
        for line in stream:
            if line[:1] == '#' or line[:5] == 'track':
                continue
            chrom, start, end = line.split('\t', 3)[:3]
            start = int(start)
            end = int(end)
            if chrom in 'XY':
                for a, b in PAR[chrom]:
                    if a <= start < b or a <= end < b:
                        break
                else:
                    intervals[chrom].append((start, end))
            else:
                intervals[chrom].append((start, end))
    return intervals


def get_coverage(bamfile, intervals):
    chromosomes = ['4', '15', 'X', 'Y']

    with pysam.Samfile(bamfile, 'rb') as bstream:
        chromosome_coverage = Counter()
        for chrom in chromosomes:
            clen = n = 0
            for a, b in intervals[chrom]:
                clen += b - a
                for rec in bstream.fetch(chrom, a, b):
                    if rec.is_unmapped or rec.is_duplicate or \
                       rec.is_secondary or \
                       rec.mapping_quality < 20:
                        continue
                    n += 1
            chromosome_coverage[chrom] = n / clen
    return chromosome_coverage

if __name__ == "__main__":
    import argparse
    desc = """Sex assignment from chromosome coverage

    Analyze the coverage of X and Y chromosomes relative to the autosomal
    chromosomes and decide if the sample is consistent with Male or
    Female sex. Chromosomes 4 and 15 are chosen for the comparision,
    since they have similar number of capture base positions
    (2561961 and 2339757 respectively) as chromosome X (2340535).
    In computing the coverage for X and Y chromosomes, intervals that
    overlap with the pseudoautosomal regions (PAR) are excluded.
    
    The procedure is as follows
    1. For each chromosome compute the coverage as the number of reads mapping
       to the capture intervals on that chromosome divided by the number of
       capture bases on that chromosome
    
    2. Compute the ratio of the coverage of the X and Y chromosomes to that of
       chromosomes 4 and 15. If the ratio for X is close to 1, and Y is close
       to 0 then the sample is Female. If the ratio of X and Y are both
       close to 0.5, then the sample is Male

    The output contains the following values:

    1. Name of sample
    2. Coverage of Chromosome 4
    3. Coverage of Chromosome 15
    4. Coverage of Chromosome X
    5. Coverage of Chromosome Y
    6. Ratio of coverage of X to 4
    7. Ratio of coverage of X to 15
    8. Ratio of coverage of Y to 4
    9. Ratio of coverage of Y to 15
    10. Distance to ideal female
    11. Distance to ideal male
    12. Assigned sex (M = Male, F = Female)
    """
    
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-b', '--bedfile', dest='bedfile', type=str,
                        required=True,
                        help='BED file of capture regions')
    parser.add_argument('-i', '--bamfile', dest='bamfile', type=str,
                        required=True,
                        help='BAM file')
    parser.add_argument('-o', '--outfile', dest='outfile', type=str,
                        required=True,
                        help='File for results')
    args = parser.parse_args()

    ivls = get_intervals(args.bedfile)
    cvg = get_coverage(args.bamfile, ivls)
    c4 = cvg['4']
    c15 = cvg['15']
    x = cvg['X']
    y = cvg['Y']
    x4 = x / c4
    x15 = x / c15
    y4 = y / c4
    y15 = y / c15
    dfemale = (math.sqrt((x15 - 1)**2 + (y15 ** 2)) +
               math.sqrt((x4 - 1)**2 + (y4 ** 2))) / 2
    dmale = (math.sqrt((x15 - 0.5)**2 + ((y15 - 0.5) ** 2)) +
             math.sqrt((x4 - 0.5)**2 + ((y4 - 0.5)** 2))) / 2

    if dmale < dfemale:
        sex = 'M'
    else:
        sex = 'F'
    sname = os.path.basename(args.bamfile).split('.')[0]
    fields = ['#sample', 'Chrom_4', 'Chrom_15', 'Chrom_X', 'Chrom_Y',
              'Chrom_X/4', 'Chrom_X/15', 'Chrom_Y/4', 'Chrom_Y/15',
              'Dist_Female', 'Dist_Male', 'Sex']
    with open(args.outfile, 'w') as ostream:
        ostream.write('\t'.join(fields))
        ostream.write('\n')
        output = [sname, '%g' % c4, '%g' % c15, '%g' % x, '%g' % y,
                  '%g' % x4, '%g' % x15, '%g' % y4, '%g' % y15,
                  '%g' % dfemale, '%g' % dmale, sex]
        ostream.write('\t'.join(output))
        ostream.write('\n')
