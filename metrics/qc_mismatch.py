from __future__ import division
import pysam
import re
from bedcheck import *

mdparse = re.compile(r'[0-9]+|[TCGNA^]+')


def mismatch(bamfile, vcf='', samplename='', bedfile=''):
    if vcf:
        variants = get_variants(vcf, samplename)
    else:
        variants = set()
    if bedfile:
        bedregion = BedFileChecker(bedfile)
    with pysam.Samfile(bamfile, 'rb') as stream:
        nr = 0
        dm = [0, 0, 0, 0]
        chroms = []
        for d in stream.header['SQ']:
            chroms.append(d['SN'])        
        for rec in stream:
            chrom = chroms[rec.tid]
            if bedfile:
                readpos = rec.get_reference_positions()
                if not readpos:
                    continue
                if (not bedregion.inbed(chrom, readpos[0])) and \
                     (not bedregion.inbed(chrom, readpos[-1])):
                    continue
            if not valid_record(rec):
                continue
            nr += 1
            MD = rec.get_tag('MD')
            if MD.isdigit():
                dm[0] += 1
                continue
            v = variants_from_md(rec, MD, stream.references[rec.tid])
            nm = len(v)
            for el in v:
                if el in variants:
                    nm -= 1
            if nm > 3:
                nm = 3
            dm[nm] += 1

    dm = [el/nr for el in dm]
    return dm


def get_variants(vcf, samplename):

    with open(vcf, 'rt') as stream:
        while 1:
            line = stream.readline()
            if not line:
                raise ValueError('File has no #CHROM record')
            if line[:6] == '#CHROM':
                samples = line.rstrip().split('\t')[9:]
                try:
                    sampidx = samples.index(samplename)
                except IndexError:
                    raise ValueError('File does not have sample')
                else:
                    break
        variants = set()
        while 1:
            line = stream.readline()
            if not line:
                return variants
            args = line.split('\t')
            if 'PASS' not in args[6].split(';'):
                continue
            pos = int(args[1])
            ref = args[3]
            alt = args[4]
            if ',' not in alt:
                for c, p, r, a in normalize(args[0], pos, ref, alt):
                    if len(r) == len(a) == 1:
                        variants.add((c, p, r, a))
            else:
                for aa in args[4].split(','):
                    for c, p, r, a in normalize(args[0], pos, ref, aa):
                        if len(r) == len(a) == 1:
                            variants.add((c, p, r, a))
                

def normalize(chrom, pos, ref, alt):
    """Normalize a variant so that it is parsimonious and left-aligned
    """

    if 1 == len(ref) == len(alt):
        return [(chrom, pos, ref, alt)]
    else:
        ref, alt = right_trim(ref, alt)
        p, r, a = left_trim(pos, ref, alt)
        rlen = len(r)
        alen = len(a)
        if rlen == alen == 1:
            return [(chrom, p, r, a)]
        if (rlen == 1 or alen == 1) and (r[0] == a[0]):
            return [(chrom, p, r, a)]
        if rlen == alen:
            return primitives(chrom, p, r, a)
        else:
            return [(chrom, p, r, a)]


def right_trim(r, a):
    """Given two strings `r` and `a` remove any common suffix
    between the characters
    """
    while r[-1] == a[-1]:
        if len(r) > 1 and len(a) > 1:
            r = r[:-1]
            a = a[:-1]
        else:
            break
    return r, a


def left_trim(pos, r, a):
    """Given two strings `r` and `a` remove any common prefix
    between the characters. The two strings passed to this
    should have *previously been trimmed* using right_trim
    """
    if len(r) == len(a):
        clen = 1
    else:
        clen = 2
    while len(r) > 1 and len(a) > 1:
        if r[:clen] == a[:clen]:
            r = r[1:]
            a = a[1:]
            pos += 1
        else:
            break
    return pos, r, a


def primitives(chrom, pos, ref, alt):
    """Decompose a variant into it's primitives
    """
    if len(ref) != len(alt):
        return [(chrom, pos, ref, alt)]
    p = []
    for idx, (r1, a1) in enumerate(zip(ref, alt)):
        if r1 != a1:
            p.append((chrom, pos + idx, r1, a1))
    return p


def valid_record(rec):
    if rec.is_unmapped:
        return 0
    c = rec.cigar
    if len(c) != 1: # No indels, hard or soft-clipped bases
        return 0
    if c[0][0] != 0:
        return 0
    if rec.mapping_quality < 20:
        return 0
    if rec.is_duplicate:
        return 0
    if rec.is_secondary:
        return 0
    if rec.is_qcfail:
        return 0
    return rec.is_proper_pair


def variants_from_md(rec, MD, chrom):
    rstart = rec.reference_start + 1
    v = []
    seq = rec.query_alignment_sequence
    sstart = 0
    for item in mdparse.findall(MD):
        if item.isdigit():
            sstart += int(item)
        else:
            for rn in item:
                v.append((chrom, rstart + sstart, rn, seq[sstart]))
                sstart += 1
    return v


if __name__ == "__main__":
    import sys
    import argparse
    desc = 'Mismatch read wise ' +\
    'defined as any substitution in the aligned' +\
    ' portion of a read compared with the reference genome'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-i', '--inputBAM', dest='bam', type=str,
                        help='Input BAM files. You may provide more than '+\
                        'one bam file by comma separated.')
    parser.add_argument('-s', '--samplename', dest='samplename', default='', type=str,
                        help='Sample Name')
    parser.add_argument('-v', '--vcffile', dest='vcffile', default='', type=str,
                        help='vcffile')            
    parser.add_argument('-o', '--out', dest='outfile', type=str,
                        help='Output tsv file')
    parser.add_argument('-r', '--bedfile', dest='bedfile', default='', type=str,
                        help='bed file with region of interest')
    args = parser.parse_args()
    
    d = mismatch(args.bam, args.vcffile, args.samplename, args.bedfile)
    tmp =  [args.samplename]
    out = open(args.outfile, 'w')
    for key in 0, 1, 2, 3:
        tmp.append(str(d[key]))
        #print '%d\t%g' % (key, d[key])
    out.write('#sample\t0_mismatch\t1_mismatch\t2_mismatch\t3_mismatch\n')
    out.write('\t'.join(tmp) + '\n')

