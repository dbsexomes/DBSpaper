from __future__ import division
import re
import pysam
from collections import defaultdict
from itertools import izip
from bedcheck import *



MD = re.compile(r'\d+|[TCGNA]+').findall


def normvar(pos, ref, alt, len=len):
    if len(ref) == len(alt) == 1:
        return [(pos, ref, alt)]
    elif (len(ref) == 1) or (len(alt) == 1):
        return [(pos, ref, alt)]

    r, a = right_trim(ref, alt)
    p, r, a = left_trim(pos, r, a)
    rlen = len(r)
    alen = len(a)
    if rlen == 1 or alen == 1:
        return [(p, r, a)]
    elif rlen == alen:
        return [(pos + idx, r1, a1) for idx, (r1, a1) in enumerate(zip(r, a))]
    else:
        return [(p, r, a)]


def right_trim(r, a):
    while r[-1] == a[-1]:
        if len(r) > 1 and len(a) > 1:
            r = r[:-1]
            a = a[:-1]
        else:
            break
    return r, a


def left_trim(pos, r, a):
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


def load_vcf(vcf, sampleid=None):
    variants = {}
    sampidx = 9
    with open(vcf, 'rt') as stream:
        for line in stream:
            if line[:1] != '#':
                args = line.split('\t')
                if 'PASS' not in args[6]:
                    continue
                chrom = args[0]
                for v in get_variants(args, sampidx):
                    try:
                        variants[chrom].add(v)
                    except KeyError:
                        s = variants[chrom] = set()
                        s.add(v)
            elif line[:6] == '#CHROM':
                if sampleid is not None:
                    sampidx = line.rstrip().split('\t').index(sampleid)
    return variants


def get_variants(args, sampidx):
    format = args[8].split(':')
    gtidx = format.index('GT')
    sample = args[sampidx].split(':')
    variants = []
    gt = sample[gtidx]
    if gt in ('0/0', './.', '.', '0'):
        return variants
    try:
        gqidx = format.index('GQ')
    except ValueError:
        gq = 30
    else:
        try:
            gq = int(sample[gqidx])
        except IndexError:
            gq = 0
    if gq < 30:
        return variants

    alts = args[4].split(',')
    alleles = []
    a1 = gt[0]
    if a1 != '0':
        alleles.append(alts[int(a1) - 1])
    try:
        a2 = gt[2]
    except IndexError:
        pass
    else:
        if a2 != '0':
            alleles.append(alts[int(a2) - 1])
    pos = int(args[1])
    ref = args[3]
    for aa in alleles:
        for p, r, a in normvar(pos, ref, aa):
            variants.append((p, r, a))
    return variants


class BaseChange(object):

    __slots__ = ['bamfile', 'mqmin', 'bqmin', 'calledvariants',
                 'basecounts', 'changecounts', 'dicounts', 'bedregion', 'bedfile']

    def __init__(self, bamfile, invcf=None, sampleid=None, mqmin=20, bqmin=20, bedfile=''):
        self.bamfile = bamfile
        self.mqmin = mqmin
        self.bqmin = bqmin
        self.bedfile = bedfile
        if invcf:
            print 'Loading VCF ...'
            self.calledvariants = load_vcf(invcf, sampleid)
            nvar = 0
            for key in self.calledvariants:
                nvar += len(self.calledvariants[key])
            print '... finished loading %d variants' % nvar
        else:
            chroms = []
            for c in range(23): chroms.append(str(c))
            chroms.append('X')
            chroms.append('Y')
            self.calledvariants = {}
            for c in chroms:
                self.calledvariants[c] = set()
    
        self.basecounts = dict(zip('ACGTN', (0, 0, 0, 0, 0)))
        cc = self.changecounts = {}
        for n1 in 'ACGTN':
            for n2 in 'ACGTN':
                if n1 != n2:
                    cc[n1 + n2] = 0
        self.dicounts = defaultdict(int)
        if self.bedfile:
            self.bedregion = BedFileChecker(self.bedfile)
            
    def compute(self):
        with pysam.Samfile(self.bamfile, 'rb') as bstream:
            chroms = []
            for d in bstream.header['SQ']:
                chroms.append(d['SN'])
            mqmin = self.mqmin
            bqmin = self.bqmin
            cc = self.changecounts
            dc = self.dicounts
            cv = self.calledvariants
            bc = self.basecounts
            ntot = nproc = 0
            for rec in bstream:
                ntot += 1
                chrom = chroms[rec.tid]
                if chrom not in cv:
                    continue
                a, b = variants_in_read(rec, bc, mqmin, bqmin)
                nproc += a
                if not (a and b):
                    continue
                for p, r, a, prev, nxt in b:
                    if self.bedfile and not self.bedregion.inbed(chrom, p):
                        continue
                    if (p, r, a) not in cv[chrom]:
                        cc[r + a] += 1
                        dc[prev + r + '>' + prev + a] += 1
                        dc[r + nxt + '>' + a + nxt] += 1
        return ntot, nproc

def variants_in_read(bamrecord, counter, mqmin=20, bqmin=20):
    """Generate a list of all variants in a read.  Returns a
    tuple of 2 items.  The first item in either 0 or 1 indicating
    respectively whether the read was processed or not.  The
    second item is a list of 4-tuples, with the tuple items
    being chrom, pos, ref and alt.  This function only processes
    reads with single nucleotide variants. Reads with indels and
    soft/hard clipped bases are ignored."""
    
    if not bamrecord.cigar:
        return 0, ()

    for op, count in bamrecord.cigar:
        if op:
            return 0, ()

    if (bamrecord.is_duplicate or bamrecord.is_unmapped or
        bamrecord.is_qcfail or bamrecord.is_secondary or
        bamrecord.mapq < mqmin):
        return 0, ()
    if not bamrecord.is_proper_pair:
        return 0, ()
    
    try:
        md = bamrecord.opt('MD')
    except KeyError:
        return 0, ()
    if md.isdigit(): # read is reference
        for c, q in izip(bamrecord.seq, bamrecord.query_qualities):
            if q >= bqmin:
                counter[c] += 1
        return 1, ()

    qref = get_ref(bamrecord.seq, md)
    qlen = len(bamrecord.seq) - 1
    start = bamrecord.pos
    variants = []
    for idx, (r, a, q) in enumerate(izip(qref, bamrecord.seq,
                                         bamrecord.query_qualities)):
        if q >= bqmin:
            counter[r] += 1
            if r != a:
                if idx > 0:
                    prev = qref[idx - 1]
                else:
                    prev = 'N'
                if idx < qlen:
                    next = qref[idx + 1]
                else:
                    next = 'N'
                variants.append((start + idx, r, a, prev, next))

    return 1, variants


def get_ref(seq, md):
    ref = list(seq)
    offset = 0
    for item in MD(md):
        if item.isdigit():
            offset += int(item)
        else:
            ref[offset] = item
            offset += 1
    return ''.join(ref)


if __name__ == "__main__":
    import argparse
    import sys
    import os

    
    desc = """Compute number and fraction of non-variant base changes"""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-v', '--invcf', dest='invcf', type=str,
                        default='',
                        help='Name of VCF file with called variants')
    parser.add_argument('-s', '--sampleid', dest='sampleid', type=str,
                        default='',
                        help='Sample name whose variants are to be analyzed')
    parser.add_argument('-q', '--mqmin', dest='mqmin', type=float,
                        default=20.0,
                        help='Minimum mapping quality of a read to be analyzed')
    parser.add_argument('-b', '--bqmin', dest='bqmin', type=float,
                        default=20.0,
                        help='Minimum quality of a base to be considered')
    parser.add_argument('-f', '--bamfile', dest='bamfile', type=str,
                        help='Name of BAM file')
    parser.add_argument('-r', '--bedfile', dest='bedfile', type=str,
                        help='Name of bed file')
    parser.add_argument('-o', '--outfile', dest='outfile', type=str,
                        default='-',
                        help='Name of file for results')
    
    inputs = parser.parse_args()

    bc = BaseChange(inputs.bamfile, inputs.invcf or None,
                    inputs.sampleid or None, inputs.mqmin,
                    inputs.bqmin, inputs.bedfile)

    ntot, nproc = bc.compute()
    sys.stderr.write('Total number of alignments in file = %d\n' % ntot)
    sys.stderr.write('Number of alignments processed = %d\n' % nproc)

    #header = '#Substitution Frequencies for %s with High Quality Base ' + \
    #         'Count=A=%d,C=%d,T=%d,G=%d\n'
    if inputs.outfile == '-':
        ostream = sys.stdout
    else:
        ostream = open(inputs.outfile, 'w')
    basename = os.path.basename(inputs.bamfile)
    #ostream.write(header % (basename, bc.basecounts['A'], bc.basecounts['C'],
    #                        bc.basecounts['G'], bc.basecounts['T']))
    bsub = basename + '_Substitution'
    hdr = ['#sample']
    fvals = [basename]
    for key in sorted(bc.changecounts.keys()):
        if 'N' in key:
            continue
        keycount = bc.changecounts[key]
        keyfrac = keycount / bc.basecounts[key[0]]
        hdr.append(key[0] + '>' + key[1])
        fvals.append(str(keyfrac))
        
        #ostream.write('%s\t%s>%s\t%d\t%s\n' % (bsub, key[0], key[1],
        #                                       keycount, keyfrac))
    bsubs = basename + '_Subsequent'
    bsuba = basename + '_Antecedent'
    for key in sorted(bc.dicounts.keys()):
        if 'N' in key:
            continue
        keycount = bc.dicounts[key]
        hdr.append(key)
        if key[0] == key[3]:
            keyfrac = keycount / bc.basecounts[key[1]]
            bsub = bsuba
        elif key[1] == key[4]:
            keyfrac = keycount / bc.basecounts[key[0]]
            bsub = bsubs
        fvals.append(str(keyfrac))
        #ostream.write('%s\t%s\t%d\t%s\n' % (bsub, key, keycount, keyfrac))
    ostream.write('\t'.join(hdr) + '\n')
    ostream.write('\t'.join(fvals) + '\n')
    if ostream != sys.stdout:
        ostream.close()

                                               
    

