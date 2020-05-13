#!python3
# This code is part of the Genome Commons Navigator (GCN) distribution and
# governed by its license. Please see the LICENSE file that should have been
# included as part of this package.
#
import tempfile
import os
import multiprocessing
import sys
from collections import defaultdict

import pysam

import anyopen3


def gcov(bamfile, nproc, outname, mqmin=1, bqmin=1, bedfile=None):

    if bedfile is None:
        ivls = [(bamfile, part, mqmin, bqmin) for part in
                partition_intervals(get_intervals(bamfile), nproc)]
    else:
        ivls = [(bamfile, part, mqmin, bqmin) for part in
                partition_intervals(load_intervals(bedfile), nproc)]

    p = multiprocessing.Pool(nproc)
    r = p.map(coverage, ivls)
    if outname == '-':
        ostream = sys.stdout
    else:
        ostream = open(outname, 'wt')
    for fs in r:
        with open(fs, 'rt') as stream:
            for line in stream:
                ostream.write(line)
        os.unlink(fs)
    if ostream is not sys.stdout:
        ostream.close()

    sys.exit(0)


def get_intervals(bamfile):
    intervals = []
    with pysam.Samfile(bamfile, 'rb') as bamstream:
        sq = bamstream.header['SQ']
        for item in sq:
            intervals.append((item['SN'], 0, item['LN']))
    return intervals


def load_intervals(bedfile):
    with anyopen3.openfile(bedfile, 'rb') as stream:
        intervals = []
        for line in stream:
            try:
                chrom, start, end = line.split(b'\t', 3)[:3]
                #chrom = chrom.replace(b'chr', b'')
                intervals.append((chrom.decode(), int(start), int(end)))
            except ValueError:
                pass
    sys.stderr.write('Loaded %d intervals\n' % len(intervals))
    return intervals


def partition_intervals(intervals, nproc):
    bases_per_proc, remainder = divmod(countbases(intervals), nproc)
    partitions = [[] for _ in range(nproc)]
    nbases = 0
    curpart = 0
    for ivl in intervals:
        ilen = ivl[2] - ivl[1]
        ntot = nbases + ilen
        if ntot <= bases_per_proc:
            partitions[curpart].append(ivl)
            nbases = ntot
        else:
            excess = ntot - bases_per_proc
            jvl = ivl[1] + ilen - excess
            partitions[curpart].append([ivl[0], ivl[1], jvl])
            while excess:
                curpart += 1
                if curpart >= nproc:
                    curpart = nproc - 1
                if excess < bases_per_proc:
                    kvl = jvl + excess
                    partitions[curpart].append([ivl[0], jvl, kvl])
                    excess = 0
                    nbases = excess
                else:
                    kvl = jvl + bases_per_proc
                    partitions[curpart].append([ivl[0], jvl, kvl])
                    excess -= bases_per_proc
                    jvl = kvl
    return partitions


def countbases(intervals):
    nb = 0
    for _, start, end in intervals:
        nb += end - start
    sys.stderr.write('Total number of bases = %d\n' % nb)
    return nb


def coverage(ivls):
    bamfile, intervals, mqmin, bqmin = ivls
    fd, tmpfilename = tempfile.mkstemp(dir=os.getcwd())

    with pysam.Samfile(str(bamfile), 'rb') as bamstream:
        chromosomes = set([el['SN'] for el in bamstream.header['SQ']])
        for chrom, start, end in intervals:
            end = int(end)
            start = int(start)
            bam_has_chrom = False
            if chrom in chromosomes:
                bam_has_chrom = True
            else:
                if chrom == 'M':
                    if 'MT' in chromosomes:
                        chrom = 'MT'
                        bam_has_chrom = True
                elif chrom == 'MT':
                    if 'M' in chromosomes:
                        chrom = 'M'
                        bam_has_chrom = True
            if not bam_has_chrom:
                s = '%s\t%d\t%d\t0\n' % (chrom, start, end)
                os.write(fd, s.encode('ascii'))
            else:
                for i in range(start, end, 10000):
                    istart = i
                    iend = min(i + 10000, end)
                    counts = [0] * (iend - istart)
                    for arec in bamstream.fetch(chrom, istart, iend):
                        if (not (arec.is_duplicate or
                                arec.is_qcfail or
                                arec.is_unmapped or
                                arec.is_secondary or
                                (arec.mapping_quality < mqmin))) and arec.is_proper_pair:
                            for pos in arec.positions:
                                if istart <= pos < iend:
                                    counts[pos - istart] += 1
                    report(chrom, istart, counts, fd)
    os.close(fd)
    return tmpfilename


def report(chrom, start, counts, fd):
    fmt = '%s\t%d\t%d\n'
    for idx, value in enumerate(counts):
        os.write(fd, (fmt % (chrom, start + idx, value)).encode('ascii'))


if __name__ == "__main__":
    import argparse

    desc = 'Compute coverage per position for intervals in a BED file'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-o', '--output', dest='outname', type=str,
                        default='-', help='file to write results to')
    parser.add_argument('-n', '--cpus', dest='cpus', type=int, default=1,
                        help='Number of cpus to use')
    parser.add_argument('-b', '--bedfile', dest='bedfile', type=str,
                        help='BED file of intervals for which coverage '
                             'is to be computed')
    parser.add_argument('-q', '--mqmin', dest='minmq', type=int, default=1,
                        help='Minimum Mapping quality')
    parser.add_argument('-Q', '--bqmin', dest='minbq', type=int, default=1,
                        help='Minimum base quality')
    parser.add_argument('bamfile', nargs='+')

    args = parser.parse_args()
    inbed = args.bedfile or None
    gcov(args.bamfile[0], args.cpus, args.outname, args.minmq, args.minbq,
         inbed)

