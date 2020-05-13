from bisect import bisect
import sys


class BedFileChecker():

    def __init__(self, bedfile):
        self.bedfile = bedfile 
        self.ld_bed_file = self.load_bed()

    def load_bed(self):
        d = {}
        for line in open(self.bedfile):
            if 'track' in line:
                continue
            chrom, p1, p2 = line.split(None, 3)[:3]
            if 'chr' in chrom:
                chrom = chrom[3:]

            try:
                pos = d[chrom]
            except KeyError:
                pos = []
                d[chrom] = pos
            p1 = int(p1)
            p2 = int(p2) - 1 #bed file half open
            pos.append((p1, p2))

        for key in d:
            v = d[key]
            v.sort()
            d[key] = zip(*v)
        return d

    def inbed(self, chrom, pos):
        pos = int(pos)
        bedflag = False
        if chrom not in self.ld_bed_file:
            return False
        else:
            s, e = self.ld_bed_file[chrom]
            idx = bisect(s, pos)
            if idx == 0:
                return False
            else:
                if pos <= e[idx-1]:
                    return True

        return bedflag



if __name__ in "__main__":
    bedfile = sys.argv[1]
    chrom = sys.argv[2]
    pos = sys.argv[3]
    b = BedFileChecker(bedfile)
    print b.inbed(chrom, pos)
    


