"""
.. module:: ostats
    :platform: Unix, Windows, MacOSX
    :synopsis: Compute basic statistics for a series of integers in
               an online fashion

Compute basic statistics for a series of integers in an online fashion

"""
from __future__ import division
from collections import Counter
import math


class OStats(object):
    __slots__ = ['_c', '_mean', '_var', '_n', '_ndp']

    """Compute mean, median and mode on streaming (online) data

    A typical usage for doing bed coverage woule be as follows::

    >>> from bed import readfile
    >>> import pysam
    >>> bam = pysam.Samfile('somefile.bam', 'rb')
    >>> gstats = Ostats()
    >>> for rec in readfile('somefile.bed'):
    ...     reads = bam.fetch(rec.chrom, rec.start, rec.end)
    ...     counts = [0] * (rec.end - rec.start)
    ...     for pos in rec.positions:
    ...         if rec.start <= pos < rec.end:
    ...             counts[pos - rec.start] += 1
    ...     cs = Ostats(counts)
    ...     print rec.chrom, rec.start, rec.end, cs.mean, c.median, cs.mode
    ...     gstats.update(cs)
    >>>  print 'Global Stats:', gstats.mean, gstats.median, gstats.mode
    """

    def __init__(self, data=None):
        self._c = Counter()
        self._mean = self._var = 0.0
        self._n = self._ndp = 0
        if data:
            self.update(data)

    def update(self, d):
        """Add data for computation

        Args:
            d : The data to be added. Can be a single integer,
                a sequence (list/tuple/...) of integers, a dictionary
                of integer:count, or an instance of this class
        """
        if isinstance(d, self.__class__):
            self._c.update(d._c)
        else:
            self._c.update(d)
        self._update_mean(d)


    def _update_mean(self, d):
        if isinstance(d, self.__class__):
            self._update_mean_dict(d._c.items())
        elif isinstance(d, dict):
            self._update_mean_dict(d.items())
        elif isinstance(d[0], int):
            self._update_mean_seq(d)
        else:
            self._update_mean_dict(d)

    def _update_mean_dict(self, d):
        for key, value  in d:
            temp = self._n + value
            delta = key - self._mean
            R = delta * value / temp
            self._mean += R
            self._var += self._n * delta * R
            self._n = temp
        self._ndp += len(d)

    def _update_mean_seq(self, d):
        for el in d:
            self._n += 1
            delta = el - self._mean
            self._mean += delta / self._n
            self._var += delta * (el - self._mean)

    @property
    def nobs(self):
        """Number of observations"""
        return self._n

    @property
    def mean(self):
        """Mean of the observations"""
        return self._mean

    @property
    def var(self):
        """Variance of the observations"""
        if self._ndp > 1 and self._n:
            return (self._var / self._n) * (self._ndp / (self._ndp - 1))
        else:
            try:
                return self._var / (self._n - 1)
            except ZeroDivisionError:
                return 0.0

    @property
    def sdev(self):
        """Standard Deviation"""
        return math.sqrt(self.var)

    @property
    def median(self):
        """Median of the observations"""
        d = self._n / 2
        n = 0
        for i, j in sorted(self._c.items()):
            n += j
            if n > d:
                return i
        return i

    @property
    def mode(self):
        """Compute the mode(s) for the integers"""
        mc = None
        modes = []
        for i, j in self._c.most_common():
            if mc is None:
                mc = j
                modes.append(i)
            elif j == mc:
                modes.append(i)
            else:
                break
        return modes

    @property
    def minobs(self):
        return min(self._c.keys())

    @property
    def maxobs(self):
        return max(self._c.keys())

    @property
    def mincount(self):
        return min(self._c.values())

    @property
    def maxcount(self):
        return self._c.most_common(1)[0]

    def binned_mode(self, binsize=10):
        d = Counter()
        for key, value in self._c.items():
            d[binsize * (key // binsize)] += value
        return d.most_common(1)[0][0] + binsize/2

    def binned_counts(self, binvalues, asfraction=False):
        d = dict([(v, 0) for v in binvalues])
        for k, v in self._c.items():
            for b in binvalues:
                if k >= b:
                    d[b] += v
                else:
                    break

        if asfraction:
            denom = 1.0 / self.nobs
            for v in d:
                d[v] = (d[v], d[v] * denom)
        return sorted(d.items())

    def cumulative_fraction(self, f):
        denom = 1.0 / self._n
        nf = 0.0
        v = []
        for i, j in sorted(self._c.items()):
            nf += j * denom
            v.append((i, nf))
            if nf >= f:
                break
        return v

if __name__ == "__main__":
    c = OStats()
    for i in range(10):
        d = OStats(range(i, i + 10))
        print d.nobs, d.mean, d.median, d.sdev
        c.update(d)
    print c.nobs, c.mean, c.median, c.sdev
    v = []
    for i in range(10):
        v.extend(range(i, i+10))
    v.sort()
    n = len(v)
    mean = sum(v)/len(v)
    mid, rem = divmod(len(v), 2)
    if rem:
        median = v[mid]
    else:
        median = (v[mid - 1] + v[mid]) / 2

    var = 0
    for el in v:
        var += (el - mean) ** 2
    print n, mean, median, math.sqrt(var / (n - 1))


