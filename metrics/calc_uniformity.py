import numpy as np
import multiprocessing as mp
import gzip
#import matplotlib.pyplot as plt
import glob
import re
import os
import sys
from scipy.stats import linregress
from scipy.interpolate import interp1d

def get_stats(filename):
    with open(filename, 'rt') as stream:
        a = []
        for line in stream:
            if line[:1] == '#': continue
            a.append(float(line.split('\t')[2]))
        a = np.array(a)
        mean = a.mean()
        q1, median, q3 = np.percentile(a, (25, 50, 75))
        numpoints = len(a)
        fractions = []
        sname = os.path.basename(filename).split('.')[0].split('_')[0]
        for i in (1, 2, 3, 4, 5, 6, 7, 8, 9):
            ix = mean * i / 10.0
            low = mean - ix
            hgh = mean + ix
            tmp = a[np.where(a >= low)[0]]
            inrange = len(np.where(tmp < hgh)[0])
            fractions.append(inrange * 1.0 / numpoints)
    return sname, mean, median, q1, q3, fractions

files = glob.glob(sys.argv[1] + '*basecov.tsv')

p = mp.Pool(20)
cov_stats = p.map(get_stats, files)

sdata = []
for name, r in zip(files, cov_stats):
    bname = os.path.basename(name)
    sdata.append(r)
x = list(range(10, 100, 10))
hdr = '#sample mean median q1 q3 10 20 30 40 50 60 70 80 90 lr0 lr1 perc_mean_fr0.5 iqr_med'.split()
print('\t'.join(hdr))

for item in sdata:
    t = list(item[:5])
    for p in item[5]:
        t.append(p)
    lr = linregress(x, item[-1])
    t.extend([str(lr[0]), str(lr[1])])
    inter = interp1d(item[-1], x)
    t.append(str(inter(0.5)))
    t.append(str((float(item[4]) - float(item[3]))/float(item[2])))
    t = [str(el) for el in t]
    print('\t'.join(t))

