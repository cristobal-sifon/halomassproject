#!/usr/bin/env python
import glob
import pylab
import readfile
import scipy
import stattools

import locale
locale.setlocale(locale.LC_NUMERIC, 'C')

from astLib import astCoords
from astro import cosmology, members

no = 6

path = '/Users/cristobal/Documents/act/sifon+2012/paper/ref_review/'
# relevant files
filename = path + 'catalogs/'
ls = sorted(glob.glob(filename + '*.cat'))
cluster = ls[no].split('/')[-1][:-4]
filename += cluster + '.cat'
dynfile = path + 'dynamics.dat'

# read redshift
zo = readfile.table(dynfile, cols=2, include=cluster)
print zo

# read BCG coordinates
bcg = open(filename).readline().split()[1:]
bcg = [float(x) for x in bcg]

data = readfile.table(filename,
                      cols=(0,1,2,3),
                      dtype=(int,float,float,float))
obj = data[0]
ra = data[1]
dec = data[2]
z = data[3]
Ngal = len(obj)

dist = astCoords.calcAngSepDeg(bcg[0], bcg[1], ra, dec)
r = 2 * 1e3 * scipy.array(map(cosmology.dProj, zo*scipy.ones(Ngal), dist))

m = members.sgapper(r, z, maingap=500, min_binsize=250, converge=False,
                    plot_output=True, full_output=True)
m, binsize, nit, incut = m
binloc = [sum(binsize[:i+1]) for i in xrange(len(binsize))]

c = 299792.458
zo = stattools.Cbi(z[m])
v = c * (z-zo) / (1+zo)
s = c * stattools.Sbi(z[m]) / (1+zo)

print cluster
print zo
print s

exit()

for i in xrange(-5000, 5000, 500):
  pylab.axhline(i, ls=':', color='k')
pylab.plot(r, abs(v), 'wo')
pylab.plot(r[m], abs(v[m]), 'ko')
# cheating
pylab.plot(r[abs(v) < 2500], abs(v[abs(v) < 2500]), 'rx')
pylab.axhline(0, ls='--', color='k')
for bin in binloc:
  pylab.axvline(bin, ls=':')
pylab.ylim(0, 5000)
pylab.show()