#!/usr/bin/env python
import pylab
import readfile

path = '/Users/cristobal/Documents/act/sifon+2012/paper/ref_review/'
filename = path + 'dynamics.dat'
output = 'ngal-hist.png'

data = readfile.table(filename, cols=(0,1,2), dtype=(str,int,float))
clusters = [cl[4:9] for cl in data[0]]
Ngal = data[1]
z = data[2]
Ncl = len(clusters)
index = pylab.arange(Ncl, dtype=int)
median = pylab.median(Ngal)

fig = pylab.axes([0.1, 0.13, 0.86, 0.83])
pylab.bar(index, Ngal, fc='b', align='center')
pylab.axhline(median, ls='--', color='r')
pylab.annotate('median = %d' %int(round(median, 0)),
               xy=(11,median), color='r', va='bottom')
fig.set_xticks(index)
fig.set_xticklabels(clusters, rotation=40, fontsize='10')
pylab.ylabel('Ngal')
pylab.xlim(-1, Ncl)
pylab.ylim(0, 100)
pylab.savefig(output, format=output[-3:])
print 'Saved to', output