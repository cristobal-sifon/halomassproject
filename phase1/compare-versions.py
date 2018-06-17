#!/usr/bin/env python
import pylab
import readfile
import scipy

sgf = 'returns/LC1Halo_SGF'
sgf1 = readfile.table(sgf+'_1')
sgf2 = readfile.table(sgf+'_2')

mratio = sgf1[4] / sgf2[4]

bins = scipy.arange(0, 5.01, 0.05)
pylab.hist(mratio, bins=bins)
pylab.show()