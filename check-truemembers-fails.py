#!/usr/bin/env python
import pylab
import readfile
import scipy

def main():
  data = readfile.table('phase1/ret/HMRC_Halo_HOD1_SGF-true')
  Nh = len(data[0])
  print 'Total %d clusters' %Nh
  failed = scipy.arange(Nh)[data[-1] == -1]
  print '%d clusters failed' %len(failed)
  catalog = 'phase1/cat/phase1_HOD-true.txt'
  for line in open(catalog):
    line = line.split()
    if len(line) == 2:
      id = int(line[0])
      n = int(line[1])
      ra = []
      dec = []
      v = []
      mabs = []
      color = []
    else:
      if id in failed:
        ra.append(float(line[0]))
        dec.append(float(line[1]))
        v.append(float(line[2]))
        mabs.append(float(line[3]))
        color.append(float(line[4]))
    # once we've read all members from a halo
    # note that this will only happen for clusters that failed
    if len(ra) == n:
      
    