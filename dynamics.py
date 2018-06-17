#!/usr/bin/env python
import glob
import os
import pylab
import readfile
import stattools
import scipy
import sys
from astro import cosmology
from astro.clusters import conversions, members, scalings
from itertools import count, izip

import locale
locale.setlocale(locale.LC_NUMERIC, 'C')

from astLib import astCoords

cosmology.Omega_M = 0.27
cosmology.Omega_L = 0.73
cosmology.h = 0.7

c = 299792.458
deg2rad = scipy.pi/180

def main(version='1'):
  out = open(output, 'w')
  head = '# halo-id  ra  dec  vres  m200  vdisp200  r200  N200 '
  head += ' M200err(stat)  M200err(syst)  Ngal'
  print >>out, head
  out.close()
  mout = open(mfile, 'w')
  head = '# halo-id  galaxy-id  within-r200?'
  print >>mout, head
  mout.close()
  
  par = parameters(version)
  min_Nperbin, min_binsize, maingap, maxgap = par
  masses = [run(center, min_Nperbin, min_binsize, maingap, maxgap,
                version=version) for center in centers]
  print 'Saved to', output
  masses = scipy.array(masses)
  # show histogram
  #bins = scipy.arange(12, 16.01, 0.1)
  #pylab.title(output)
  #pylab.hist(scipy.log10(masses[masses > 0]), bins=bins,
             #histtype='step', color='r')
  #pylab.xlabel('log M')
  #pylab.ylabel('N')
  #pylab.show()
  return

def run(center, min_Nperbin=15, min_binsize=250,
        maingap=500, maxgap=1000, sigma=False, version='1'):
  #center = centers[i]
  id = center[0]
  xo = center[1]
  yo = center[2]
  zo = center[3]

  j = close(xo, yo, zo, ra, dec, z)
  r = astCoords.calcAngSepDeg(xo, yo, ra[j], dec[j])
  r = scipy.array([cosmology.dProj(zo, ri, input_unit='deg', unit='Mpc') \
                   for ri in r])
  if '2' in version:
      r *= (1 + zo)
  r *= 1e3
  outplot = get_plotname(version, id)
  #output = 'phase%s/plots/v%s/rv/halo-%d.png' %(phase, version, id)
  m = members.sgapper(r, z[j], zo=zo, min_Nperbin=min_Nperbin,
                      maingap=maingap, maxgap=maxgap, sigma=sigma,
                      verbose=False, debug=False, plot_output=outplot,
                      full_output=True)
  m, binsize, nit, incut = m
  Nm = len(m)

  zo = stattools.Cbi(z[j[m]])
  s, m200, r200 = mass(z[j[m]], zo)
  mo = m[r[m] <= r200]
  n200 = -1
  it = []
  while n200 != len(mo):
    try:
      n200 = len(mo)
      s, m200, r200 = mass(z[j[mo]], zo)
      #zo = stattools.Cbi(z[j[mo]])
      mo = m[r[m] <= r200]
      # should maybe find a more sophisticated solution by, say, averaging
      # all items from one of these "cycles"
      if r200 in it:
        break
      it.append(r200)
    except ZeroDivisionError:
      print 'broke'
      break

  if n200 > 15:
    stats = (stattools.Cbi, stattools.Sbi)
    z_err, s_err = stattools.bootstrap(stats, z[j[incut]],
                                       n_obj=Nm, n_samples=1000)
    s_err = c * s_err / (1+zo)
  else:
    z_err = scipy.std(z) / scipy.sqrt(Nm)
    s_err = c/(1+zo) * z_err / scipy.sqrt(2)
  m200, m200_err = scalings.sigma(s, zo, s_err, z_err,
                                  separate_errors=True)
  r200 = conversions.rsph(m200, zo)

  # to make it their format
  m200_err = [scipy.log10(1+me/m200) for me in m200_err]

  try:
    #if n200 > 
    line = '%4d  %.3f  %.3f  %6d  %.1e  %4d  %.1f  %4d  %.2f  %.2f  %4d' \
           %(id, xo*deg2rad, yo*deg2rad, round(c*zo, 0), m200,
             round(s, 0), r200/1e3, n200, m200_err[0], m200_err[1], Nm)
  except TypeError:
    line = '%4d  %.3f  %.3f  %6d     -1      -1   -1    -1    -1    -1  %4d' \
           %(id, xo*deg2rad, yo*deg2rad, round(c*zo, 0), Nm)
  print line, '%4.1f' %max(r[m]/1e3)
  out = open(output, 'a')
  print >>out, line
  out.close()

  write_members(id, galid[j[m]], galid[j[mo]])
  return m200

def close(xo, yo, zo, ra1, dec1, z1):
  """
  conditions: within 1 deg and 5,000 km/s
  """
  j = scipy.arange(len(ra), dtype=int)[(abs(ra1 - xo) < 1) & \
                                       (abs(dec1 - yo) < 1) & \
                                       (c * abs(z1-zo) / (1+zo) < 5e3) & \
                                       (z1 > 0)]
  return j

def dProj(z, r):
  return cosmology.dProj(z, r, input_unit='deg', unit='kpc')

def mass(z, zo):
  s = c * stattools.Sbi(z) / (1+zo)
  # using z=0 as the cosmological redshift
  if '--z=0' in sys.argv:
    zo = 0
  m = scalings.sigma(s, zo)
  r = conversions.rsph(m, zo, ref='200c', unit='kpc')
  return s, m, r

def parameters(version):
  if version == '1':
    return 15, 250, 500, 1000
  elif version == '2':
    return 10, 150, 300, 500
  return

def write_members(id, members, mem200):
    def print_line(mi):
        if mi in mem200:
            print >>out, '%4d  %6d  1' %(id, mi)
        else:
            print >>out, '%4d  %6d  0' %(id, mi)
        return
    out = open(mfile, 'a')
    map(print_line, members)
    out.close()
    return

def get_output(galaxyfile, version='1'):
    output = galaxyfile.replace('cat/', 'ret/')
    output = output.split('/')
    nameparts = output[-1].split('_')
    if '1' in phase:
        output[-1] = 'HMRC_Halo_%s_%s_SGF%s.txt' \
                     %(nameparts[1], nameparts[0], version)
    else:
        if phase == '2h':
            output[-1] = 'HMRC_Halo_HOD2_SW_SG%s.txt' %version
        elif phase == '2s':
            output[-1] = 'HMRC_Halo_SAM2_SW_SG%s.txt' %version
        elif phase == '2d':
            output[-1] = 'HMRC_Halo_SAM2_DN_SG%s.txt' %version
    return '/'.join(output)

def get_plotname(version, id):
    if phase in ('2s', '2h'):
        path = 'phase2/shallow/'
    elif phase == '2d':
        path = 'phase2/deep/'
    else:
        path = 'phase1/'
    output = path + 'plots/v%s/rv/halo-%d.png' %(version, id)
    return output

version = '1'
if '-v' in sys.argv:
  i = sys.argv.index('-v')
  version = sys.argv[i+1]

# get files
# which phase?
if '--shallow' in sys.argv:
    if '--hod' in sys.argv:
        phase = '2h'
        catnames = 'phase2/shallow/cat/SW_HOD2_Galaxies.txt'
    else:
        phase = '2s'
        catnames = 'phase2/shallow/cat/SW_SAM2_Galaxies.txt'
elif '--deep' in sys.argv:
    phase = '2d'
    catnames = 'phase2/deep/cat/DN_SAM2_Galaxies.txt'
elif '--1b' in sys.argv:
    phase = '1b'
    catnames = 'phase1/cat/SW_HOD1b_Galaxies.txt'
else:
    phase = '1'
    catnames = 'phase1/cat/SW_HOD1_Galaxies.txt'
galaxyfiles = sorted(glob.glob(catnames))
halofiles = [i.replace('Galaxies', 'Haloes') for i in galaxyfiles]
outputs = [get_output(i, version=version) for i in galaxyfiles]
mfiles = [i.replace('Halo', 'Galaxies') for i in outputs]
#print outputs
#exit()

for gf, cf, output, mfile in izip(galaxyfiles, halofiles, outputs, mfiles):
  print '  Will save to', output
  data = readfile.table(gf, dtype=(int,float,float,float,float,float))
  galid = data[0]
  ra = data[1] / deg2rad
  dec = data[2] / deg2rad
  z = data[3] / c
  centers = readfile.table(cf, dtype=(int,float,float,float))
  centers[1] /= deg2rad
  centers[2] /= deg2rad
  centers[3] /= c
  centers = scipy.transpose(centers)
  main(version=version)
