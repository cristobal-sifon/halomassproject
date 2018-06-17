#!/usr/bin/env python
import pylab
import readfile
import scipy

from matplotlib import rcParams
rcParams['font.size'] = 16

import locale
locale.setlocale(locale.LC_NUMERIC, 'C')

from astLib import astCoords

# my modules
import pytools
import stattools
from astro import cosmology

c = 299792.458

def main():
  #output_gal = make_output_gal()
  output_halo = make_output_halo()

  filename = 'phase1/cat/phase1_HOD-true.txt'
  file = open(filename)
  while True:
    try:
      id, center, ngal, n200, z, s, m, r = run_halo(file)
      line_halo = '%4d  %.5f  %.5f  %5d  %8.1e  %4d' \
                  %(id, center[0]*scipy.pi/180, center[1]*scipy.pi/180,
                    round(c*z, 0), m[0], round(s, 0))
      line_halo += '  %4.1f  %3d  %5.2f  %5.2f  %4d' \
                   %(min(r,r/1e3), ngal, m[1][0], m[1][1], n200)
      hout = open(output_halo, 'a')
      print >>hout, line_halo
      hout.close()
      #print line_halo
    except IndexError:
      break
  return

def run_halo(file):
  halo = [int(x) for x in file.readline().split()]
  ra = scipy.zeros(halo[1])
  dec = scipy.zeros(halo[1])
  z = scipy.zeros(halo[1])
  mag = scipy.zeros(halo[1])
  for i in xrange(halo[1]):
    line = file.readline().split()
    ra[i] = float(line[0])
    dec[i] = float(line[1])
    z[i] = float(line[2])
    mag[i] = float(line[3])
  ra *= 180/scipy.pi
  dec *= 180/scipy.pi
  z /= c
  zo = stattools.Cbi(z)
  bcg = scipy.argmin(mag)
  center = (ra[bcg], dec[bcg])
  dist = astCoords.calcAngSepDeg(center[0], center[1], ra, dec)
  rproj = 1e3 * scipy.array(map(cosmology.dProj, zo*scipy.ones(len(z)), dist))
  failed =  (63, 86, 202, 276, 396, 401, 558,
             625, 653, 655, 665, 676, 693, 836)
  if halo[0] in failed:
    print ''
    print 'Halo %d' %halo[0]
    verbose = True
  else:
    verbose = False
  zo, s, m200, r200, n200 = pytools.M200(rproj, z, errors=True,
                                         converge=False, verbose=verbose)
  outpath = 'phase1/plots/true-members/rv/'
  if n200 == -1:
    outpath += 'failed/'
  plot_rv(halo[0], rproj, z, zo, n200, m200, r200, outpath)
  return halo[0], center, halo[1], n200, zo, s, m200, r200

def make_output_gal():
  output = 'phase1/ret/HMRC_Galaxies_HOD1_SGF-true'
  out = open(output, 'w')
  print >>out, 'ID  RA  Dec  Vres  M200  Vdisp  R200  N200  Merr(stat)',
  print >>out, ' Merr(syst)  Ntot'
  out.close()
  return output

def make_output_halo():
  output = 'phase1/ret/HMRC_Halo_HOD1_SGF-true'
  out = open(output, 'w')
  print >>out, 'Halo  RA  Dec  Vres  M200  sigma200  r200  Ngal ',
  print >>out, ' Merr_stat  Merr_syst  N200'
  out.close()
  return output

def plot_rv(halo, rproj, z, zo, n200, m200, r200, outpath,
            ylim=(-5000,5000)):
  v = c * (z-zo) / (1+zo)
  dm = [m200[0]*(10**me-1) for me in m200[1]]
  dm = scipy.hypot(*dm)
  output = outpath + 'halo-%d.png' %halo
  fig = pylab.axes([0.1, 0.1, 0.85, 0.85])
  pylab.axhline(0, ls='-', color='k')
  pylab.plot(rproj, v, 'b.')
  msg = 'Halo %d\n$z=%.3f$\n' %(halo, zo)
  if n200 == -1:
    msg += 'Failed'
  else:
    pylab.axvline(r200, ls='--', color='k')
    msg += r'$M_{200}=%.1f\pm%.1f\times10^{14}M_\odot$' \
           %(m200[0]/1e14, dm/1e14)
  pylab.annotate(msg, xy=(0.55,0.95), xycoords='axes fraction', va='top')
  for y in xrange(ylim[0]+1000, ylim[1], 1000):
    pylab.axhline(y, ls=':', color='k')
  pylab.ylim(*ylim)
  pylab.savefig(output, format=output[-3:])
  pylab.close()
  return

main()
  