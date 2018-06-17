#!/usr/bin/env python
import glob
import os
import pylab
import readfile
import scipy
import sys
from astLib import astCoords
from astro import cosmology

"""
Design MOS-like observations of clusters.

Follow the code created for the project with Nick.

Three setups: random, "ACT-south", and "ACT-eq".

Two galaxy samples: red and random.

"""

def main():
    
    return

def major_axis(x, y, ndist=10, mask=1):
    dist = scipy.array([scipy.hypot(x-i,y-j) \
                        for i, j in itertools.izip(x, y)])
    shape = dist.shape
    imax = scipy.zeros(ndist, dtype=int)
    jmax = scipy.zeros(ndist, dtype=int)
    # dummy run
    for i in xrange((mask-1)*ndist):
        ii, jj = scipy.unravel_index(scipy.argmax(dist), shape)
        dist[ii] = scipy.zeros(shape[0])
        dist[jj] = scipy.zeros(shape[0])
        dist[:,ii] = scipy.zeros(shape[1])
        dist[:,jj] = scipy.zeros(shape[1])
    for i in xrange(ndist):
        imax[i], jmax[i] = scipy.unravel_index(scipy.argmax(dist), shape)
        dist[imax[i]] = scipy.zeros(shape[0])
        dist[jmax[i]] = scipy.zeros(shape[0])
        dist[:,imax[i]] = scipy.zeros(shape[1])
        dist[:,jmax[i]] = scipy.zeros(shape[1])
    #print imax, jmax
    slopes = [(y[i]-y[j])/(x[i]-x[j]) \
              for i, j in itertools.izip(imax, jmax)]
    zeros = [y[i] - m*x[i] for i, m in itertools.izip(imax, slopes)]
    m = scipy.median(slopes)
    n = scipy.median(zeros)
    return imax, jmax, m, n

def make_mask(X, Y, v, fov=5, slit_size=None, ngals=40,
               width=1., fill=False):
    ntot = len(X)
    imax, jmax, m, n = major_axis(X, X, mask=1)
    observed = observe_mask(X, Y, m, n, fov=fov,
                            slit_size=slit_size,
                            ngals=ngals, width=width)
    return observed

def observe_mask(X, Y, m, n, fov=5,
                 slit_size=None, ngals=40, width=1):
    """
    can either give the slit size, or the number of objects. If the number
    of objects is given, the slit size parameter is ignored, and instead the
    slit size is calculated to have that number of objects in one observation.
    """
    ngal = len(X)
    j = scipy.arange(ngal, dtype=int)
    dist = scipy.hypot(X, Y)
    # the BCG is always the first element (bcg = 0)
    bcg = j[dist == 0][0]
    # rotate the coordinates to the mask direction
    theta = scipy.arctan(m)
    Xrot = X*scipy.cos(theta) - Y*scipy.sin(theta)
    Yrot = X*scipy.sin(theta) + Y*scipy.cos(theta)
    pylab.plot(Xrot, Yrot, 'k.', mew=2)
    # bin the galaxies in Xrot
    in_fov = j[(abs(Xrot) < fov[0]/2.) & (abs(Yrot) < fov[1]/2.)]
    Xrot = Xrot[in_fov]
    Yrot = Yrot[in_fov]
    n_in_fov = len(in_fov)
    j = scipy.arange(n_in_fov, dtype=int)
    if ngals:
        slit_size = fov[0]/float(ngals)
        Xbins = scipy.linspace(-fov[0]/2., fov[0]/2., ngals+1)
    else:
        Xbins = scipy.arange(-fov[0]/2., fov[0]/2., slit_size/60.)
    Xbins += slit_size / 2.
    Xbinned = scipy.digitize(Xrot, Xbins)
    # observe one galaxy per Xbin with a Gaussian probability with width
    # defined below (in arcmin), and excluding the bin containing the
    # BCG, which is necessarily observed (we always do!).
    # This width is such that I preferentially observe galaxies near the
    # center of the image.
    observed = [j[Xbinned == i][scipy.argmin(abs(Yrot[Xbinned == i] -
                                                 random.normal(0, width)))] \
                for i in xrange(len(Xbins)-1) \
                if (len(Yrot[Xbinned == i]) > 0) & (i != Xbinned[0])]
    observed = scipy.append(bcg, observed)
    # plot (just once, to show -- and make sure!)
    #for x in Xbins:
        #pylab.axvline(x, ls='-', color='0.7')
    #pylab.plot(0, 0, 'o', ms=8, mfc='orange', mec='orange')
    #pylab.plot(Xrot[observed], Yrot[observed], 'rx', ms=6, mew=2)
    ## field of view
    #pylab.plot([-fov[0]/2., -fov[0]/2.], [-fov[1]/2., fov[1]/2.], 'b-', lw=2)
    #pylab.plot([fov[0]/2., fov[0]/2.], [-fov[1]/2., fov[1]/2.], 'b-', lw=2)
    #pylab.plot([-fov[0]/2., fov[0]/2.], [-fov[1]/2., -fov[1]/2.], 'b-', lw=2)
    #pylab.plot([-fov[0]/2., fov[0]/2.], [fov[1]/2., fov[1]/2.], 'b-', lw=2)
    #pylab.xlabel('rotated x (arcmin)')
    #pylab.ylabel('rotated y (arcmin)')
    #pylab.xlim(-6, 6)
    #pylab.ylim(-6, 6)
    #output = 'plots/mask1_sample.png'
    #pylab.savefig(output, format=output[-3:])
    #pylab.close()
    #print 'Saved to', output
    #exit()
    return in_fov[observed]

def observe_mos(hostids, bcgs, j, x, y, v, zobs=0.5, fov=(10,5),
                slit_size=None, ngals=40, width=1., fill=False, nmin=100,
                mbm=True, iterate=True):
    ngr = len(bcgs)
    # doing the last one separately should be faster since the code
    # doesn't have to check the condition below every time
    output = 'outputs_%s/' \
             %(''.join(str(datetime.today()).split()[0].split('-')))
    if fov is None:
        output += 'dynamics-fullobs.out'
    else:
        output += 'equator'
        output += '-nmin_%d-zobs_%.2f-fov_%.1f_%.1f' \
                  %(nmin, zobs, fov[0], fov[1])
        if slit_size is None:
            output += '-slitsize_None'
        else:
            output += '-slitsize_%d' %slit_size
        if ngals is None:
            output += '-ngals_None'
        else:
            output += '-ngals_%d' %ngals
        output += '-width_%.1f-fill_%s.out' %(width, fill)
    if mbm:
        output = output.replace('.out', '_MBM10.out')
    if not iterate:
        output = output.replace('.out', '_noiter.out')
    print 'Will save to', output
    hdr = '# hostid  ntot  nobs  n200'
    hdr += '  s200  s200_err  m200  m200_err  r200  r200_err'
    out = open(output, 'w')
    print >>out, hdr
    aux = [one_group(out, hostids[i],
                     x[bcgs[i]:bcgs[i+1]], y[bcgs[i]:bcgs[i+1]],
                     v[bcgs[i]:bcgs[i+1]], zobs=zobs, fov=fov,
                     slit_size=slit_size, ngals=ngals, width=width,
                     fill=fill, mbm=mbm, iterate=iterate) \
           for i in j[:-1]]
    # the last halo, in case it is the last one on the list
    i = j[-1]
    if i == ngr - 1:
        one_group(out, hostids[i],
                  x[bcgs[i]:], y[bcgs[i]:], v[bcgs[i]:],
                  zobs=zobs, fov=fov, slit_size=slit_size, ngals=ngals,
                  width=width, fill=fill, mbm=mbm, iterate=iterate)
    else:
        one_group(out, hostids[i],
                  x[bcgs[i]:bcgs[i+1]], y[bcgs[i]:bcgs[i+1]],
                  v[bcgs[i]:bcgs[i+1]], zobs=zobs, fov=fov,
                  slit_size=slit_size, ngals=ngals, width=width,
                  fill=fill, mbm=mbm, iterate=iterate)
    out.close()
    print 'Saved to', output
    return

def observe_random(hostids, bcgs, j, x, y, v, zobs=0.5,
                   fov=5, ngals=60, nmin=100, mbm=True, iterate=True):
    ngr = len(bcgs)
    #nmin = min(
    output = 'outputs_%s/south'\
             %(''.join(str(datetime.today()).split()[0].split('-')))
    output += '-nmin_%d-zobs_%.2f-fov_%.1f-ngals_%d.out' \
              %(nmin, zobs, fov, ngals)
    if mbm:
        output = output.replace('.out', '_MBM10.out')
    if not iterate:
        output = output.replace('.out', '_noiter.out')
    print 'Will save to', output
    hdr = '# hostid  ntot  nobs  n200'
    hdr += '  s200  s200_err  m200  m200_err  r200  r200_err'
    #out = None
    out = open(output, 'w')
    print >>out, hdr
    aux = [one_group(out, hostids[i],
                     x[bcgs[i]:bcgs[i+1]], y[bcgs[i]:bcgs[i+1]],
                     v[bcgs[i]:bcgs[i+1]], zobs=zobs, fov=fov,
                     ngals=ngals, rnd=True, mbm=mbm, iterate=iterate) \
           for i in j[:-1]]
    # the last halo, in case it is the last one on the list
    i = j[-1]
    if i == ngr - 1:
        one_group(out, hostids[i],
                  x[bcgs[i]:], y[bcgs[i]:], v[bcgs[i]:],
                  zobs=zobs, fov=fov, ngals=ngals, rnd=True,
                  mbm=mbm, iterate=iterate)
    else:
        one_group(out, hostids[i],
                  x[bcgs[i]:bcgs[i+1]], y[bcgs[i]:bcgs[i+1]],
                  v[bcgs[i]:bcgs[i+1]], zobs=zobs, fov=fov,
                  ngals=ngals, rnd=True, mbm=mbm, iterate=iterate)
    out.close()
    print 'Saved to', output
    return

def one_group(out, hostid, x, y, v, zobs=0.5, fov=10,
              slit_size=None, ngals=40, width=1.,
              rnd=False, fill=False, mbm=True, iterate=True):
    # relative positions -- is the BCG actually the first object in the list?
    x -= x[0]
    y -= y[0]
    r = scipy.hypot(x, y)
    if fov is None:
        obs = scipy.arange(len(x), dtype=int)
    else:
        # coordinates in arcmin
        X = scipy.array([cosmology.dProj(zobs, xi, input_unit='Mpc',
                                         unit='arcmin') \
                         for xi in x])
        Y = scipy.array([cosmology.dProj(zobs, yi, input_unit='Mpc',
                                         unit='arcmin') \
                         for yi in y])
        if rnd:
            size = min(ngals,len(x))
            j = scipy.arange(len(x))
            obs = random.permutation(j[(X < fov) & (Y < fov)])[:size]
        else:
            obs = make_mask(X, Y, v, fov=fov,
                            slit_size=slit_size, ngals=ngals,
                            width=width, fill=fill)
    dynamics(hostid, r, v, obs, out, mbm=mbm, iterate=iterate)
    return

main()
