#!/usr/bin/env python
import itertools
import lnr
import os
import pylab
import readfile
import scipy
import stattools
from astro import cosmology
from astro.clusters import scalings

c = 299792.458

def main():
    xo = 0.12
    dx = 0.42
    yo = 0.14
    dy = 0.83
    pylab.figure(figsize=(7,4))
    pylab.axes([xo, yo, dx, dy], xscale='log', yscale='log')
    run('HOD2_SW')
    pylab.ylabel(r'$S_{200}\,(\mathrm{km\,s^{-1}})$')
    pylab.xlabel(r'$E(z)M_{200}\,(M_\odot)$')
    pylab.axes([xo+dx+0.02, yo, dx, dy], xscale='log', yscale='log',
               yticklabels=[])
    run('SAM2_SW')
    #pylab.xlabel(r'$S_{200}\,(\mathrm{km\,s^{-1}})$')
    pylab.xlabel(r'$E(z)M_{200}\,(M_\odot)$')
    output = 'plots/internal-scalings.png'
    pylab.savefig(output, format=output[-3:])
    pylab.close()
    print 'Saved to', output
    return

def run(catalog='HOD2_SW'):
    print catalog
    halofile = os.path.join('Phase_2_input_catalogues',
                            '{0}_input_halo_data_M200c.txt'.format(catalog))
    truehalos = readfile.dict(halofile)
    truehalos['logm200'] = truehalos.pop('m200')
    truehalos['m200'] = 10**truehalos['logm200']
    truehalos['z'] = truehalos['v_res'] / c
    truehalos['m200hz'] = cosmology.E(truehalos['z']) * truehalos['m200']
    Nhalo = len(truehalos['halo-id'])
    galfile = os.path.join('Phase_2_input_catalogues',
                            '{0}_input_galaxy_data.txt'.format(catalog))
    truegals = readfile.dict(galfile)
    Ngals = len(truegals)
    s = [stattools.Sbi(truegals['v_res'][truegals['halo-id'] == i])
         for i in xrange(1, Nhalo+1)]
    truehalos['sigma'] = scipy.array(s) / (1+truehalos['z'])
    # write this out
    output = '{0}_input_halo_data_M200c_sigma.txt'.format(catalog)
    out = open(output, 'w')
    print >>out, '# halo-id  m200  ra  dec  z  vdisp'
    for i in xrange(Nhalo):
        print >>out, '%d  %.3e  %.9f  %.9f  %.6f  %7.2f' \
                     %(truehalos['halo-id'][i], truehalos['m200'][i],
                       truehalos['ra'][i], truehalos['dec'][i],
                       truehalos['z'][i], truehalos['sigma'][i])
    out.close()
    print 'Saved to', output

    # fit scaling relation
    j = scipy.arange(Nhalo)
    #j = j[truehalos['m200hz'] >= 1e14]
    Mo = int(round(scipy.median(truehalos['m200'][j]), 0))
    So = int(round(scipy.median(truehalos['sigma'][j]), 0))
    s = scipy.logspace(2, 3.3, 100)
    t = scipy.logspace(13, 15.5, 100)
    #Amle, Bmle, Smle = lnr.mle(truehalos['sigma']/So, truehalos['m200'],
                               #logify=True)
    Amle, Bmle, Smle = lnr.mle(truehalos['m200hz'][j]/Mo,
                               truehalos['sigma'][j],
                               logify=True)
    print Amle, Bmle, Smle
    #A, B, Sint = lnr.mcmc(truehalos['sigma']/So, truehalos['m200'],
                          #nsteps=2000, logify=True, output=(50,16,84))
    A, B, Sint = lnr.mcmc(truehalos['m200hz'][j]/Mo,
                          truehalos['sigma'][j],
                          nsteps=2000, logify=True, output=(50,16,84))
    print A
    print B
    print Sint
    #pylab.plot(truehalos['sigma'], truehalos['m200'], 'k.')
    pylab.plot(truehalos['m200hz'], truehalos['sigma'], 'k.')
    if len(j) < Nhalo:
        pylab.plot(truehalos['m200hz'][j], truehalos['sigma'][j], 'b.')
    #bestfit = '$S_0=%d$\n$A=%.3f$\n$B=%.3f$' %(So, A, B)
    #bestfit = r'$S_0=%d\,\mathrm{km\,s^{-1}}$' %(So)
    #bestfit += '\n'
    #bestfit += r'$A=%.3f_{-%.4f}^{+%.4f}$' %(A[0], A[0]-A[1], A[2]-A[0])
    #bestfit += '\n'
    bestfit = r'$B=%.3f_{-%.3f}^{+%.3f}$' %(B[0], B[0]-B[1], B[2]-B[0])
    pylab.plot(t, 10**A[0] * (t/Mo)**B[0], 'r-', lw=2)#, label=bestfit)
    #Slabel = r'$\sigma_{M|S}=%.3f_{-%.3f}^{+%.3f}$' \
             #%(Sint[0], Sint[0]-Sint[1], Sint[2]-Sint[0])
    Slabel = r'$\sigma_{S|M}=%.3f_{-%.3f}^{+%.3f}$' \
             %(Sint[0], Sint[0]-Sint[1], Sint[2]-Sint[0])
    #pylab.plot(t, (1-Sint[0]) * 10**A[0] * (t/Mo)**B[0], 'r--')
    #pylab.plot(t, (1+Sint[0]) * 10**A[0] * (t/Mo)**B[0], 'r--')#,
               #label=Slabel)
    pylab.plot(scalings.sigma(s, 0, scaling='evrard08'), s, 'c-', lw=2,
               label='Evrard+08')
    pylab.plot(scalings.sigma(s, 0, scaling='munari13dm'), s, '-',
               color='orange', lw=2, label='Munari+13 DM')
    pylab.plot(scalings.sigma(s, 0, scaling='munari13sub'), s, '--',
               color='orange', lw=2, label='Munari+13 sub')
    pylab.plot(scalings.sigma(s, 0, scaling='munari13gal'), s, ':',
               color='orange', lw=2, label='Munari+13 gal')
    if catalog[-2:] == 'SW':
        depth = 'shallow'
    elif catalog[-2:] == 'DN':
        depth = 'deep'
    msg = '{0} {1}'.format(catalog[:4], depth)
    pylab.annotate(msg, xy=(0.04,0.96), xycoords='axes fraction',
                   ha='left', va='top')
    l = pylab.legend(loc='lower right')
    l.get_frame().set_alpha(0)
    pylab.ylim(80, 2000)
    pylab.xlim(1e13, 2e15)
    return

main()
