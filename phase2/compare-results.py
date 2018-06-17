#!/usr/bin/env python
import glob
import itertools
import os
import pylab
import readfile
import scipy
import stattools

c = 299792.458

def main():
    xo = 0.09
    dx = 0.44
    bins = scipy.linspace(0, 2200, 50)
    pylab.figure(figsize=(9,9))
    axes = [pylab.axes([xo, xo+dx+0.02, dx, dx], xticklabels=[])]
    run(bins, catalog='HOD2_SW', version='SG1')
    plot_massimg('HOD2', 'SG1', max(bins))
    pylab.ylabel(r'$S_{\rm Sif\'on}\,(\mathrm{km\,s^{-1}})$')
    axes.append(pylab.axes([xo+dx+0.02, xo+dx+0.02, dx, dx],
                           xticklabels=[], yticklabels=[]))
    run(bins, catalog='SAM2_SW', version='SG1')
    plot_massimg('SAM2', 'SG1', max(bins))
    l = pylab.legend(loc='upper left', numpoints=1)
    l.get_frame().set_alpha(0)
    axes.append(pylab.axes([xo, xo, dx, dx]))
    run(bins, catalog='HOD2_SW', version='SG2')
    plot_massimg('HOD2', 'SG2', max(bins))
    pylab.xlabel(r'$S_{\rm True}\,(\mathrm{km\,s^{-1}})$')
    pylab.ylabel(r'$S_{\rm Sif\'on}\,(\mathrm{km\,s^{-1}})$')
    axes.append(pylab.axes([xo+dx+0.02, xo, dx, dx], yticklabels=[]))
    run(bins, catalog='SAM2_SW', version='SG2')
    plot_massimg('SAM2', 'SG2', max(bins))
    pylab.xlabel(r'$S_{\rm True}\,(\mathrm{km\,s^{-1}})$')
    #axes[-1].get_xticklabels(0).set_visible(False)
    for ax in axes:
        ax.set_xlim(0, 2200)
        ax.set_ylim(0, 2200)
    output = 'plots/compare-sigma.png'
    pylab.savefig(output, format=output[-3:])
    pylab.close()
    print 'Saved to', output
    return

def run(bins, catalog='HOD2_SW', version='SG2'):
    if catalog[-2:] == 'SW':
        folder = 'shallow'
    elif catalog[-2:] == 'DN':
        folder = 'deep'
    else:
        print 'ERROR: catalog does not exist'
        exit()
    myfile = 'HMRC_Halo_{0}_{1}.txt'.format(catalog, version)
    myfile = os.path.join(folder, 'ret', myfile)
    mydata = readfile.dict(myfile)
    Nmine = len(mydata['halo-id'])
    halofile = os.path.join('Phase_2_input_catalogues',
                            '{0}_input_halo_data_M200c.txt'.format(catalog))
    truehalos = readfile.dict(halofile)
    truehalos['z'] = truehalos['v_res'] / c
    Nhalo = len(truehalos['halo-id'])
    galfile = os.path.join('Phase_2_input_catalogues',
                            '{0}_input_galaxy_data.txt'.format(catalog))
    truegals = readfile.dict(galfile)
    Ngals = len(truegals)
    s = [stattools.Sbi(truegals['v_res'][truegals['halo-id'] == i])
         for i in xrange(1, Nhalo+1)]
    truehalos['sigma'] = scipy.array(s) / (1+truehalos['z'])

    pylab.plot(truehalos['sigma'], mydata['vdisp200'], 'k.')
    pylab.hist(truehalos['sigma'], bins, histtype='step',
               color='r', weights=10*scipy.ones(Nhalo), zorder=10,
               label=r'$S_{\rm True}$')
    pylab.plot(bins, bins, 'b-')
    pylab.hist(mydata['vdisp200'], bins, histtype='step',
               color='c', weights=10*scipy.ones(Nhalo), zorder=10,
               label=r'$S_{\rm Sif\'on}$')
    pylab.annotate('{0} {1}\n{2}'.format(catalog[:4], folder, version),
                   xy=(0.60,0.95), xycoords='axes fraction',
                   ha='left', va='top')
    return

def plot_massimg(catalog, version, smax):
    extent = (0.55*smax, 0.99*smax, 0.10*smax, 0.54*smax)
    img = 'plots/paper2/{0}_{1}.png'.format(catalog, version)
    img = pylab.imread(img)
    pylab.imshow(img, extent=extent)
    return

main()
