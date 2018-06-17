#!/usr/bin/env python
import glob
import os
import pylab
import readfile
import scipy
from astLib import astCoords
from astro import cosmology
from astro.clusters import redsequence

act_path = '/Users/cristobal/Documents/act/'

def main(survey='s82'):
    phot_path = os.path.join(act_path, 'catalogs-photometry/')
    if survey == 'south':
        spec_path = os.path.join(act_path, 'sifon+2012/paper/members/')
    elif survey in ('dr8', 's82'):
        spec_path = os.path.join(act_path, 'equator/members/')
    ls = sorted(glob.glob(os.path.join(spec_path, '*.cat')))
    for i in ls:
        cluster = i.split('/')[-1]
        if '_finder' in cluster or '_merged' in cluster:
            cluster = '_'.join(cluster.split('_')[:-1])
        else:
            cluster = cluster[:-4]
        try:
            run(cluster, survey, spec_path, phot_path)
        except TypeError:
            print '%s: TypeError' %cluster
            pass
    return

def run(cluster, survey, spec_path, phot_path):
    center = get_center(cluster, survey)
    if center is None:
        return
    zo = get_z(cluster, survey)
    print cluster, zo, center
    if survey == 'south':
        spec = os.path.join(spec_path, cluster + '.cat')
        phot = os.path.join(phot_path, cluster + '_merged.cat')
        phot = readfile.table(phot, cols=(1,2,6,7,9,10))
    elif survey == 's82':
        spec = os.path.join(spec_path, cluster + '.cat')
        phot = os.path.join(phot_path, cluster + '_finder.dat')
        try:
            phot = readfile.table(phot, cols=(1,2,7,8,9,10,35))
        except IOError:
            return
        #phot = readfile.table(phot, cols=(1,2,16,17,18,19))
        zphot = phot[-1]
    spec = readfile.table(spec, cols=(1,2))
    print '%d members' %len(spec[0])
    rmag = phot[2]
    color = phot[2] - phot[4]
    cerr = scipy.hypot(phot[3], phot[5])
    #if survey == 's82':
        #j = scipy.arange(len(rmag))[(zphot > 0.4) & (zphot < 0.6)]
    #pylab.errorbar(rmag[j], color[j], yerr=cerr[j], fmt='r.')
    #pylab.show()
    j = scipy.arange(len(rmag))[(cerr < 0.2)]
    rs = redsequence.fit(rmag[j], color[j], cerr[j], pivot=20,
                         plot_output='cmr_%s.png' %cluster,
                         mag_label='m_r', color_label='r-i')
    rsg = rs[0]
    dphot, rphot = distances(phot[0][j[rsg]], phot[1][j[rsg]], center, zo)
    dspec, rspec = distances(spec[0], spec[1], center, zo)

    dbins = scipy.arange(0, 10.01, 0.5)
    rbins = scipy.arange(0, 2.01, 0.10)
    pylab.figure(figsize=(5,4))
    pylab.axes([0.12, 0.12, 0.85, 0.85])
    pylab.hist(rphot, rbins, color='r',
               histtype='step', lw=2, label='red sequence')
    pylab.hist(rspec, rbins, color='b',
               histtype='step', lw=2, label='spec')
    pylab.xlabel('distance from BCG (Mpc)')
    pylab.ylabel('N')
    pylab.legend(loc='upper right')
    pylab.xlim(0, 2.2)
    output = 'distance_%s.png' %cluster
    pylab.savefig(output, format='png')
    pylab.close()
    print 'Saved to', output
    return

def distances(ra, dec, center, z):
    d = astCoords.calcAngSepDeg(ra, dec, center[0], center[1])
    r = scipy.array([cosmology.dProj(z, di) for di in d])
    return 60*d, r

def get_center(cluster, survey):
    try:
        if survey == 'south':
            spec = os.path.join(act_path, 'sifon+2012/paper/catalogs/')
            spec = os.path.join(spec, cluster + '.cat')
            file = open(spec)
            center = file.readline().split()[1:]
            file.close()
            return [float(i) for i in center]
        else:
            bcgs = readfile.table(os.path.join(act_path, 'equator/bcgs.dat'),
                                  cols=(2,3), include=cluster)
            return bcgs
    except IndexError:
        return

def get_z(cluster, survey):
    if survey == 'south':
        dyn = os.path.join(act_path,
                           'sifon+2012/paper/ref_review/dynamics_final.dat')
        dyn = readfile.table(dyn, cols=2, include=cluster)
    else:
        dyn = os.path.join(act_path, 'equator/dynamics/dynamics.out')
        dyn = readfile.table(dyn, cols=5, include=cluster)
    return dyn

main()
