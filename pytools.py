import scipy
import stattools
from astro.clusters import conversions, scalings

def M200(r, z, ntot=0, errors=False, converge=True, verbose=False):
  """
  ntot is to be used in the bootstrap resampling
  """
  c = 299792.458
  def _mass(z, zo):
    s = c * stattools.Sbi(z) / (1 + zo)
    m = scalings.sigma(s, zo)
    r200 = conversions.rsph(m, zo, ref='200c', unit='kpc')
    return s, m, r200

  zo = stattools.Cbi(z)
  s, m, r200 = _mass(z, zo)
  if verbose:
    print 'iteration 0:'
    print '  sigma = %d km/s' %int(s)
    print '  m200 = %.1e' %m
    print '  r200 = %d' %int(r200)
    print '  Ngal = %d' %len(r)
  j = scipy.arange(len(r), dtype=int)
  if converge:
    jo = j[r < r200]
    n200 = 0
    it = []
    while len(jo) != n200:
      try:
        n200 = len(jo)
        s, m, r200 = _mass(z[jo], zo)
        if verbose:
          print 'iteration %d:' %(len(it)+1)
          try:
            print '  sigma = %d km/s' %int(s)
            print '  m200 = %.1e' %m
            print '  r200 = %d' %int(r200)
            print '  n200 = %d' %n200
          except ValueError:
            print '  Measurement failed'
        jo = j[r < r200]
        if r200 in it:
          break
        it.append(r200)
      except ZeroDivisionError:
        break
  else:
    n200 = len(r)
  if n200 == 0:
    return zo, -1, (-1, (-1,-1)), -1, -1
  if errors:
    if n200 > 15:
      stats = (stattools.Cbi, stattools.Sbi)
      z_err, s_err = stattools.bootstrap(stats, z, ntot)
      s_err *= c / (1 + zo)
    else:
      z_err = scipy.std(z) / scipy.sqrt(n200)
      s_err = c/(1+zo) * z_err / scipy.sqrt(2)
    m = scalings.sigma(s, zo, ds=s_err, dz=z_err,
                       scaling='evrard08', separate_errors=True)
    r = conversions.rsph(m[0], zo, ref='200c', unit='kpc')
    merr = [scipy.log10(1+me/m[0]) for me in m[1]]
    m = (m[0], merr)
  else:
    m = scalings.sigma(s, zo)
    r = conversions.rsph(m, zo, ref='200c', unit='kpc')
  return zo, s, m, r, n200
