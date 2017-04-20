import scipy
from scipy import linspace, stats, optimize
from scipy import pi,sqrt,exp
from scipy.special import erf

from pylab import plot,show

def pdf(x):
    return 1/sqrt(2*pi) * exp(-x**2/2)

def cdf(x):
    return (1 + erf(x/sqrt(2))) / 2

def skew(x, *args):
    e, w, a = args
    t = (x-e) / w
    return 2 / w * pdf(t) * cdf(a*t)

def double(x, *args):
    m1, m2, s1, s2, k1, k2 = args
    ret = k1*scipy.stats.norm.pdf(x, loc=m1 ,scale=s1)
    ret += k2*scipy.stats.norm.pdf(x, loc=m2 ,scale=s2)
    return ret

