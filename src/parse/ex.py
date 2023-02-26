import sys
import os
sys.path.append(os.getcwd() + '/..')
from base import Parameters

test = r"""
icntl = 0
maxcas = 10
maxbch = 1
emin(1) = 0.001
emin(2) = 1.0E-10
dmax(2) = 20.0
emin(12) = 0.001
emin(13) = 0.001
emin(14) = 0.001
emin(15) = 0.001
emin(16) = 0.001
emin(17) = 0.001
emin(18) = 0.001
emin(19) = 0.00001
esmin = 0.000001
igamma = 1
itall = 1
igchk = 1
file(6) = shuttle-block_validation.out
file(7) = /Users/jeffchancellor/Applications/phits/data/xsdir.jnd
file(20) = /Users/jeffchancellor/Applications/phits/data/xsdir.jnd
ides = 0


nedisp = 1
nspred = 1
nlost = 200
istdev = 1
irqmd = 1
ismm = 1
icrhi = 1
ndedx = 3
mdbatima = 100"""

print(Parameters.parse(test).pretty())
