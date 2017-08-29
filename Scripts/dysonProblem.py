import numpy as np
import matplotlib.pyplot as plt
import os

from read1D import read1D
from DysonAndFit import DysonAndFit


# file = 'data/EandZvsP/' + 'a=1.000000 p=0.000000 mu=-1.020000 tmax=60.000000 dt=0.005000 MCsecs=600.000000 DMCsecs=86400.000000 date=2017-08-25 00:39:01 numcores=4 unique=0.566519.txt'
# file = 'data/EandZvsP/' + 'a=1.000000 p=1.000000 mu=-0.620000 tmax=80.000000 dt=0.005000 MCsecs=600.000000 DMCsecs=86400.000000 date=2017-08-25 00:39:03 numcores=4 unique=0.025428.txt'
file = 'data/EandZvsP/' + 'a=1.000000 p=1.730000 mu=-0.050000 tmax=170.000000 dt=0.005000 MCsecs=600.000000 DMCsecs=86400.000000 date=2017-08-25 00:39:02 numcores=4 unique=0.344755.txt'
# file = 'data/EandZvsP/' + 'a=1.000000 p=1.820000 mu=-0.018000 tmax=200.000000 dt=0.005000 MCsecs=600.000000 DMCsecs=86400.000000 date=2017-08-25 02:57:42 numcores=4 unique=0.605376.txt'
# file = 'data/EandZvsP/' + 'a=1.000000 p=1.850000 mu=-0.013000 tmax=200.000000 dt=0.005000 MCsecs=600.000000 DMCsecs=86400.000000 date=2017-08-25 03:24:33 numcores=4 unique=0.713220.txt'




t, S, parameters = read1D(file)

a = parameters['a']
p = parameters['p']
mu = parameters['mu']
unique = parameters['unique']

E0, Z0 = DysonAndFit (t, S, a, p, mu, unique, False)

print(E0, Z0)