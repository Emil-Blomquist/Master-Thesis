import numpy as np
import matplotlib.pyplot as plt
from decimal import Decimal

from scipy.special import erf


a = 1
mu = -1.1
dt = 0.02

def S (ts):
  return a*(np.pi*ts)**-0.5 * np.exp((mu - 1)*ts)

def mS (ts):
  t0 = ts - 0.5*dt
  t1 = ts + 0.5*dt
  p = (1 - mu)**0.5

  return a/dt/p * (erf(p*t1**0.5) - erf(p*t0**0.5))


n = 8
ts = np.linspace(0, (n - 1)*dt, n) + 0.5*dt

# print(["{%.2f}".format(i) for i in ts])

# for t in ts:
  # print(t, np.log(t))

  # print(np.round((mS(t) - S(t))*10000)/1000)

# plt.plot(ts, S(ts))
# plt.plot(ts, mS(ts))
# plt.show()


Ds = mS(ts) - S(ts) 

for D in Ds:
  # m = np.floor(np.log(D)/np.log(10))
  # d = D/10**m
  # print(np.round(d*100)/100 * 10**m)
  print('%.2E' % Decimal(D))