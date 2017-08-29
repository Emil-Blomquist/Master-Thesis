import numpy as np
import matplotlib.pyplot as plt

from Dyson import Dyson
from fit2exponential import fit2exponential

def DysonAndFit (t, S, a, p, mu, unique, plot = False):
  # increase mu in order to make dSE_w less sharply peaked
  dMus = np.linspace(-0.1, 0, 20)
  deltaMu = dMus[1] - dMus[0]

  DMus = []
  E0s = []
  Z0s = []


  for dMu in dMus:

    G = Dyson(S*np.exp(dMu*t), t, a, p, mu + dMu)

    # plt.plot(t, G)
    # plt.show()


    E0, Z0 = fit2exponential(G, t, p, mu + dMu, True, False)

    if E0 is not False:
      DMus += [dMu]
      E0s += [E0]
      Z0s += [Z0]

      print('\t', str(dMu) + ': ', E0, Z0)


  #  to litte data -> cant get information about E0 and Z0
  if len(DMus) < 2:
    return False, False


  # select the latest possible mu for which Z0 is not increasing to fast
  dZdMus = np.diff(Z0s)/deltaMu
  chosen = len(Z0s) - 1
  for i, dZdMu in enumerate(dZdMus):
    if dZdMu > 1e-4:
      chosen = i
      break

  f, axarr = plt.subplots(4, sharex=True, figsize=(20, 10))
  axarr[0].plot(DMus, E0s)
  axarr[0].plot(DMus[chosen], E0s[chosen], 'o', color='green')
  axarr[1].plot(DMus, Z0s)
  axarr[1].plot(DMus[chosen], Z0s[chosen], 'o', color='green')
  axarr[2].plot(DMus, Z0s)
  axarr[2].set_ylim([Z0s[chosen] - 0.0001, Z0s[chosen] + 0.0001])
  axarr[2].plot(DMus[chosen], Z0s[chosen], 'o', color='green')
  axarr[3].plot(np.array(DMus[1:]) - 0.5*deltaMu, np.diff(Z0s)/deltaMu)

  if chosen != len(Z0s) - 1:
    axarr[3].plot(DMus[chosen] + 0.5*deltaMu, (np.diff(Z0s)/deltaMu)[chosen], 'o', color='green')

  axarr[3].set_xlabel(r'$ \delta \mu $')
  axarr[3].set_ylabel(r'$ \partial_\mu Z_0 $')
  axarr[3].set_xlabel(r'$ \Delta \mu $')
  axarr[0].set_ylabel(r'$ E_0 $')
  axarr[1].set_ylabel(r'$ Z_0 $')
  axarr[2].set_ylabel(r'$ Z_0 $')
  axarr[3].set_ylabel(r'$ \partial_\mu Z_0 $')

  plt.tight_layout()
  plt.savefig('plots/DysonAndFit/' + r'a={0}_p={1}_u={2}'.format(a, p, unique) + '.pdf')
  
  if plot:
    plt.show()

  # clear figure
  f.clf()
  plt.close()
  del f

  return E0s[chosen], Z0s[chosen]