import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import json

from read2D import read2D
from read1D import read1D
from DysonAndFit import DysonAndFit


path = 'data/EandZvsA'
refinedDataPath = 'data/refined/EandZvsA.txt'

data = {}
if os.path.isfile(refinedDataPath):
  # file exist, parse content
  with open(refinedDataPath) as dataFile:    
    data = json.load(dataFile)
    data = {float(k1):{float(k2):v2 for k2, v2 in v1.items()} for k1, v1 in data.items()}


# Fixed momentum simulations
for root, dirs, files in os.walk(os.getcwd() + '/' + path):
  # filter files to only accept data files (.txt)
  files = [path + '/' + file for file in files if (file.endswith(".txt") and 'pmax' not in file)]

  for file in files:
    t, S, parameters = read1D(file)

    a = parameters['a']
    p = parameters['p']
    mu = parameters['mu']
    unique = parameters['unique']

    if ((a in data) and (unique in data[a])):
      continue
    
    print(a)
    E0, Z0 = DysonAndFit(t, S, a, p, mu, unique)
    print(E0, Z0)
    print('---------')

    if a in data:
      data[a][unique] = [E0, Z0]
    else:
      data[a] = {
        unique: [E0, Z0]
      }

    with open(refinedDataPath, 'w') as dataFile:
      json.dump(data, dataFile)


a_s = []
avgE0s = []
avgZ0s = []
stdE0s = []
stdZ0s = []

# calculate mean and deviation
for a, val in sorted(data.items()):
  E0sAtA = np.array([])
  Z0sAtA = np.array([])

  for k, v in val.items():
    if v[0] is not False:
      E0sAtA = np.append(E0sAtA, v[0])
      Z0sAtA = np.append(Z0sAtA, v[1])

  # more than one item in order to calculate statistics
  if E0sAtA.size > 1:
    avgE0 = np.average(E0sAtA)
    stdE0 = np.sqrt(np.average((E0sAtA - avgE0)**2))

    avgZ0 = np.average(Z0sAtA)
    stdZ0 = np.sqrt(np.average((Z0sAtA - avgZ0)**2))

    a_s += [a]
    avgE0s += [avgE0]
    avgZ0s += [avgZ0]
    stdE0s += [stdE0]
    stdZ0s += [stdZ0]



print(stdE0s)
print(stdZ0s)



f, axarr = plt.subplots(2, 1, sharex=True, figsize=(10, 10))


a = np.linspace(0, 6, 1000)
E0first = -a
E0second = -a - 2/3*(1/8 - 1/(3*np.pi))*a**2
axarr[0].plot(a, E0first, '--', color='green', label=r'$E_0^{(1)}$')
axarr[0].plot(a, E0second, '--', color='red', label=r'$E_0^{(2)}$')

a = np.linspace(0, 1, 1000)
Z0first = 1 - a/2
axarr[1].plot(a, Z0first, '--', color='green', label=r'$Z_0^{(1)}$')



# axarr[0].errorbar(a_s, avgE0s, yerr=stdE0s, marker='x', ls=":")
# axarr[1].errorbar(a_s, avgZ0s, yerr=stdZ0s, marker='x', ls=":")
axarr[0].plot(a_s, avgE0s, marker='.', ls=":", label='$E_0^\mathrm{DMC} $')
axarr[1].plot(a_s, avgZ0s, marker='.', ls=":", label='$Z_0^\mathrm{DMC} $')

f.suptitle(r'$ p = {0} $'.format(0), fontsize=18)
# axarr[0].set_ylabel('$ E_0 $', fontsize=18)
# axarr[1].set_ylabel('$ Z_0 $', fontsize=18)
# axarr[0].set_xlabel(r'$ \alpha $', fontsize=18)
axarr[1].set_xlabel(r'$ \alpha $', fontsize=18)

axarr[0].autoscale(True, 'both', True)
axarr[0].margins(0.02, 0.03)
axarr[1].autoscale(True, 'both', True)
axarr[1].margins(0.02, 0.03)


axarr[0].legend(loc=1, fontsize=18)
axarr[1].legend(loc=1, fontsize=18)

# plt.tight_layout(rect=[-0.01, 0, 1, 0.96])
# plt.savefig('plots/EandZvsA.pdf')
# plt.show()




for i in range(0, len(a_s)):
  print(r'{0:g} & {1:.4f} & {2:.4f} \\'.format(a_s[i], avgE0s[i], avgZ0s[i]))