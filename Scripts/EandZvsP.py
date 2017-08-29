import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import json

from read2D import read2D
from read1D import read1D
from DysonAndFit import DysonAndFit


path = 'data/EandZvsP'
refinedDataPath = 'data/refined/EandZvsP.txt'

data = {}
if os.path.isfile(refinedDataPath):
  # file exist, parse content
  with open(refinedDataPath) as dataFile:    
    data = json.load(dataFile)
    data = {float(k1):{float(k2):v2 for k2, v2 in v1.items()} for k1, v1 in data.items()}

# Varying momentum simulations
for root, dirs, files in os.walk(os.getcwd() + '/' + path):
  # filter files to only accept data files (.txt)
  files = [path + '/' + file for file in files if (file.endswith(".txt") and ' p=' not in file)]

  for file in files:
    Ts, Ps, S, parameters = read2D(file)

    a = parameters['a']
    mu = parameters['mu']
    unique = parameters['unique']

    for i, p in enumerate(Ps[:, 0]):
      if ((p in data) and (unique in data[p])):
        continue

      print(p)
      E0, Z0 = DysonAndFit(Ts[i,:], S[i,:], a, p, mu, unique)
      print(E0, Z0)
      print('---------')

      if p in data:
        data[p][unique] = [E0, Z0]
      else:
        data[p] = {
          unique: [E0, Z0]
        }

      with open(refinedDataPath, 'w') as dataFile:
        json.dump(data, dataFile)


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

    if ((p in data) and (unique in data[p])):
      continue
    
    print(p)
    E0, Z0 = DysonAndFit(t, S, a, p, mu, unique)
    print(E0, Z0)
    print('---------')

    if p in data:
      data[p][unique] = [E0, Z0]
    else:
      data[p] = {
        unique: [E0, Z0]
      }

    with open(refinedDataPath, 'w') as dataFile:
      json.dump(data, dataFile)


ps = []
avgE0s = []
avgZ0s = []
stdE0s = []
stdZ0s = []

# calculate mean and deviation
for p, val in sorted(data.items()):
  E0sAtP = np.array([])
  Z0sAtP = np.array([])

  for k, v in val.items():
    if v[0] is not False:
      E0sAtP = np.append(E0sAtP, v[0])
      Z0sAtP = np.append(Z0sAtP, v[1])

  # more than one item in order to calculate statistics
  if E0sAtP.size > 1:
    avgE0 = np.average(E0sAtP)
    stdE0 = np.sqrt(np.average((E0sAtP - avgE0)**2))

    avgZ0 = np.average(Z0sAtP)
    stdZ0 = np.sqrt(np.average((Z0sAtP - avgZ0)**2))

    ps += [p]
    avgE0s += [avgE0]
    avgZ0s += [avgZ0]
    stdE0s += [stdE0]
    stdZ0s += [stdZ0]



print(stdE0s)
print(stdZ0s)



f, axarr = plt.subplots(1, 2, sharex=True, figsize=(12, 5))

# axarr[0].errorbar(ps, avgE0s, yerr=stdE0s, marker='x', ls=":")
# axarr[1].errorbar(ps, avgZ0s, yerr=stdZ0s, marker='x', ls=":")
axarr[0].plot(ps, avgE0s, marker='.', ls=":")
axarr[1].plot(ps, avgZ0s, marker='.', ls=":")

f.suptitle(r'$ \alpha = {0} $'.format(1), fontsize=18)
axarr[0].set_ylabel('$ E_0 $', fontsize=18)
axarr[1].set_ylabel('$ Z_0 $', fontsize=18)
axarr[0].set_xlabel('$ p $', fontsize=18)
axarr[1].set_xlabel('$ p $', fontsize=18)

axarr[0].autoscale(True, 'both', True)
axarr[0].margins(0.02, 0.03)
axarr[1].autoscale(True, 'both', True)
axarr[1].margins(0.02, 0.03)

plt.tight_layout(rect=[-0.01, 0, 1, 0.96])
plt.savefig('plots/EandZvsP.pdf')
plt.show()

