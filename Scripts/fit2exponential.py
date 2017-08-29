import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sig
from scipy.odr import *

import numbers

import warnings
warnings.filterwarnings(action="ignore", module="matplotlib", message="^tight_layout")

def fit2exponential (G, t, p, mu, Dyson = False, plot = False):


  dt = t[1] - t[0]
  length = t[-1] + 0.5*dt
  minFittingLength = length*0.2


  # find first zero in G and cut
  makeFinite = np.where(np.real(G) < 1e-9)# & np.abs(np.imag(np.log(G))) > 1e-2)
  if makeFinite[0].size:
    end = makeFinite[0][0]
  else:
    end = G.size - 1

  # removed zeros which will be troublesome when using log later
  t = t[:end]
  G = np.real(G[:end])

  
  # remove oscillations from Dyson equation 
  upperBound = end - 1
  if Dyson: 
    oscillations = np.abs(np.diff(np.log(G)))
    oscDiff = np.diff(oscillations)

    for i in range(0, oscDiff.size):
      if oscillations[i] > 1e-2 and oscDiff[i] > 1e-5:
        upperBound = i - 1
        break




  Z = 1
  E = 0
  lowerBound = 0
  fittingData = True
  
  m = np.log(Z)
  k = mu - E
  y = np.log(G)

  while fittingData:

    # Define a function to fit the data with.
    def lin_fit(params, t):
      m, k = params
      return m + k*t

    # Create a model for fitting.
    lin_model = Model(lin_fit)

    # Create a RealData object using our initiated data from above.
    data = RealData(t[lowerBound:upperBound], y[lowerBound:upperBound])

    # Set up ODR with the model and data.
    odr = ODR(data, lin_model, beta0=[m, k])

    # Run the regression.
    out = odr.run()

    m = out.beta[0]
    k = out.beta[1]

    # print and display while fitting
    if False:
      Z = np.exp(m)
      E = mu - k
      g_fit = Z*np.exp(-(E - mu)*t)

      diff = y[lowerBound:upperBound] - lin_fit([m, k], t[lowerBound:upperBound])
      signs = np.sign(diff)
      rmse = np.mean(diff**2)**0.5
      cumsum = np.cumsum(diff)*dt

      print(t[lowerBound], 'rmse=' + str(rmse), 'cumsum=' + str(np.max(np.cumsum(diff)*dt)))

      if True:
        f, axarr = plt.subplots(4, sharex=True, figsize=(20, 10))

        axarr[0].set_title(r'$p = {0}, \, \mu = {1}, \, E_0 = {2}, \, Z_0 = {3}$'.format(p, mu, E, Z))

        axarr[0].plot(t, G)
        axarr[0].plot(t[lowerBound], G[lowerBound], 'o', color='green')
        axarr[0].plot(t, Z*np.exp(-(E - mu)*t), '--', color='red')
        axarr[0].plot(t[upperBound], G[upperBound], 'o', color='green')
        axarr[0].set_xlabel(r'$\tau$')
        axarr[0].set_ylabel(r'$G$')

        axarr[1].semilogy(t, G)
        axarr[1].semilogy(t[lowerBound], G[lowerBound], 'o', color='green')
        axarr[1].semilogy(t, Z*np.exp(-(E - mu)*t), '--', color='red')
        axarr[1].semilogy(t[upperBound], G[upperBound], 'o', color='green')
        axarr[1].set_xlabel(r'$\tau$')
        axarr[1].set_ylabel(r'$\log G$')

        axarr[2].plot(t[lowerBound:upperBound], G[lowerBound:upperBound] - g_fit[lowerBound:upperBound])
        axarr[2].plot(t[lowerBound:upperBound], np.zeros(upperBound - lowerBound), color='red')
        axarr[2].set_xlabel(r'$\tau$')
        axarr[2].set_ylabel(r'$G - G_{\mathrm{fit}}$')

        axarr[3].plot(t[lowerBound:upperBound], diff)
        axarr[3].plot(t[lowerBound:upperBound], np.zeros(upperBound - lowerBound), color='red')
        axarr[3].plot(t[lowerBound:upperBound], cumsum)
        axarr[3].set_xlabel(r'$\tau$')
        axarr[3].set_ylabel(r'$ \log G - \log G_{\mathrm{fit}}$')

        plt.tight_layout()
        plt.show()



    diff = y[lowerBound:upperBound] - lin_fit([m, k], t[lowerBound:upperBound])
    signs = np.sign(diff)
    cumsum = np.cumsum(diff)*dt

    # if the fit is extremely good
    if np.max(cumsum) < 1e-6:
      break



    # otherwise find intersection
    for i, s in enumerate(signs):

      if s != signs[0]:
        if i > 1:
          lowerBound += i;
          break
        else:
          fittingData = False
          break

      if i == signs.size - 1:
        fittingData = False


  Z = np.exp(m)
  E = mu - k

  g_fit = Z*np.exp(-(E - mu)*t)

  cumsum = np.cumsum(diff)*dt

  if plot:

    f, axarr = plt.subplots(4, sharex=True, figsize=(20, 10))

    axarr[0].set_title(r'$p = {0}, \, \mu = {1}, \, E_0 = {2}, \, Z_0 = {3}$'.format(p, mu, E, Z))

    axarr[0].plot(t, G)
    axarr[0].plot(t[lowerBound], G[lowerBound], 'o', color='green')
    axarr[0].plot(t, Z*np.exp(-(E - mu)*t), '--', color='red')
    axarr[0].plot(t[upperBound], G[upperBound], 'o', color='green')
    axarr[0].set_xlabel(r'$\tau$')
    axarr[0].set_ylabel(r'$G$')

    axarr[1].semilogy(t, G)
    axarr[1].semilogy(t[lowerBound], G[lowerBound], 'o', color='green')
    axarr[1].semilogy(t, Z*np.exp(-(E - mu)*t), '--', color='red')
    axarr[1].semilogy(t[upperBound], G[upperBound], 'o', color='green')
    axarr[1].set_xlabel(r'$\tau$')
    axarr[1].set_ylabel(r'$\log G$')

    axarr[2].plot(t[lowerBound:upperBound], G[lowerBound:upperBound] - g_fit[lowerBound:upperBound])
    axarr[2].plot(t[lowerBound:upperBound], np.zeros(upperBound - lowerBound), color='red')
    axarr[2].set_xlabel(r'$\tau$')
    axarr[2].set_ylabel(r'$G - G_{\mathrm{fit}}$')

    axarr[3].plot(t[lowerBound:upperBound], diff)
    axarr[3].plot(t[lowerBound:upperBound], np.zeros(upperBound - lowerBound), color='red')
    axarr[3].plot(t[lowerBound:upperBound], cumsum)
    axarr[3].set_xlabel(r'$\tau$')
    axarr[3].set_ylabel(r'$ \log G - \log G_{\mathrm{fit}}$')

    plt.tight_layout()
    plt.savefig('plots/plot.pdf')
    plt.show()

  maxCumsum = np.max(cumsum)
  if minFittingLength < (upperBound - lowerBound)*dt or maxCumsum < 0.001:
    return E, Z
  else:
    print('percentage of length:' + str((upperBound - lowerBound)*dt/length), 'maxCumsum=' + str(maxCumsum))
    return False, False
