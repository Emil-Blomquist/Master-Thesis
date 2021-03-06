from __future__ import division
import numpy as np

def read1D (path):
  with open(path) as fp:

    t = np.array([])
    H = np.array([])
    scaleFactor = 0
    parameters = {}
    stringCharacters = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', ' ']

    for i, line in enumerate(fp):

      if i == 0:
        # fetch the parameters from the first line
        for d in line[9: -10].split(' '):
          separated = d.split('=', 1)
          if len(separated) == 2:
            key = separated[0]
            parameters[key] = separated[-1]
          else:
            # if the value is space separated
            parameters[key] += ' ' + separated[0]

        # convert string to numbers
        for key, value in parameters.items():
          if not any(char in value for char in stringCharacters):
            parameters[key] = float(parameters[key])

      elif i == 2:
        # useful quantities
        Nt = int(round(parameters['tmax']/parameters['dt']))
        dt = parameters['dt']
        t = np.linspace(0, parameters['tmax'] - dt, Nt) + 0.5*dt

        # parse line
        line = line.split('\n')[0]
        line = line.split(' ')

        H = [float(part) for part in line]

        return t, H, parameters
  
  # pi not found
  return [], [], parameters