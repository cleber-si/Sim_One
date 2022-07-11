# -*- coding: utf-8 -*-
"""Sim_Two.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1EasJjhr9O1DjzrBQhB2D4wsnG9UX-ysJ

# Libraries
"""

from scipy.optimize import fsolve
import numpy as np
from scipy.fft import fft, ifft

import matplotlib.pyplot as plt

"""# Gera Ruído"""

def ffgn(flux_, sigma = 0.0006, H = 0.500001, n = 1, force = False):
  '''
   Generates *exact* paths of Fractional Gaussian Noise by using
   circulant embedding (for 1/2<H<1) and Lowen's method (for 0<H<1/2).  
  
   Input: 
     sigma <- standard deviation of the time series
     H     <- Hurst exponent
     n     <- number of independent paths
     N     <- the length of the time series
     force <- if *1*, then works with the smallest integer power of 2,
              greater than or equal to N.
  
   Output:
     f     <- a (nxN) matrix.  Each row contains an independently generated
              path of FGN.
   
   Usage:
     f = ffgn(sigma,H,n,N,force);
  
   Written jointly by Yingchun Zhou (Jasmine), zhouyc@math.bu.edu 
   and Stilian Stoev, sstoev@umich.edu, September, 2005.
  
   Example:
  
    f = ffgn(1,0.8,1,2^16,0); plot(cumsum(f)); title('FBM, H=0.8');
  '''
  
  flux = np.copy(flux_)
  N = len(flux)

  if (H>0.5) and (H<1):
    '''
      Use the "circulant ebedding" technique.  This method
      works only in the case when 1/2 < H < 1.
    '''

    if force:
      N2 = 2**(np.floor(np.log2(N)));
      if N2 < N:
        N = 2*N2
      else:
        N  = N2

    # First step: specify the covariance
    Ns = np.linspace(0, N, N+1)
    c = sigma**2/2 * ((Ns+1)**(2*H) - 2*(Ns**(2*H)) + abs(Ns-1)**(2*H))
    v = [c[0:N], c[N+1:0:-1]]
    v = np.concatenate(v, axis=None)


    # Second step: calculate Fourier transform of c
    g = (fft(v)).real
    if min(g) < 0:
      raise Exception (' Some of the g(k) are negative!')
    g = abs(g)


    # Third step: generate {z(1),...,z(2*N)}
    c1 = np.sqrt(2)*np.random.normal()
    c2 = np.sqrt(2)*np.random.normal()

    z = np.zeros(N+1)
    z[0], z[-1] = c1, c2

    a = np.random.normal(size = (n, N-1))
    b = np.random.normal(size = (n, N-1))
    z1 = [complex(a[0][i], b[0][i]) for i in range(len(a[0]))]

    z = [complex(z[i], 0) for i in range(len(z))]

    for i in range(len(z)):
      if i != 0 and i != len(z)-1:
        z[i] = z1[i-1]

    y = np.copy(z)
    z1_flip = z1[::-1]
    y_sub = [z1_flip[i].conjugate() for i in range(len(z1_flip))]
    
    y = np.concatenate([y, y_sub], axis=None)
    y = y * (np.ones(n)*np.sqrt(g))


    # Fourth step: calculate the stationary process f
    f = (fft(y)).real / np.sqrt(4*N)
    f = f[0:N]

  elif H == 0.5:
    f = np.random.normal(size = (n,N))

  elif (H<0.5) and (H>0):
    f = np.random.normal(size = (n,N))

    G1 = np.random.normal(size = (n,N-1))
    G2 = np.random.normal(size = (n,N-1))
    G1_complex = [complex(0, G1[0][i]) for i in range(len(G1[0]))]
    G2_complex = [complex(0, G2[0][i]) for i in range(len(G2[0]))]

    G = [(G1_complex[i] + G2_complex[i])/np.sqrt(2) for i in range(len(G1[0]))]

    GN = np.random.normal(size = (n,1))[0][0]
    G0 = np.zeros(n)

    H2 = 2*H

    R=(1-(np.linspace(1,N-1, N-1)/N)**H2)
    R = np.concatenate([[1], R, [0], R[::-1]])

    S = np.ones(n)*(abs(fft(R, 2*N))**0.5)

    G_flip = G[::-1]
    G_sub = [G_flip[i].conjugate() for i in range(len(G_flip))]
    X = np.concatenate([[complex(0,0)], G, GN, G_sub], axis=None)
    X *= S

    x = ifft(X,2*N)
    y = np.sqrt(N)*((x[0:N] - x[0]*np.ones(N))).real
    f = sigma*N**H*np.concatenate([[y[0]], np.diff(y)], axis=None)

  else:
    raise Exception(' The value of the Hurst parameter H must be in (0,1).')

  return f

"""# Variables"""

''' Variables Declaration '''
# Solar parameters
M_sun = 1.989 * 10**30                          # Sun's mass in kg
M_earth = 5.972 * 10**24
M_jup = 1.899 * 10**27

# Input Parameters
P_orb = 4.23                                      # Orbital period in days
ai = 90                                         # Inclination angle of star (ranged over 0 and 90 degrees)
omega = (np.pi/180)*(0)                       # Aargument of pericenter in true orbital plane (rad)
excent = 0                                   # Eccentricity
tp = 2449797.77

np_ = 0.5                                         # Planet mass fraction (in Jupter mass unit)
M_p = np_ * M_jup                               # Planet mass

ns = 1.11                                      # Stellar mass fraction (in Sun mass unit)
M_star = ns * M_sun                             # Stellar mass

"""# Execution"""

# Eq Torres
K = 203.29 * ((M_p * np.sin(np.pi/180*ai))/M_jup) * ((M_sun/(M_star + M_p))**(2/3)) * ((P_orb)**(-1/3)) * (1/np.sqrt(1-excent**2))

Q = P_orb                  			  # Quarter duration (in days)
cadence = 0.020833      			  # Kepler's cadence (in days)
n = int(Q/cadence)					    # Number of points based on cadence and quarter lenght
times = np.linspace(0, n, n+1)		# Creates time array with lenght equals to n+1
times *= cadence						      # Ajust the time array dimension to be in days

'''
cadence = 8
times = np.linspace(1, 100, 100)
range_ = P_orb/len(times)
times *= range_
'''


def transcendental(Et, excent, eqLeft):
  return Et - excent*np.sin(Et) - eqLeft

EtS = []

for time in times:
  eqLeft = (2*np.pi/P_orb) * (time-tp)                                       # Mass function
  res = fsolve(lambda x: transcendental(x, excent, eqLeft), x0=1)
  if time == 20:
    print(res[0])
  #print(res[0], '\n')
  EtS.append(res[0])

ni_t = np.arccos((np.cos(EtS) - excent) / (1 - excent*np.cos(EtS)))

#times = times + 1
#times = P_orb/times

RV = K * (np.cos(omega + ni_t) + excent*np.cos(omega))

tam = 16

def plotar(t, y, yerr=14, save = False):
  plt.rcParams["figure.figsize"] = (16,8)
  plt.plot(t, y, '.k')
  plt.errorbar(t, y, yerr=yerr, fmt=".k", capsize=4)

  # Curve Fitting
  model5 = np.poly1d(np.polyfit(t, y, 5))
  polyline = np.linspace(0, t[-1], 100)
  plt.plot(polyline, model5(polyline), color='red')

  plt.xlabel('Time (days)', fontsize = tam+5)
  plt.xticks(fontsize = tam)

  plt.ylabel('Radial Velocity (m/s)', fontsize = tam+5)
  plt.yticks(fontsize = tam)

  plt.title('Radial Velocity Reconstruction of 51 Peg', fontsize = tam+10)

  folder = '/content/drive/MyDrive/_Pesquisa/Astrofísica Observacional'
  path = folder + '/RV.pdf'
  if save:
    plt.savefig(path, bbox_inches = 'tight')

  plt.show()


def plotar_fase(f, y, yerr=14, tam = 16, save = False):
  plt.plot(f, y, '.k')
  plt.errorbar(f, y, yerr=yerr, fmt=".k", capsize=3)

  # Curve Fitting
  model5 = np.poly1d(np.polyfit(f, y, 5))
  polyline = np.linspace(0, 1, 100)
  plt.plot(polyline, model5(polyline), color='red')

  plt.xlabel('Orbital Phase $\phi$', fontsize = tam+5)
  plt.ylabel('Radial Velocity (m/s)', fontsize = tam+5)

  plt.title('Phase Reconstruction of 51 Peg', fontsize = tam+10)
  folder = '/content/drive/MyDrive/_Pesquisa/Astrofísica Observacional'
  path = folder + '/phaseRV.pdf'
  if save:
    plt.savefig(path, bbox_inches = 'tight')

  plt.show()

"""# Noise"""

''' Noisy time series '''

stdnoise = 11
noise = ffgn(RV, sigma = stdnoise)

RV_noise = RV + noise

# Remove periodically poins (2 in each 3 points)
rem = 6
RV_noise_red = RV_noise[0::rem]
times_red = times[0::rem]

plotar(times_red, RV_noise_red, save = True)

fase = ((times_red - np.mean(times_red)) % P_orb)/P_orb
plotar_fase(fase, RV_noise_red, save = True)