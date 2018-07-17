#!/usr/bin/python
import math, sys, os, glob

from sympy import *
from pylab import *

import numpy as npy
import scipy as spy
import scipy.interpolate
import matplotlib as mp
import matplotlib.pyplot as pl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib.ticker import MultipleLocator

ev = []
ev2=[]
coupling = npy.arange(-1, 1, 0.05)

for g_val in coupling:
    H1 = matrix([[2-g_val , -g_val/2.,  -g_val/2., -g_val/2., -g_val/2.,     0],
                         [-g_val/2.,   4-g_val,  -g_val/2., -g_val/2.,    0., -g_val/2.],
                         [-g_val/2., -g_val/2.,    6-g_val,     0, -g_val/2., -g_val/2.],
                 [-g_val/2., -g_val/2.,      0,   6-g_val, -g_val/2., -g_val/2.],
                 [-g_val/2.,     0,  -g_val/2., -g_val/2.,   8-g_val, -g_val/2.],
                 [0    , -g_val/2.,  -g_val/2., -g_val/2., -g_val/2.,  10-g_val]])
    
    u1, v1 = linalg.eig(H1)
    ev.append(min(u1))
    ev2.append(min(u1)-2+g_val)
        


fig, ax = pl.subplots()
#pl.plot(coupling, ev)
pl.plot(coupling, ev2)

pl.xlim(-1.1, 1.1)
ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(0.5))

pl.ylim(-0.6, 0.05)
ax.yaxis.set_minor_locator(MultipleLocator(0.05))
ax.yaxis.set_major_locator(MultipleLocator(0.2))

pl.xlabel('coupling g')
pl.ylabel('Correlation energy')


pl.savefig('pairing_model_Hamiltonian', format='PDF')                  #save as pdf

    


    
