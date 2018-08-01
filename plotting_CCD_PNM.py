import numpy as np
import math, sys, os, glob
import scipy as spy
import scipy.interpolate
import matplotlib as mp
import matplotlib.pyplot as pl
from matplotlib import rc
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib.ticker import MultipleLocator

#--------------------
# RCParams
#--------------------
pl.rcParams['font.size'] = 10
pl.rcParams['font.family'] = 'Sans-Serif'
pl.rcParams['axes.labelsize'] =10
pl.rcParams['xtick.major.pad'] = 10
pl.rcParams['ytick.major.pad'] = 10

colors=["darkorange", "#009E73", "#0072B2",] # "firebrick"]          # "colorblind"

fig= pl.figure(figsize=(6,5),dpi=60)
fig.subplots_adjust(left=0.12, right=0.985, bottom=0.11, top=0.98, hspace=0.36, wspace=0)
ax=fig.add_subplot(1, 1, 1)
ax.set_facecolor('0.9')

#ecorr of hf energy
plotpar = 0

magica=14
magica2=66


i = 0
if plotpar == 0:
    for Nmax in range(1,3,1):
        datacorr = np.loadtxt ('data/correlation_energy_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '.txt')
        pl.plot( datacorr[:,0], datacorr[:,1], color = colors[i],  marker="o",\
         markersize = 4.5, linewidth=1, label ="N$_{\\rm{Max}}$=" + str(Nmax))
        i += 1
    i = 1
    #for Nmax in range(2,4,1):
    #    datacorr = np.loadtxt ('data/correlation_energy_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica2) + '.txt')
    #    pl.plot( datacorr[:,0], datacorr[:,1], color = colors[i], linestyle='--',  marker="o", markersize = 4.5, linewidth=1)
    #    i += 1
    pl.ylim( -1.5, 0)
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax.yaxis.set_major_locator(MultipleLocator(0.2))
    pl.xlim( 0.005, 0.185)
    ax.xaxis.set_minor_locator(MultipleLocator(0.005))
    ax.xaxis.set_major_locator(MultipleLocator(0.02))
    pl.tick_params(axis='both', which='major', labelsize=10)
    pl.xlabel('$\\rho$ (fm$^{-3}$)')
    pl.ylabel('E$_{\\rm{Corr}}$', fontsize=14)
    pl.legend(loc='lower left', fontsize=12,ncol=1,markerscale=1,numpoints=1,prop={'size':12})

    rc('text', usetex=True)
    pl.text(0.408, 0.1,'A=' +str(magica),
    horizontalalignment='center',
    verticalalignment='center',
    transform = ax.transAxes,
    fontsize = 10)
    pl.text(0.346, 0.102,'\\rule{0.5cm}{0.025cm}',
        horizontalalignment='center',
        verticalalignment='center',
        transform = ax.transAxes)
    pl.text(0.44,0.04, "$--$ A=" + str(magica2), horizontalalignment='right', transform=ax.transAxes, fontsize=10)

    pl.savefig('CCD_PNM_results_ecorr.pdf', format='PDF')


else:
    for Nmax in range(1,3,1):
        datahf = np.loadtxt ('data/hf_energy_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '.txt')
        datacorr = np.loadtxt ('data/correlation_energy_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '.txt')
        pl.plot( datahf[:,0], datahf[:,1]+datacorr[:,1], color = colors[i],\
          marker="o", markersize = 4.5, linewidth=1, label ="N$_{\\rm{Max}}$=" + str(Nmax))
        i += 1
    #i = 1
    #for Nmax in range(2,4,1):
    #    datahf = np.loadtxt ('data/hf_energy_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica2) + '.txt')
    #    datacorr = np.loadtxt ('data/correlation_energy_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica2) + '.txt')
    #    pl.plot( datahf[:,0], datahf[:,1]+datacorr[:,1], color = colors[i],\
    #     linestyle='--',  marker="o", markersize = 4.5, linewidth=1, label ="N$_{\\rm{Max}}$=" + str(Nmax))

    rc('text', usetex=True)
    pl.text(0.408, 0.1,'A=' +str(magica),
    horizontalalignment='center',
    verticalalignment='center',
    transform = ax.transAxes,
    fontsize = 10)
    pl.text(0.346, 0.102,'\\rule{0.5cm}{0.025cm}',
        horizontalalignment='center',
        verticalalignment='center',
        transform = ax.transAxes)
    pl.text(0.44,0.04, "$--$ A=" + str(magica2), horizontalalignment='right', transform=ax.transAxes, fontsize=10)

    pl.ylim( 0, 15)
    ax.yaxis.set_minor_locator(MultipleLocator(0.5))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    pl.xlim( 0.005, 0.185)
    ax.xaxis.set_minor_locator(MultipleLocator(0.005))
    ax.xaxis.set_major_locator(MultipleLocator(0.02))
    pl.tick_params(axis='both', which='major', labelsize=10)
    pl.ylabel('E$_{\\rm{HF}}$ + E$_{\\rm{Corr}}$', fontsize=14)
    pl.xlabel('$\\rho$ (fm$^{-3}$)')
    pl.legend(loc='lower left', fontsize=12,ncol=1,markerscale=1,numpoints=1,prop={'size':12})
    pl.savefig('CCD_PNM_results_ehf.pdf', format='PDF')
