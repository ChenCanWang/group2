import numpy as np
import math, sys, os, glob
import matplotlib.pyplot as plt

sys.setrecursionlimit(30000)

#constants
mneutron = 939.56541
mnucleon = 938.92
mass = mnucleon
hbarc = 197.33


def nmaxcomp(Nmax, lspace, cutoff):
    if (Nmax*2*math.pi)/lspace < cutoff:
        print('!!! Nmax too small !!!')
#set up single-particle basis
def spstates(Nmax):
    dim=(2*Nmax+1)**3*2
    spstates=np.zeros(shape=(dim,5))
    count=0
    for i in range(0, 2*Nmax+1):
        for j in range(0, 2*Nmax+1):
            for k in range(0, 2*Nmax+1):
                for l in (0,1):
                    spstates[count][0]=i-Nmax
                    spstates[count][1]=j-Nmax
                    spstates[count][2]=k-Nmax
                    spstates[count][3]=2*l-1
                    radius = (i-Nmax)**2+(j-Nmax)**2+(k-Nmax)**2
                    spstates[count][4] = radius
                    count += 1

    spstates=spstates[np.argsort(spstates[:,4])]
    fcutsp=magica
    holerad = int(spstates[fcutsp-1][4])
    return spstates, fcutsp, dim
#set up hole hole channels with total momentum P
#set Pmax to large number here, e.g. 60 for the hole states
def setchannelhh(Nmax,fcutsp,spstates):
    channel=np.zeros(shape=(300,4))
    chcount=0
    for px in range(-2*Nmax, 2*Nmax+1):
        for py in range(-2*Nmax, 2*Nmax+1):
            for pz in range(-2*Nmax, 2*Nmax+1):
                count=0
                for i in range(fcutsp):
                    for j in range(fcutsp):
                        if px == spstates[i][0]+ spstates[j][0] and \
                        py == spstates[i][1]+ spstates[j][1] and \
                        pz == spstates[i][2]+ spstates[j][2] and \
                        i !=j:
                            count += 1
                if count == 0:
                    continue
                else:
                    channel[chcount][0]= px
                    channel[chcount][1]= py
                    channel[chcount][2]= pz
                    channel[chcount][3]= count
                    chcount += 1
    #delete empty entries from channel for total momentum
    j=0
    for i in range(len(channel)):
        if channel[j][3] == 0:
            channel=np.delete(channel, j, 0)
        else:
            j += 1
            if j == len(channel):
                break
    return channel
#set up two-body relative momentum basis with all states in the order of P
def settbbasishh(channel,fcutsp,spstates):
    tbrel=np.zeros(shape=(int(sum(channel[:,3])),5))
    count = 0
    for l in range(len(channel)):
        for i in range(fcutsp):
            for j in range(fcutsp):
                if channel[l][0] == spstates[i][0]+ spstates[j][0] and \
                channel[l][1] == spstates[i][1]+ spstates[j][1] and \
                channel[l][2] == spstates[i][2]+ spstates[j][2] and \
                i != j:
                    tbrel[count][0] = (spstates[i][0]-spstates[j][0])/2
                    tbrel[count][1] = (spstates[i][1]-spstates[j][1])/2
                    tbrel[count][2] = (spstates[i][2]-spstates[j][2])/2
                    tbrel[count][3] = spstates[i][3]
                    tbrel[count][4] = spstates[j][3]
                    count += 1
                else:
                    continue
    return tbrel
#now set up particle states with equal total momentum
def setchannelpp(channel,fcutsp,dim,spstates):
    channelpp = np.zeros(shape=(len(channel),4))
    for l in range(len(channel)):
        count = 0
        for i in range(fcutsp, dim):
            for j in range(fcutsp, dim):
                if channel[l][0] == spstates[i][0]+ spstates[j][0] and \
                channel[l][1] == spstates[i][1]+ spstates[j][1] and \
                channel[l][2] == spstates[i][2]+ spstates[j][2] and \
                i != j:
                    count += 1
                if count == 0:
                    continue
                else:
                    channelpp[l][0]= channel[l][0]
                    channelpp[l][1]= channel[l][1]
                    channelpp[l][2]= channel[l][2]
                    channelpp[l][3]= count
    return channelpp
#set up two-body relative momentum basis for particle particle case order as given above
def settbbasispp(channel,fcutsp,dim,spstates):
    tbrelpp=np.zeros(shape=(int(sum(channelpp[:,3])),5))
    count = 0
    for l in range(len(channel)):
        for i in range(fcutsp, dim):
            for j in range(fcutsp, dim ):
                if channel[l][0] == spstates[i][0]+ spstates[j][0] and \
                channel[l][1] == spstates[i][1]+ spstates[j][1] and \
                channel[l][2] == spstates[i][2]+ spstates[j][2] and \
                i != j :
                    tbrelpp[count][0] = (spstates[i][0]-spstates[j][0])/2
                    tbrelpp[count][1] = (spstates[i][1]-spstates[j][1])/2
                    tbrelpp[count][2] = (spstates[i][2]-spstates[j][2])/2
                    tbrelpp[count][3] = spstates[i][3]
                    tbrelpp[count][4] = spstates[j][3]
                    count += 1
                else:
                    continue
    return tbrelpp
#set up particle hole channels with total momentum P,set Pmax to largest possible number
def setchannelph(Nmax,fcutsp,dim,spstates):
    channelph=np.zeros(shape=((2*(Nmax+holerad)+1)**3,4))
    chcount=0
    for px in range(-2*Nmax, 2*Nmax+1):
        for py in range(-2*Nmax, 2*Nmax+1):
            for pz in range(-2*Nmax, 2*Nmax+1):
                count=0
                for i in range(fcutsp,dim):
                    for j in range(fcutsp):
                        if px == spstates[i][0]+ spstates[j][0] and \
                        py == spstates[i][1]+ spstates[j][1] and \
                        pz == spstates[i][2]+ spstates[j][2]:
                            count += 1
                if count == 0:
                    continue
                else:
                    channelph[chcount][0]= px
                    channelph[chcount][1]= py
                    channelph[chcount][2]= pz
                    channelph[chcount][3]= count
                    chcount += 1
    #delete empty entries from channel for total momentum
    k=0
    for l in range(len(channelph)):
        if channelph[k][3] == 0:
            channelph=np.delete(channelph, k, 0)
        else:
            k += 1
            if k == len(channelph):
                break
    return channelph
# !!!  ph channel so far not needed, so not called
def settbbasisph(channelph,fcutsp,dim,spstates):
    tbrelph=np.zeros(shape=(int(sum(channelph[:,3])),5))
    count = 0
    for l in range(len(channelph)):
        for i in range(fcutsp,dim):
            for j in range(fcutsp):
                if channelph[l][0] == spstates[i][0]+ spstates[j][0] and \
                channelph[l][1] == spstates[i][1]+ spstates[j][1] and \
                channelph[l][2] == spstates[i][2]+ spstates[j][2] and \
                i != j :
                    tbrelph[count][0] = (spstates[i][0]-spstates[j][0])/2
                    tbrelph[count][1] = (spstates[i][1]-spstates[j][1])/2
                    tbrelph[count][2] = (spstates[i][2]-spstates[j][2])/2
                    tbrelph[count][3] = spstates[i][3]
                    tbrelph[count][4] = spstates[j][3]
                    count += 1
                else:
                    continue
    return tbrelph
#define interaction, here we use the Minnesota potential in momentum space (see lecture notes)
def minnesota(tbstatep,tbstatek,lspace):
    '''define minnesoat potential, here special case for spin singlet states only'''
    vrad = 200
    vsing = -91.85
    vtrip = -178
    krad = 1.487
    ksing = 0.465
    ktrip = 0.639
    factork = 2*math.pi/lspace
    qsquar = factork**2*((tbstatep[0]-tbstatek[0])**2 + (tbstatep[1]-tbstatek[1])**2 + (tbstatep[2]-tbstatek[2])**2)
    qsquarp = factork**2*((tbstatep[0]+tbstatek[0])**2 + (tbstatep[1]+tbstatek[1])**2 + (tbstatep[2]+tbstatek[2])**2)
    vint = vrad / (lspace**3) * (math.pi/krad)**(3/2)*math.exp(-qsquar/(4*krad)) \
         + vsing / (lspace**3) * (math.pi/ksing)**(3/2)*math.exp(-qsquar/(4*ksing)) \
         + vrad / (lspace**3) * (math.pi/krad)**(3/2)*math.exp(-qsquarp/(4*krad)) \
              + vsing / (lspace**3) * (math.pi/ksing)**(3/2)*math.exp(-qsquarp/(4*ksing))
    if tbstatep[3] == tbstatep[4] or tbstatek[3] == tbstatek[4]:
        return 0
    elif tbstatep[3] == tbstatek[3]:
        return 1/2*vint
    elif tbstatep[3] == tbstatek[4]:
        return -1/2*vint
#set up vhhhh for all total momenta P
def fillhhhh(channel,tbrel,lspace):
    '''fill hhhh interaction list, with each listelement corresponding to an
    array for a block of conserved total momentum P'''
    vhhhhlist=[]
    count = 0
    for l in range(len(channel)):
        vhhhh = np.zeros(shape=(int(channel[l][3]),int(channel[l][3])))
        for i in range(count,int(channel[l][3])+count):
            for j in range(count,int(channel[l][3])+count):
                vhhhh[i-count][j-count] = minnesota(tbrel[i],tbrel[j],lspace)
        vhhhhlist.append(vhhhh)
        count += int(channel[l][3])
    np.save('me_minnesota/V_hhhh_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho),vhhhhlist)
#set up vpppp for all total momenta P
def fillpppp(channelpp,tbrelpp,lspace):
    '''fill pppp interaction list, with each listelement corresponding to an
    array for a block of conserved total momentum P'''
    vpppplist=[]
    count = 0
    for l in range(len(channelpp)):
        vpppp = np.zeros(shape=(int(channelpp[l][3]),int(channelpp[l][3])))
        for i in range(count,int(channelpp[l][3])+count):
            for j in range(count,int(channelpp[l][3])+count):
                vpppp[i-count][j-count] = minnesota(tbrelpp[i],tbrelpp[j],lspace)
        vpppplist.append(vpppp)
        count += int(channelpp[l][3])
    np.save('me_minnesota/V_pppp_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho),vpppplist)
#set up vpphh for all total momenta P
def fillpphh(channel,channelpp,tbrel,tbrelpp,lspace):
    '''fill pphh interaction list, with each listelement corresponding to an
    array for a block of conserved total momentum P'''
    vpphhlist=[]
    countpp=0
    counthh=0
    for l in range(len(channel)):
        vpphh = np.zeros(shape=(int(channelpp[l][3]),int(channel[l][3])))
        for i in range(countpp,int(channelpp[l][3])+countpp):
            for j in range(counthh,int(channel[l][3])+counthh):
                vpphh[i-countpp][j-counthh] = minnesota(tbrelpp[i],tbrel[j],gridspace)
        vpphhlist.append(vpphh)
        counthh += int(channel[l][3])
        countpp += int(channelpp[l][3])
    np.save('me_minnesota/V_pphh_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho),vpphhlist)
# !!! vphph not used so far
def fillphph(channelph,tbrelph,lspace):
    '''fill phph interaction list, with each listelement corresponding to an
    array for a block of conserved total momentum P, but this P can be different
    form the hhhh, pppp, pphh total P'''
    vphphlist=[]
    count = 0
    for l in range(len(channelph)):
        vphph = np.zeros(shape=(int(channelph[l][3]),int(channelph[l][3])))
        for i in range(count,int(channelph[l][3])+count):
            for j in range(count,int(channelph[l][3])+count):
                vphph[i-count][j-count] = minnesota(tbrelph[i],tbrelph[j],gridspace)
        vphphlist.append(vphph)
        count += int(channelph[l][3])
    np.save('me_minnesota/V_phph_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho),vphphlist)
# load interaction files from file
def loadint(Nmax, magica, density):
    '''load all hhhh, pppp, and pphh interations from file'''
    vhhhhlist=np.load('me_minnesota/V_hhhh_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho) + '.npy')
    vpppplist=np.load('me_minnesota/V_pppp_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho) + '.npy')
    vpphhlist=np.load('me_minnesota/V_pphh_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho) + '.npy')
    #vphphlist=np.load('me_minnesota/V_phph_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho) + '.npy')
    return vhhhhlist, vpppplist, vpphhlist
#set up fock matrices for the pp and hh case
#use sign 0, and 1
def fockhh(tbstatep, postotalP, sign, lspace, channel, magica, hbarc, mass):
    ksp = np.zeros(shape=(1,4))
    factork = 2*math.pi/lspace
    ksp[0][0] = channel[postotalP][0]/2 + (-1)**sign * tbstatep[0]
    ksp[0][1] = channel[postotalP][1]/2 + (-1)**sign * tbstatep[1]
    ksp[0][2] = channel[postotalP][2]/2 + (-1)**sign * tbstatep[2]
    ksp[0][3] = tbstatep[3+sign]
    sumint = 0
    for i in range(magica):
        relp = np.zeros(shape=(1,5))
        relp[0][0] = (ksp[0][0] - spstates[i][0])/2
        relp[0][1] = (ksp[0][1] - spstates[i][1])/2
        relp[0][2] = (ksp[0][2] - spstates[i][2])/2
        relp[0][3] = ksp[0][3]
        relp[0][4] = spstates[i][3]
        sumint += minnesota(relp[0],relp[0],lspace)
    fhh = (ksp[0][0]**2+ksp[0][1]**2+ksp[0][2]**2)*(factork**2*hbarc**2/(2*mass)) \
        + sumint
    return fhh
def fockpp(tbstatep, postotalP, sign, lspace, channelpp, magica, hbarc, mass):
    ksp = np.zeros(shape=(1,4))
    factork = 2*math.pi/lspace
    ksp[0][0] = channelpp[postotalP][0]/2 + (-1)**sign * tbstatep[0]
    ksp[0][1] = channelpp[postotalP][1]/2 + (-1)**sign * tbstatep[1]
    ksp[0][2] = channelpp[postotalP][2]/2 + (-1)**sign * tbstatep[2]
    ksp[0][3] = tbstatep[3+sign]
    sumint = 0
    for i in range(magica):
        relp = np.zeros(shape=(1,5))
        relp[0][0] = (ksp[0][0] - spstates[i][0])/2
        relp[0][1] = (ksp[0][1] - spstates[i][1])/2
        relp[0][2] = (ksp[0][2] - spstates[i][2])/2
        relp[0][3] = ksp[0][3]
        relp[0][4] = spstates[i][3]
        sumint += minnesota(relp[0],relp[0],lspace)
    fpp = (ksp[0][0]**2+ksp[0][1]**2+ksp[0][2]**2)*(factork**2*hbarc**2/(2*mass)) \
        + sumint
    return fpp
#possibility to store fock matrices to file with loadfock and save in fill funcitons
def fillfockhh(tbrel,lspace,channel,magica,hbarc,mass):
    fhhlist=[]
    count = 0
    for l in range(len(channel)):
        fhh = np.zeros(shape=(int(channel[l][3]),3))
        for i in range(count,int(channel[l][3])+count):
            fhh[i-count][0] = i - count
            fhh[i-count][1] = fockhh(tbrel[i],l,0,lspace,channel,magica,hbarc,mass)
            fhh[i-count][2] = fockhh(tbrel[i],l,1,lspace,channel,magica,hbarc,mass)
        fhhlist.append(fhh)
        count += int(channel[l][3])
    #np.save('me_minnesota/f_hh_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho),fhhlist)
    return fhhlist
def fillfockpp(tbrelpp,lspace,channelpp,magica,hbarc,mass):
    fpplist=[]
    count = 0
    for l in range(len(channelpp)):
        fpp = np.zeros(shape=(int(channelpp[l][3]),3))
        for i in range(count,int(channelpp[l][3])+count):
            fpp[i-count][0] = i - count
            fpp[i-count][1] = fockpp(tbrelpp[i],l,0,lspace,channelpp,magica,hbarc,mass)
            fpp[i-count][2] = fockpp(tbrelpp[i],l,1,lspace,channelpp,magica,hbarc,mass)
        fpplist.append(fpp)
        count += int(channelpp[l][3])
    #np.save('me_minnesota/f_pp_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho),fpplist)
    return fpplist
def loadfock(Nmax, magica, density):
    '''load fock matrices for pp and hh from file'''
    fhhlist=np.load('me_minnesota/f_hh_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho) + '.npy')
    fpplist=np.load('me_minnesota/f_pp_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho) + '.npy')
    return fhhlist, fpplist
#calculate HF energy
def hfe(spstates,tbrel,lspace,magica,hbarc,mass):
    hfkin = 0
    hfint = 0
    faktlspace = 2*math.pi/lspace
    for i in range(magica):
        hfkin += (spstates[i][0]**2 + spstates[i][1]**2 + spstates[i][2]**2) \
            *(faktlspace**2*hbarc**2/(2*mass))
    for i in range(len(tbrel)):
            hfint += 1/2*minnesota(tbrel[i],tbrel[i],lspace)

    return (hfkin + hfint) / magica
#start initial T2 guess from MBPT2
def filltcoef(channel,channelpp,vpphhlist,fhhlist,fpplist):
    tstartlist = []
    countpp=0
    counthh=0
    for l in range(len(channel)):
        tstart = np.zeros(shape=(int(channelpp[l][3]),int(channel[l][3])))
        for i in range(countpp,int(channelpp[l][3])+countpp):
            for j in range(counthh,int(channel[l][3])+counthh):
                tstart[i-countpp][j-counthh] = vpphhlist[l][i-countpp][j-counthh] \
                / (fhhlist[l][j-counthh][1]+fhhlist[l][j-counthh][2]
                -fpplist[l][i-countpp][1]-fpplist[l][i-countpp][2])
        tstartlist.append(tstart)
        counthh += int(channel[l][3])
        countpp += int(channelpp[l][3])
    return tstartlist
# set up the Hamiltonian
def fillhbar(channel,channelpp,vpphhlist,tlist,fpplist,fhhlist,vpppplist,vhhhhlist):
    hbarlist=[]
    countpp=0
    counthh=0
    for l in range(len(channel)):
        hbar = np.zeros(shape=(int(channelpp[l][3]),int(channel[l][3])))
        hbar = vpphhlist[l] \
             + np.einsum('ai,a->ai', tlist[l], fpplist[l][:,2]) \
             + np.einsum('ai,a->ai', tlist[l], fpplist[l][:,1]) \
             - np.einsum('ai,i->ai', tlist[l], fhhlist[l][:,2]) \
             - np.einsum('ai,i->ai', tlist[l], fhhlist[l][:,1]) \
             + 1/2*np.einsum('ac,ci->ai', vpppplist[l], tlist[l]) \
             + 1/2*np.einsum('ki,ak->ai', vhhhhlist[l], tlist[l])
        hbarlist.append(hbar)
    return hbarlist
# fockmatrix in pphh format, saved to file, then loaded
def fillfockmatrix(channel,channelpp,fhhlist,fpplist):
    fockmatrixlist=[]
    counthh = 0
    countpp = 0
    for l in range(len(channel)):
        fpphhmatrix = np.zeros(shape=(int(channelpp[l][3]),int(channel[l][3])))
        for i in range(countpp,int(channelpp[l][3])):
            for j in range(counthh,int(channel[l][3])):
                fpphhmatrix[i-countpp][j-counthh] = \
                    fhhlist[l][j-counthh][1]+fhhlist[l][j-counthh][2] \
                    -fpplist[l][i-countpp][1]-fpplist[l][i-countpp][2]
        fockmatrixlist.append(fpphhmatrix)
    np.save('me_minnesota/fockmatrixlist_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho),fockmatrixlist)
def loadfockmatrix(Nmax,magica,lspace):
    fockmatrixlist=np.load('me_minnesota/fockmatrixlist_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '_dens_' + str(rho) + '.npy')
    return fockmatrixlist
#iterate and update T2 list
def iterative(channel,vpphhlist,tlist,magica,channelpp,fpplist,fhhlist,vpppplist,vhhhhlist,fockmatrixlist):
    ecorr = 0
    ecorrnew = 0
    mix = 1
    for l in range(len(channel)):
        ecorr += 1/4.0*np.einsum('ia,ai->', vpphhlist[l].T, tlist[l]) / magica

        tlist[l] = tlist[l] + mix * fillhbar(channel,channelpp,vpphhlist,tlist,fpplist,fhhlist,vpppplist,vhhhhlist)[l] \
        / fockmatrixlist[l]

    for l in range(len(channel)):
        ecorrnew += 1/4.0*np.einsum('ia,ai->', vpphhlist[l].T, tlist[l]) / magica
    return ecorrnew, abs(ecorrnew-ecorr), tlist


#if save=1 -> create and save to file if save=0 only load from file
save = 1
rholist = np.linspace(0.01,0.18,18)
#set input paramters
magica = 14
Nmax=2

spstates, fcutsp, dim = spstates(Nmax)

hfelist = []
ecorrlist = []
f = open('data/correlation_energy_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '.txt','w')
g = open('data/hf_energy_minnesota_Nmax_' + str(Nmax) + '_A_' + str(magica) + '.txt','w')
for rho in rholist:
    print('density ', round(rho,3))
    rho = round(rho,3)
    #rho = 0.08

    gridspace = (magica / rho)**(1/3)
    cutoff = (2*0.456)**(1/2)

    nmaxcomp(Nmax, gridspace, cutoff)

    channel = setchannelhh(Nmax,fcutsp,spstates)
    tbrel = settbbasishh(channel,fcutsp,spstates)
    channelpp = setchannelpp(channel,fcutsp,dim,spstates)
    tbrelpp = settbbasispp(channel,fcutsp,dim,spstates)
    if save == 1:
        # save interaction matrices to file
        fillhhhh(channel,tbrel,gridspace)
        fillpppp(channelpp,tbrelpp,gridspace)
        fillpphh(channel,channelpp,tbrel,tbrelpp,gridspace)
        #fillphph(Nmax,magica,rho)

    vhhhhlist, vpppplist, vpphhlist = loadint(Nmax,magica,rho)

    fhhlist = fillfockhh(tbrel,gridspace,channel,magica,hbarc,mass)
    fpplist = fillfockpp(tbrelpp,gridspace,channelpp,magica,hbarc,mass)

    hfenergy = hfe(spstates,tbrel,gridspace,magica,hbarc,mass)
    print('Hartree-Fock Energy ', hfenergy)

    tstartlist = filltcoef(channel,channelpp,vpphhlist,fhhlist,fpplist)

    if save == 1:
        # save fock matrix to file
        fillfockmatrix(channel,channelpp,fhhlist,fpplist)
    #load fockmatrix form file
    fockmatrixlist = loadfockmatrix(Nmax,magica,gridspace)

    #iteration parameters
    itermax=500
    eps=1e-9
    titer = tstartlist

    for iter in range(itermax):
    	result=iterative(channel,vpphhlist,titer,magica,channelpp,fpplist,fhhlist,vpppplist,vhhhhlist,fockmatrixlist)
    	ecorr=result[0]
    	deltecorr=result[1]
    	if iter < 10:
    		print('step ', iter, '   correlation energy  ', round(ecorr,6), '  delta E  ', round(deltecorr,6))
    	elif iter > 9:
    		print('step ', iter, '  correlation energy  ', round(ecorr,6), '  delta E  ', round(deltecorr,6))
    	elif iter > 99:
    		print('step ', iter, ' correlation energy  ', round(ecorr,6), '  delta E  ', round(deltecorr,6))

    	titer=result[2]

    	if deltecorr < eps:
    		break
    	if iter == itermax-1 and deltecorr > eps:
    		print('!!! no converged results after ', iter + 1, ' iteration steps !!!')

    print('---> Total E/N = ', hfenergy + ecorr, '\n')
    print(rho, ecorr, file = f)
    print(rho, hfenergy, file = g)



f.close()
g.close()
