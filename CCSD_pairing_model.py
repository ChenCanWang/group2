import numpy as np
#from sympy import *
import math, sys, os, glob
import matplotlib.pyplot as plt

#give input values
#where lines means degenerate states
particle = 4
line = 4
gval=0.7
delta = 1.0
#coupling = np.linspace(-1,1,11)


# set the particles to pairs
pair = particle/2
if(particle%2 != 0):
	print('Odd number of particles.')
	exit(0)

if(pair>line):
	print('Wrong configuration.')
	exit(0)

if(pair==1):
	states=[ i for i in range(line)]
elif(pair==2):
	states = [(i,j) for i in range(line)
	                for j in range(i,line) if j!=i]
elif(pair==3):
	states = [(i,j,k) for i in range(line)
	                  for j in range(i,line)
	                  for k in range(j,line)
	                  if j!=i and k!=j]
elif(pair==4):
	states = [(i,j,k,l) for i in range(line)
	                    for j in range(i,line)
	                    for k in range(j,line)
	                    for l in range(k,line)
	                    if j!=i and k!=j and l!=k]

elif(pair==5):
	states = [(i,j,k,l,m) for i in range(line)
	                      for j in range(i,line)
	                      for k in range(j,line)
	                      for l in range(k,line)
	                      for m in range(l,line)
	                      if j!=i and k!=j and l!=k and m!=l]

config_num  = len(states)
#print('There are %d states: \n'  %config_num)
#print(states)

# evaluate the sp operator
def h_sp(A):
	if(pair==1):
		return A
	else:
		return sum(A)

# evaluate the pair term
def h_pr(A,B):
	if(pair==1):
		if(A==B):
			return 2
		else:
			return 1
	else:
		if(len(A)!=len(B)):
			print('Something wrong.')
			exit(0)
		elif(pair!=1):
			n = len(A)
		s1 = set(A)
		s2 = set(B)
		m = len(s1&s2)

		if(n==m):
			return 2
		elif(n - m ==1):
			return 1
		else:
			return 0

# initialize a unit matrix with define dimension
H = np.eye(config_num)

Hc= [ ['123' for i in range(config_num)]
      for i in range(config_num)]

# The form of the matrix...
#print('The form of the matrix is: \n')
for i in range(config_num):
	for j in range(config_num):
		Hc[i][j]='-'+str(int(h_pr(states[i],states[j])))+'g'
		if(i==j):
			Hc[i][j] = Hc[i][j] \
                   + '+' \
			       + str(int(2*delta*h_sp(states[i])))   \
			       +'d'

#for i in range(len(states)):
#	print(Hc[i])



egv = [] # as the eigenvalue of H
crt = [] # as the correlation energy

# find the minimum of the eigenvalues
#for g in coupling:
g = gval
if g == gval:
# construction of the Hamiltonian
	for i in range(config_num):
		for j in range(config_num):
			H[i,j]=-g*h_pr(states[i],states[j])
			if(i==j):
				H[i,j] = H[i,j] + 2*delta*h_sp(states[i])

    # get the eigenvalues of H
	u, v = np.linalg.eig(H)
	egv.append(min(u))

	if(pair==1):
	    crt.append(min(u) + 2*g)
	else:
		crt.append(min(u)- sum(states[0])*2*delta + 2*g)
	print('\n', 'correlation energy by FCI     ', round(crt[0],4), '\n')

#####----set new part for CCD of pairing model----#####

#set fermi level, start counting states at 0
def flevel(number):
    return number-1

#set valence space defined by total number of lines, start counting at 0
def exspace(levels):
    return (2*levels)-1

fermi=flevel(particle)
vspace=exspace(line)

#fockt matrix, single particle part contributes only for diagonal
def fpq(p,q,gval):
	if p != q:
		return 0
	elif p in range(particle):
		return math.floor((p)/2.0)*delta - gval
	else:
		return math.floor((p)/2.0)*delta

#define two-body interaction and check of a,b and i,j are paired
#imply antisymmety for different cases
def vint(a,b,i,j,gval):
	if math.floor(a/2.0) == math.floor(b/2.0) and math.floor(i/2.0) == math.floor(j/2.0) and a!=b and i!=j:
		if a < b and i < j:
			return -gval
		if a < b and i > j:
			return gval
		if a > b and i < j:
			return gval
		if a > b and i > j:
			return -gval
	else:
		return 0

#differ pppp pphh (transpose to hhpp) and V_hhhh
#start counting a and b from Fermi level, so substract fermilevel in array
v_pphh = np.zeros(shape=(vspace-fermi, vspace-fermi, fermi+1, fermi+1))
v_pppp = np.zeros(shape=(vspace-fermi, vspace-fermi, vspace-fermi, vspace-fermi))
v_hhhh = np.zeros(shape=(fermi+1, fermi+1, fermi+1, fermi+1))

for matrix in (v_pphh, v_pppp, v_hhhh):
	for a in range(np.size(matrix,0)):
		for b in range(np.size(matrix,1)):
			for i in range(np.size(matrix,2)):
				for j in range(np.size(matrix,3)):
					matrix[a,b,i,j]=vint(a,b,i,j,gval)
v_hhpp = v_pphh.T

#defining the two fock matrices
f_hh=np.zeros(shape=(fermi+1,fermi+1))
f_pp=np.zeros(shape=(vspace-fermi,vspace-fermi))

for p in range(np.size(f_hh,0)):
	for q in range(np.size(f_hh,1)):
		f_hh[p,q] = fpq(p,q,gval)

for p in range(np.size(f_pp,0)):
	for q in range(np.size(f_pp,1)):
		f_pp[p,q] = fpq(p+particle,q+particle,gval)

#define T2 array for starting point, e.g. PT
tstart= np.zeros(shape=(vspace-fermi,vspace-fermi,fermi+1,fermi+1))
#set up Hamiltonian
hbar= np.zeros(shape=(vspace-fermi,vspace-fermi,fermi+1,fermi+1))

#starting point from PT for t0
for a in range(np.size(tstart,0)):
	   for b in range(np.size(tstart,1)):
		   for i in range(np.size(tstart,2)):
			   for j in range(np.size(tstart,3)):
				   tstart[a,b,i,j]= -v_pphh[a,b,i,j] / (f_pp[a,a]+f_pp[b,b]-f_hh[i,i]-f_hh[j,j])

#initialize division matrix in recursive formel for t
fsum= np.zeros(shape=(vspace-fermi,vspace-fermi,fermi+1,fermi+1))
for a in range(np.size(fsum,0)):
		  for b in range(np.size(fsum,1)):
			  for i in range(np.size(fsum,2)):
				  for j in range(np.size(fsum,3)):
					  fsum[a,b,i,j]= f_pp[a,a]+f_pp[b,b]-f_hh[i,i]-f_hh[j,j]

def tfunc(t):
	#intermediate chis for diagrams 7,8,9, and 10
	intermedchi_alid= 1/2.0*np.einsum('klcd,acik->alid', v_hhpp, t)
	intermedchi_li= 1/2.0*np.einsum('klcd,cdik->li', v_hhpp, t)
	intermedchi_ad= 1/2.0*np.einsum('klcd,ackl->ad', v_hhpp, t)
	intermedchi_klij= 1/4.0*np.einsum('klcd,cdij->klij', v_hhpp, t)

	hbar= v_pphh \
	 	+ np.einsum('bc,acij->abij', f_pp, t) - np.einsum('ac,bcij->abij', f_pp, t) \
		- np.einsum('kj,abik->abij', f_hh, t) + np.einsum('ki,abjk->abij', f_hh, t) \
		+ 1/2.0*np.einsum('abcd,cdij->abij', v_pppp, t) \
		+ 1/2.0*np.einsum('klij,abkl->abij', v_hhhh, t) \
		+ np.einsum('alid,dblj->abij', intermedchi_alid, t) - np.einsum('blid,dalj->abij', intermedchi_alid, t) \
		- np.einsum('aljd,dbli->abij', intermedchi_alid, t) + np.einsum('bljd,dali->abij', intermedchi_alid, t) \
		+ np.einsum('li,ablj->abij', intermedchi_li, t) - np.einsum('lj,abli->abij', intermedchi_li, t) \
		+ np.einsum('ad,dbij->abij', intermedchi_ad, t) - np.einsum('bd,daij->abij', intermedchi_ad, t) \
		+ np.einsum('klij,abkl->abij', intermedchi_klij, t)

	t = t - hbar / fsum
	return t

#iterative function to return new correlation energy
def iterative(t):
	ecorr= 1/4.0*np.einsum('ijab,abij->', v_hhpp, t)
	tnew=tfunc(t)
	ecorrnew= 1/4.0*np.einsum('ijab,abij->', v_hhpp, tnew)
	return (ecorrnew, abs(ecorr-ecorrnew), tnew)

#set maximal steps of iterations and convergence criterium
itermax=100
eps=0.00001
titer=tstart

for iter in range(itermax):
	t=titer
	result=iterative(t)
	ecorr=result[0]
	deltecorr=result[1]
	if iter < 10:
		print('step ', iter, '   correlation energy  ', round(ecorr,4), '  delta E  ', round(deltecorr,4))
	elif iter > 9:
		print('step ', iter, '  correlation energy  ', round(ecorr,4), '  delta E  ', round(deltecorr,4))
	elif iter > 99:
		print('step ', iter, ' correlation energy  ', round(ecorr,4), '  delta E  ', round(deltecorr,4))
	titer=result[2]
	if deltecorr < eps:
		break
	if iter == itermax-1 and deltecorr > eps:
		print('!!! no convergeged results after ', iter + 1, ' iteration steps !!!')

print('difference of FCI and CCD  ' ,round(abs(ecorr-crt[0]),10))

plt.tick_params(width=20,direction='in',
	            labelsize='large')
plt.figure(figsize=(6,5),dpi=60)

#plt.plot(coupling,crt,'r-o',lw=2.5)
#plt.plot(coupling,crt,'r-o',lw=2.5)
#plt.plot((-1,1),(0,0),'b--',lw=1.0)
#plt.xlim(-1,1)
#plt.ylim(-5,5)
#plt.xlabel('coupling stength $g$ ',fontsize=12)
#plt.ylabel('Correlation energy',fontsize=12)
#plt.savefig('comp_FCI_CCD_pairing.pdf',format='pdf')
