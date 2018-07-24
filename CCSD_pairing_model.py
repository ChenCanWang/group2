import numpy as np
#from sympy import *
import math, sys, os, glob
import matplotlib.pyplot as plt



#give input values
particle = 6
line = 6
gval=0.0
delta = 1.0



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

config_num  = len(states)

print('There are %d states: \n'  %config_num)
print(states)

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
print('The form of the matrix is: \n')
for i in range(config_num):
	for j in range(config_num):
		Hc[i][j]='-'+str(int(h_pr(states[i],states[j])))+'g'
		if(i==j):
			Hc[i][j] = Hc[i][j]
#                   + '+' \
#			       + str(int(2*delta*h_sp(states[i])))   \
#			       +'d'


#for i in range(len(states)):
#	print(Hc[i])

print('\n')

coupling = np.linspace(-1,1,11)

egv = [] # as the eigenvalue of H
crt = [] # as the correlation energy



# find the minimum of the eigenvalues
for g in coupling:
# construction of the Hamiltonian
	for i in range(config_num):
		for j in range(config_num):
			H[i,j]=-g*h_pr(states[i],states[j])
			if(i==j):
				H[i,j] = H[i,j] + delta*h_sp(states[i])


    # get the eigenvalues of H
	u, v = np.linalg.eig(H)
	egv.append(min(u))

	if(pair==1):
	    crt.append(min(u) + 2*g)
	else:
		crt.append(min(u)- sum(states[0])* delta + 2*g)

#set fermi level
def flevel(number):
    return number-1

#set valence space defined by number of lines
def exspace(levels):
    return (2*levels)-1

fermi=flevel(particle)
vspace=exspace(line)

#single particle part contributes only for diagonal
def fpq(p,q,gval):
	if p != q:
		return 0
	elif p in range(particle):
		return math.ceil((p-fermi)/2.0)*delta - gval
	else:
		return math.ceil((p-fermi)/2.0)*delta

#define two-body interaction and check of a,b and i,j are paired
#differ V_pppp V_pphh (transpose to V_hhpp) and V_hhhh
#start counting a and b from Fermi level, so substract number of particles in array
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


print(f_pp)
print(f_hh)

#define T2 array
tfactor= np.zeros(shape=(vspace-fermi,vspace-fermi,fermi+1,fermi+1))


#set up Hamiltonian
#def hbar(t,i,j,a,b,gval):
#    return tbint(a,b,i,j,gval) +  t[i,j,a,b]*(fpq(b,g)+fpq(a,g)-fpq(i,g)-fpq(j,g))

#print(hbar(1,0,0,1,2))
