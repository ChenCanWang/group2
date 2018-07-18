import numpy as np 
#from sympy import *
import matplotlib.pyplot as plt




particle = 4
# Maybe we deal with only the 
# two-pair problem first.


line = 5


delta = 1.0



# set the particles to pairs
pair = particle/2
if(particle%2 != 0):
	print('Odd number of particles.')
	exit(0)


# Get the coefficients of binomial coefficients 
def fact(n):
	if n == 0: 
	    return 1
	else:
		return n*fact(n-1) 
     

def Cmn(n,m):
    if n>=m: 
	    return fact(n)/(fact(n-m)*fact(m))
    else:
     	print('Wrong input.')
     	exit(0)
# Just use Cmn to decide the number of configurations.

config_num  = Cmn(line,pair)

states = [(i,j) for i in range(line) for  
           j in range(i,line) if i!=j]

print('There are %d states: \n'  %config_num)
print(states)

print('\n')

# evaluate the sp operator 
def h_sp((a,b)):
	return a+b



# evaluate the pair term 
def h_pr((a,b),(c,d)):
	h = 0.0
	if(a!=c and a!=d  and
	   b!=c and b!=d  ):
	  return h
	elif(a==c and b ==d ):
	  return h+2
	else:
	  return h+1



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
			Hc[i][j] = Hc[i][j] + '+' \
			       + str(int(2*delta*h_sp(states[i])))   \
			       +'d'


for i in range(len(states)):
	print(Hc[i])

print('\n')

# Now consider the  


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
	crt.append(min(u)-delta+2*g)

for i in range(len(crt)): 		
	 print(crt[i])



plt.tick_params(width=20,direction='in',
	            labelsize='large')
plt.figure(figsize=(6,5),dpi=60)

plt.plot(coupling,crt,'r-o',lw=2.5)
plt.plot((-1,1),(0,0),'g--',lw=1.0)
plt.xlim(-1,1)
plt.ylim(-5,5)
plt.xlabel('coupling stength $g$ ',fontsize=12)
plt.ylabel('Correlation energy',fontsize=12)
plt.savefig('result.pdf',format='pdf')



