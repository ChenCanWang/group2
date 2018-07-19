import numpy as np


line = 8
particle = 4

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

minbit  = 0
for i in range(0,particle):
	minbit += 2**i 

dim = Cmn(line,particle)
states=[]
i = minbit
while 1:
	if bin(i).count('1')==particle:
		states.append(i) 		
	if len(states) == Cmn(line,particle): 
		break  
	i=i+1

print(states)
def trans(n):
  s = []
  s = list(bin(n))
  s.reverse()
  del(s[-1])
  del(s[-1])
  while len(s)<line:
   	 s.append('0')  
  return s

k = 0 
pair = []
fullpair = []

for i in range(dim):
	for j in range(line/2):
		k = 2*j
		if(trans(states[i])[k]=='1' and trans(states[i])[k+1]=='1'):
			pair.append(i)

for i in range(dim):
 	if pair.count(i) == particle/2:
 		fullpair.append(i)

print(fullpair)
H = np.eye(dim)
g=1.0

for i in range(dim):
	for j in range(dim):
		if bin(states[i]&states[j]).count('1') - particle == 0:
			if i in fullpair:
				H[i,j] = -g/2* particle/2
		if 	bin(states[i]&states[j]).count('1') - particle + 2 == 0:
			H[i,j]=0
		if 	bin(states[i]&states[j]).count('1') - particle + 4 == 0:
			if i in pair and j in pair:
				H[i,j] = -g/2
#uncomplete 
print(H)			





	

