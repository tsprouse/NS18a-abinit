fspnum=open('nucleispnumbers.dat','r')
ftbme=open('nucleitwobody.dat','r')
w=0.0 # oscillator freq. (don't know value, so I can't do this yet!)
m=0.0 # mass of nucleon (unknown units, so I can't do this yet!)

import kinematic
import numpy as np
from util import hat

spnum = []
invspnum = dict()
i=0

for l in fspnum:
    tup= (int(l.split()[1]),int(l.split()[2]),int(l.split()[3]),int(l.split()[4]),int(l.split()[5]), 2*int(l.split()[1])+int(l.split()[2]))
    spnum.append(tup)
    invspnum[tup]=i
    i=i+1




# for key in invspnum:
    # print invspnum[key],key

# define a basis for a given N=N1+N2, for Mj1+Mj2=0
class deuteron:
    tbme=None
    tbhash=list()
    dim=0

    def __init__(self, nm):
        for s1 in spnum:
            for s2 in spnum:
                log=True # boolean algebra var.
                log = log and (s1[3]+s2[3]==0)  # Check that Mj=0 (fix direction)
                log = log and (s1[4]+s2[4]==0)  # Check that Mt=0 (fix to deuteron nucleus)
                log = log and (spnum.index(s1)<=spnum.index(s2)) # Fix phase of antisymmetrized states
                log = log and (s1[5]+s2[5]<=nm) # Check N1+N2<=Nmax
                log = log and (s1[1]+s2[1])%2 == nm%2 # check that g=Nmax%2
                if log: self.tbhash.append(tuple([spnum.index(s1),spnum.index(s2)])) # If all above are true, add to 2-body basis
        print "deuteron 2-b basis initialized to space of dimension: ",len(self.tbhash), ' for N2,max= ',nm # Print diagnostic info
        self.dim = len(self.tbhash)
        self.tbme = [[0.0 for x in xrange(self.dim)] for y in xrange(self.dim)] # allocate Hamiltonian
        for l in ftbme:
            a,b,c,d,en = l.split()[0:5]
            b1 = (int(a)-1, int(b)-1) # bra
            b2 = (int(c)-1, int(d)-1) # ket
            if b1 in self.tbhash and b2 in self.tbhash:
                tbme = float(en)
                self.tbme[self.tbhash.index(b1)][self.tbhash.index(b2)] = tbme
        ftbme.close()

        for x in xrange(self.dim):
            for y in xrange(self.dim):
                self.tbme[x][y]+=self.kin2bme(x,y)

    def kin2bme(self,x,y):
        x1=spnum[self.tbhash[x][0]]
        x2=spnum[self.tbhash[x][1]]
        y1=spnum[self.tbhash[y][0]]
        y2=spnum[self.tbhash[y][1]]

        res = 0.0

        if x1==y1 and y2[1:5]==x2[1:5]:
            res+=1/hat(x2[1])*kinematic.tkrme(x2[0],x2[1],y2[0],y2[1])
        if x2==y2 and y1[1:5]==x1[1:5]:
            res+=1/hat(x1[1])*kinematic.tkrme(x1[0],x1[1],y1[0],y1[1])

        return res

    def solve(self):
        np.savetxt('matrix.txt', self.tbme)
        print(self.tbhash)
        eigval,eigvec = np.linalg.eigh(np.array(self.tbme))
        for e in eigval:
            print e


b=deuteron(2)
b.solve()
