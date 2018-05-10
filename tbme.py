fspnum=open('nucleispnumbers.dat','r')
ftbme=open('nucleitwobody.dat','r')
w=0.0 # oscillator freq. (don't know value, so I can't do this yet!)
m=0.0 # mass of nucleon (unknown units, so I can't do this yet!)

import kinematic
import numpy as np

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
                log=1 # boolean algebra var.
                if s1[3]+s2[3]!=0: log*=0 # Check that Mj=0 (fix direction)
                if s1[4]+s2[4]!=0: log*=0 # Check that Mt=0 (fix to deuteron nucleus)
                if spnum.index(s1)>=spnum.index(s2): log*=0 # Fix phase of antisymmetrized states
                if s1[5]+s2[5]>nm: log*=0 # Check N1+N2<=Nmax
                if log==1: self.tbhash.append(tuple([spnum.index(s1),spnum.index(s2)])) # If all above are true, add to 2-body basis
        print "deuteron 2-b basis initialized to space of dimension: ",len(self.tbhash), ' for N2,max= ',nm # Print diagnostic info
        self.dim = len(self.tbhash)
        self.tbme = [[0.0 for x in xrange(self.dim)] for x in xrange(self.dim)] # allocate Hamiltonian
        for l in ftbme:
            b1 = tuple([int(l.split()[0]),int(l.split()[1])]) # bra
            b2 = tuple([int(l.split()[2]),int(l.split()[3])]) # ket
            if b1 in self.tbhash and b2 in self.tbhash:
                tbme = float(l.split()[4])
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

        if x1==y1 and y2[2:5]==x2[2:5]:
            res+=kinematic.tke(x2[0],x2[1],y2[0],y2[1])
        if x2==y2 and y1[2:5]==x1[2:5]:
            res+=kinematic.tke(x1[0],x1[1],y1[0],y1[1])

        return res

    def solve(self):
        eigval,eigvec = np.linalg.eigh(np.array(self.tbme))
        for e in eigval:
            print e


b=deuteron(3)
b.solve()
