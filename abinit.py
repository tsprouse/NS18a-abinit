fspnum='nucleispnumbers.dat'
ftbme='nucleitwobody.dat'

import kinematic
import numpy as np


class sp_state:
    def __init__(self):
        self.n = 0
        self.N = 0
        self.j = 0
        self.m = 0
        self.l = 0
        self.tz = 0

    def read(self,l):
        self.n = int(l.split()[1])
        self.l = int(l.split()[2])
        self.j = float(l.split()[3])/2.0
        self.m = float(l.split()[4])/2.0
        self.tz = float(l.split()[5])/2.0
        self.N = 2*self.n+self.l



class spnum:

    def __init__(self):
        self.spec = list()
        f = open(fspnum)
        for l in f:
            b=sp_state()
            b.read(l)
            self.spec.append(b)

        f.close()


class nb_state:
    def __init__(self):
        self.null = False
        self.phase  = 1
        self.spec = list()
    def a(self,ind):
        if len(self.spec) == 0 or self.null:
            self.null = True
        elif ind in self.spec:
            i = self.spec.index(ind)
            self.phase = self.phase*((-1)**(i))
            self.spec.remove(ind)
        else:
            self.null = True
    def adag(self,ind):
        if len(self.spec)== 0 and not self.null:
            self.spec.append(ind)
        elif ind in self.spec or self.null:
            self.null = True
        else:
            for i in xrange(len(self.spec)):
                if ind < self.spec[i]:
                    self.spec.insert(i,ind)
                    self.phase = self.phase*((-1)**i)
                    return
            self.spec.append(ind)
            self.phase = self.phase*((-1)**(len(self.spec)))
    def fix(self):
        self.null = False
        self.phase = 1
        self.save = tuple(self.spec)
        self.rep = set(self.save)
    def reset(self):
        self.null = False
        self.phase = 1
        self.spec = list(self.save)

def ip(bra,ket):
    if bra.null or ket.null:
        return 0.0
    elif bra.spec == ket.spec:
        return 1.*bra.phase*ket.phase
    else:
        return 0.0

class nb_basis:
    def __init__(self, sp):
        self.spec = list()
        self.sp = sp
    def build(self,z,a,nmax,m):
        self.z = z
        self.a = a
        self.nmax = nmax
        self.m = m
        g = gen(len(self.sp.spec),a)
        pf = True
        while pf:
            if self.check(g.vec):
                nb = nb_state()
                for x in g.vec: nb.adag(x)
                nb.fix()
                self.spec.append(nb)
            pf = g.nxt()

    def check(self,tup):
        m = 0
        z = 0
        n = 0
        for x in tup:
            m+=self.sp.spec[x].m
            n+=self.sp.spec[x].N
            if self.sp.spec[x].tz>0:
                z+=1
        if m==self.m and n == self.nmax and z==self.z:
            return True
        else:
            return False

class gen:
    def __init__(self,splen,a):
        self.splen = splen
        self.a = a
        self.vec = [0 for x in xrange(a)]
    def nxt(self):
        for x in xrange(self.a):
            if self.vec[x] < self.splen-1:
                self.vec[x] = self.vec[x]+1
                if x == self.a-1: print self.vec[x]
                return True
            else:
                self.vec[x] = 0
        return False


class tbint:
    def __init__(self):
        self.spec=dict()
        f = open(ftbme)
        for l in f:
            bra1 = int(l.split()[0])-1
            bra2 = int(l.split()[1])-1
            ket1 = int(l.split()[2])-1
            ket2 = int(l.split()[3])-1
            val  = float(l.split()[4])
            self.spec[(bra1,bra2,ket1,ket2)] = val

class hamiltonian:

    def __init__(self,bas,v):
        self.basis = bas
        self.v = v

        self.spec = [[0.0 for x in xrange(len(self.basis.spec))] for x in xrange(len(self.basis.spec))]

        for xx in xrange(len(self.basis.spec)):
            for yy in xrange(len(self.basis.spec)):
                x = self.basis.spec[xx]
                y = self.basis.spec[yy]
                diff = x.rep.difference(y.rep)
                if len(diff)<=2:
                    for a in x.rep:
                        for b in x.rep:
                            for g in y.rep:
                                for d in y.rep:
                                    if (a,b,g,d) in self.v.spec:
                                        y.a(g)
                                        y.a(d)
                                        y.adag(b)
                                        y.adag(a)
                                        self.spec[xx][yy]+=ip(y,x)*0.25*self.v.spec[(a,b,g,d)]
                                        y.reset()
        for xx in xrange(len(self.basis.spec)):
            for yy in xrange(len(self.basis.spec)):
                x = self.basis.spec[xx]
                y = self.basis.spec[yy]
                diff = x.rep.difference(y.rep)
                if len(diff)<=1:
                    for a in x.rep:
                        for b in y.rep:
                            y.a(b)
                            y.adag(a)
                            self.spec[xx][yy]+=ip(y,x)*KinMatEl(self.basis.sp[a],self.basis.sp[b])
                            y.reset()

def KinMatEl(bra,ket):
    if bra.j==ket.j and bra.m==ket.m and bra.tz==ket.tz:
        return kinematic.tke()


this = spnum()
that = nb_basis(this)
that.build(2,4,3,0)
these = tbint()
those = hamiltonian(that,these)

















#end
