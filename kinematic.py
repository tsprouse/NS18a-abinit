from util import phase, hat


def rme_bdag(n1, l1, n2, l2):
    N1 = 2*n1+l1
    N2 = 2*n2+l2
    if N1-1==N2 and l1-1==l2:
        res = hat(l2+1)*((l2+1.)*(N2+l2+3.))**0.5
    elif N1-1==N2 and l1+1==l2:
        res = hat(l2-1)*((l2*(N2-l2+2.)))**0.5
    else:
        res = 0.0
    return res

def rme_b(n1,l1,n2,l2):
    N1 = 2*n1+l1
    N2 = 2*n2+l2
    if N1+1==N2 and l1-1==l2:
        res = hat(l2+1)*((l2+1.)*(N2-l2))**0.5
    elif N1+1==N2 and l1+1==l2:
        res = hat(l2-1)*((l2*(N2+l2+1.)))**0.5
    else:
        res = 0.0
    return res

def tke(n1,l1,n2,l2):
    res = 2.5 / hat(l1)**2
    bdbd = 0.
    bbd  = 0.
    bb   = 0.
    bdb  = 0.

    # QMN = n l j mj tz
    # Only loop over n and l

    if l1 != l2:
        return 0.

    for dn in [-1, 1]:
        n = n2 + dn
        for dl in [-1, 1]:
            l = l2 + dl
            bdbd += phase(l-l1)*rme_bdag(n1,l1,n,l)*rme_bdag(n,l,n2,l2)
            bbd  += phase(l-l1)*rme_b(n1,l1,n,l)*rme_bdag(n,l,n2,l2)
            bb   += phase(l-l1)*rme_b(n1,l1,n,l)*rme_b(n,l,n2,l2)
            bdb  += phase(l-l1)*rme_bdag(n1,l1,n,l)*rme_b(n,l,n2,l2)

    res = res*(bdbd + bb - bbd - bdb)

    return res
