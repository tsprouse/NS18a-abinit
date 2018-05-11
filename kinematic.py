from util import phase, hat


def rme_bdag(n1, l1, n2, l2):
    if abs(l2-l1) != 1:
        return 0.
    N1 = 2*n1+l1
    N2 = 2*n2+l2
    if N1==N2+1 and l1==l2+1:
        res =  ((l2+1.)*(N2+l2+3.))**0.5
    elif N1==N2+1 and l1==l2-1:
        res =  ((l2*(N2-l2+2.)))**0.5
    else:
        res = 0.0
    return res

def rme_b(n1,l1,n2,l2):
    if abs(l2-l1) != 1:
        return 0.
    N1 = 2*n1+l1
    N2 = 2*n2+l2
    if N1==N2-1 and l1==l2+1:
        res =  ((l2+1.)*(N2-l2))**0.5
    elif N1==N2-1 and l1==l2-1:
        res =  ((l2*(N2+l2+1.)))**0.5
    else:
        res = 0.0
    return res

def tkrme(n1,l1,n2,l2):
    res = -2.5 / hat(l1)
    bdbd = 0.
    bbd  = 0.
    bb   = 0.
    bdb  = 0.

    # QMN = n l j mj tz
    # Only loop over n and l

    if l1 != l2:
        return 0.

    for dn in [-2, -1, 0, 1, 2]:
        n = n2 + dn
        if n < 0:
            continue
        for dl in [-1, 1]:
            l = l2 + dl
            if l < 0:
                continue
            #bdbd += phase(l-l1)*rme_bdag(n1,l1,n,l)*rme_bdag(n,l,n2,l2)
            bdbd += rme_bdag(n1,l1,n,l)*rme_bdag(n,l,n2,l2)
            bbd  += rme_b(n1,l1,n,l)*rme_bdag(n,l,n2,l2)
            bb   += rme_b(n1,l1,n,l)*rme_b(n,l,n2,l2)
            bdb  += rme_bdag(n1,l1,n,l)*rme_b(n,l,n2,l2)

    res = res*(bdbd + bb - bbd - bdb)

    return res
