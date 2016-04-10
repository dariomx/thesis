from datetime import datetime
from numpy import loadtxt, sign
from scipy.linalg import norm
from scipy.io import mmread, mminfo

def load(L_file, ac_file=None, fv_file=None):
    if L_file.endswith(".mat"):
        L  = loadtxt(L_file)
    elif L_file.endswith(".mtx") or L_file.endswith(".mtz.gz"):
        print "loading mm format:  " + str(mminfo(L_file))
        start = datetime.now()
        L  = mmread(L_file).toarray()
        end = datetime.now()
        print "loaded in " + str(end-start)
    ac = None if ac_file is None else loadtxt(ac_file)[1]
    fv = None if fv_file is None else loadtxt(fv_file)
    return L, ac, fv

# invert y if the signs are opposite as x
def invsign(y, x):
    return -y if (sign(x) == -sign(y)).all() else y

def relres(L, ac, fv):
    return norm(L*fv - ac*fv) / norm(L)

def relerr(x, y):
    return norm(x - invsign(y, x)) / norm(x)

def cmp(L, ac, fv, cac, cfv):
    # for i in xrange(len(fv)):
    #     args = (i, fv[i], cfv[i], abs(fv[i] - cfv[i]))
    #     print "fv[%03d]: %+.16f - %+.16f = %.16f" % args
    print "relres (calc): %.16f" % relres(L, cac, cfv)
    if ac is not None and fv is not None:
        print "relres (inpt): %.16f" % relres(L, ac, fv)
        args = (relerr(ac, cac), relerr(fv, cfv))
        print "relerr:  ac = %.16f,  fv = %.16f" % args
    
def test(argv, calc):
    if len(argv) == 4:
        L, ac, fv = load(argv[1], argv[2], argv[3])
    else:
        L, ac, fv = load(argv[1])  
    start = datetime.now()
    cac, cfv = calc(L)
    end = datetime.now()
    print "calc took %s" % (end - start)
    cmp(L, ac, fv, cac, cfv)
