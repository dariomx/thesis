import sys, slepc4py
slepc4py.init(sys.argv)

from petsc4py import PETSc
from slepc4py import SLEPc
import numpy as np
from test_util import conv_mat

def fiedler_lobpcg(L):
    Lr = conv_mat(L, "csr")
    ptrs = (Lr.indptr, Lr.indices, Lr.data)
    Lp = PETSc.Mat().createAIJ(size=L.shape, csr=ptrs)
    P = SLEPc.EPS(); prob.create()
    P.setOperators(Lp)
    P.setProblemType(SLEPc.EPS.ProblemType.HEP)
    P.setWhichEigenpairs(SLEPc.EPS.Which.SMALLEST_MAGNITUDE)    
    P.setType(SLEPc.EPS.Type.LOBPCG)
    P.setFromOptions()
    P.setDimenstions(2)
    
    
