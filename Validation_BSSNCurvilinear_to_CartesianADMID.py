import NRPy_param_funcs as par
import reference_metric as rfm
import indexedexp as ixp
import numpy as np 
from BSSN_SF.BSSNCurvilinear_to_CartesianADMID import * 
from BSSN_SF.CartesianADMID_to_BSSNCurvilinearID import * 
from sympy import lambdify

DIM = 3
par.set_parval_from_str("grid::DIM",DIM)

par.set_parval_from_str("reference_metric::CoordSystem","Spherical")
rfm.reference_metric() 
par.set_parval_from_str("indexedexp::symmetry_axes","12")

Cartxyz = ixp.declarerank1("Cartxyz") 
betaCartU = ixp.zerorank1() 

BCartU    = ixp.zerorank1() 
alphaCart = np.random.rand()+1

gammaCartDD = ixp.zerorank2() 
KCartDD     = ixp.zerorank2()

N = 10

for i in range(DIM):
    for j in range(DIM):
        if j >= i : 
            gammaCartDD[i][j] =  sp.Rational(np.random.randint(1,N),np.random.randint(1,N))
            KCartDD[i][j] =  sp.Rational(np.random.randint(1,N),np.random.randint(1,N))
        else :
            gammaCartDD[i][j] =  gammaCartDD[j][i] 
            KCartDD[i][j] = KCartDD[j][i] 


for i in range(DIM):
    gammaCartDD[i][i] += 1
    betaCartU[i] =  sp.Rational(np.random.randint(1,N),np.random.randint(1,N))
    BCartU[i] =  sp.Rational(np.random.randint(1,N),np.random.randint(1,N))

print "Input gamma"
for i in range(DIM):
    print gammaCartDD[:][i]

print "Input K "
for i in range(DIM):
    print KCartDD[:][i]

print "Input beta " 
for i in range(DIM):
    print betaCartU[i]

print "Input B  " 
for i in range(DIM):
    print BCartU[i]

print "Input alpha " 
print alphaCart

cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU = Convert_Cartesian_ADM_to_BSSN_curvilinear(Cartxyz, gammaCartDD,KCartDD,alphaCart,betaCartU,BCartU)

print "--------------------------------------------------"

gammaCartDD_out,KCartDD_out,alphaCart_out,betaCartU_out,BCartU_out = Convert_BSSN_curvilinear_Cartesian_ADM(Cartxyz,cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU)


print "Ouyput gamma"
for i in range(DIM):
    print gammaCartDD_out[:][i]

print "Ouyput K "
for i in range(DIM):
    print KCartDD_out[:][i]

print "Output beta " 
for i in range(DIM):
    print betaCartU_out[i]

print "Output B  " 
for i in range(DIM):
    print BCartU_out[i]

print "Output alpha " 
print alphaCart_out

