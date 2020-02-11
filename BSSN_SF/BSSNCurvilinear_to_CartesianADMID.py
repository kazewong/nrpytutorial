# This module converts Cartesian ADM initial data to BSSN
# curvilinear initial data in terms of the variables used
# in BSSN_RHSs.py

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

import sympy as sp
import NRPy_param_funcs as par
from outputC import *
import indexedexp as ixp
import reference_metric as rfm
import BSSN_SF.BSSN_RHSs as bssn # needed for parameters

def Convert_BSSN_curvilinear_Cartesian_ADM(Cartxyz,cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU):

    
    DIM = 3


    Jac_dUCart_dDrfmUD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            Jac_dUCart_dDrfmUD[i][j] = sp.diff(rfm.xxCart[i], rfm.xx[j])

    Jac_dUrfm_dDCartUD, dummyDET = ixp.generic_matrix_inverter3x3(Jac_dUCart_dDrfmUD)  

    gammabarDD = ixp.zerorank2()
    gammaDD = ixp.zerorank2()
    AbarDD = ixp.zerorank2()
    betaU = ixp.zerorank1()
    BU = ixp.zerorank1()
    for i in range(DIM):
        betaU[i] = vetU[i]*rfm.ReU[i]
        BU[i] = betU[i]*rfm.ReU[i]  
        for j in range(DIM):
            gammabarDD[i][j] = hDD[i][j]*rfm.ReDD[i][j] + rfm.ghatDD[i][j]
            AbarDD[i][j] = aDD[i][j] * rfm.ReDD[i][j]
            

    if par.parval_from_str("ConformalFactor") == "W":
        gammaCartDET = (cf ** -6) 
    else:
        print("Error ConformalFactor type = \"" + par.parval_from_str("ConformalFactor") + "\" unknown.")
        exit(1)

    ACartbarDD = ixp.zerorank2()
    gammaCartbarDD = ixp.zerorank2()
    betaCartU = ixp.zerorank1()
    BCartU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            betaCartU[i] += Jac_dUCart_dDrfmUD[i][j] * betaU[j]
            BCartU[i]    += Jac_dUCart_dDrfmUD[i][j] * BU[j]
            for k in range(DIM):
                for l in range(DIM):
                    gammaCartbarDD[i][j] += Jac_dUrfm_dDCartUD[k][i] *Jac_dUrfm_dDCartUD[l][j] * gammabarDD[k][l]
                    ACartbarDD[i][j] += Jac_dUrfm_dDCartUD[k][i] *Jac_dUrfm_dDCartUD[l][j] * AbarDD[k][l]
    
    KCartDD = ixp.zerorank2()
    gammaCartDD = ixp.zerorank2()
    for i in range(DIM):
        betaCartU[i] = sp.simplify(betaCartU[i])
        BCartU[i] = sp.simplify(BCartU[i])
        for j in range(DIM):
            gammaCartDD[i][j] = gammaCartDET ** ( sp.Rational(1,3))*gammaCartbarDD[i][j]
            KCartDD[i][j] = gammaCartDET**(sp.Rational(1,3))*ACartbarDD[i][j] + sp.Rational(1,3)*gammaCartDD[i][j] * trK
            gammaCartDD[i][j] = sp.simplify(gammaCartDD[i][j])
            KCartDD[i][j] = sp.simplify(KCartDD[i][j])

    alphaCart = alpha

    return gammaCartDD,KCartDD,alphaCart,betaCartU,BCartU 


