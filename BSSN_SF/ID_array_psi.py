# This module sets up arbitrary initial data to be input later froma file in terms of
# the variables used in BSSN_RHSs.py
# 
# **Input variables needed for spacetime evolution**:
# * Conformal factor, psi_in
# * Desired coordinate system
# * Desired initial lapse $\alpha_in$ and shift $\beta^i_in$
# 
# **Transformation to curvilinear coordinates**:
# * Once the above variables have been set in Cartesian coordinates, we will apply the appropriate coordinate transformations and tensor rescalings ([described in the BSSN NRPy+ tutorial module](Tutorial-BSSNCurvilinear.ipynb))

# Step P0: Load needed modules
import numpy as np
import sympy as sp
import NRPy_param_funcs as par
from outputC import *
import indexedexp as ixp
import BSSN_SF.CartesianADMID_to_BSSNCurvilinearID as ctob
import BSSN_SF.BSSN_ID_function_string as bIDf

thismodule = "ID_array_psi"
psi_in = par.Cparameters("REAL", thismodule, ["psi_in"], "")
alpha_in = par.Cparameters("REAL", thismodule, ["alpha_in"], "")

def ID_array_psi():
    global Cartxyz, gammaCartDD, KCartDD, alphaCart, betaCartU, BCartU
    
    # Step 0: Set spatial dimension (must be 3 for BSSN)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)
    Cartxyz = ixp.declarerank1("Cartxyz")

    # Step 1: Set psi, the conformal factor:
    psi = psi_in 

    # Step 2: Set all needed ADM variables in Cartesian coordinates
    gammaCartDD = ixp.zerorank2()
    KCartDD     = ixp.zerorank2() # K_{ij} = 0 for these initial data
    for i in range(DIM):
        gammaCartDD[i][i] = psi**4

    alphaCart = alpha_in
    betaCartU = ixp.zerorank1() # We generally choose \beta^i = 0 for these initial data
    BCartU    = ixp.zerorank1() # We generally choose B^i = 0 for these initial data
    
    cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU = \
        ctob.Convert_Cartesian_ADM_to_BSSN_curvilinear(Cartxyz, gammaCartDD,KCartDD,alphaCart,betaCartU,BCartU)
    
    global returnfunction
    returnfunction = bIDf.BSSN_ID_function_string(cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU)