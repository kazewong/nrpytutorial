from outputC import *
import NRPy_param_funcs as par
import grid as gri
import loop as lp
import indexedexp as ixp
import finite_difference as fin
import reference_metric as rfm

DIM = 3
par.set_parval_from_str("grid::DIM",DIM)
par.set_parval_from_str("reference_metric::CoordSystem","Spherical")
rfm.reference_metric() # Create ReU, ReDD needed for rescaling B-L initial data, generating BSSN RHSs, etc.
par.set_parval_from_str("indexedexp::symmetry_axes","12")


from BSSN_SF.BSSNCurvilinear_to_CartesianADMID import * 
import BSSN_SF.BSSN_RHSs as bssnrhs 
import BSSN_SF.ADM_ID_function_string as ADIF

bssnrhs.BSSN_RHSs()

Cartxyz = ixp.declarerank1("Cartxyz") 

cf       = bssnrhs.cf
hDD      = bssnrhs.hDD
lambdaU  = None
aDD      = bssnrhs.aDD
trK      = bssnrhs.trK
alpha    = bssnrhs.alpha
vetU     = bssnrhs.vetU
betU     = bssnrhs.betU

gammaCartDD,KCartDD,alphaCart,betaCartU,BCartU = Convert_BSSN_curvilinear_Cartesian_ADM(Cartxyz,cf,hDD,lambdaU,aDD,trK,alpha,vetU,betU)
returnfunction = ADIF.ADM_ID_function_string(gammaCartDD,KCartDD,alphaCart,betaCartU,BCartU)

with open("BSSN_SF/ID_array_ADM.h","w") as file:
        file.write(returnfunction)
