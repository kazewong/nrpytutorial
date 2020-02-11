# This module sets up an initial data function meant to
# be called in a pointwise manner at all gridpoints.

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

from outputC import *

def ADM_ID_function_string(gammaCartDD,KCartDD,alphaCart,betaCartU,BCartU):
    returnstring = "void BSSN_ID(FILE *out2D, REAL time, REAL xx0,REAL xx1,REAL xx2,REAL Cartxyz0,REAL Cartxyz1,REAL Cartxyz2,\n"
    returnstring += "\tREAL *hDD00,REAL *hDD01,REAL *hDD02,REAL *hDD11,REAL *hDD12,REAL *hDD22,\n"
    returnstring += "\tREAL *aDD00,REAL *aDD01,REAL *aDD02,REAL *aDD11,REAL *aDD12,REAL *aDD22,\n"
    returnstring += "\tREAL *trK,\n"
    returnstring += "\tREAL *lambdaU0,REAL *lambdaU1,REAL *lambdaU2,\n"
    returnstring += "\tREAL *vetU0,REAL *vetU1,REAL *vetU2,\n"
    returnstring += "\tREAL *betU0,REAL *betU1,REAL *betU2,\n"
    returnstring += "\tREAL *alpha,REAL *cf) {\n"
    returnstring += "\tdouble gammaCartDD00;\n "
    returnstring += "\tdouble gammaCartDD01;\n "
    returnstring += "\tdouble gammaCartDD02;\n "
    returnstring += "\tdouble gammaCartDD11;\n "
    returnstring += "\tdouble gammaCartDD12;\n "
    returnstring += "\tdouble gammaCartDD22;\n "
    returnstring += "\tdouble KCartDD00;\n "
    returnstring += "\tdouble KCartDD01;\n "
    returnstring += "\tdouble KCartDD02;\n "
    returnstring += "\tdouble KCartDD11;\n "
    returnstring += "\tdouble KCartDD12;\n "
    returnstring += "\tdouble KCartDD22;\n "
    returnstring += "\tdouble betaCartU0;\n "
    returnstring += "\tdouble betaCartU1;\n "
    returnstring += "\tdouble betaCartU2;\n "
    returnstring += "\tdouble BCartU0;\n "
    returnstring += "\tdouble BCartU1;\n "
    returnstring += "\tdouble BCartU2;\n "
    returnstring += "\tdouble alpha;\n "
    returnstring += outputC([gammaCartDD[0][0], gammaCartDD[0][1], gammaCartDD[0][2], gammaCartDD[1][1], gammaCartDD[1][2], gammaCartDD[2][2],
                             KCartDD[0][0], KCartDD[0][1], KCartDD[0][2], KCartDD[1][1], KCartDD[1][2], KCartDD[2][2],
                             betaCartU[0], betaCartU[1], betaCartU[2],
                             BCartU[0], BCartU[1], BCartU[2],
                             alphaCart],
                            ["gammaCartDD00", "gammaCartDD01", "gammaCartDD02", "gammaCartDD11", "gammaCartDD12", "gammaCartDD22",
                             "KCartDD00", "KCartDD01", "KCartDD02", "KCartDD11", "KCartDD12", "KCartDD22",
                             "betaCartU0", "betaCartU1", "betaCartU2",
                             "BCartU0", "BCartU1", "BCartU2",
                             "alpha"], filename="returnstring",
                            params="preindent=1,CSE_enable=True,outCverbose=False",  # outCverbose=False to prevent
                            # enormous output files.
                            prestring="", poststring="")
    returnstring += 'fprintf(out2D,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e  \\n", \n '
    returnstring += '       time, Cartxyz0, gammaCartDD00, gammaCartDD01, gammaCartDD02, gammaCartDD11, gammaCartDD12, gammaCartDD22,  \n '
    returnstring += '       KCartDD00, KCartDD01, KCartDD02, KCartDD11, KCartDD12, KCartDD22,  \n '
    returnstring += '       betaCartU0, betaCartU1, betaCartU2, \n '
    returnstring += '       BCartU0, BCartU1, BCartU2, \n '
    returnstring += '       alpha) \n '
    returnstring += "}\n"
    return returnstring
