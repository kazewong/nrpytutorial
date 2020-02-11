# This module will declare and set many of the quantities that are useful for numerical relativity.
# This includes the Ricci tensor, metric tensor, and Christoffel Symbols

# Author: Zachariah B. Etienne
#         zachetie **at** gmail **dot* com

# Step P1: import all needed modules from NRPy+:
import NRPy_param_funcs as par
import indexedexp as ixp
import grid as gri
import finite_difference as fin
import reference_metric as rfm
from outputC import *
import sympy as syp

# Calculates the useful tensors and tensor-like quantities in the BSSN metric.
# Step P2: Initialize BSSN_RHS parameters
thismodule = __name__
par.initialize_param(par.glb_param("char", thismodule, "ConformalFactor", "W"))
par.initialize_param(par.glb_param("char", thismodule, "ShiftEvolutionOption", "GammaDriving2ndOrder_NoCovariant"))
par.initialize_param(par.glb_param("bool", thismodule, "detgbarOverdetghat_equals_one", "True"))

PI = par.Cparameters("#define", thismodule, "M_PI","")
br_on = par.Cparameters("REAL", thismodule, "br_on",1.) # Switch on/off back-reaction
pot1_on = par.Cparameters("REAL", thismodule, "pot1_on",1.) # Switch on/off quadratic potential
pot2_on = par.Cparameters("REAL", thismodule, "pot2_on",1.) # Switch on/off self-interacting potential
wavespeed = par.Cparameters("REAL", thismodule, "wavespeed",1.) # scalar wave speed
scalarmass = par.Cparameters("REAL", thismodule, "scalarmass",1.) # scalar mass
fa = par.Cparameters("REAL", thismodule, "fa",0.5) # axion decay constant

def BSSN_RHSs():
    # Step 1: Set up reference metric
    rfm.reference_metric()

    # Step 2: Set spatial dimension (must be 3 for BSSN)
    DIM = 3
    par.set_parval_from_str("grid::DIM",DIM)

    # Step 3: Register all needed *evolved* gridfunctions.
    # Step 3a: Register indexed quantities, using ixp.register_... functions
    global hDD # Needed for setting detgammabar = detgammahat constraint, etc.
    hDD = ixp.register_gridfunctions_for_single_rank2("EVOL","hDD", "sym01")
    global aDD,vetU,betU # Needed for BSSN constraints, etc.
    aDD = ixp.register_gridfunctions_for_single_rank2("EVOL","aDD", "sym01")
    lambdaU = ixp.register_gridfunctions_for_single_rank1("EVOL","lambdaU")
    vetU = ixp.register_gridfunctions_for_single_rank1("EVOL","vetU")
    betU = ixp.register_gridfunctions_for_single_rank1("EVOL","betU")
    # Step 3b: Register scalar quantities, using gri.register_gridfunctions()
    global trK,cf,alpha # Needed as global for BSSN matter source terms, etc.
    trK, cf, alpha = gri.register_gridfunctions("EVOL",["trK","cf","alpha"])
    Rtrace = gri.register_gridfunctions("EVOL",["Rtrace"])

    global uu, vv # Needed as global for BSSN matter source terms, etc.
    uu, vv = gri.register_gridfunctions("EVOL",["uu","vv"])
  
    # Step 4: Register all *auxiliary* gridfunctions.
    # Step 4a: Register indexed quantities, using ixp.register_... functions
    #RbarDD = ixp.register_gridfunctions_for_single_rank2("EVOL","RbarDD", "sym01")
    # Step 4b: Register scalar quantities, using gri.register_gridfunctions()
    if par.parval_from_str("BSSN_SF.BSSN_RHSs::detgbarOverdetghat_equals_one") == "False":
        detgbarOverdetghat = gri.register_gridfunctions("AUX", ["detgbarOverdetghat"])
        detgbarOverdetghatInitial = gri.register_gridfunctions("AUX", ["detgbarOverdetghatInitial"])
        print("Error: detgbarOverdetghat_equals_one=\"False\" is not fully implemented yet.")
        exit(1)
    else:
        detgbarOverdetghat = sp.sympify(1)

    # Step 5a: Define \varepsilon_{ij} = epsDD[i][j]
    epsDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            epsDD[i][j] = hDD[i][j]*rfm.ReDD[i][j]

    # Step 5b: Define epsDD_dD[i][j][k]
    hDD_dD = ixp.declarerank3("hDD_dD","sym01")
    epsDD_dD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                epsDD_dD[i][j][k] = hDD_dD[i][j][k]*rfm.ReDD[i][j] + hDD[i][j]*rfm.ReDDdD[i][j][k]

    # Step 5c: Define epsDD_dDD[i][j][k][l]
    hDD_dDD = ixp.declarerank4("hDD_dDD","sym01_sym23")
    epsDD_dDD = ixp.zerorank4()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    epsDD_dDD[i][j][k][l] = hDD_dDD[i][j][k][l]*rfm.ReDD[i][j] + \
                                            hDD_dD[i][j][k]*rfm.ReDDdD[i][j][l] + \
                                            hDD_dD[i][j][l]*rfm.ReDDdD[i][j][k] + \
                                            hDD[i][j]*rfm.ReDDdDD[i][j][k][l]

    # Step 5d: Define DhatgammabarDDdD[i][j][l] = \bar{\gamma}_{ij;\hat{l}}
    # \bar{\gamma}_{ij;\hat{l}} = \varepsilon_{i j,l}
    #                           - \hat{\Gamma}^m_{i l} \varepsilon_{m j}
    #                           - \hat{\Gamma}^m_{j l} \varepsilon_{i m}
    global gammabarDD_dHatD # Needed for BSSN constraints, etc.
    gammabarDD_dHatD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for l in range(DIM):
                gammabarDD_dHatD[i][j][l] = epsDD_dD[i][j][l]
                for m in range(DIM):
                    gammabarDD_dHatD[i][j][l] += - rfm.GammahatUDD[m][i][l]*epsDD[m][j] \
                                                 - rfm.GammahatUDD[m][j][l]*epsDD[i][m]

    # Step 5e: Define \bar{\gamma}_{ij;\hat{l},k} = DhatgammabarDD_dHatD_dD[i][j][l][k]:
    #        \bar{\gamma}_{ij;\hat{l},k} = \varepsilon_{ij,lk}
    #                                      - \hat{\Gamma}^m_{i l,k} \varepsilon_{m j}
    #                                      - \hat{\Gamma}^m_{i l} \varepsilon_{m j,k}
    #                                      - \hat{\Gamma}^m_{j l,k} \varepsilon_{i m}
    #                                      - \hat{\Gamma}^m_{j l} \varepsilon_{i m,k}
    gammabarDD_dHatD_dD = ixp.zerorank4()
    for i in range(DIM):
        for j in range(DIM):
            for l in range(DIM):
                for k in range(DIM):
                    gammabarDD_dHatD_dD[i][j][l][k] = epsDD_dDD[i][j][l][k]
                    for m in range(DIM):
                        gammabarDD_dHatD_dD[i][j][l][k] += -rfm.GammahatUDDdD[m][i][l][k]*epsDD[m][j]  \
                                                           -rfm.GammahatUDD[m][i][l]*epsDD_dD[m][j][k] \
                                                           -rfm.GammahatUDDdD[m][j][l][k]*epsDD[i][m]  \
                                                           -rfm.GammahatUDD[m][j][l]*epsDD_dD[i][m][k]

    # Step 5f: Define \bar{\gamma}_{ij;\hat{l}\hat{k}} = DhatgammabarDD_dHatDD[i][j][l][k]
    #          \bar{\gamma}_{ij;\hat{l}\hat{k}} = \partial_k \hat{D}_{l} \varepsilon_{i j}
    #                                           - \hat{\Gamma}^m_{lk} \left(\hat{D}_{m} \varepsilon_{i j}\right)
    #                                           - \hat{\Gamma}^m_{ik} \left(\hat{D}_{l} \varepsilon_{m j}\right)
    #                                           - \hat{\Gamma}^m_{jk} \left(\hat{D}_{l} \varepsilon_{i m}\right)
    gammabarDD_dHatDD = ixp.zerorank4()
    for i in range(DIM):
        for j in range(DIM):
            for l in range(DIM):
                for k in range(DIM):
                    gammabarDD_dHatDD[i][j][l][k] = gammabarDD_dHatD_dD[i][j][l][k]
                    for m in range(DIM):
                        gammabarDD_dHatDD[i][j][l][k] += - rfm.GammahatUDD[m][l][k]*gammabarDD_dHatD[i][j][m] \
                                                         - rfm.GammahatUDD[m][i][k]*gammabarDD_dHatD[m][j][l] \
                                                         - rfm.GammahatUDD[m][j][k]*gammabarDD_dHatD[i][m][l]

    # Step 5g: Compute \bar{\gamma}_{ij} and its inverse (using built-in function ixp.symm_matrix_inverter3x3()):
    global gammabarDD # Needed as global for BSSN matter source terms, etc.
    gammabarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            gammabarDD[i][j] = hDD[i][j]*rfm.ReDD[i][j] + rfm.ghatDD[i][j]
    global gammabarUU # Needed as global for BSSN matter source terms, etc.
    gammabarUU, dummydet = ixp.symm_matrix_inverter3x3(gammabarDD)

    # Step 5h: Add the first term to RbarDD:
    #         - \frac{1}{2} \bar{\gamma}^{k l} \hat{D}_{k} \hat{D}_{l} \bar{\gamma}_{i j}
    global RbarDD # Needed as global for BSSN constraints, etc.
    RbarDD = ixp.zerorank2()
    RbarDDpiece = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    RbarDD[i][j] += -sp.Rational(1,2) * gammabarUU[k][l]*gammabarDD_dHatDD[i][j][l][k]
                    RbarDDpiece[i][j] += -sp.Rational(1,2) * gammabarUU[k][l]*gammabarDD_dHatDD[i][j][l][k]

    # Step 6a: Second term of RhatDD: compute \hat{D}_{j} \bar{\Lambda}^{k} = LambarU_dHatD[k][j]
    lambdaU_dD = ixp.declarerank2("lambdaU_dD","nosym")
    LambarU_dHatD = ixp.zerorank2()
    for j in range(DIM):
        for k in range(DIM):
            LambarU_dHatD[k][j] = lambdaU_dD[k][j]*rfm.ReU[k] + lambdaU[k]*rfm.ReUdD[k][j]
            for m in range(DIM):
                LambarU_dHatD[k][j] += rfm.GammahatUDD[k][m][j]*lambdaU[m]*rfm.ReU[m]

    # Step 6b: Add the second term to the Ricci tensor
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                RbarDD[i][j] += sp.Rational(1,2) * (gammabarDD[k][i]*LambarU_dHatD[k][j] + \
                                                    gammabarDD[k][j]*LambarU_dHatD[k][i])

    # Step 7a: Define \bar{\gamma}_{ij,k} = gammabarDDdD[i][j][k]
    #          = h_{ij,k} \text{ReDD[i][j]} + h_{ij} \text{ReDDdD[i][j][k]} + \hat{\gamma}_{ij,k}.
    gammabarDD_dD = ixp.zerorank3()
    hDD_dupD = ixp.declarerank3("hDD_dupD","sym01") # Needed for \bar{\gamma}_{ij} RHS
    gammabarDD_dupD = ixp.zerorank3()  # Needed for \bar{\gamma}_{ij} RHS
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gammabarDD_dD[i][j][k] = hDD_dD[i][j][k]*rfm.ReDD[i][j] + hDD[i][j]*rfm.ReDDdD[i][j][k] \
                                       + rfm.ghatDDdD[i][j][k]
                gammabarDD_dupD[i][j][k] = hDD_dupD[i][j][k]*rfm.ReDD[i][j] + hDD[i][j]*rfm.ReDDdD[i][j][k] \
                                         + rfm.ghatDDdD[i][j][k]

    # Step 7b: Define barred Christoffel symbol \bar{\Gamma}^{i}_{jk} = GammabarUDD[i][j][k]
    GammabarUDD = ixp.zerorank3()
    for i in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                for m in range(DIM):
                    GammabarUDD[i][k][l] += (sp.Rational(1,2))*gammabarUU[i][m]* \
                                            (gammabarDD_dD[m][k][l] + gammabarDD_dD[m][l][k] - gammabarDD_dD[k][l][m])

    # Step 7c: Define \Delta^i_{jk} = \bar{\Gamma}^i_{jk} - \hat{\Gamma}^i_{jk} = DGammaUDD[i][j][k]
    global DGammaUDD # Needed for BSSN constraints, etc.
    DGammaUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                DGammaUDD[i][j][k] = GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k]

    # Step 7d: Define \Delta^i = \bar{\gamma}^{jk} \Delta^i_{jk}
    DGammaU = ixp.zerorank1()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                DGammaU[i] += gammabarUU[j][k] * DGammaUDD[i][j][k]

    # Step 7e: Define \Delta_{ijk} = \bar{\gamma}_{im} \Delta^m_{jk}
    DGammaDDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for m in range(DIM):
                    DGammaDDD[i][j][k] += gammabarDD[i][m] * DGammaUDD[m][j][k]

    # Step 7e: Add third term to Ricci tensor:
    #      \Delta^{k} \Delta_{(i j) k} = 1/2 \Delta^{k} (\Delta_{i j k} + \Delta_{j i k})
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                RbarDD[i][j] += sp.Rational(1,2) * DGammaU[k] * (DGammaDDD[i][j][k] + DGammaDDD[j][i][k])

    # Step 7f: Add remaining terms to Ricci tensor:
    # \bar{\gamma}^{k l} (\Delta^{m}_{k i} \Delta_{j m l}
    #                   + \Delta^{m}_{k j} \Delta_{i m l}
    #                   + \Delta^{m}_{i k} \Delta_{m j l})
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    for m in range(DIM):
                        RbarDD[i][j] += gammabarUU[k][l] * (DGammaUDD[m][k][i]*DGammaDDD[j][m][l] +
                                                            DGammaUDD[m][k][j]*DGammaDDD[i][m][l] +
                                                            DGammaUDD[m][i][k]*DGammaDDD[m][j][l])

    # Step 8a: Define \beta^i and \beta^i_{,k} in terms of rescaled quantity vetU[i] and vetU_dD[i][j]:
    global betaU # Needed as global for BSSN matter source terms, etc.
    betaU = ixp.zerorank1()
    for i in range(DIM):
        betaU[i] = vetU[i] * rfm.ReU[i]

    vetU_dD = ixp.declarerank2("vetU_dD", "nosym")
    vetU_dupD = ixp.declarerank2("vetU_dupD", "nosym")  # Needed for \beta^i RHS
    vetU_dDD = ixp.declarerank3("vetU_dDD", "sym12")  # Needed for \bar{\Lambda}^i RHS
    betaU_dD = ixp.zerorank2()
    betaU_dupD = ixp.zerorank2()  # Needed for \beta^i RHS
    betaU_dDD = ixp.zerorank3()  # Needed for \bar{\Lambda}^i RHS
    for i in range(DIM):
        for j in range(DIM):
            betaU_dD[i][j] = vetU_dD[i][j] * rfm.ReU[i] + vetU[i] * rfm.ReUdD[i][j]
            betaU_dupD[i][j] = vetU_dupD[i][j] * rfm.ReU[i] + vetU[i] * rfm.ReUdD[i][j]  # Needed for \beta^i RHS
            for k in range(DIM):
                # Needed for \bar{\Lambda}^i RHS:
                betaU_dDD[i][j][k] = vetU_dDD[i][j][k] * rfm.ReU[i] + vetU_dD[i][j] * rfm.ReUdD[i][k] + \
                                     vetU_dD[i][k] * rfm.ReUdD[i][j] + vetU[i] * rfm.ReUdDD[i][j][k]

    # Step 8b: First term of \partial_t \bar{\gamma}_{i j} right-hand side:
    # \beta^k \bar{\gamma}_{ij,k} + \beta^k_{,i} \bar{\gamma}_{kj} + \beta^k_{,j} \bar{\gamma}_{ik}
    gammabar_rhsDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                gammabar_rhsDD[i][j] += betaU[k] * gammabarDD_dupD[i][j][k] + betaU_dD[k][i] * gammabarDD[k][j] \
                                        + betaU_dD[k][j] * gammabarDD[i][k]

    # Step 8c: Define \bar{A}_{ij} = a_{ij} \text{ReDD[i][j]} = AbarDD[i][j], and its contraction trAbar = \bar{A}^k_k
    global AbarDD # Needed for BSSN constraints, etc.
    AbarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            AbarDD[i][j] = aDD[i][j] * rfm.ReDD[i][j]
    trAbar = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            trAbar += gammabarUU[i][j] * AbarDD[i][j]

    # Step 8d: Define detgammabar, detgammabar_dD, and detgammabar_dDD (needed for \partial_t \bar{\Lambda}^i below)
    global detgammabar # Needed for BSSN constraints, etc.
    detgammabar = detgbarOverdetghat * rfm.detgammahat
    global detgammabar_dD # Needed for BSSN constraints, etc.
    detgammabar_dD = ixp.zerorank1()
    if par.parval_from_str("BSSN_SF.BSSN_RHSs::detgbarOverdetghat_equals_one") == "False":
        detgbarOverdetghat_dD = ixp.declarerank1("detgbarOverdetghat_dD")
    else:
        detgbarOverdetghat_dD = ixp.zerorank1()
    for i in range(DIM):
        detgammabar_dD[i] = detgbarOverdetghat_dD[i] * rfm.detgammahat + detgbarOverdetghat * rfm.detgammahatdD[i]

    detgammabar_dDD = ixp.zerorank2()
    if par.parval_from_str("BSSN_SF.BSSN_RHSs::detgbarOverdetghat_equals_one") == "False":
        detgbarOverdetghat_dDD = ixp.declarerank2("detgbarOverdetghat_dDD", "sym01")
    else:
        detgbarOverdetghat_dDD = ixp.zerorank2()

    for i in range(DIM):
        for j in range(DIM):
            detgammabar_dDD[i][j] = detgbarOverdetghat_dDD[i][j] * rfm.detgammahat + \
                                    detgbarOverdetghat_dD[i] * rfm.detgammahatdD[j] + \
                                    detgbarOverdetghat_dD[j] * rfm.detgammahatdD[i] + \
                                    detgbarOverdetghat * rfm.detgammahatdDD[i][j]

    # Step 8e: Compute the contraction \bar{D}_k \beta^k = \beta^k_{,k} + \frac{\beta^k \bar{\gamma}_{,k}}{2 \bar{\gamma}}
    Dbarbetacontraction = sp.sympify(0)
    for k in range(DIM):
        Dbarbetacontraction += betaU_dD[k][k] + betaU[k] * detgammabar_dD[k] / (2 * detgammabar)

    # Step 8f: Second term of \partial_t \bar{\gamma}_{i j} right-hand side:
    # \frac{2}{3} \bar{\gamma}_{i j} \left (\alpha \bar{A}_{k}^{k} - \bar{D}_{k} \beta^{k}\right )
    for i in range(DIM):
        for j in range(DIM):
            gammabar_rhsDD[i][j] += sp.Rational(2, 3) * gammabarDD[i][j] * (alpha * trAbar - Dbarbetacontraction)

    # Step 8g: Third term of \partial_t \bar{\gamma}_{i j} right-hand side:
    # -2 \alpha \bar{A}_{ij}
    for i in range(DIM):
        for j in range(DIM):
            gammabar_rhsDD[i][j] += -2 * alpha * AbarDD[i][j]

    # Step 9a: First term of \partial_t \bar{A}_{i j}:
    # \beta^k \partial_k \bar{A}_{ij} + \partial_i \beta^k \bar{A}_{kj} + \partial_j \beta^k \bar{A}_{ik}

    # First define aDD_dupD:
    AbarDD_dupD = ixp.zerorank3()
    aDD_dupD = ixp.declarerank3("aDD_dupD","sym01")
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                AbarDD_dupD[i][j][k] += aDD_dupD[i][j][k]*rfm.ReDD[i][j] + aDD[i][j]*rfm.ReDDdD[i][j][k]

    Abar_rhsDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                Abar_rhsDD[i][j] += betaU[k]*AbarDD_dupD[i][j][k] + betaU_dD[k][i]*AbarDD[k][j] \
                                                                  + betaU_dD[k][j]*AbarDD[i][k]

    # Step 9b: Second term of \partial_t \bar{A}_{i j}:
    # - (2/3) \bar{A}_{i j} \bar{D}_{k} \beta^{k} - 2 \alpha \bar{A}_{i k} {\bar{A}^{k}}_{j} + \alpha \bar{A}_{i j} K
    AbarUD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                AbarUD[i][j] += gammabarUU[i][k]*AbarDD[k][j]

    for i in range(DIM):
        for j in range(DIM):
            Abar_rhsDD[i][j] += -sp.Rational(2,3)*AbarDD[i][j]*Dbarbetacontraction + alpha*AbarDD[i][j]*trK
            for k in range(DIM):
                Abar_rhsDD[i][j] += -2*alpha * AbarDD[i][k]*AbarUD[k][j]

    # Step 9c: Define partial derivatives of \phi in terms of evolved quantity "cf":
    cf_dD = ixp.declarerank1("cf_dD")
    cf_dupD = ixp.declarerank1("cf_dupD")  # Needed for \partial_t \phi next.
    cf_dDD = ixp.declarerank2("cf_dDD", "sym01")
    global phi_dD # Needed for BSSN constraints, etc.
    phi_dD = ixp.zerorank1()
    phi_dupD = ixp.zerorank1()
    phi_dDD = ixp.zerorank2()
    global exp_m4phi # Needed as global for BSSN matter source terms, etc.
    exp_m4phi = sp.sympify(0)
    if par.parval_from_str("BSSN_SF.BSSN_RHSs::ConformalFactor") == "phi":
        for i in range(DIM):
            phi_dD[i] = cf_dD[i]
            phi_dupD[i] = cf_dupD[i]
            for j in range(DIM):
                phi_dDD[i][j] = cf_dDD[i][j]
        exp_m4phi = sp.exp(-4 * cf)
    elif par.parval_from_str("BSSN_SF.BSSN_RHSs::ConformalFactor") == "W":
        # \partial_i W = \partial_i (e^{-2 phi}) = -2 e^{-2 phi} \partial_i phi
        # -> \partial_i phi = -\partial_i cf / (2 cf)
        for i in range(DIM):
            phi_dD[i] = - cf_dD[i] / (2 * cf)
            phi_dupD[i] = - cf_dupD[i] / (2 * cf)
            for j in range(DIM):
                # \partial_j \partial_i phi = - \partial_j [\partial_i cf / (2 cf)]
                #                           = - cf_{,ij} / (2 cf) + \partial_i cf \partial_j cf / (2 cf^2)
                phi_dDD[i][j] = (- cf_dDD[i][j] + cf_dD[i] * cf_dD[j] / cf) / (2 * cf)
        exp_m4phi = cf * cf
    elif par.parval_from_str("BSSN_SF.BSSN_RHSs::ConformalFactor") == "chi":
        # \partial_i chi = \partial_i (e^{-4 phi}) = -4 e^{-4 phi} \partial_i phi
        # -> \partial_i phi = -\partial_i cf / (4 cf)
        for i in range(DIM):
            phi_dD[i] = - cf_dD[i] / (4 * cf)
            phi_dupD[i] = - cf_dupD[i] / (4 * cf)
            for j in range(DIM):
                # \partial_j \partial_i phi = - \partial_j [\partial_i cf / (4 cf)]
                #                           = - cf_{,ij} / (4 cf) + \partial_i cf \partial_j cf / (4 cf^2)
                phi_dDD[i][j] = (- cf_dDD[i][j] + cf_dD[i] * cf_dD[j] / cf) / (4 * cf)
        exp_m4phi = cf
    else:
        print("Error: ConformalFactor == " + par.parval_from_str("BSSN_RHSs::ConformalFactor") + " unsupported!")
        exit(1)

    # Step 9d: Define phi_dBarD = phi_dD (since phi is a scalar) and phi_dBarDD (covariant derivative)
    #          \bar{D}_i \bar{D}_j \phi = \phi_{;\bar{i}\bar{j}} = \bar{D}_i \phi_{,j}
    #                                   = \phi_{,ij} - \bar{\Gamma}^k_{ij} \phi_{,k}
    global phi_dBarD # Needed for BSSN constraints, etc.
    phi_dBarD = phi_dD
    global phi_dBarDD # Needed for BSSN constraints, etc.
    phi_dBarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            phi_dBarDD[i][j] = phi_dDD[i][j]
            for k in range(DIM):
                phi_dBarDD[i][j] += - GammabarUDD[k][i][j] * phi_dD[k]

    # Step 9e: Define first and second derivatives of \alpha, as well as
    #         \bar{D}_i \bar{D}_j \alpha, which is defined just like phi
    alpha_dD = ixp.declarerank1("alpha_dD")
    alpha_dDD = ixp.declarerank2("alpha_dDD", "sym01")
    alpha_dBarD = alpha_dD
    alpha_dBarDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            alpha_dBarDD[i][j] = alpha_dDD[i][j]
            for k in range(DIM):
                alpha_dBarDD[i][j] += - GammabarUDD[k][i][j] * alpha_dD[k]

    # Step 9f: Define the terms in curly braces:
    curlybrackettermsDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            curlybrackettermsDD[i][j] = -2 * alpha * phi_dBarDD[i][j] + 4 * alpha * phi_dBarD[i] * phi_dBarD[j] \
                                        + 2 * alpha_dBarD[i] * phi_dBarD[j] \
                                        + 2 * alpha_dBarD[j] * phi_dBarD[i] \
                                        - alpha_dBarDD[i][j] + alpha * RbarDD[i][j]

    # Step 9g: Compute the trace:
    curlybracketterms_trace = sp.sympify(0)
    for i in range(DIM):
        for j in range(DIM):
            curlybracketterms_trace += gammabarUU[i][j] * curlybrackettermsDD[i][j]

    # Step 9h: Third and final term of Abar_rhsDD[i][j]:
    for i in range(DIM):
        for j in range(DIM):
            Abar_rhsDD[i][j] += exp_m4phi * (curlybrackettermsDD[i][j] - \
                                             sp.Rational(1, 3) * gammabarDD[i][j] * curlybracketterms_trace)

    # Step 10: right-hand side of conformal factor variable "cf". Supported
    #          options include: cf=phi, cf=W=e^(-2*phi) (default), and cf=chi=e^(-4*phi)
    # \partial_t phi = \left[\beta^k \partial_k \phi \right] <- TERM 1
    #                  + \frac{1}{6} \left (\bar{D}_{k} \beta^{k} - \alpha K \right ) <- TERM 2
    global cf_rhs
    cf_rhs = sp.Rational(1,6) * (Dbarbetacontraction - alpha*trK) # Term 2
    for k in range(DIM):
        cf_rhs += betaU[k]*phi_dupD[k] # Term 1

    # Next multiply to convert phi_rhs to cf_rhs.
    if par.parval_from_str("BSSN_SF.BSSN_RHSs::ConformalFactor") == "phi":
        pass # do nothing; cf_rhs = phi_rhs
    elif par.parval_from_str("BSSN_SF.BSSN_RHSs::ConformalFactor") == "W":
        cf_rhs *= -2*cf # cf_rhs = -2*cf*phi_rhs
    elif par.parval_from_str("BSSN_SF.BSSN_RHSs::ConformalFactor") == "chi":
        cf_rhs *= -4*cf # cf_rhs = -4*cf*phi_rhs
    else:
        print("Error: ConformalFactor == "+par.parval_from_str("BSSN_RHSs::ConformalFactor")+" unsupported!")
        exit(1)

    # Step 11: right-hand side of trK (trace of extrinsic curvature):
    # \partial_t K = \beta^k \partial_k K <- TERM 1
    #           + \frac{1}{3} \alpha K^{2} <- TERM 2
    #           + \alpha \bar{A}_{i j} \bar{A}^{i j} <- TERM 3
    #           - - e^{-4 \phi} (\bar{D}_{i} \bar{D}^{i} \alpha + 2 \bar{D}^{i} \alpha \bar{D}_{i} \phi ) <- TERM 4
    global trK_rhs
    trK_rhs = sp.Rational(1,3)*alpha*trK*trK # Term 2
    trK_dupD = ixp.declarerank1("trK_dupD")
    for k in range(DIM):
        trK_rhs += betaU[k]*trK_dupD[k] # Term 1
    for i in range(DIM):
        for j in range(DIM):
            trK_rhs += -exp_m4phi*gammabarUU[i][j]*(alpha_dBarDD[i][j] + 2*alpha_dBarD[j]*phi_dBarD[i]) # Term 4
    global AbarUU # Needed for BSSN constraints, etc.
    AbarUU = ixp.zerorank2() # Needed also for \partial_t \bar{\Lambda}^i
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for l in range(DIM):
                    AbarUU[i][j] += gammabarUU[i][k]*gammabarUU[j][l]*AbarDD[k][l]
    for i in range(DIM):
        for j in range(DIM):
            trK_rhs += alpha*AbarDD[i][j]*AbarUU[i][j] # Term 3

    # Step 12: right-hand side of \partial_t \bar{\Lambda}^i:
    # \partial_t \bar{\Lambda}^i = \beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k <- TERM 1
    #                            + \bar{\gamma}^{j k} \hat{D}_{j} \hat{D}_{k} \beta^{i} <- TERM 2
    #                            + \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j} <- TERM 3
    #                            + \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j} <- TERM 4
    #                            - 2 \bar{A}^{i j} (\partial_{j} \alpha - 6 \partial_{j} \phi) <- TERM 5
    #                            + 2 \alpha \bar{A}^{j k} \Delta_{j k}^{i} <- TERM 6
    #                            - \frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K <- TERM 7

    # Step 12a: Term 1 of \partial_t \bar{\Lambda}^i: \beta^k \partial_k \bar{\Lambda}^i - \partial_k \beta^i \bar{\Lambda}^k
    # First we declare \bar{\Lambda}^i and \bar{\Lambda}^i_{,j} in terms of \lambda^i and \lambda^i_{,j}
    LambarU = ixp.zerorank1()
    LambarU_dupD = ixp.zerorank2()
    lambdaU_dupD = ixp.declarerank2("lambdaU_dupD","nosym")
    for i in range(DIM):
        LambarU[i] = lambdaU[i]*rfm.ReU[i]
        for j in range(DIM):
            LambarU_dupD[i][j] = lambdaU_dupD[i][j]*rfm.ReU[i] + lambdaU[i]*rfm.ReUdD[i][j]

    Lambar_rhsU = ixp.zerorank1()
    for i in range(DIM):
        for k in range(DIM):
            Lambar_rhsU[i] += betaU[k]*LambarU_dupD[i][k] - betaU_dD[i][k]*LambarU[k] # Term 1

    # Step 12b: Term 2 of \partial_t \bar{\Lambda}^i = \bar{\gamma}^{jk} (Term 2a + Term 2b + Term 2c)
    # Term 2a: \bar{\gamma}^{jk} \beta^i_{,kj}
    Term2aUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                Term2aUDD[i][j][k] += betaU_dDD[i][k][j]
    # Term 2b: \hat{\Gamma}^i_{mk,j} \beta^m + \hat{\Gamma}^i_{mk} \beta^m_{,j}
    #          + \hat{\Gamma}^i_{dj}\beta^d_{,k} - \hat{\Gamma}^d_{kj} \beta^i_{,d}
    Term2bUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for m in range(DIM):
                    Term2bUDD[i][j][k] += rfm.GammahatUDDdD[i][m][k][j]*betaU[m]    \
                                          + rfm.GammahatUDD[i][m][k]*betaU_dD[m][j] \
                                          + rfm.GammahatUDD[i][m][j]*betaU_dD[m][k] \
                                          - rfm.GammahatUDD[m][k][j]*betaU_dD[i][m]
    # Term 2c: \hat{\Gamma}^i_{dj}\hat{\Gamma}^d_{mk} \beta^m - \hat{\Gamma}^d_{kj} \hat{\Gamma}^i_{md} \beta^m
    Term2cUDD = ixp.zerorank3()
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                for m in range(DIM):
                    for d in range(DIM):
                        Term2cUDD[i][j][k] += ( rfm.GammahatUDD[i][d][j]*rfm.GammahatUDD[d][m][k] \
                                               -rfm.GammahatUDD[d][k][j]*rfm.GammahatUDD[i][m][d])*betaU[m]

    Lambar_rhsUpieceU = ixp.zerorank1()

    # Put it all together to get Term 2:
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                Lambar_rhsU[i] += gammabarUU[j][k] * (Term2aUDD[i][j][k] + Term2bUDD[i][j][k] + Term2cUDD[i][j][k])
                Lambar_rhsUpieceU[i] += gammabarUU[j][k] * (Term2aUDD[i][j][k] + Term2bUDD[i][j][k] + Term2cUDD[i][j][k])

    # Step 12c: Term 3 of \partial_t \bar{\Lambda}^i:
    #    \frac{2}{3} \Delta^{i} \bar{D}_{j} \beta^{j}
    for i in range(DIM):
        Lambar_rhsU[i] += sp.Rational(2,3)*DGammaU[i]*Dbarbetacontraction # Term 3

    # Step 12d: Term 4 of \partial_t \bar{\Lambda}^i:
    #           \frac{1}{3} \bar{D}^{i} \bar{D}_{j} \beta^{j}
    Dbarbetacontraction_dBarD = ixp.zerorank1()
    for k in range(DIM):
        for m in range(DIM):
            Dbarbetacontraction_dBarD[m] += betaU_dDD[k][k][m] + \
                                            (betaU_dD[k][m]*detgammabar_dD[k] +
                                             betaU[k]*detgammabar_dDD[k][m])/(2*detgammabar) \
                                            -betaU[k]*detgammabar_dD[k]*detgammabar_dD[m]/(2*detgammabar*detgammabar)
    for i in range(DIM):
        for m in range(DIM):
            Lambar_rhsU[i] += sp.Rational(1,3)*gammabarUU[i][m]*Dbarbetacontraction_dBarD[m]

    # Step 12e: Term 5 of \partial_t \bar{\Lambda}^i:
    #           - 2 \bar{A}^{i j} (\partial_{j} \alpha - 6 \alpha \partial_{j} \phi)
    for i in range(DIM):
        for j in range(DIM):
            Lambar_rhsU[i] += -2*AbarUU[i][j]*(alpha_dD[j] - 6*alpha*phi_dD[j])

    # Step 12f: Term 6 of \partial_t \bar{\Lambda}^i:
    #           2 \alpha \bar{A}^{j k} \Delta^{i}_{j k}
    for i in range(DIM):
        for j in range(DIM):
            for k in range(DIM):
                Lambar_rhsU[i] += 2*alpha*AbarUU[j][k]*DGammaUDD[i][j][k]

    # Step 12g: Term 7 of \partial_t \bar{\Lambda}^i:
    #           -\frac{4}{3} \alpha \bar{\gamma}^{i j} \partial_{j} K
    trK_dD = ixp.declarerank1("trK_dD")
    for i in range(DIM):
        for j in range(DIM):
            Lambar_rhsU[i] += -sp.Rational(4,3)*alpha*gammabarUU[i][j]*trK_dD[j]

    # Step 13: \partial_t \alpha = \beta^i \alpha_{,i} - 2*\alpha*K
    global alpha_rhs
    alpha_rhs = -2*alpha*trK
    alpha_dupD = ixp.declarerank1("alpha_dupD")
    for i in range(DIM):
        alpha_rhs += betaU[i]*alpha_dupD[i]

    # Step 14: Set \partial_t \beta^i
    beta_rhsU = ixp.zerorank1()
    if par.parval_from_str("BSSN_SF.BSSN_RHSs::ShiftEvolutionOption") == "GammaDriving2ndOrder_NoCovariant":
        # Step 14 Option 1: \partial_t \beta^i = \beta^j \beta^i_{,j} + B^i
        # First define BU, in terms of rescaled variable \bet^i
        BU = ixp.zerorank1()
        for i in range(DIM):
            BU[i] = betU[i]*rfm.ReU[i]

        # Then compute right-hand side:
        for i in range(DIM):
            beta_rhsU[i] += BU[i]
            for j in range(DIM):
                beta_rhsU[i] += betaU[j]*betaU_dupD[i][j]
    if par.parval_from_str("BSSN_SF.BSSN_RHSs::ShiftEvolutionOption") == "GammaDriving2ndOrder_Covariant":
        # Step 14 Option 2: \partial_t \beta^i = \left[\beta^j \bar{D}_j \beta^i\right] + B^{i}
        # First define BU, in terms of rescaled variable \bet^i
        BU = ixp.zerorank1()
        for i in range(DIM):
            BU[i] = betU[i]*rfm.ReU[i]

        # Then compute right-hand side:
        # Term 1: \beta^j \beta^i_{,j}
        for i in range(DIM):
            for j in range(DIM):
                beta_rhsU[i] += betaU[j]*betaU_dupD[i][j]

        # Term 2: \beta^j \bar{\Gamma}^i_{mj} \beta^m
        for i in range(DIM):
            for j in range(DIM):
                for m in range(DIM):
                    beta_rhsU[i] += betaU[j]*GammabarUDD[i][m][j]*betaU[m]
        # Term 3: B^i
        for i in range(DIM):
            beta_rhsU[i] += BU[i]
        
    # Step 15: Evaluate \partial_t B^i

    # Step 15a: Declare free parameter eta and B_rhsU
    eta = par.Cparameters("REAL",thismodule,["eta"],1.)
    B_rhsU = ixp.zerorank1()

    if par.parval_from_str("BSSN_SF.BSSN_RHSs::ShiftEvolutionOption") == "GammaDriving2ndOrder_NoCovariant":
        # Step 15: Non-covariant option:
        #  \partial_t B^i = \beta^j \partial_j B^i
        #                 + \frac{3}{4} \partial_{0} \bar{\Lambda}^{i} - \eta B^{i}
        # Step 15b: Define BU_dupD, in terms of derivative of rescaled variable \bet^i
        BU_dupD = ixp.zerorank2()
        betU_dupD = ixp.declarerank2("betU_dupD","nosym")
        for i in range(DIM):
            for j in range(DIM):
                BU_dupD[i][j] = betU_dupD[i][j]*rfm.ReU[i] + betU[i]*rfm.ReUdD[i][j]

        # Step 15c: Compute \partial_0 \bar{\Lambda}^i = (\partial_t - \beta^i \partial_i) \bar{\Lambda}^j 
        Lambar_partial0 = ixp.zerorank1()
        for i in range(DIM):
            Lambar_partial0[i] = Lambar_rhsU[i]
        for i in range(DIM):
            for j in range(DIM):
                Lambar_partial0[j] += -betaU[i]*LambarU_dupD[j][i]

        # Step 15d: Evaluate RHS of B^i:
        for i in range(DIM):
            B_rhsU[i] += sp.Rational(3,4)*Lambar_partial0[i] - eta*BU[i]
            for j in range(DIM):
                B_rhsU[i] += betaU[j]*BU_dupD[i][j]

    if par.parval_from_str("BSSN_SF.BSSN_RHSs::ShiftEvolutionOption") == "GammaDriving2ndOrder_Covariant":
        # Step 15: Covariant option:
        #  \partial_t B^i = \beta^j \bar{D}_j B^i
        #               + \frac{3}{4} ( \partial_t \bar{\Lambda}^{i} - \beta^j \bar{D}_j \bar{\Lambda}^{i} ) 
        #               - \eta B^{i}
        #                 = \beta^j B^i_{,j} + \beta^j \bar{\Gamma}^i_{mj} B^m
        #               + \frac{3}{4}[ \partial_t \bar{\Lambda}^{i} 
        #                            - \beta^j (\bar{\Lambda}^i_{,j} + \bar{\Gamma}^i_{mj} \bar{\Lambda}^m)] 
        #               - \eta B^{i}
        # Term 1, part a: First compute B^i_{,j} using upwinded derivative
        BU_dupD = ixp.zerorank2()
        betU_dupD = ixp.declarerank2("betU_dupD","nosym")
        for i in range(DIM):
            for j in range(DIM):
                BU_dupD[i][j] = betU_dupD[i][j]*rfm.ReU[i] + betU[i]*rfm.ReUdD[i][j]
        # Term 1: \beta^j B^i_{,j}
        for i in range(DIM):
            for j in range(DIM):
                B_rhsU[i] += betaU[j]*BU_dupD[i][j]
        # Term 2: \beta^j \bar{\Gamma}^i_{mj} B^m
        for i in range(DIM):
            for j in range(DIM):
                for m in range(DIM):
                    B_rhsU[i] += betaU[j]*GammabarUDD[i][m][j]*BU[m]
        # Term 3: \frac{3}{4}\partial_t \bar{\Lambda}^{i}
        for i in range(DIM):
            B_rhsU[i] += sp.Rational(3,4)*Lambar_rhsU[i]
        # Term 4: -\frac{3}{4}\beta^j \bar{\Lambda}^i_{,j}
        for i in range(DIM):
            for j in range(DIM):
                B_rhsU[i] += -sp.Rational(3,4)*betaU[j]*LambarU_dupD[i][j]
        # Term 5: -\frac{3}{4}\beta^j \bar{\Gamma}^i_{mj} \bar{\Lambda}^m
        for i in range(DIM):
            for j in range(DIM):
                for m in range(DIM):
                    B_rhsU[i] += -sp.Rational(3,4)*betaU[j]*GammabarUDD[i][m][j]*LambarU[m]
        # Term 6: - \eta B^i
        for i in range(DIM):
            B_rhsU[i] += -eta*BU[i]          
            
            
   ## Step 16: Add scalar field     
   # Step 16a: Define physical metric 
    gammaUU = ixp.zerorank2()
    gammaDD = ixp.zerorank2()
    
    for i in range(DIM):
        for j in range(DIM):
            gammaUU[i][j] = gammabarUU[i][j] * exp_m4phi
            gammaDD[i][j] = gammabarDD[i][j] / exp_m4phi
            
    # Step 16b: Declare the rank-1 indexed expressions \partial_{i} u and \partial_{i} v
    uu_dD = ixp.declarerank1("uu_dD")
    vv_dD = ixp.declarerank1("vv_dD")
    
    # Step 16c: Declare the rank-2 indexed expression \partial_{ij} u,
    #          which is symmetric about interchange of indices i and j
    #          Derivative variables like these must have an underscore
    #          in them, so the finite difference module can parse the
    #          variable name properly.
    uu_dDD = ixp.declarerank2("uu_dDD","sym01")

    # Step 16d: Specify RHSs as global variables,
    #         to enable access outside this
    #         function (e.g., for C code output)
    global uu_rhs,vv_rhs
    
    # Define potential and derivative
    # if (pot1_on == 1): # quadratic potential
    #    Vofu = scalarmass*scalarmass*uu*uu / 2.
    #    dVdu = scalarmass*scalarmass*uu
    # if (pot2_on == 1): # self interacting potential
    #    Vofu = fa*fa*scalarmass*scalarmass*(1. - syp.cos(uu/fa))
    #    dVdu = fa*scalarmass*scalarmass*syp.sin(uu/fa)
    Vofu = pot1_on * scalarmass*scalarmass*uu*uu / 2. + pot2_on * fa*fa*scalarmass*scalarmass*(1. - syp.cos(uu/fa))
    dVdu = pot1_on * scalarmass*scalarmass*uu + pot2_on * fa*scalarmass*scalarmass*syp.sin(uu/fa)
    # Step 17: Define right-hand sides for the evolution.
    # Step 17a: uu_rhs = \beta^i \partial_i u - \alpha*vv:
    uu_rhs = - alpha*vv
    
    uu_dupD = ixp.declarerank1("uu_dupD")
    for i in range(DIM):
        uu_rhs += betaU[i]*uu_dupD[i] # Term 1
    
#    for i in range(DIM):
#        uu_rhs += betaU[i]*uu_dD[i]
    # Step 17b: The right-hand side of the \partial_t v equation
    #          is given by:
    #    Part 1:  
    #       \alpha*trK*vv + dV/du 
    #    Part 2:  
    #       \beta^i \partial_i v
    #    Part 3:
    #       - \alpha \gamma^{ij} \partial_i \partial_j u 
    #       - \alpha \gamma^{ij} \partial_j u \partial_i \alpha 
    #       - 2* \gamma^{ij} \partial_j \phi \partial_i u
    #    Part 4:
    #       + \alpha * \gamma^{ij} * \bar{\Gamma}^k_{ij} * \partial_k u
    
    # PART 1: 
    vv_rhs = alpha*trK*vv + alpha*dVdu
    # PART 2: 
    for i in range(DIM):  
        vv_rhs += betaU[i]*vv_dD[i]
    # PART 3:     
        for j in range(DIM):
            vv_rhs += - alpha*gammaUU[i][j]*uu_dDD[j][i] - gammaUU[i][j]*uu_dD[j]*alpha_dD[i] - 2.0*alpha*gammaUU[i][j]*phi_dD[j]*uu_dD[i]
    # PART 4: 
            for k in range(DIM):
                vv_rhs +=  alpha*gammaUU[i][j]*GammabarUDD[k][i][j]*uu_dD[k] 
    
    
    ### Step 17: Adding backreaction terms
    
    global rho
    global SU
    SDD = ixp.zerorank2()
    SD  = ixp.zerorank1()
    SU  = ixp.zerorank1()
    S   = sp.sympify(0)
    rho = sp.sympify(0)
    
    # Step 17a: Define some useful quantities:
    #     vv2 = v*v and psi2 = partial_i u * \partial^i u
    vv2 = vv*vv
    psi2 = 0
    for i in range(DIM):
        for j in range(DIM):
            psi2 += gammaUU[i][j]*uu_dD[i]*uu_dD[j]
    
    # Step 17b: Calculate components of EM Tensor
    #     S_ij = \partial_i uu * partial_j uu + 0.5 * gamma_{ij} * (vv2 - psi2) - \gamma_{ij} V(u)
    for i in range(DIM):
        for j in range(DIM):
            SDD[i][j] = uu_dD[i]*uu_dD[j] + 0.5*gammaDD[i][j] * (vv2 - psi2) - gammaDD[i][j]*Vofu

    #     S = Tr_S_ij = exp_m4phi*gammabarUU[i][j]*SDD[i][j]
    #    for i in range(DIM):
    #        for j in range(DIM):
    #            S += gammaUU[i][j]*SDD[i][j]
    S = 1.5*vv2 - 0.5*psi2 - 3.*Vofu
    
    # S_i
    for i in range(DIM):
        SD[i] = uu_dD[i]*vv
        
    # S^i
    for i in range(DIM):
        for j in range(DIM):
            SU[i] += gammaUU[i][j]*SD[j]
            
    # rho 
    rho = 0.5 * (vv2 + psi2) + Vofu        
            
    # Step 18: Add matter source term to BSSN RHS:
    # Step 18a: Abar_RHS
    S_TFDD = ixp.zerorank2()
    sTFDD = ixp.zerorank2()
    for i in range(DIM):
        for j in range(DIM):
            S_TFDD[i][j] = SDD[i][j] - sp.Rational(1,3)*gammaDD[i][j]*S
            Abar_rhsDD[i][j] += (-8*PI*alpha * exp_m4phi * S_TFDD[i][j]) * br_on

    # Step 18b: \partial_t K 
    trK_rhs += (4*PI*alpha*(rho + S)) * br_on
    
    # Step 18c: Lambar_RHS
    for i in range(DIM):
        for j in range(DIM):
            Lambar_rhsU[i] += (-16*PI*alpha*gammabarUU[i][j]*SD[j]) * br_on
            
    # Step 19: Rescale the RHS quantities so that the evolved
    #          variables are smooth across coord singularities
    global h_rhsDD,a_rhsDD,lambda_rhsU,vet_rhsU,bet_rhsU
    h_rhsDD     = ixp.zerorank2()
    a_rhsDD     = ixp.zerorank2()
    lambda_rhsU = ixp.zerorank1()
    vet_rhsU    = ixp.zerorank1()
    bet_rhsU    = ixp.zerorank1()
    for i in range(DIM):
        lambda_rhsU[i] = Lambar_rhsU[i] / rfm.ReU[i]
        vet_rhsU[i]    =   beta_rhsU[i] / rfm.ReU[i]
        bet_rhsU[i]    =      B_rhsU[i] / rfm.ReU[i]
        for j in range(DIM):
            h_rhsDD[i][j] = gammabar_rhsDD[i][j] / rfm.ReDD[i][j]
            a_rhsDD[i][j] =     Abar_rhsDD[i][j] / rfm.ReDD[i][j]
