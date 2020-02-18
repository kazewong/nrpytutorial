import numpy as np
import matplotlib.pyplot as plt 

data = np.genfromtxt("BSSN_SF-output2D/quad_pot_2d_t-00000000.txt")

N1 = int(len(data[:,0])/4)

time = data[:N1,0]
x = data[:N1,1]

gammaCartDD00 = data[:N1,2]
gammaCartDD01 = data[:N1,3]
gammaCartDD02 = data[:N1,4]
gammaCartDD11 = data[:N1,5]
gammaCartDD12 = data[:N1,6]
gammaCartDD22 = data[:N1,7]
KCartDD00   = data[:N1,8]
KCartDD01   = data[:N1,9]
KCartDD02   = data[:N1,10]
KCartDD11   = data[:N1,11]
KCartDD12   = data[:N1,12]
KCartDD22   = data[:N1,13]
betaCartU0  = data[:N1,14]
betaCartU1  = data[:N1,15]
betaCartU2  = data[:N1,16]
BCartU0  = data[:N1,17]
BCartU1  = data[:N1,18]
BCartU2  = data[:N1,19]
alphaCart = data[:N1,20]

uu = data[:N1,21]
vv = data[:N1,22]

dx = x[1] - x[0]
order = 2 
dgammaCartDD00dr = np.gradient(gammaCartDD00,dx,edge_order = order)
dgammaCartDD01dr = np.gradient(gammaCartDD01,dx,edge_order = order)
dgammaCartDD02dr = np.gradient(gammaCartDD02,dx,edge_order = order)
dgammaCartDD11dr = np.gradient(gammaCartDD11,dx,edge_order = order)
dgammaCartDD12dr = np.gradient(gammaCartDD12,dx,edge_order = order)
dgammaCartDD22dr = np.gradient(gammaCartDD22,dx,edge_order = order)


d2gammaCartDD00dr2 = np.gradient(dgammaCartDD00dr,dx,edge_order = order)
d2gammaCartDD01dr2 = np.gradient(dgammaCartDD01dr,dx,edge_order = order)
d2gammaCartDD02dr2 = np.gradient(dgammaCartDD02dr,dx,edge_order = order)
d2gammaCartDD11dr2 = np.gradient(dgammaCartDD11dr,dx,edge_order = order)
d2gammaCartDD12dr2 = np.gradient(dgammaCartDD12dr,dx,edge_order = order)
d2gammaCartDD22dr2 = np.gradient(dgammaCartDD22dr,dx,edge_order = order)

gammaDD = np.zeros((N1,3,3))
KDD     = np.zeros((N1,3,3))

for i in range(N1):
    gammaDD[i] =  np.array([[gammaCartDD00[i],gammaCartDD01[i],gammaCartDD02[i]],
                            [gammaCartDD01[i],gammaCartDD11[i],gammaCartDD12[i]],
                            [gammaCartDD02[i],gammaCartDD12[i],gammaCartDD22[i]]])
    
    KDD[i]     =  np.array([[KCartDD00[i],KCartDD01[i],KCartDD02[i]],
                            [KCartDD01[i],KCartDD11[i],KCartDD12[i]],
                            [KCartDD02[i],KCartDD12[i],KCartDD22[i]]])


gammaUU = np.linalg.inv(gammaDD)

TrK = np.zeros(N1)

for i in range(N1):
    for k in range(3):
        for d in range(3):
            TrK[i] += KDD[i][k,d]*gammaUU[i][k,d]
        
KijKij = np.zeros(N1)
for i in range(N1):
    for k in range(3):
        for d in range(3):
            for l in range(3):
                for m in range(3):
                    KijKij[i] += KDD[i][k,d]*KDD[i][l,m]*gammaUU[i][k,l]*gammaUU[i][d,m]


uu_dD = np.zeros((N1,3))

vv2 = vv*vv
psi2 = 0
for i in range(DIM):
    for j in range(DIM):
        psi2 += gammaUU[i][j]*uu_dD[i]*uu_dD[j]
    
    # rho 
rho = 0.5 * (vv2 + psi2) + Vofu        
     


print(gammaDD[0,:,:])
    





