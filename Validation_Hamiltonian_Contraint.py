import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("BSSN_SF-output2D/quad_pot_2d_t-00000000.txt")

DIM = 3
N1 = int(len(data[:, 0])/4)

time = data[:N1, 0]
x = data[:N1, 1]
y = data[:N1, 2]
z = data[:N1, 3]

r = np.sqrt(x**2 + y**2 + z**2)

gammaCartDD00 = data[:N1, 4]
gammaCartDD01 = data[:N1, 5]
gammaCartDD02 = data[:N1, 6]
gammaCartDD11 = data[:N1, 7]
gammaCartDD12 = data[:N1, 8]
gammaCartDD22 = data[:N1, 9]
KCartDD00 = data[:N1, 10]
KCartDD01 = data[:N1, 11]
KCartDD02 = data[:N1, 12]
KCartDD11 = data[:N1, 13]
KCartDD12 = data[:N1, 14]
KCartDD22 = data[:N1, 15]
betaCartU0 = data[:N1, 16]
betaCartU1 = data[:N1, 17]
betaCartU2 = data[:N1, 18]
BCartU0 = data[:N1, 19]
BCartU1 = data[:N1, 20]
BCartU2 = data[:N1, 21]
alphaCart = data[:N1, 22]
hamreal = data[:N1,23]

uu = data[:N1, 23]
vv = data[:N1, 24]

dr = r[1] - r[0]
order = 2
dgammaCartDD00dr = np.gradient(gammaCartDD00, dr, edge_order=order)
dgammaCartDD01dr = np.gradient(gammaCartDD01, dr, edge_order=order)
dgammaCartDD02dr = np.gradient(gammaCartDD02, dr, edge_order=order)
dgammaCartDD11dr = np.gradient(gammaCartDD11, dr, edge_order=order)
dgammaCartDD12dr = np.gradient(gammaCartDD12, dr, edge_order=order)
dgammaCartDD22dr = np.gradient(gammaCartDD22, dr, edge_order=order)


d2gammaCartDD00dr2 = np.gradient(dgammaCartDD00dr, dr, edge_order=order)
d2gammaCartDD01dr2 = np.gradient(dgammaCartDD01dr, dr, edge_order=order)
d2gammaCartDD02dr2 = np.gradient(dgammaCartDD02dr, dr, edge_order=order)
d2gammaCartDD11dr2 = np.gradient(dgammaCartDD11dr, dr, edge_order=order)
d2gammaCartDD12dr2 = np.gradient(dgammaCartDD12dr, dr, edge_order=order)
d2gammaCartDD22dr2 = np.gradient(dgammaCartDD22dr, dr, edge_order=order)

gammaDD = np.zeros((N1, 3, 3))
dgammaDDdr = np.zeros(( 3, 3))
dgammaDDdx = np.zeros((N1, 3, 3, 3))
KDD = np.zeros((N1, 3, 3))

for i in range(N1):
    gammaDD[i] = np.array([[gammaCartDD00[i],
                            gammaCartDD01[i],
                            gammaCartDD02[i]],
                           [gammaCartDD01[i],
                            gammaCartDD11[i],
                            gammaCartDD12[i]],
                           [gammaCartDD02[i],
                            gammaCartDD12[i],
                            gammaCartDD22[i]]])

    dgammaDDdr = np.array([[dgammaCartDD00dr[i],
                            dgammaCartDD01dr[i],
                            dgammaCartDD02dr[i]],
                           [dgammaCartDD01dr[i],
                            dgammaCartDD11dr[i],
                            dgammaCartDD12dr[i]],
                           [dgammaCartDD02dr[i],
                            dgammaCartDD12dr[i],
                            dgammaCartDD22dr[i]]])

    dgammaDDdx[i][0] = x[i]/r[i] * dgammaDDdr
    dgammaDDdx[i][1] = y[i]/r[i] * dgammaDDdr
    dgammaDDdx[i][2] = z[i]/r[i] * dgammaDDdr

    KDD[i] = np.array([[KCartDD00[i], KCartDD01[i], KCartDD02[i]],
                       [KCartDD01[i], KCartDD11[i], KCartDD12[i]],
                       [KCartDD02[i], KCartDD12[i], KCartDD22[i]]])


gammaUU = np.linalg.inv(gammaDD)

TrK = np.zeros(N1)

for i in range(N1):
    for k in range(DIM):
        for d in range(DIM):
            TrK[i] += KDD[i][k, d] * gammaUU[i][k, d]

KijKij = np.zeros(N1)
for i in range(N1):
    for k in range(DIM):
        for d in range(DIM):
            for l in range(DIM):
                for m in range(DIM):
                    KijKij[i] += KDD[i][k, d] * KDD[i][l, m] * \
                        gammaUU[i][k, l] * gammaUU[i][d, m]


uu_dD = np.zeros((N1, 3))
uu_dr = np.gradient(dgammaCartDD22dr, dr, edge_order=order)
uu_dD[:, 0] = x / r * uu_dr
uu_dD[:, 1] = y / r * uu_dr
uu_dD[:, 2] = z / r * uu_dr

Vofu = 0

vv2 = vv * vv
psi2 = 0
for ind in range(N1):
    for j in range(DIM):
        for k in range(DIM):
            psi2 += gammaUU[ind][j, k] * uu_dD[ind][j] * uu_dD[ind][k]

rho = 0.5 * (vv2 + psi2) + Vofu

ChristoffelUDD = np.zeros((N1,DIM,DIM,DIM))
for ind in range(N1):
    for i in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                for m in range(DIM):
                    ChristoffelUDD[ind][i,k,l] += 0.5  * gammaUU[ind][i,m] * ( dgammaDDdx[ind][l,m,k]  + dgammaDDdx[ind][k,m,l] - dgammaDDdx[ind][m,k,l])

dChristoffelUDDdx = np.zeros((N1,DIM,DIM,DIM,DIM))
dChristoffelUDDdr = np.zeros((N1,DIM,DIM,DIM))
for i in range(DIM):
    for k in range(DIM):
        for l in range(DIM):
            dChristoffelUDDdr[:][i,k,l] = np.gradient(ChristoffelUDD[:][i,k,l], dr, edge_order=order)

for ind in range(N1):
    for i in range(DIM):
        for k in range(DIM):
            for l in range(DIM):
                dChristoffelUDDdx[ind][0,i,k,l] += x[ind]/r[ind] * dChristoffelUDDdr[ind][i,k,l] 
                dChristoffelUDDdx[ind][1,i,k,l] += y[ind]/r[ind] * dChristoffelUDDdr[ind][i,k,l] 
                dChristoffelUDDdx[ind][2,i,k,l] += z[ind]/r[ind] * dChristoffelUDDdr[ind][i,k,l]  

RiemannUDDD = np.zeros((N1,DIM,DIM,DIM,DIM))
for ind in range(N1):
    for l in range(DIM):
        for i in range(DIM):
            for j in range(DIM):
                for k in range(DIM):
                    RiemannUDDD[ind][l,i,j,k] += dChristoffelUDDdx[ind][j,l,i,k] - dChristoffelUDDdx[ind][k,l,i,j]  
                    for s in range(DIM):
                        RiemannUDDD[ind][l,i,j,k] += ChristoffelUDD[ind][l,j,s]*ChristoffelUDD[ind][s,i,k] - ChristoffelUDD[ind][l,k,s] * ChristoffelUDD[ind][s,i,j]

RicciTensorDD = np.zeros((N1,DIM,DIM))

for ind in range(N1):
    for l in range(DIM):
        for i in range(DIM):
            for j in range(DIM):
                RicciTensorDD[ind][l,i] += RiemannUDDD[ind][j,i,j,k]

RicciScalar = np.zeros(N1)
for ind in range(N1):
    for l in range(DIM):
        for i in range(DIM):
            RicciScalar[ind] += RicciTensorDD[ind][l,i]*gammaUU[ind][l][i] 

HamConstraint = RicciScalar + KijKij - TrK**2 - 16*np.pi*rho

#plt.plot(rho, label =  "rho ",linestyle = ":" )
#plt.plot(HamConstraint,label = "Calculated here ")
#plt.plot(hamreal, label =  "Calculated internally " )
plt.plot(RicciScalar, label =  "Ricci ",linestyle = "--" )
plt.legend()
plt.show()
