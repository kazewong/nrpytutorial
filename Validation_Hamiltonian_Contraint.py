import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("BSSN_SF-output2D/quad_pot_2d_t-00000000.txt")

DIM = 3
N1 = int(len(data[:, 0]) / 4)

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

    KDD[i] = np.array([[KCartDD00[i], KCartDD01[i], KCartDD02[i]],
                       [KCartDD01[i], KCartDD11[i], KCartDD12[i]],
                       [KCartDD02[i], KCartDD12[i], KCartDD22[i]]])


gammaUU = np.linalg.inv(gammaDD)

TrK = np.zeros(N1)

for i in range(N1):
    for k in range(3):
        for d in range(3):
            TrK[i] += KDD[i][k, d] * gammaUU[i][k, d]

KijKij = np.zeros(N1)
for i in range(N1):
    for k in range(3):
        for d in range(3):
            for l in range(3):
                for m in range(3):
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
for i in range(DIM):
    for j in range(DIM):
        for k in range(DIM):
            psi2 += gammaUU[i][j, k] * uu_dD[i][j] * uu_dD[i][k]

rho = 0.5 * (vv2 + psi2) + Vofu

RicciScalar = (
    2 * (
        dgammaCartDD00dr * dgammaCartDD12dr - 2 * d2gammaCartDD12dr2 * gammaCartDD00) * gammaCartDD12**3 + gammaCartDD12**2 * (
            dgammaCartDD12dr**2 * gammaCartDD00 + 3 * dgammaCartDD11dr * dgammaCartDD22dr * gammaCartDD00 - 4 * dgammaCartDD02dr * dgammaCartDD12dr * gammaCartDD01 + (
                -4 * dgammaCartDD01dr * dgammaCartDD12dr + 8 * d2gammaCartDD12dr2 * gammaCartDD01) * gammaCartDD02 - dgammaCartDD00dr * dgammaCartDD22dr * gammaCartDD11 + 2 * d2gammaCartDD22dr2 * gammaCartDD00 * gammaCartDD11 + (
                    -(
                        dgammaCartDD00dr * dgammaCartDD11dr) + 2 * d2gammaCartDD11dr2 * gammaCartDD00) * gammaCartDD22) + gammaCartDD11**2 * (
                            dgammaCartDD22dr**2 * gammaCartDD00 - 2 * dgammaCartDD02dr * dgammaCartDD22dr * gammaCartDD02 + 2 * d2gammaCartDD22dr2 * gammaCartDD02**2 + (
                                dgammaCartDD00dr * dgammaCartDD22dr - 2 * d2gammaCartDD22dr2 * gammaCartDD00) * gammaCartDD22) + gammaCartDD11 * (
                                    -(
                                        dgammaCartDD22dr**2 * gammaCartDD01**2) + 2 * dgammaCartDD12dr * dgammaCartDD22dr * gammaCartDD01 * gammaCartDD02 + (
                                            -3 * dgammaCartDD12dr**2 + 2 * dgammaCartDD11dr * dgammaCartDD22dr) * gammaCartDD02**2 + (
                                                3 * dgammaCartDD12dr**2 * gammaCartDD00 - dgammaCartDD11dr * dgammaCartDD22dr * gammaCartDD00 - 2 * dgammaCartDD01dr * dgammaCartDD22dr * gammaCartDD01 + 2 * d2gammaCartDD22dr2 * gammaCartDD01**2 - 2 * dgammaCartDD02dr * dgammaCartDD11dr * gammaCartDD02 + 2 * d2gammaCartDD11dr2 * gammaCartDD02**2) * gammaCartDD22 + (
                                                    dgammaCartDD00dr * dgammaCartDD11dr - 2 * d2gammaCartDD11dr2 * gammaCartDD00) * gammaCartDD22**2) + gammaCartDD22 * (
                                                        gammaCartDD01**2 * (
                                                            -3 * dgammaCartDD12dr**2 + 2 * dgammaCartDD11dr * dgammaCartDD22dr + 2 * d2gammaCartDD11dr2 * gammaCartDD22) + 2 * dgammaCartDD11dr * gammaCartDD01 * (
                                                                dgammaCartDD12dr * gammaCartDD02 - dgammaCartDD01dr * gammaCartDD22) + dgammaCartDD11dr**2 * (
                                                                    -gammaCartDD02**2 + gammaCartDD00 * gammaCartDD22)) + 2 * gammaCartDD12 * (
                                                                        -2 * dgammaCartDD12dr * dgammaCartDD22dr * gammaCartDD00 * gammaCartDD11 + gammaCartDD02**2 * (
                                                                            dgammaCartDD11dr * dgammaCartDD12dr - 2 * d2gammaCartDD12dr2 * gammaCartDD11) - 2 * dgammaCartDD11dr * dgammaCartDD12dr * gammaCartDD00 * gammaCartDD22 - dgammaCartDD00dr * dgammaCartDD12dr * gammaCartDD11 * gammaCartDD22 + 2 * d2gammaCartDD12dr2 * gammaCartDD00 * gammaCartDD11 * gammaCartDD22 + gammaCartDD01**2 * (
                                                                                dgammaCartDD12dr * dgammaCartDD22dr - 2 * d2gammaCartDD12dr2 * gammaCartDD22) + gammaCartDD01 * (
                                                                                    dgammaCartDD02dr * dgammaCartDD22dr * gammaCartDD11 + (
                                                                                        dgammaCartDD02dr * dgammaCartDD11dr + 2 * dgammaCartDD01dr * dgammaCartDD12dr) * gammaCartDD22) + gammaCartDD02 * (
                                                                                            (dgammaCartDD12dr**2 - 3 * dgammaCartDD11dr * dgammaCartDD22dr) * gammaCartDD01 + (
                                                                                                2 * dgammaCartDD02dr * dgammaCartDD12dr + dgammaCartDD01dr * dgammaCartDD22dr - 2 * d2gammaCartDD22dr2 * gammaCartDD01) * gammaCartDD11 + (
                                                                                                    dgammaCartDD01dr * dgammaCartDD11dr - 2 * d2gammaCartDD11dr2 * gammaCartDD01) * gammaCartDD22))) / (
                                                                                                        2. * (
                                                                                                            gammaCartDD02**2 * gammaCartDD11 - 2 * gammaCartDD01 * gammaCartDD02 * gammaCartDD12 + gammaCartDD01**2 * gammaCartDD22 + gammaCartDD00 * (
                                                                                                                gammaCartDD12**2 - gammaCartDD11 * gammaCartDD22))**2)



HamConstraint = RicciScalar + KijKij - TrK**2 - 16*np.pi*rho

plt.plot(HamConstraint)
plt.show()
