import numpy as np
import matplotlib.pyplot as plt 

data = np.genfromtxt("BSSN_SF-output2D/quad_pot_2d_t-00000050.txt")

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


R = (2*(dgammaCartDD00dr*dgammaCartDD12dr - 2*d2gammaCartDD12dr2*gammaCartDD00)*gammaCartDD12**3 + gammaCartDD12**2*(dgammaCartDD12dr**2*gammaCartDD00 + 3*dgammaCartDD11dr*dgammaCartDD22dr*gammaCartDD00 - 4*dgammaCartDD02dr*dgammaCartDD12dr*gammaCartDD01 + (-4*dgammaCartDD01dr*dgammaCartDD12dr + 8*d2gammaCartDD12dr2*gammaCartDD01)*gammaCartDD02 - dgammaCartDD00dr*dgammaCartDD22dr*gammaCartDD11 + 2*d2gammaCartDD22dr2*gammaCartDD00*gammaCartDD11 + (-(dgammaCartDD00dr*dgammaCartDD11dr) + 2*d2gammaCartDD11dr2*gammaCartDD00)*gammaCartDD22) + gammaCartDD11**2*(dgammaCartDD22dr**2*gammaCartDD00 - 2*dgammaCartDD02dr*dgammaCartDD22dr*gammaCartDD02 + 2*d2gammaCartDD22dr2*gammaCartDD02**2 + (dgammaCartDD00dr*dgammaCartDD22dr - 2*d2gammaCartDD22dr2*gammaCartDD00)*gammaCartDD22) + gammaCartDD11*(-(dgammaCartDD22dr**2*gammaCartDD01**2) + 2*dgammaCartDD12dr*dgammaCartDD22dr*gammaCartDD01*gammaCartDD02 + (-3*dgammaCartDD12dr**2 + 2*dgammaCartDD11dr*dgammaCartDD22dr)*gammaCartDD02**2 + (3*dgammaCartDD12dr**2*gammaCartDD00 - dgammaCartDD11dr*dgammaCartDD22dr*gammaCartDD00 - 2*dgammaCartDD01dr*dgammaCartDD22dr*gammaCartDD01 + 2*d2gammaCartDD22dr2*gammaCartDD01**2 - 2*dgammaCartDD02dr*dgammaCartDD11dr*gammaCartDD02 + 2*d2gammaCartDD11dr2*gammaCartDD02**2)*gammaCartDD22 + (dgammaCartDD00dr*dgammaCartDD11dr - 2*d2gammaCartDD11dr2*gammaCartDD00)*gammaCartDD22**2) + gammaCartDD22*(gammaCartDD01**2*(-3*dgammaCartDD12dr**2 + 2*dgammaCartDD11dr*dgammaCartDD22dr + 2*d2gammaCartDD11dr2*gammaCartDD22) + 2*dgammaCartDD11dr*gammaCartDD01*(dgammaCartDD12dr*gammaCartDD02 - dgammaCartDD01dr*gammaCartDD22) + dgammaCartDD11dr**2*(-gammaCartDD02**2 + gammaCartDD00*gammaCartDD22)) + 2*gammaCartDD12*(-2*dgammaCartDD12dr*dgammaCartDD22dr*gammaCartDD00*gammaCartDD11 + gammaCartDD02**2*(dgammaCartDD11dr*dgammaCartDD12dr - 2*d2gammaCartDD12dr2*gammaCartDD11) - 2*dgammaCartDD11dr*dgammaCartDD12dr*gammaCartDD00*gammaCartDD22 - dgammaCartDD00dr*dgammaCartDD12dr*gammaCartDD11*gammaCartDD22 + 2*d2gammaCartDD12dr2*gammaCartDD00*gammaCartDD11*gammaCartDD22 + gammaCartDD01**2*(dgammaCartDD12dr*dgammaCartDD22dr - 2*d2gammaCartDD12dr2*gammaCartDD22) + gammaCartDD01*(dgammaCartDD02dr*dgammaCartDD22dr*gammaCartDD11 + (dgammaCartDD02dr*dgammaCartDD11dr + 2*dgammaCartDD01dr*dgammaCartDD12dr)*gammaCartDD22) + gammaCartDD02*((dgammaCartDD12dr**2 - 3*dgammaCartDD11dr*dgammaCartDD22dr)*gammaCartDD01 + (2*dgammaCartDD02dr*dgammaCartDD12dr + dgammaCartDD01dr*dgammaCartDD22dr - 2*d2gammaCartDD22dr2*gammaCartDD01)*gammaCartDD11 + (dgammaCartDD01dr*dgammaCartDD11dr - 2*d2gammaCartDD11dr2*gammaCartDD01)*gammaCartDD22)))/(2.*(gammaCartDD02**2*gammaCartDD11 - 2*gammaCartDD01*gammaCartDD02*gammaCartDD12 + gammaCartDD01**2*gammaCartDD22 + gammaCartDD00*(gammaCartDD12**2 - gammaCartDD11*gammaCartDD22))**2)


gammaDD = np.array([[gammaCartDD00,gammaCartDD01,gammaCartDD02],
                    [gammaCartDD01,gammaCartDD11,gammaCartDD12],
                    [gammaCartDD02,gammaCartDD12,gammaCartDD22]]

gammaUU = np.linalg.inv(gammaDD)



