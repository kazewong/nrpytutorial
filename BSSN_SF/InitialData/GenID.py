import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d

phi0 = 0.0006
r0 = 5
delta = 1
q = 2

def phi(r,phi0,r0,delta,q):
  phi  = phi0*r**3*np.exp(-((r-r0)/delta)**q)
  return phi

def dphidr(r,phi0,r0,delta,q):
	a = 3*phi0*np.exp(-((r-r0)/delta)**q)*r**2
	b = (phi0*np.exp(-((r-r0)/delta)**q)*q*r**3*((r-r0)/delta)**(q-1))/delta
	return a-b


r_axis = np.arange(1e-10,100003*0.01,0.01)
#dphi_interp = interp1d(r_axis,np.gradient(phi(r_axis,phi0,r0,delta,q),r_axis),bounds_error=False,fill_value='extrapolate')

def RHS(r,y):
	return [y[1],-np.pi*y[0]*dphidr(r,phi0,r0,delta,q)**2-2*y[1]/r]

y = solve_ivp(RHS,[1e-10,1001],np.array([1.,0.]),t_eval=r_axis).y

tag = 'CC1'

np.savetxt('alpha'+tag+'.csv',np.ones(r_axis.size))
np.savetxt('phi'+tag+'.csv',phi(r_axis,phi0,r0,delta,q))
np.savetxt('psi'+tag+'.csv',y[0])
np.savetxt('Pi'+tag+'.csv',np.zeros(r_axis.size))
