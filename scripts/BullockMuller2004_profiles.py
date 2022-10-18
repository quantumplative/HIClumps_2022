import pynbody 
import numpy as np 
import matplotlib.pyplot as plt 

Lz = 1 #LambdaZ corresponds to 0.1 solar 
mp = 1.67262e-24 #g
Cv = 13 # concentration used in Muller&Bullock 2004 for Mvir ~ 10^12 Msol halos 
tf = 8 # Gyr formation time of halo with Mvir ~ 10^12 Msol
t8 = tf/8
z  = 0 #MB04
DeltaV = 360/(1+z)
Omegam = 0.3 
Mvir   = 1e12 #Msol 
mue    = 0.62 # Jess's value 
#mue = 1.18  #MB04
nH_fac = 1.16 #Werk 2014  
Tcool = 2e4 #K
Thalo = 1.3*1e6 #K

def Rvir(Mvir, z): 
	return 206 * (DeltaV*Omegam/97.2)**(-1/3) * (Mvir/1e12)**(1/3) * (1 + z)**(-1) #kpc

def Rs(Mvir, z): 
	return Rvir(Mvir, z)/Cv #kpc

def Rc(Mvir, z): 
	return 157*T_c(Mvir, z)**(-1/8)*(Lz*t8)**(1/3) #kpc 

def Vvir(Mvir, z): 
	return 144 * (DeltaV*Omegam/97.2)**(1/6)  * (Mvir/1e12)**(1/3) * (1 + z)**(1/2) * 1e5 #cm/s

def Vmax(Mvir, z):
	return (Vvir(Mvir, z)/1e5/0.468)**(1/1.1) #km/s equation 7 of MB04

def n_c(Mvir, z): 
	return 6.1*1e-5*T_c(Mvir, z)**2*(Lz*t8)**(-1) #cm^-3

def T6(T):  
	return T/1e6 #unitless 

def T_c(Mvir, z): 
	#return T6(1e6*(Vvir(Mvir, z)/1e5/163)**2) # divide by 1e5 bc need Vvir in km/s
	return T6(1e6*(Vmax(Mvir, z)/163)**2)

def rho_c(Mvir, z): 
	return mue*mp*n_c(Mvir, z) #g/cm^3

def P_c(Mvir, z): 
	#return (Vvir(Mvir, z)**2/2)*rho_c(Mvir, z) #g s^-2 cm^-1
	return ((Vmax(Mvir, z)*1e5)**2/2)*rho_c(Mvir, z) #g s^-2 cm^-1

def adiabatic_function_form(x, Mvir, z):
	Cc = Rc(Mvir, z)/Rs(Mvir, z)
	return (1 + (3.7/x)*np.log(1+x) - (3.7/Cc)*np.log(1+Cc))

def adiabatic_density_profile(r, Mvir, z):
	x = r/Rs(Mvir, z)
	return rho_c(Mvir, z)*(adiabatic_function_form(x, Mvir, z))**(3/2)

def adiabatic_pressure_profile(r, Mvir, z): 
	x = r/Rs(Mvir, z) 
	return P_c(Mvir, z)*(adiabatic_function_form(x, Mvir, z))**(5/2) #g s^-2 cm^-1 

def werk2014_fig12():
	print(r'T$_c$ = ', T_c(Mvir, z))
	rvir   = Rvir(Mvir, z)
	rbins  = np.linspace(1, rvir, 100)	
	rho = adiabatic_density_profile(rbins, 1e12, 0)
	ne  = rho / (mp*mue)
	ne_cold = ne *(Thalo/Tcool) *2.7 #extra factor of 2.7 for highly ionized (Jess says)
	fig = plt.figure(figsize=(12,10))
	plt.plot(rbins/rvir, nH_fac*ne, color='maroon', lw=2, label='Hot')
	plt.plot(rbins/rvir, nH_fac*ne_cold, color='b', lw=2, label='Cold') 
	plt.loglog()
	plt.xlabel(r'r/r$_{\rm vir}$', fontsize=16)
	plt.ylabel(r'n$_{\rm e}$ [cm$^{-3}$]', fontsize=16)
	plt.tick_params(top=True, bottom=True, left=True, right=True, direction='in',labelsize=16, which='both')
	plt.legend(fontsize=14)
	plt.ylim(1e-5, 1e0)
	plt.xlim(1e-2, 1)
	plt.savefig('recreate_werk2014_fig12.pdf')
	plt.show()	
