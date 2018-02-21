#nodes / dim
nd=1

#box
L=0.03*nd

#particle number
NPBH=1e2#1*nd**3
NCDM=2*(192*nd)**3

eps=1.0*NPBH/NCDM

fpbh=1e-3#1e-5 
omega_c=0.26*(1.-fpbh)
omega_b=0.05
omega_p=0.26*fpbh
omega_m=omega_c+omega_b+omega_p

M=2.775*10.**11.*omega_m*L**3. #Total mass in box

print("NCDM: ",NCDM)
print("NPBH: ",NPBH)
print("fpbh: ",fpbh)

print("")
print("Box Mass (Msun/h):          ",M)
print("CDM Particle Mass (Msun/h): ",M*omega_c/omega_m/NCDM)
print("PBH Particle Mass (Msun/h): ",M*omega_p/omega_m/NPBH)

print("")
print("m/f^2/M:    ",omega_p/omega_m/NPBH/fpbh**2)


#print("Total box (Mpc/h): ",L)
#print("Min k (h/Mpc): ",2*3.1415926535/L)
#print("Max k (h/Mpc): ",(NCDM/2)**(1./3.)*2*3.1415926535/L)
