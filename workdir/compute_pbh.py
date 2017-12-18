#nodes / dim
nd=1

#box
L=0.1*nd

#particle number
NPBH=10000*nd**3
NCDM=2*(2*128*nd)**3

fpbh=0.001
omega_c=0.27*(1.-fpbh)
omega_b=0.05
omega_p=0.27*fpbh
omega_m=omega_c+omega_b+omega_p

M=2.775*10.**11.*omega_m*L**3. #Total mass in box

print("CDM Particle Mass (Msun/h): ",M*omega_c/omega_m/NCDM)
print("PBH Particle Mass (Msun/h): ",M*omega_p/omega_m/NPBH)

print("Total box (Mpc/h): ",L)

print("Min k (h/Mpc): ",2*3.1415926535/L)
print("Max k (h/Mpc): ",(NCDM/2)**(1./3.)*2*3.1415926535/L)
