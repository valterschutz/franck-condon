from scipy.constants import h,c,eV

nu_e = 125.69*100
x_ee__nu_e = 0.764*100
x_ee = x_ee__nu_e / nu_e
D_ee = nu_e/(4*x_ee)
# print(f"D_ee: {D_ee}")
def E_vib_e(v):
    return h*c*((v+1/2)*nu_e-(v+1/2)**2*x_ee__nu_e)

nu_g = 214.50*100
x_eg__nu_g = 0.614*100
x_eg = x_eg__nu_g / nu_g
D_eg = nu_g/(4*x_eg)
# print(f"D_eg: {D_eg}")
def E_vib_g(v):
    return h*c*((v+1/2)*nu_g-(v+1/2)**2*x_eg__nu_g)

E_max_e = h*c*D_ee
E_max_g = h*c*D_eg
print(f"E_max: {E_max_g/eV} eV")
