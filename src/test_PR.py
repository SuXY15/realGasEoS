from utils import *
from copy import deepcopy

mechs = ['mech/nDodecane_AlphaGP.yaml', 'mech/nDodecane_temp.yaml', 'mech/nDodecane_Reitz.yaml', 'mech/nDodecane_Reitz.yaml']
names = ['nDodecane_PR_ALphaGP', 'nDodecane_PR', 'nDodecane_RK', 'nDodecane_IG']
lines = ['o', 'd', 'v', '.']
colors = ['r', 'b', 'g', 'c']

# settings for C12
fluid = "C12"
X = {'c12h26':1}
T_step = 10
D_step = 20
T_lo, T_hi = 500, 800
P_arr = np.array([1.817]) * 1e6

# # settings for o2
# fluid = "oxygen"
# X = {'o2':1}
# T_step = 20
# D_step = 40
# T_lo, T_hi = 60, 400
# P_arr = np.array([5]) * 1e6

## get adaptive TP list and NIST data
TPD_arr = []
for P in P_arr:
    TPD_arr += get_TPD_under_P(fluid, P, T_lo, T_hi, T_step, D_step)
TPD_arr = np.array(TPD_arr)

## =========================
## Density
plt.figure()
plt.plot(TPD_arr[:,0], TPD_arr[:,2], 'ks', label="NIST", fillstyle='none')

for k,name in enumerate(names):
    gas = ct.Solution(mechs[k], name)
    TPD_calc = deepcopy(TPD_arr)
    
    t0 = time.time()
    for i,(T,P,_) in enumerate(TPD_calc):
        gas.TPX = T,P, X
        TPD_calc[i,2] = gas.density
    print("Density cost of %-20s = %.5f s"%(name, time.time()-t0))
    plt.plot(TPD_calc[:,0], TPD_calc[:,2], colors[k]+lines[k], label=name, alpha=0.8, fillstyle='none')
plt.xlabel("Temperature [K]")
plt.ylabel("Density [kg/m^3]")
plt.legend()

## =========================
## Cp_mass
TPV_arr = deepcopy(TPD_arr)
for i,(T,P,_) in enumerate(TPD_calc):
    TPV_arr[i,2] = CP.PropsSI("C", "T", T, "P", P, fluid)

plt.figure()
plt.plot(TPV_arr[:,0], TPV_arr[:,2], 'ks', label="NIST", fillstyle='none')

for k,name in enumerate(names):
    gas = ct.Solution(mechs[k], name)
    TPV_calc = deepcopy(TPV_arr)
    
    t0 = time.time()
    for i,(T,P,_) in enumerate(TPD_calc):
        gas.TPX = T, P, X
        TPV_calc[i,2] = gas.cp_mass
    print("Cp_mass cost of %-20s = %.5f s"%(name, time.time()-t0))
    plt.plot(TPV_calc[:,0], TPV_calc[:,2], colors[k]+lines[k], label=name, alpha=0.8, fillstyle='none')
    # print("Cp_mass", TPV_calc[:,2])
    
plt.ylim([0,6000])
plt.xlabel("Temperature [K]")
plt.ylabel("Cp_mass [J/kg/K]")
plt.legend()

plt.show()
