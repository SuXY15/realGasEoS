from utils import *
from copy import deepcopy

mechs = ['mech/nDodecane_AlphaGP.yaml', 'mech/nDodecane_temp.yaml', 'mech/nDodecane_Reitz.yaml', 'mech/nDodecane_Reitz.yaml']
names = ['nDodecane_PR_ALphaGP', 'nDodecane_PR', 'nDodecane_RK', 'nDodecane_IG']
lines = ['o', 'd', 'v', '.']
colors = ['r', 'b', 'g', 'c']

## settings for C12
#fluid = "C12"
#X = {'c12h26':1}
#T_step = 10
#D_step = 20
#T_lo, T_hi = 500, 800
#P_arr = np.array([1.817, 0.4]) * 1e6

## settings for o2
fluid = "oxygen"
X = {'o2':1}
T_step = 20
D_step = 40
T_lo, T_hi = 60, 400
P_arr = np.array([5]) * 1e6

## get adaptive TP list and NIST data
TPD_arr = []
for P in P_arr:
    TPD_arr += get_TPD_under_P(fluid, P, T_lo, T_hi, T_step, D_step)
TPD_arr = np.array(TPD_arr)
plt.plot(TPD_arr[:,0], TPD_arr[:,2], 'ks', label="NIST", fillstyle='none')

for k,name in enumerate(names):
    gas = ct.Solution(mechs[k], name)
    gas.X = X
    TPD_calc = deepcopy(TPD_arr)
    
    t0 = time.time()
    for i,(T,P,_) in enumerate(TPD_calc):
        gas.TP = T,P
        TPD_calc[i,2] = gas.density
    
    print("cost of", name , time.time() - t0)
    
    plt.plot(TPD_calc[:,0], TPD_calc[:,2], colors[k]+lines[k], label=name, alpha=0.8, fillstyle='none')
    
    
    


plt.xlabel("Temperature [K]")
plt.ylabel("Density [kg/m^3]")
plt.legend()
plt.show()
