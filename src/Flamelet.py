from utils import *

phi_arr = np.linspace(0.5, 2.0, 31)

## compute for real gas
real_gas = ct.Solution('mech/nDodecane_Reitz.yaml', 'nDodecane_RK')

real_T = np.zeros(phi_arr.shape)
for i,phi in enumerate(phi_arr):
    print("Runing for Real Gas: phi=%.3f"%phi)
    real_gas.TP = 300, 50*ct.one_atm
    real_gas.set_equivalence_ratio(phi, 'c12h26', 'o2:1,n2:3.76')
    real_gas.equilibrate('HP')
    real_T[i] = real_gas.T

    
## compute for ideal gas
ideal_gas = ct.Solution('mech/nDodecane_Reitz.yaml', 'nDodecane_IG')

ideal_T = np.zeros(phi_arr.shape)
for i,phi in enumerate(phi_arr):
    print("Runing for Ideal Gas: phi=%.3f"%phi)
    ideal_gas.TP = 300, 50*ct.one_atm
    ideal_gas.set_equivalence_ratio(phi, 'c12h26', 'o2:1,n2:3.76')
    ideal_gas.equilibrate('HP')
    ideal_T[i] = ideal_gas.T  

plt.plot(phi_arr, real_T, label='Real Gas RK', lw=2)
plt.plot(phi_arr, ideal_T, label='Ideal Gas', lw=2)
plt.xlabel('Equivalence ratio, $\phi$')
plt.ylabel('Temperature [K]');
plt.legend()
plt.show()
