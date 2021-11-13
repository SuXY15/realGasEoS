from utils import *

# settings for C12
fluid = "C12"
name = "c12h26"
M = 170.33 / 1000
Pcrit = CP.PropsSI(fluid, 'pcrit')
Tcrit = CP.PropsSI(fluid, 'Tcrit')
a, b = PR(Tcrit, Pcrit)
T_step = 20
D_step = 40
T_lo, T_hi = 500, 800
P_arr = np.array([0.5, 1, 1.817, 3, 5, 10, 20, 30, 40, 50]) * 1e6
# P_arr = np.linspace(0.5, 50, 20) * 1e6


# settings for O2
#fluid = "oxygen"
#name = "o2"
#M = 32 / 1000
#Pcrit = CP.PropsSI(fluid, 'pcrit')
#Tcrit = CP.PropsSI(fluid, 'Tcrit')
#a, b = PR(Tcrit, Pcrit)
#T_step = 20
#D_step = 40
#T_lo, T_hi = 60, 400
#P_arr = np.array([1, 3, 5]) * 1e6

# load omega
AS = CP.AbstractState("HEOS", fluid)
omega = AS.acentric_factor()

print("Critical Properties: \r\nTc: %f \r\nPc: %f"%(Tcrit, Pcrit))
print("YAML \r\na: %.4e \r\nb: %.4e"%PR(Tcrit, Pcrit, R_u=R*1e6))
print("acentric-factor:", omega)
TPD_arr = []
for P in P_arr:
    TPD_arr += get_TPD_under_P(fluid, P, T_lo, T_hi, T_step, D_step)

TPD_arr = np.array(TPD_arr)

T = TPD_arr[:,0]
P = TPD_arr[:,1]
D = TPD_arr[:,2]
V = M / D
Alpha = (R*T/(V-b) - P) / a * (V*(V+b) + b*(V-b))
PR_Alpha  = PR_alpha(T, P, Tcrit, Pcrit, omega)

# ===============
# # show results
plt.figure()
plt.plot(T, D, 's', label="NIST Desity", alpha=0.5)

# plt.figure()
# plt.plot(T, V, 'p', label="NIST Desity", alpha=0.5)

plt.figure()
plt.plot(T, Alpha, 's', label="NIST alpha", alpha=0.5)
plt.plot(T, PR_Alpha, 'o', label="NIST alpha", alpha=0.5)

# ===============
# # save results
Data = np.zeros(TPD_arr.shape)
Data[:,0] = T / Tcrit
Data[:,1] = P / Pcrit
Data[:,2] = Alpha
np.savetxt("mech/Alpha/%s.csv"%name, Data, delimiter=', ')

plt.show()
