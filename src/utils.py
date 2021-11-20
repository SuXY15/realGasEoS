import os, sys, time
import numpy as np
import cantera as ct
import CoolProp.CoolProp as CP
import matplotlib.pyplot as plt
print('Running Cantera version: ' + ct.__version__)

plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['figure.autolayout'] = True
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['font.family'] = 'serif'

R = 8.31446
def PR(T_c, P_c, R_u=R):
    a = 0.45724 * R_u**2 * T_c**2 / P_c
    b = 0.07780 * R_u * T_c / P_c
    return a, b

def PR_alpha(T, P, T_c, P_c, omega):
    a = 0.45724 * R**2 * T_c**2 / P_c
    b = 0.07780 * R * T_c / P_c
    kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2
    alpha = ( 1 + kappa * (1-np.sqrt(T/T_c)) )**2
    return alpha

def PR_dalphadT(T, P, T_c, P_c, omega):
    a = 0.45724 * R**2 * T_c**2 / P_c
    b = 0.07780 * R * T_c / P_c
    kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2
    dalphadT = kappa*kappa/T_c - (kappa*kappa+kappa)/np.sqrt(T_c) * T**(-1/2)
    return dalphadT

def PR_d2alphadT2(T, P, T_c, P_c, omega):
    a = 0.45724 * R**2 * T_c**2 / P_c
    b = 0.07780 * R * T_c / P_c
    kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2
    d2alphadT2 = 0.5*(kappa*kappa+kappa)/np.sqrt(T_c) * T**(-3/2)
    return d2alphadT2

def get_TPD_under_P(fluid,P, T_lo, T_hi, T_step=20, D_step=40):
    TPD_arr = []
    T = T_lo
    D = CP.PropsSI("D", "T", T, "P", P, fluid)
    T_old = T
    D_old = D
    count = 0
    alpha = 0.50
    while T <= T_hi:
        D = CP.PropsSI("D", "T", T, "P", P, fluid)
        if abs(D_old - D) < D_step or count > 5:
            TPD_arr.append([T, P, D])
            T_old = T
            D_old = D
            T += T_step
            count = 0
        else:
            T = T - alpha*(T - T_old)
            count += 1
    return TPD_arr
