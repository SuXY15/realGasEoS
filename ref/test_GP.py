from utils import *
from copy import deepcopy

## settings for C12
fluid = "C12"
name = 'c12h26'
T_step = 20
D_step = 40
T_lo, T_hi = 500, 800
P_arr = np.array([1.817, 3, 5]) * 1e6

## settings for O2
#fluid = "oxygen"
#name = "o2"
#T_step = 20
#D_step = 40
#T_lo, T_hi = 60, 400
#P_arr = np.array([5]) * 1e6

# load data
Data = np.loadtxt("mech/Alpha/%s.csv"%name, delimiter=',')
X = Data[:,:2]
y = Data[:,2]
N = len(y)
dim = 2

x1lim = [np.min(X[:,0]), np.max(X[:,0])]
x2lim = [np.min(X[:,1]), np.max(X[:,1])]

## define covariance
def k0(X1, X2, kernel_size=[20,20]):
    cov = np.zeros((len(X1), len(X2)))
    for i in range(len(X1)):
        for j in range(len(X2)):
            cov[i,j] = np.exp(-np.sum((X1[i] - X2[j])**2/kernel_size))
    return cov

## GP predict

# random queries
M = 100
Xnew = np.random.rand(M, dim)
Xnew[:,0] = Xnew[:,0] * (x1lim[1] - x1lim[0]) + x1lim[0]
Xnew[:,1] = Xnew[:,1] * (x2lim[1] - x2lim[0]) + x2lim[0]
Xnew = Xnew[np.argsort(Xnew[:,1]),:]
Xnew = Xnew[np.argsort(Xnew[:,0]),:]

## original queries
#Pc = CP.PropsSI(fluid, 'pcrit')
#Tc = CP.PropsSI(fluid, 'Tcrit')
#TPD_arr = []
#for P in P_arr:
    #TPD_arr += get_TPD_under_P(fluid, P, T_lo, T_hi, T_step, D_step)
#TPD_arr = np.array(TPD_arr)

#Xnew = np.zeros((len(TPD_arr), 2))
#Xnew[:,0] = TPD_arr[:,0] / Tc + np.random.rand(len(TPD_arr[:,1]))*0.05
#Xnew[:,1] = TPD_arr[:,1] / Pc + np.random.rand(len(TPD_arr[:,1]))*0.05

Kxx = k0(Xnew, Xnew)
Kx = k0(X, Xnew)
K  = k0(X, X)
Ki = np.linalg.inv(K + 1e-4*np.diag(np.ones(N)))

mu = Kx.T @ Ki @ y
cov = Kxx - Kx.T @ Ki @ Kx
sigma = np.diag(cov)

plt.plot(X[:,0], y, 'rv')
plt.plot(Xnew[:,0], mu, 'gp')
plt.plot(Xnew[:,0], mu+3*sigma, 'g.', ms=0.5)
plt.plot(Xnew[:,0], mu-3*sigma, 'g.', ms=0.5)

plt.figure()
plt.plot(X[:,0]+X[:,1], y, 'rv')
plt.plot(Xnew[:,0]+Xnew[:,1], mu, 'gp')

plt.legend(["Groundtruth ", "Prediction"])
plt.show()
