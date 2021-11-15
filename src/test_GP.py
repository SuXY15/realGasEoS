from utils import *
from copy import deepcopy

## settings for C12
fluid = "C12"
name = 'c12h26'

# # settings for O2
# fluid = "oxygen"
# name = "o2"

# ================================================
# load data
data = np.loadtxt("mech/Alpha/%s.csv"%name, delimiter=',')
dim = data.shape[1]-1
X = data[:,:dim]
y = data[:,dim]
N = len(y)

para = np.loadtxt('mech/Alpha/%s_para.csv'%name, delimiter=',')
𝛾 = para[0,:dim] # kernel size
σ = para[0,dim]  # kernel multiplier
θ = para[1,:]    # basis function's parameters

# ================================================
# define covariance
def k0(X1, X2, 𝛾=𝛾, σ=σ):
    cov = np.zeros((len(X1), len(X2)))
    for i in range(len(X1)):
        for j in range(len(X2)):
            cov[i,j] = σ**2 * np.exp(-np.sum((X1[i] - X2[j])**2/ 2 / 𝛾**2))
    return cov

# ================================================
# random queries
M = 100
x1lim = [np.min(X[:,0]), np.max(X[:,0])]
x2lim = [np.min(X[:,1]), np.max(X[:,1])]
Xnew = np.random.rand(M, dim)
Xnew[:,0] = Xnew[:,0] * (x1lim[1] - x1lim[0]) + x1lim[0]
Xnew[:,1] = Xnew[:,1] * (x2lim[1] - x2lim[0]) + x2lim[0]
Xnew = Xnew[np.argsort(Xnew[:,1]),:]
Xnew = Xnew[np.argsort(Xnew[:,0]),:]

# ================================================
# GP predict
Kxx = k0(Xnew, Xnew)
Kx = k0(X, Xnew)
K  = k0(X, X)
Ki = np.linalg.inv(K + 1e-4*np.diag(np.ones(N)))

mu = Kx.T @ Ki @ y
cov = Kxx - Kx.T @ Ki @ Kx
sigma = np.diag(cov)

# ================================================
# GP gradient
m = Ki @ y
Tc = CP.PropsSI(fluid, 'Tcrit')
dalphadT = [(-1/Tc * (Xnew[i,0] - X[:,0]) / 𝛾[0]**2 * Kx[:,i]) @ m for i in range(M)]
d2alphadT2 = [(((Xnew[i,0] - X[:,0])**2/𝛾[0]**2 - 1)/Tc**2/𝛾[0]**2 * Kx[:,i]) @ m for i in range(M)]

# ================================================
# PR gradient
AS = CP.AbstractState("HEOS", fluid)
Tc = CP.PropsSI(fluid, 'Tcrit')
Pc = CP.PropsSI(fluid, 'pcrit')
omega = AS.acentric_factor()
dalphadT_PR = [PR_dalphadT(Xnew[i,0]*Tc, Xnew[i,1]*Pc, Tc, Pc, omega) for i in range(M)]
d2alphadT2_PR = [PR_d2alphadT2(Xnew[i,0]*Tc, Xnew[i,1]*Pc, Tc, Pc, omega) for i in range(M)]

# 3d plot of alpha
fig, ax = plt.subplots(subplot_kw={"projection":"3d"})
ax.set_title("Alpha")
ax.scatter(X[:,0], X[:,1], y, 'r')
ax.scatter(Xnew[:,0], Xnew[:,1], mu, 'g')
ax.set_xlabel("Tr")
ax.set_ylabel("Pr")
plt.legend(["Groundtruth ", "Prediction"])
fig.savefig("figs/PythonAlphaGP_alpha_%s.png"%fluid)

# 3d plot of dalphadT
fig, ax = plt.subplots(subplot_kw={"projection":"3d"})
ax.set_title("dalphadT")
ax.scatter(Xnew[:,0], Xnew[:,1], dalphadT_PR, 'r')
ax.scatter(Xnew[:,0], Xnew[:,1], dalphadT, 'g')
ax.set_xlabel("Tr")
ax.set_ylabel("Pr")
plt.legend(["PR", "AlphaGP"])
fig.savefig("figs/PythonAlphaGP_dalphadT_%s.png"%fluid)

# 3d plot of d2alphadT2
fig, ax = plt.subplots(subplot_kw={"projection":"3d"})
ax.set_title("d2alphadT2")
ax.scatter(Xnew[:,0], Xnew[:,1], d2alphadT2_PR, 'r')
ax.scatter(Xnew[:,0], Xnew[:,1], d2alphadT2, 'g')
ax.set_xlabel("Tr")
ax.set_ylabel("Pr")
plt.legend(["PR", "AlphaGP"])
fig.savefig("figs/PythonAlphaGP_d2alphadT2_%s.png"%fluid)

plt.show()
