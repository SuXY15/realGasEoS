from utils import *
from copy import deepcopy
np.random.seed(0)

## settings for C12
fluid = "C12"
name = 'c12h26'

# # settings for O2
#fluid = "oxygen"
#name = "o2"

# ================================================
# load data
data = np.loadtxt("mech/Alpha/%s.csv"%name, delimiter=',')
dim = data.shape[1]-1
X = data[:,:dim]
y = data[:,dim]
N = len(y)
print(N)

Ntrain = int(N*0.7)
Ntest = N - Ntrain
idx = np.random.permutation(N)
Xall = deepcopy(X)
yall = deepcopy(y)
X = Xall[idx[:Ntrain], :]
y = yall[idx[:Ntrain]]
Xtest = Xall[idx[Ntrain:], :]
ytest = yall[idx[Ntrain:]]

para = np.loadtxt('mech/Alpha/%s_para.csv'%name, delimiter=',')
𝛾 = para[0,:dim]*5 # kernel size
σ = para[0,dim]  # kernel multiplier
θ = para[1,:]    # basis function's parameters

AS = CP.AbstractState("HEOS", fluid)
Tc = CP.PropsSI(fluid, 'Tcrit')
Pc = CP.PropsSI(fluid, 'pcrit')
omega = AS.acentric_factor()

print(𝛾, σ)
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
M = 200

deltaTr = 0.5/1000
Xnew = np.empty([M,dim])
Xnew_l = np.empty([M,dim]) # left
Xnew_r = np.empty([M,dim]) # right 

Xnew[:,0] = np.linspace(np.min(X[:,0]), np.max(X[:,0]),M)
Xnew[:,1] = np.ones(M) * 3

Xnew_l[:,0] = Xnew[:,0]-deltaTr
Xnew_r[:,0] = Xnew[:,0]+deltaTr

Xnew_l[:,1] = Xnew[:,1]
Xnew_r[:,1] = Xnew[:,1]

#x1lim = [np.min(X[:,0]), np.max(X[:,0])]
#x2lim = [np.min(X[:,1]), np.max(X[:,1])]
#Xnew = np.random.rand(M, dim)
#Xnew[:,0] = Xnew[:,0] * (x1lim[1] - x1lim[0]) + x1lim[0]
#Xnew[:,1] = Xnew[:,1] * (x2lim[1] - x2lim[0]) + x2lim[0]
#Xnew = Xnew[np.argsort(Xnew[:,1]),:]
#Xnew = Xnew[np.argsort(Xnew[:,0]),:]

# ================================================
# GP predict
Kxx = k0(Xnew, Xnew)
Kx = k0(X, Xnew)
K  = k0(X, X)
Ki = np.linalg.inv(K + 1e-8*np.diag(np.ones(len(X))))

mu = Kx.T @ Ki @ y
cov = Kxx - Kx.T @ Ki @ Kx
sigma = np.diag(cov)

# GP predict right
Kxx1 = k0(Xnew_r, Xnew_r)
Kx1 = k0(X, Xnew_r)
mu1 = Kx1.T @ Ki @ y
cov1 = Kxx1 - Kx1.T @ Ki @ Kx1
sigma1 = np.diag(cov1)

# GP predict left
Kxx2 = k0(Xnew_l, Xnew_l)
Kx2 = k0(X, Xnew_l)
mu2 = Kx2.T @ Ki @ y 
cov2 = Kxx2 - Kx2.T @ Ki @ Kx2
sigma2 = np.diag(cov2)

# ================================================
# GP gradient, analytical
m = Ki @ y
Tc = CP.PropsSI(fluid, 'Tcrit')
dalphadT = [(-1/Tc * (Xnew[i,0] - X[:,0]) / 𝛾[0]**2 * Kx[:,i]) @ m for i in range(M)]
#print(dalphadT)
d2alphadT2 = [(((Xnew[i,0] - X[:,0])**2/𝛾[0]**2 - 1)/Tc**2/𝛾[0]**2 * Kx[:,i]) @ m for i in range(M)]

# ================================================
# GP gradient, numerical
dalphadTnu = [(mu1[i]-mu2[i])/2/deltaTr/Tc for i in range(M)]
d2alphadT2nu = [(mu1[i]-2*mu[i]+mu2[i])/deltaTr**2/Tc**2 for i in range(M)]
#print(dalphadT)
#print(d2alphadT2)

# ================================================
# PR gradient
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

# 2d plot of dalphadT
fig=plt.figure()
plt.plot(Xnew[:,0],dalphadT,'gp')
plt.plot(Xnew[:,0],dalphadTnu,'r--')
plt.plot(Xnew[:,0],dalphadT_PR,'k--')
plt.xlabel("Tr")
plt.ylabel("dalphadT")
plt.legend(["GP", "GPnumerical","PR"])

# 2d plot of da2lphadT2
fig=plt.figure()
plt.plot(Xnew[:,0],d2alphadT2,'gp')
plt.plot(Xnew[:,0],d2alphadT2nu,'r--')
plt.plot(Xnew[:,0],d2alphadT2_PR,'k--')
plt.xlabel("Tr")
plt.ylabel("d2alphadT2")
plt.legend(["GP","GPnumerical","PR"])

# 3d plot of d2alphadT2
fig, ax = plt.subplots(subplot_kw={"projection":"3d"})
ax.set_title("d2alphadT2")
ax.scatter(Xnew[:,0], Xnew[:,1], d2alphadT2_PR, 'r')
ax.scatter(Xnew[:,0], Xnew[:,1], d2alphadT2, 'g')
ax.set_xlabel("Tr")
ax.set_ylabel("Pr")
plt.legend(["PR", "AlphaGP"])
fig.savefig("figs/PythonAlphaGP_d2alphadT2_%s.png"%fluid)

# ================================================
# Extra validation
ypred_train = k0(X, X).T @ Ki @ y
ypred_test = k0(X, Xtest).T @ Ki @ y

fig = plt.figure()
plt.plot(y, ypred_train, 'rs', fillstyle="none", label='Training Data')
plt.plot(ytest, ypred_test, 'gv', fillstyle="none", label='Test Data')
plt.xlabel("Groundtruth")
plt.ylabel("Prediction")
plt.legend()
fig.savefig("figs/PythonAlphaGP_Validation_%s.png"%fluid)

plt.show()
