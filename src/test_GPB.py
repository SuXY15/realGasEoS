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

# define basis function
def f(x, θ):
    return θ[0] + θ[1]*x[:,0] + θ[2]*x[:,1]

def dfdθ(x, θ):
    if len(x.shape)<=1:
        return np.array([1, x[0], x[1]])
    else:
        N = len(x)
        return np.vstack([np.ones(N), x[:,0], x[:,1]]).T

K0  = k0(X, X)
df = dfdθ(X,θ)
Hm = [np.outer(df[i], df[j]) * K0[i,j] for i in range(len(X)) for j in range(len(X))]
H = np.mean(Hm, axis=0)

def h(X1, X, θ):
    K0x = k0(X1, X)
    df = dfdθ(X1,θ)
    return df.T * np.mean(K0x, axis=1)

def k(X1, X2, X, H, θ):
    return k0(X1,X2) + h(X1,X,θ).T @ np.linalg.inv(H) @ h(X2,X,θ)


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
#print(Xnewfront)

# ================================================
# GP predict
K = k(X, X, X, H, θ)
Ki = np.linalg.inv(K + 1e-8*np.diag(np.ones(len(X))))
Kxx = k(Xnew, Xnew, X, H, θ)
Kx = k(X, Xnew, X, H, θ)

mu = f(Xnew,θ) + Kx.T @ Ki @ (y-f(X,θ))
cov = Kxx - Kx.T @ Ki @ Kx
sigma = np.diag(cov)

Ypred = f(X,θ) + K.T @ Ki @ (y-f(X,θ))
# GP predict behind
Kxx1 = k(Xnew_r, Xnew_r, X, H, θ)
Kx1 = k(X, Xnew_r, X, H, θ)

mu1 = f(Xnew_r,θ) + Kx1.T @ Ki @ (y-f(X,θ))
cov1 = Kxx1 - Kx1.T @ Ki @ Kx1
sigma1 = np.diag(cov1)
#print(mu1)
# GP predict front
Kxx2 = k(Xnew_l, Xnew_l, X, H, θ)
Kx2 = k(X, Xnew_l, X, H, θ)

mu2 = f(Xnew_l,θ) + Kx2.T @ Ki @ (y-f(X,θ))
cov2 = Kxx2 - Kx2.T @ Ki @ Kx2
sigma2 = np.diag(cov2)

# ================================================
# GP gradient, numerical
dalphadT = [(mu1[i]-mu2[i])/2/deltaTr/Tc for i in range(M)]
d2alphadT2 = [(mu1[i]-2 * mu[i]+mu2[i])/deltaTr**2/Tc**2 for i in range(M)]
#print(dalphadT)
#print(d2alphadT2)

# ================================================
# PR gradient
dalphadT_PR = [PR_dalphadT(Xnew[i,0]*Tc, Xnew[i,1]*Pc, Tc, Pc, omega) for i in range(M)]
d2alphadT2_PR = [PR_d2alphadT2(Xnew[i,0]*Tc, Xnew[i,1]*Pc, Tc, Pc, omega) for i in range(M)]

# ================================================
# 2d plot
plt.plot(y, y, 'k--')
plt.plot(y, Ypred, 'gp', fillstyle='none')
plt.xlabel("Groundtruth")
plt.ylabel("Prediction")

# 2d plot of dalphadT
fig=plt.figure()
plt.plot(Xnew[:,0],dalphadT,'gp')
plt.plot(Xnew[:,0],dalphadT_PR,'k--')
plt.xlabel("Tr")
plt.ylabel("dalphadT")
plt.legend(["f+GP","PR"])

# 2d plot of da2lphadT2
fig=plt.figure()
plt.plot(Xnew[:,0],d2alphadT2,'gp')
plt.plot(Xnew[:,0],d2alphadT2_PR,'k--')
plt.xlabel("Tr")
plt.ylabel("d2alphadT2")
plt.legend(["f+GP","PR"])

# 3d plot
fig, ax = plt.subplots(subplot_kw={"projection":"3d"})
ax.scatter(X[:,0], X[:,1], y, 'r')
ax.scatter(Xnew[:,0], Xnew[:,1], mu, 'g')
ax.set_xlabel("Tr")
ax.set_ylabel("Pr")
ax.set_zlabel("Alpha")
plt.legend(["Groundtruth ", "Prediction"])

# ================================================
# Extra validation
ypred_train = f(X,θ) + k(X, X, X, H, θ).T @ Ki @ (y-f(X,θ))
ypred_test = f(Xtest,θ) + k(X, Xtest, X, H, θ).T @ Ki @ (y-f(X,θ))

fig = plt.figure()
plt.plot(y, ypred_train, 'rs', fillstyle="none", label='Training Data')
plt.plot(ytest, ypred_test, 'gv', fillstyle="none", label='Test Data')
plt.xlabel("Groundtruth")
plt.ylabel("Prediction")
plt.legend()
fig.savefig("figs/PythonAlphaGPB_Validation_%s.png"%fluid)

plt.show()
