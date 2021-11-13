from utils import *
from copy import deepcopy

# name = "o2"
name = 'c12h26'

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
Hm = [np.outer(df[i], df[j]) * K0[i,j] for i in range(N) for j in range(N)]
H = np.mean(Hm, axis=0)

def h(X1, X, θ):
    K0x = k0(X1, X)
    df = dfdθ(X1,θ)
    return df.T * np.mean(K0x, axis=1)

def k(X1, X2, X, H, θ):
    return k0(X1,X2) + h(X1,X,θ).T @ H @ h(X2,X,θ)

# ================================================
# random queries
M = 200
x1lim = [np.min(X[:,0]), np.max(X[:,0])]
x2lim = [np.min(X[:,1]), np.max(X[:,1])]
Xnew = np.random.rand(M, dim)
Xnew[:,0] = Xnew[:,0] * (x1lim[1] - x1lim[0]) + x1lim[0]
Xnew[:,1] = Xnew[:,1] * (x2lim[1] - x2lim[0]) + x2lim[0]
Xnew = Xnew[np.argsort(Xnew[:,1]),:]
Xnew = Xnew[np.argsort(Xnew[:,0]),:]

K = k(X, X, X, H, θ)
Ki = np.linalg.inv(K + 1e-4*np.diag(np.ones(N)))
Kxx = k(Xnew, Xnew, X, H, θ)
Kx = k(X, Xnew, X, H, θ)

mu = f(Xnew,θ) + Kx.T @ Ki @ (y-f(X,θ))
cov = Kxx - Kx.T @ Ki @ Kx
sigma = np.diag(cov)

Ypred = f(X,θ) + K.T @ Ki @ (y-f(X,θ))

# ================================================
# 2d plot
plt.plot(y, y, 'k--')
plt.plot(y, Ypred, 'gp', fillstyle='none')
plt.xlabel("Groundtruth")
plt.ylabel("Prediction")

# 3d plot
fig, ax = plt.subplots(subplot_kw={"projection":"3d"})
ax.scatter(X[:,0], X[:,1], y, 'r')
ax.scatter(Xnew[:,0], Xnew[:,1], mu, 'g')
ax.set_xlabel("Tr")
ax.set_ylabel("Pr")
ax.set_zlabel("Alpha")
plt.legend(["Groundtruth ", "Prediction"])

plt.show()
