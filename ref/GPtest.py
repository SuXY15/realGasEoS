import numpy as np
import matplotlib.pyplot as plt
np.random.seed(0)

N = 20
dim = 1
xlim = [-1, 1]

##generate data
X = np.random.rand(N, dim)*(xlim[1]-xlim[0]) + xlim[0]
f0 = 0
g = np.random.rand(dim)
H = np.random.rand(dim,dim)
f = np.array([f0 + X[i].T@g + X[i] @ H @ X[i].T for i in range(N)])

## add noise
sigma_noise = 0.1
y = f + np.random.normal(0, sigma_noise, size=(N))

## define covariance
def k0(X1, X2, kernel_size=1):
    cov = np.zeros((len(X1), len(X2)))
    for i in range(len(X1)):
        for j in range(len(X2)):
            cov[i,j] = np.exp(-np.sum((X1[i] - X2[j])**2)/kernel_size)
    return cov

## GP predict
M = 100
Xnew = sorted(np.random.rand(M, dim)*(xlim[1]-xlim[0]) + xlim[0])
Kxx = k0(Xnew, Xnew)
Kx = k0(X, Xnew)
K  = k0(X, X)
Ki = np.linalg.inv(K + 1e-1*np.diag(np.ones(N)))

mu = Kx.T @ Ki @ y
cov = Kxx - Kx.T @ Ki @ Kx
sigma = np.diag(cov)

plt.plot(X, f, 'ko')
plt.plot(X, y, 'rv')
plt.plot(Xnew, mu, 'g-')
plt.plot(Xnew, mu+3*sigma, 'g--')
plt.plot(Xnew, mu-3*sigma, 'g--')

plt.legend(["Groundtruth ", "Observation", "Preduction"])
plt.show()
