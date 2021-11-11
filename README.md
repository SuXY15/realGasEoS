### Real Gas Effect in Reacting Flow

>  This code repository is used to evaluate the advantages of the proposed real gas equation of state (EoS)

#### 1. Code Environment

```shell
# Cantera used in /opt/cantera_libs/cantera, one can modified code in `src` and build & install by
conda activate cantera
python3 /usr/bin/scons build
python3 /usr/bin/scons install
```

```shell
# The code in this repo can be run by
source /usr/local/bin/setup_cantera
python3 src/xxx.py
```



#### 2. Basic Theory

##### 1.0 Generalized EoS

Helmholtz Free Energy: $ F \equiv U-TS =  A(T,v,n)$, where $U$ is the internal energy of the system, $T$ is the absolute temperature, $S$ is the entropy of the system.

+ 1st law of thermodynamics: $dU = \delta Q + \delta W$
+ 2nd law of thermodynamics: $\delta Q = TdS$

+ pressure work: $\delta W = -pdV$

Thus we have:

+ $ dU = TdS-pdV $

+ => $dU = d(TS) - SdT - pdV$
+ => $d(U-TS) = -SdT - pdV$
+ => $dF = -SdT - pdV$

Therefore,
$$
\begin{align*}
S &= -\left( \frac{\partial F}{\partial T} \right)_V \\
p &= -\left( \frac{\partial F}{\partial V} \right)_T \\
U &= F + TS = F - T \left( \frac{\partial F}{\partial T} \right)_V \\
H &= U + pV = F - T \left( \frac{\partial F}{\partial T} \right)_V - V \left( \frac{\partial F}{\partial V} \right)_T \\
C_p &= \left( \frac{\partial H}{\partial T} \right)_p = T \left[ \left( \frac{\partial^2 F}{\partial V\partial T} \right)_{V,T}^2 / \left( \frac{\partial^2 F}{\partial V^2} \right)_T - \left( \frac{\partial^2 F}{\partial T^2} \right)_V\right] \\
C_V &= \left( \frac{\partial U}{\partial T} \right)_V = -T\left( \frac{\partial^2 F}{\partial T^2} \right)_V
\end{align*}
$$

##### 1.1 Ideal Gas

$$
p = \frac{nRT}{V_m}
$$

where $V_m$ is the molar volume of the gas .

##### 1.2 Van der Waals (vdW) real gas equation _1783_

$$
p=\frac{RT}{V_m-b} - \frac{a}{V_m^2}
$$

where $b$ is the volume that is occupied by one mole of the molecules. vdW is a `cubic` type EoS.

##### 1.4 Redlich-Kwong (RK) real gas equation _1948_

$$
p = \frac{RT}{V_m - b} - \frac{a}{\sqrt{T} V_m (V_m+b)}
$$

where $a$ is a constant that corrects for attractive potential of molecules and $b$ is a constant that corrects for volume. The constants are different depending on which gas is being analyzed and can be calculated from the critical point data of corresponding gas:
$$
\begin{align*}
a &= \frac{1}{9(\sqrt[3]{2}-1)} \frac{R^2 T_c^{2.5}}{P_c} \\
b &= \frac{\sqrt[3]{2}-1}{3} \frac{RT_c}{P_c}
\end{align*}
$$
Where $T_c$ is the temperature at the critical point and $P_c$ is the pressure at the critical point. Noted that $a$ and $b$ are solved analytically from the thermodynamic criteria for the critical point:
$$
\left( \frac{\partial P}{\partial V} \right)_T = 0, \left( \frac{\partial^2 P}{\partial V^2} \right)_T = 0.
$$
For more details, please refer to original paper or wikipedia of Redlich Kwong equation of state.

##### 1.5 Soave modified Redlich-Kwong (SRK) real gas equation _1972_

In 1966, Barner noted that the acentric factor $\omega$ can be used to improve RK EoS. SRK is proposed by Soave as:
$$
p = \frac{RT}{V_m - b} - \frac{a \alpha}{V_m (V_m+b)}
$$
where
$$
\begin{align*}
\alpha &= \left [ 1+(0.480+1.574\omega-0.176\omega^2)(1-\sqrt{T_r}) \right]^2 \\

a &= \frac{1}{9(\sqrt[3]{2}-1)} \frac{R^2 T_c^{2}}{P_c} \\
b &= \frac{\sqrt[3]{2}-1}{3} \frac{RT_c}{P_c}
\end{align*}
$$
where $T_r=T/T_c$ is the reduced temperature of the compound. Noted that $\frac{1}{\sqrt{T}}$ term in RK is vanished and its effect is put into $a$ in SRK.

##### 1.6 Peng-Robinson (PR) real gas equation _1976_

The Peng-Robinson EoS further modified the RK EoS by modifying the attractive term, giving
$$
p = \frac{R T}{V_m-b} - \frac{a\alpha}{V_m(V_m+b) + b(V_m-b)}
$$
where
$$
\begin{align*}
\alpha &= \left [ 1+(0.37464+1.54226\omega-0.26992\omega^2)(1-\sqrt{T_r}) \right]^2 \\

a &= 0.457235 \frac{R^2 T_c^{2}}{P_c} \\
b &= 0.077796 \frac{RT_c}{P_c}
\end{align*}
$$

##### 1.7 PR-RK real gas equation _2005_

The limited accuracy of the two-parameter cubic EoS is caused by the inconsistency in the density dependence of its equation rather than the empirical nature of model constants $a$ and $b$. Many researchers proposed three parameters EoS.  In 2005, Cismondi and Mollerup have  proposed a three-parameter cubic EoS (hereafter RK-PR EoS) by introducing an additional parameter, $\delta_1$. 
$$
p = \frac{R T}{V_m-b} - \frac{a\alpha}{(V_m+\delta_1b)(V_m+\delta_2b)}
$$
where $\delta_2$ is not independent and $\delta_2 = (1-\delta_1)(1+\delta_1)$. And it will degenerate to SRK when $\delta_1=1, \delta_2=0$; and degenerate to PR when  $\delta_1 = 1+\sqrt2, \delta_2 = 1-\sqrt2 $.

##### 1.7 Mixture properties

When handling mixture, to obtain the aggregate parameters $a_m$, $b_m$ of the cubic EoS, mixing rules from the van der Waals one-fluid theory is applied:
$$
\begin{align*}
a_m &= \sum_i^{N_s} \sum_j^{N_s} X_i X_j a_{ij} \\
b_m &= \sum_i^{N_s} X_i b_i
\end{align*}
$$
where $X_i$ is the mole fraction of the $i^{th}$ components of the mixture,  $X_j$ is the mole fraction of the $j^{th}$ components of the mixture. Besides, $a_{ij}$ is the attractive term between a molecule of species $i$ and species $j$, which is usually calculated by $a_{ij} = (a_i a_j)^{1/2}$  or  $a_{ij}=(a_i a_j)^{1/2}(1-\delta_{ij})$.  $\delta_{ij}$

is an empirically determined binary interaction coefficient characterizing the binary formed by component i and component j.

#### 2. The calibration model of alpha

The $\alpha$ is previously calculated by Sovae-type functions like

+ In SRK: $\alpha = \left [ 1+(0.480+1.574\omega-0.176\omega^2)(1-\sqrt{T_r}) \right]^2$
+ In PR: $\alpha = \left [ 1+(0.37464+1.54226\omega-0.26992\omega^2)(1-\sqrt{T_r}) \right]^2 $

But those are empirical approximations and would lead to notable effects in particular thermodynamically sensitive simulations. Here **a Gaussian Process (GP) based semi-parametric method** is used to utilize the information of experimental data as most as possible.

##### 2.1 The GP model

Here we're trying to fit a function that can map the inputs including $T_r$ to the ouput $\alpha$. Using the concept of Gaussian Process (GP) model, given the observations $(\boldsymbol X,y)$, when applying GP to predict the target $y_*$ at new input $\boldsymbol X_*$, one have
$$
\begin{align*}
\mu(y_*) & = \boldsymbol K_*^T (\boldsymbol K+\sigma^2_n\boldsymbol I)^{-1} y \\
\sigma^2(y_*) &= \boldsymbol K_{**} -\boldsymbol K_*^T (\boldsymbol K+\sigma^2_n\boldsymbol I)^{-1}\boldsymbol K_*
\end{align*}
$$
where $\sigma^2_n$ is the covariance of estimated Gaussian noise from observed data $(\boldsymbol X,y)$, $\boldsymbol K$ related terms are the covariance matrices as
$$
\begin{align*}
\boldsymbol K  &= cov(\boldsymbol X, \boldsymbol X) \\
\boldsymbol K_* &= cov(\boldsymbol X, \boldsymbol X_*) \\
\boldsymbol K_{**} &= cov(\boldsymbol X_*, \boldsymbol X_*) 
\end{align*}
$$

##### 2.2 Kernel functions

The covariance are usually obtained via kernel functions, which are very important for GP's data training. The most famous ones are:

+ Squared exponential (or say Radial Basis Function, RBF):
  $$
  cov(x_i, x_j) = \exp(-\frac{d(x_i,x_j)^2}{2l^2})
  $$
  where $d(x_i,x_j)$ is the Euclidean distance of $x_i$ and $x_j$,  $l$ is the length scale of kernel size.

+ Matern:
  $$
  k(x_i, x_j) = \frac{1}{\Gamma(\nu)2^{\nu-1}}\Bigg(\frac{\sqrt{2\nu}}{l} d(x_i , x_j )\Bigg)^\nu K_\nu\Bigg(\frac{\sqrt{2\nu}}{l} d(x_i , x_j )\Bigg),
  $$
  where $\Gamma$ is the Gamma function, $K_\nu$ is the modified Bessel function of order $\nu$.

##### 2.3 Basis functions

The basis function provides an important ability to model the mean of $y$ to zero.
$$
y(\boldsymbol x) = f(\boldsymbol x; \boldsymbol \theta) + b_{\boldsymbol \theta}(\boldsymbol x) + \epsilon
$$
where $\epsilon$ is the (Gaussian) stochastic error in the database, $f(\boldsymbol x; \boldsymbol \theta)$  is a basis function that can model $y(\boldsymbol x)$ as much as possible, and $b_{\boldsymbol \theta}(\boldsymbol x) $ is the bias term that can compensate the error between $f(\boldsymbol x; \boldsymbol \theta)$ and $y(\boldsymbol x)$.

To solve this optimization problem, let $f(\boldsymbol x; \boldsymbol \theta) = \theta_3 x_2 + \theta_2 x_1 + \theta_1 $ as a linear combination of the input variables $\boldsymbol x$, where the input variable are $\boldsymbol x=[x_1, x_2]$, $x_1$ and $x_2$ are the input from each dimension; and $\boldsymbol \theta = [\theta_1, \theta_2, \theta_3]$ are corresponding parameters in the linear combination function $f$.

The corresponding kernel function then be formed as:
$$
k(x_i,x_j) = k0(x_i,x_j) - h_\theta(x_i)^T \boldsymbol H_\theta h_\theta(x_j)
$$
where $h_\theta(x)$ is the 1st order global sensitivity obtained from:
$$
h_\theta(x) = \int_\boldsymbol X \nabla_\theta f (\boldsymbol x) k_0(x,\xi)d\xi
$$
the 2nd order global sensitivity, i.e., Hessian matrix is
$$
\boldsymbol H_\theta = \int_\boldsymbol X \int_\boldsymbol X \left (\nabla_\theta f(\xi)\right) \left (\nabla_\theta f(\zeta)\right)^T k_0( \xi, \zeta) d \xi d \zeta
$$
The parameter gradient of $f$ shows 
$$
\nabla_\theta f (\boldsymbol x) = \frac{\partial }{\partial \boldsymbol \theta}f(\boldsymbol x; \boldsymbol \theta) = \left[ \begin{align*} 1 \\ x_1 \\ x_2 \end{align*} \right]
$$
The matrix constructed by gradients' outer product of randomly samples $ \xi, \zeta$
$$
\left (\nabla_\theta f(\boldsymbol \xi)\right) \left (\nabla_\theta f(\boldsymbol \zeta)\right)^T = \left[ \begin{matrix} 1 \\ \xi_1 \\ \xi_2 \end{matrix} \right] \left[ \begin{matrix} 1 & \zeta_1 & \zeta_2 \end{matrix} \right] = \left[ \begin{matrix} 1 & \zeta_1 & \zeta_2 \\ \xi_1 & \xi_1\zeta_1 & \xi_1\zeta_2 \\ \xi_2 & \xi_2\zeta_1 & \xi_2\zeta_2 \end{matrix} \right]
$$
The kernel function $k_0$ is:
$$
k_0( \xi, \zeta) = \exp(-\sum_i\frac{(\xi_i-\zeta_i)^2}{\gamma_i}) = \exp(-\frac{(\xi_1-\zeta_1)^2}{\gamma_1})\exp(-\frac{(\xi_2-\zeta_2)^2}{\gamma_2})
$$
where $\gamma_1, \gamma_2$ are the parameters of kernel function in dimension $x_1$ and $x_2$.

The corresponding distribution of the evaluated output $y_*$ at new  $x_*$ should be:
$$
\begin{align*}

\mu(y_*) & = f(x_*;\theta) + \boldsymbol K_*^T \left(\boldsymbol K+\frac{\boldsymbol K_{**}}{\nu} \boldsymbol I\right)^{-1} \left[y-f(\boldsymbol X;\theta)\right] \\
\sigma(y_*) &= \nu + \boldsymbol K_{**} - \nu\boldsymbol K_*^T \left(\boldsymbol K+\frac{\boldsymbol K_{**}}{\nu} \boldsymbol I\right)^{-1} \boldsymbol K_*

\end{align*}
$$
where $\nu$ is a measurement of possible noises,  its maximum likelihood estimation (MLE)  $\hat \nu$ is
$$
\hat \nu = \frac{1}{N} [y-f(\boldsymbol X;\theta)]^T (\boldsymbol K+g\boldsymbol I)^{-1} [y-f(\boldsymbol X;\theta)]
$$
where $g$ is the original possible noise level, which still should be estimated by optimization program.



##### 2.4 Partial derivative

The partial derivative with respect to temperature shows
$$
\begin{align*}

\mu(y_*) & = \boldsymbol K_*^T (\boldsymbol K+\sigma^2_n\boldsymbol I)^{-1} y \\
\frac{\partial\mu(y_\star)}{\partial T}&=(\frac{\partial \boldsymbol K_*}{\partial T})^T(\boldsymbol K+\sigma^2_n\boldsymbol I)^{-1} y\\
k_{ij}&=exp(-\frac{(X_1(i,1)-X_*(j,1))^2}{\gamma_1}--\frac{(X_1(i,2)-X_*(j,2))^2}{\gamma_2})\\
\frac{\partial k_{ij}}{\partial T}&=\frac{2X_1(i,1)-2X_*(j,1)}{T_r\gamma_1}k_{ij}\\
\frac{\partial^2k_{ij}}{\partial T^2}&=\frac{-2}{T_r\gamma_1}k_{ij}+\frac{(2X_1(i,1)-2X_*(j,1))^2}{T_r^2\gamma_1^2}k_{ij}\\
h_\theta

\end{align*}
$$


##### 2.5 Validation



#### 3. Appendix

Using the chain rule from Table 5, Page 4 of [Partial derivatives of thermodynamic state properties for dynamic simulation](https://link.springer.com/article/10.1007/s12665-013-2394-z),
$$
\begin{align*}
C_p = \left( \frac{\partial H}{\partial T} \right)_p = \left( \frac{\partial H}{\partial T} \right)_V - \left( \frac{\partial H}{\partial V} \right)_T \left( \frac{\partial p}{\partial T} \right)_V \left( \frac{\partial p}{\partial V} \right)_T^{-1} 
\end{align*}
$$
where
$$
\begin{align*}
\left( \frac{\partial H}{\partial T} \right)_V
	&= \left( \frac{\partial F}{\partial T} \right)_V - \left( \frac{\partial F}{\partial T} \right)_V - T \left( \frac{\partial^2 F}{\partial T^2} \right)_V  - V \left( \frac{\partial^2 F}{\partial V\partial T} \right)_{T,V} \\ 
	&= - T \left( \frac{\partial^2 F}{\partial T^2} \right)_V  - V \left( \frac{\partial^2 F}{\partial V\partial T} \right)_{T,V} \\

\left( \frac{\partial H}{\partial V} \right)_T
	&=  \left( \frac{\partial F}{\partial V} \right)_T - T \left( \frac{\partial^2 F}{\partial T\partial V} \right)_{V,T} - \left( \frac{\partial F}{\partial V} \right)_T - V \left( \frac{\partial^2 F}{\partial V^2} \right)_T \\
	&= - T \left( \frac{\partial^2 F}{\partial T\partial V} \right)_{V,T} - V \left( \frac{\partial^2 F}{\partial V^2} \right)_T \\

\left( \frac{\partial p}{\partial T} \right)_V 
	&= -\left( \frac{\partial^2 F}{\partial V\partial T} \right)_{T,V} \\

\left( \frac{\partial p}{\partial V} \right)_T 
	&= -\left( \frac{\partial^2 F}{\partial V^2} \right)_T
\end{align*}
$$
therefore,
$$
\begin{align*}
C_p 
&= \left( \frac{\partial H}{\partial T} \right)_V - \left( \frac{\partial H}{\partial V} \right)_T \left( \frac{\partial p}{\partial T} \right)_V \left( \frac{\partial p}{\partial V} \right)_T^{-1} \\
&= - T \left( \frac{\partial^2 F}{\partial T^2} \right)_V  - V \left( \frac{\partial^2 F}{\partial V\partial T} \right)_{T,V} + \left[ T \left( \frac{\partial^2 F}{\partial T\partial V} \right)_{V,T} + V \left( \frac{\partial^2 F}{\partial V^2} \right)_T \right]\left( \frac{\partial^2 F}{\partial V\partial T} \right)_{T,V}/\left( \frac{\partial^2 F}{\partial V^2} \right)_T \\
&= - T \left( \frac{\partial^2 F}{\partial T^2} \right)_V  +T \left( \frac{\partial^2 F}{\partial T\partial V} \right)_{V,T}  \left( \frac{\partial^2 F}{\partial V\partial T} \right)_{T,V}/\left( \frac{\partial^2 F}{\partial V^2} \right)_T \\
&= T \left[ \left( \frac{\partial^2 F}{\partial V\partial T} \right)_{V,T}^2 / \left( \frac{\partial^2 F}{\partial V^2} \right)_T - \left( \frac{\partial^2 F}{\partial T^2} \right)_V\right]
\end{align*}
$$



