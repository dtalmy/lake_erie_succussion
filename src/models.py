import numpy as np

# global parameters
k,l = 1,2

def get_r(T,a=0.2,b=0.06):
    return np.exp(b*T)

def dyn_temp(t,tmin=0,tmax=25,period=365):
    return tmin + (tmax-tmin)/2*(np.sin((t-150)/period*2*np.pi)+1)

def corbor(cor,n):
    return cor**n

def get_params(temp,cor=1,bor=1,Tdyn=False):
    mu1max = 2.6
    mu2max = 2.6
    alpha1 = 0.12
    alpha2 = 0.12
    mp = 0.15
    mz = 0.15
    F = 5.0
    w1 = 1.0
    w2 = 1.0
    tau1 = 10
    tau2 = 10
    eps = 0.1
    Sin = 100.0
    NtoPinorg = 40
    tauN = 180.0
    tauP = 5.0
    Qmin = 5
    Qmax = 80
    rho = 3.0
    return (mu1max,alpha1,mp,mz,F,w1,tau1,eps,Sin,NtoPinorg,tauN,tauP,Qmin,Qmax,temp,Tdyn,cor,bor,rho)

# parallel model
def parallel(u,t,ps):
    mu1max,alpha1,mp,mz,F,w1,tau1,\
    w2,Sin,temp,Tdyn,mu2max,alpha2,rho = ps[0],ps[1],ps[2], \
    ps[3],ps[4],ps[5],ps[6],ps[7],ps[8],ps[9],ps[10],ps[11],ps[12],ps[13]
    if Tdyn > 0:
        r = get_r(dyn_temp(t))
    else:
        r = get_r(temp)
    eps = 0.1
    mu2max = mu2max*r
    mu1max = mu1max*r
    alpha2 = alpha2*r
    alpha1 = alpha1*r
    F = F*r
    tau2 = tau1*r
    tau1 = tau1*r
    mp = mp
    mz = mz
    Sin = Sin
    Nn,P1n,P2n,Z1n,Z2n = u[0],u[1],u[2],u[3],u[4]
    alphaN1, alphaN2 = alpha1,alpha2
    mu1 = mu1max*Nn/(Nn+mu1max/alphaN1)
    mu2 = mu2max*Nn/(Nn+mu2max/alphaN2)
    g1 = F*w1*P1n/(1+F*w1*P1n/tau1)
    g2 = F*w2*P2n/(1+F*w2*P2n/tau2)
    dNndt = Sin - mu1*P1n/16.0 - mu2*P2n/16.0 - rho*Nn
    dP1ndt = mu1*P1n - mp*P1n**k  - g1*Z1n
    dP2ndt = mu2*P2n - mp*P2n**k  - g2*Z2n
    dZ1ndt = eps*g1*Z1n  - mz*Z1n**l
    dZ2ndt = eps*g2*Z2n  - mz*Z2n**l
    return np.r_[[dNndt,dP1ndt,dP2ndt,dZ1ndt,dZ2ndt]]

# diamond model
def diamond(u,t,ps):
    mu1max,alpha1,mp,mz,F,w1,tau1,\
    w2,Sin,temp,Tdyn,mu2max,alpha2,rho = ps[0],ps[1],ps[2], \
    ps[3],ps[4],ps[5],ps[6],ps[7],ps[8],ps[9],ps[10],ps[11],ps[12],ps[13]
    if Tdyn > 0:
        r = get_r(dyn_temp(t))
    else:
        r = get_r(temp)
    eps = 0.1
    mu2max = mu2max*r
    mu1max = mu1max*r
    alpha2 = alpha2*r
    alpha1 = alpha1*r
    F = F*r
    tau2 = tau1*r
    tau1 = tau1*r
    mp = mp
    mz = mz
    Sin = Sin
    Nn,P1n,P2n,Zn = u[0],u[1],u[2],u[3]
    alphaN1, alphaN2 = alpha1,alpha2
    mu1 = mu1max*Nn/(Nn+mu1max/alphaN1)
    mu2 = mu2max*Nn/(Nn+mu2max/alphaN2)
    g1 = F*w1*P1n/(1+F*w1*P1n/tau1)
    g2 = F*w2*P2n/(1+F*w2*P2n/tau2)
    g1 = F*w1*P1n/(1+F*(w1*P1n/tau1 + w2*P2n/tau2))
    g2 = F*w2*P2n/(1+F*(w1*P1n/tau1 + w2*P2n/tau2))
    dNndt = Sin - mu1*P1n/16.0 - mu2*P2n/16.0 - rho*Nn
    dP1ndt = mu1*P1n - mp*P1n**k  - g1*Zn
    dP2ndt = mu2*P2n - mp*P2n**k  - g2*Zn
    dZndt = eps*(g1+g2)*Zn  - mz*Zn**l
    return np.r_[[dNndt,dP1ndt,dP2ndt,dZndt]]

