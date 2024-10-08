import numpy as np
import pandas as pd
import pylab as plt
import seaborn as sns
import scipy.integrate as si

################################################
# local functions for dataframe queries
################################################

# global parameters
k,l = 1,2

# diamond model
def diamond(u,t,mu1max,mu2max,alpha1,alpha2,mp,mz,F,w1,w2,tau1,tau2,eps,Sin):
    N,P1,P2,Z = u[0],u[1],u[2],u[3]
    dNdt = Sin - mu1max*N/(N+mu1max/alpha1)*P1 -mu2max*N/(N+mu2max/alpha2)*P2 \
            + 0.1*(mp*(P1**k +P2**k ) + mz*Z**l \
            + (1-eps)*F*(w1*P1 + w2*P2)/(1+F*(w1*P1/tau1 + w2*P2/tau2))*Z) \
            - 3.0*N
    dP1dt = mu1max*N/(N+mu1max/alpha1)*P1 - mp*P1**k  \
            - F*w1*P1/(1+F*(w1*P1/tau1 + w2*P2/tau2))*Z
    dP2dt = mu2max*N/(N+mu2max/alpha2)*P2 - mp*P2**k  \
            - F*w2*P2/(1+F*(w1*P1/tau1 + w2*P2/tau2))*Z
    dZdt = eps*F*(w1*P1 + w2*P2)/(1+F*(w1*P1/tau1 + w2*P2/tau2))*Z - mz*Z**l
    return np.r_[[dNdt,dP1dt,dP2dt,dZdt]]

# parallel model
def parallel(u,t,mu1max,mu2max,alpha1,alpha2,mp,mz,F,w1,w2,tau1,tau2,eps,Sin):
    N,P1,P2,Z1,Z2 = u[0],u[1],u[2],u[3],u[4]
    dNdt = Sin - mu1max*N/(N+mu1max/alpha1)*P1 -mu2max*N/(N+mu2max/alpha2)*P2 \
            + 0.1*(mp*(P1**k +P2**k ) + mz*(Z1**l+Z2**l) \
            + (1-eps)*F*(w1*P1)/(1+F*(w1*P1/tau1))*Z1 \
            + (1-eps)*F*(w2*P2)/(1+F*(w2*P2/tau2))*Z2) \
            - 3.0*N
    dP1dt = mu1max*N/(N+mu1max/alpha1)*P1 - mp*P1**k  \
            - F*w1*P1/(1+F*w1*P1/tau1)*Z1
    dP2dt = mu2max*N/(N+mu2max/alpha2)*P2 - mp*P2**k  \
            - F*w2*P2/(1+F*w2*P2/tau2)*Z2
    dZ1dt = eps*F*(w1*P1)/(1+F*(w1*P1/tau1))*Z1  - mz*Z1**l
    dZ2dt = eps*F*(w2*P2)/(1+F*(w2*P2/tau2))*Z2  - mz*Z2**l
    return np.r_[[dNdt,dP1dt,dP2dt,dZ1dt,dZ2dt]]

def get_r(a,b,T):
    return np.exp(b*T)

def get_inits(d,Ntot):
    fnut = 0.1
    felse = 1 - fnut
    if d == 'parallel':
        return np.r_[[fnut,0.25,0.25,0.25,0.25]]*Ntot
    if d == 'diamond':
        return np.r_[[fnut,felse/4.0,felse/4.0,felse/2.0]]*Ntot

def dyn_temp(t,tmin,tmax,period=365):
    return tmin + (tmax-tmin)/2*(np.sin((t-60)/period*2*np.pi)+1)

def get_params(r,cor=1,bor=1):
    bor = corbor(cor,n)
    mu1max = 1.6*r
    mu2max = 1.6*r*cor
    alpha1 = 0.08*r
    alpha2 = 0.08*cor*r
    mp = 0.15
    mz = 0.15
    F = 0.3*r
    w1 = 1.0
    w2 = 1.0*bor
    tau1 = 10*r
    tau2 = 10*r*bor
    eps = 0.1
    Sin = 5.0
    return (mu1max,mu2max,alpha1,alpha2,mp,mz,F,w1,w2,tau1,tau2,eps,Sin)

def corbor(cor,n):
    return cor**n

################################################
# make figures
################################################

# make all font Times New Roman
plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams["font.size"] = 14

# colors and line styles
nc,p1c,p2c,zc,z2c = 'bo-','g^--','g*-','ro--','ro-' # for intermittent model samples
ncp,p1cp,p2cp,zcp,z2cp = 'b','g--','g','r--','r' # for continious plots
nca,p1ca,p2ca,zca,z2ca = 'bo','g^','g*','ro','ro' # for intermittent model samples
ncb,p1cb,p2cb,zcb,z2cb = 'b-','g--','g-','r--','r-' # for intermittent model samples

# labels
p1lab,p2lab = 'Non-cyanobacteria chlorophyll','Phycocyanin'

mnl,mp1l,mp2l,mz1l,mz2l = r'$N$',r'$P_1$',r'$P_2$',r'$Z_1$',r'$Z_2$'

f1,ax1 = plt.subplots(3,1,figsize=[6,12])
f2,ax2 = plt.subplots(3,1,figsize=[6,12])
f3,ax3 = plt.subplots(1,2,figsize=[9,4])
f4,ax4 = plt.subplots(2,2,figsize=[12,12])
f5,ax5 = plt.subplots(3,1,figsize=[6,12])
f6,ax6 = plt.subplots(3,1,figsize=[6,12])
f7,ax7 = plt.subplots(3,1,figsize=[6,12])
f8,ax8 = plt.subplots(3,1,figsize=[6,12])
f9,ax9 = plt.subplots(2,2,figsize=[12,12])
f10,ax10 = plt.subplots(2,2,figsize=[12,12])
f11,ax11 = plt.subplots(2,2,figsize=[12,12])
f12,ax12 = plt.subplots()

ax3 = ax3[::-1]

fs=16

ax4[1,0].text(0.9,1.08,'Parallel model',fontsize=20,transform=ax4[1,0].transAxes)
ax4[0,0].text(0.9,1.08,'Diamond model',fontsize=20,transform=ax4[0,0].transAxes)
f4.subplots_adjust(hspace=0.4)

ax9[1,0].text(0.9,1.08,'Parallel model',fontsize=20,transform=ax9[1,0].transAxes)
ax9[0,0].text(0.9,1.08,'Diamond model',fontsize=20,transform=ax9[0,0].transAxes)
f9.subplots_adjust(hspace=0.32)

ax11[0,0].text(0.2,1.1,'High cost of resistance',fontsize=20,transform=ax11[0,0].transAxes)
ax11[0,1].text(0.2,1.1,'Low cost of resistance',fontsize=20,transform=ax11[0,1].transAxes)
ax11[0,1].text(1.1,0.3,'Diamond model',fontsize=20,transform=ax11[0,1].transAxes,rotation=-90)
ax11[1,1].text(1.1,0.3,'Parallel model',fontsize=20,transform=ax11[1,1].transAxes,rotation=-90)

for (axes,f) in zip([ax1,ax2],[f1,f2]):
    f.subplots_adjust(hspace=0.3)
    for (lab,ax,tit) in zip('abc',axes,['Nutrient','Phytoplankton','Zooplankton']):
        ax.set_xlabel('Time (days)')
        ax.set_ylabel(r'log mmol N m$^{-3}$')
        #ax.set_title(tit)
        ax.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=ax.transAxes)

for (axes,f) in zip([ax5,ax6,ax7,ax8],[f5,f6,f7,f8]):
    f.subplots_adjust(hspace=0.3)
    for (lab,ax,tit) in zip('abc',axes,['Nutrient','Phytoplankton','Zooplankton']):
        ax.set_xlabel('Time (days)')
        ax.set_ylabel(r'log mmol N m$^{-3}$')

f3.subplots_adjust(wspace=0.3)
f10.subplots_adjust(wspace=0.3)

ax3[0].set_xlabel('Time (days)')
ax3[0].set_ylabel('Temperature (Celcius)')
ax3[1].set_xlabel('Temperature (Celcius)')
ax3[1].set_ylabel('Temperature response')
for (lab,ax) in zip('ba',ax3):
    ax.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=ax.transAxes)

for (a,lab) in zip(ax4.flatten(),'abcd'):
    a.set_xlabel('Temperature (Celcius)',fontsize=fs)
    a.set_ylabel('log Nitrogen (mmol N m$^{-3}$)',fontsize=fs)
    a.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=a.transAxes)

for (a,lab) in zip(ax9.flatten(),'abcd'):
    a.set_xlabel('Time (days)',fontsize=fs)
    a.set_ylabel('Nitrogen (mmol N m$^{-3}$)',fontsize=fs)
    a.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=a.transAxes)

for (a,lab) in zip(ax10.flatten(),'abcd'):
    a.set_xlabel('Month',fontsize=fs)
    a.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=a.transAxes)
for (a,lab) in zip(ax11.flatten(),'abcd'):
    a.set_xlabel('Time (days)',fontsize=fs)
    a.set_ylabel('log Nitrogen (mmol N m$^{-3}$)',fontsize=fs)
    a.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=a.transAxes)

ax12.set_ylabel('Defense ability')
ax12.set_xlabel('Growth ability')


################################################
# viz data
################################################

import taihu_data
import viz_data as vd
import zoo_viz as zv

ax10[0,0].plot(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.nitname,50),nc,linewidth = 3,label='Nitrogen')
ax10twin = ax10[0,0].twinx()
ax10twin.plot(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.phosname,50),'bo--',linewidth = 3)
ntopdat = np.exp(vd.get_percentiles(vd.sds,vd.nitname,50)) / np.exp(vd.get_percentiles(vd.sds,vd.phosname,50))
ax10[0,1].plot(vd.get_months(vd.sds),ntopdat,label=p1lab,linewidth = 3)
ax10[1,0].plot(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.chlname,50),p1c,label=p1lab,linewidth = 3)
ax10[1,0].plot(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.pycname,50),p2c,label=p2lab,linewidth = 3)
g = sns.boxplot(x = zv.lez[zv.lez['YEAR']>2011]['Month'],\
        y = np.log(zv.lez[zv.lez['YEAR']>2011]['total_zooplankton_biomass']),\
        hue = zv.lez[zv.lez['YEAR']>2011]['Month'], ax=ax10[1,1])
ax10[1,1].set_ylabel('log Total Zooplankton Biomass\n (μg/m$^{3}$)',fontsize=fs)
leg = ax10[1,1].legend(prop={'size':12,'family': 'Times New Roman'})
new_labels = ['April', 'August']
for t, lab in zip(leg.texts, new_labels):
    t.set_text(lab)
leg = ax10[1,0].legend(prop={'size':12,'family': 'Times New Roman'},fontsize=16,loc='best')

ax10[0,0].set_ylabel('log '+vd.nitnamed,fontsize=fs)
ax10twin.set_ylabel('log '+vd.phosnamed,fontsize=fs)
ax10[1,0].set_ylabel('log Pigment (μg/L)',fontsize=fs)

f10.subplots_adjust(wspace=0.4)

################################################
# run model
################################################

# time array
ndays,delt = 1000,0.01
times = np.linspace(0,ndays,int(ndays/delt))
ndayss = 400

# frequency for symbol plotting of model
freq = 1000

# shorter time array for dyn sims
ndaysl,deltl = 365*5,0.05
timesl = np.linspace(0,ndaysl,int(ndaysl/deltl))

# nutrient array
Ntotmin,Ntotmax = 0.001,0.2
Ntots = np.linspace(Ntotmin,Ntotmax,10)
Ntot = Ntotmax

# temps
tmin,tmax = 0,25
temps = np.linspace(tmin,tmax,50)
temp = tmax # default temp

# specify what to plot on x axis in equilibrum
#xdim = Ntots
xdim = temps

# temperature parameters
a = 0.2 # eppley intercept
b = 0.061 # eppley exponent

# number years to show
nyears,nyearsl = 1,3
ncut = (ndaysl/365-nyears)*365
ncutl = (ndaysl/365-nyearsl)*365

# iteration to visualize time dependent dynamics
ndyn = -1

# default trade-off valis
ndef = 1
weakcost = 0.5
strongcost = 2.0

#################
# viz temp
#################

dtemp = dyn_temp(times[-int(365/delt):],tmin,tmax)
ax3[0].plot(times[-int(365/delt):],dtemp,lw=2)

rs = get_r(a,b,temps)
ax3[1].plot(temps,rs,lw=2)

#################
# viz trade-off function
#################

cors = np.linspace(0,1,100)
for (n,c,lab) in zip((weakcost,ndef,strongcost),('--','-.','-'),('High cost of defense','Medium cost of defense','Low cost of defense')):
    bors = corbor(cors,n)
    ax12.plot(cors,bors,ls=c,label=lab)
leg = ax12.legend(prop={'size':12,'family': 'Times New Roman'})
leg.draw_frame(False)

# trade-off options
cors = np.r_[[0.7,0.5]]
bors = np.r_[[0.5,0.7]]

ax12.plot([cors[0]],[bors[0]],'ro',label='Assumed high cost value',markersize=15)
ax12.plot([cors[1]],[bors[1]],'b*',label='Assumed low cost value',markersize=14)
leg = ax12.legend(prop={'size':12,'family': 'Times New Roman'})
leg.draw_frame(False)

#################
# diamond model
#################
Ns,P1s,P2s,Zs = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
for temp in temps:
    r = get_r(a,b,temp)
    params = get_params(r,cor=cors[0],bor=bors[0])
    inits = get_inits('diamond',Ntot)
    dtest = si.odeint(diamond,inits,times,args=params).T
    if temp == temps[ndyn]:
        dtestshow = dtest        
    Ns = np.append(Ns,dtest[0][-1])
    P1s = np.append(P1s,dtest[1][-1])
    P2s = np.append(P2s,dtest[2][-1])
    Zs = np.append(Zs,dtest[3][-1])

# show convergene to equilbrium for a special case
ax2[0].plot(times[:int(ndayss/delt)],dtestshow[0][:int(ndayss/delt)],ncp,label=mnl,lw=1.5)
ax2[1].plot(times[:int(ndayss/delt)],dtestshow[1][:int(ndayss/delt)],p1cp,label=mp1l,lw=1.5)
ax2[1].plot(times[:int(ndayss/delt)],dtestshow[2][:int(ndayss/delt)],p2cp,label=mp2l,lw=1.5)
ax2[2].plot(times[:int(ndayss/delt)],dtestshow[3][:int(ndayss/delt)],z2cp,label=mz1l,lw=1.5)

# plot equilibrium
ax4[0,0].plot(xdim,P1s,p1cp,label=mp1l,lw=3)
ax4[0,0].plot(xdim,P2s,p2cp,label=mp2l,lw=3)
ax4[0,0].plot(xdim,Ns,ncp,label=mnl,lw=3)
ax4[0,1].plot(xdim,Zs,zcp,label=mz1l,lw=3)

ax9[0,0].plot(times[:int(ndayss/delt)],dtestshow[1][:int(ndayss/delt)],p1cp,label=mp1l,lw=2)
ax9[0,0].plot(times[:int(ndayss/delt)],dtestshow[2][:int(ndayss/delt)],p2cp,label=mp2l,lw=2)
ax9[0,0].plot(times[:int(ndayss/delt)],dtestshow[0][:int(ndayss/delt)],ncp,label=mnl,lw=2)
ax9[0,1].plot(times[:int(ndayss/delt)],dtestshow[3][:int(ndayss/delt)],zcp,label=mz1l,lw=2)

tmin = 4

# plot dynamic - weak cost
u = np.r_[[Ns[-5],P1s[-5],P2s[-5],Zs[-5]]]
Ns,P1s,P2s,Zs = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
for t in timesl:
    tloc = dyn_temp(t,tmin,tmax)
    r = get_r(a,b,tloc)
    params = get_params(r,cor=cors[0],bor=bors[0])
    u = u + deltl*diamond(u,t,*params)
    Ns = np.append(Ns,u[0])
    P1s = np.append(P1s,u[1])
    P2s = np.append(P2s,u[2])
    Zs = np.append(Zs,u[3])
ax5[0].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Ns[int(-365*nyearsl/deltl)::freq]),nca,lw=1.5)
ax5[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P1s[int(-365*nyearsl/deltl)::freq]),p1ca,lw=1.5)
ax5[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P2s[int(-365*nyearsl/deltl)::freq]),p2ca,lw=1.5)
ax5[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Zs[int(-365*nyearsl/deltl)::freq]),zca,lw=1.5)
ax11[0,0].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P1s[int(-365*nyears/deltl)::freq]),p1ca,lw=2.5)
ax11[0,0].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P2s[int(-365*nyears/deltl)::freq]),p2ca,lw=2.5)
ax5[0].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Ns)[int(-365*nyearsl/deltl):],ncb,label=mnl,lw=1.5)
ax5[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P1s)[int(-365*nyearsl/deltl):],p1cb,label=mp1l,lw=1.5)
ax5[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P2s)[int(-365*nyearsl/deltl):],p2cb,label=mp2l,lw=1.5)
ax5[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Zs)[int(-365*nyearsl/deltl):],zcb,label=mz1l,lw=1.5)
ax11[0,0].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P1s)[int(-365*nyears/deltl):],p1cb,label=mp1l,lw=2.5)
ax11[0,0].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P2s)[int(-365*nyears/deltl):],p2cb,label=mp2l,lw=2.5)

# plot dynamic - strong cost
u = np.r_[[Ns[-5],P1s[-5],P2s[-5],Zs[-5]]]
Ns,P1s,P2s,Zs = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
for t in timesl:
    tloc = dyn_temp(t,tmin,tmax)
    r = get_r(a,b,tloc)
    params = get_params(r,cor=cors[1],bor=bors[1])
    u = u + deltl*diamond(u,t,*params)
    Ns = np.append(Ns,u[0])
    P1s = np.append(P1s,u[1])
    P2s = np.append(P2s,u[2])
    Zs = np.append(Zs,u[3])
ax6[0].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Ns[int(-365*nyearsl/deltl)::freq]),nca,lw=1.5)
ax6[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P1s[int(-365*nyearsl/deltl)::freq]),p1ca,lw=1.5)
ax6[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P2s[int(-365*nyearsl/deltl)::freq]),p2ca,lw=1.5)
ax6[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Zs[int(-365*nyearsl/deltl)::freq]),zca,lw=1.5)
ax11[0,1].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P1s[int(-365*nyears/deltl)::freq]),p1ca,lw=2.5)
ax11[0,1].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P2s[int(-365*nyears/deltl)::freq]),p2ca,lw=2.5)
ax6[0].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Ns)[int(-365*nyearsl/deltl):],ncb,label=mnl,lw=1.5)
ax6[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P1s)[int(-365*nyearsl/deltl):],p1cb,label=mp1l,lw=1.5)
ax6[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P2s)[int(-365*nyearsl/deltl):],p2cb,label=mp2l,lw=1.5)
ax6[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Zs)[int(-365*nyearsl/deltl):],zcb,label=mz1l,lw=1.5)
ax11[0,1].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P1s)[int(-365*nyears/deltl):],p1cb,label=mp1l,lw=2.5)
ax11[0,1].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P2s)[int(-365*nyears/deltl):],p2cb,label=mp2l,lw=2.5)


#################
# parallel model
#################
Ns,P1s,P2s,Z1s,Z2s = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
for temp in temps:
    r = get_r(a,b,temp)
    params = get_params(r,cor=cors[0],bor=bors[0])
    inits = get_inits('parallel',Ntot)
    dtest = si.odeint(parallel,inits,times,args=params).T
    if temp == temps[ndyn]:
        dtestshow = dtest
    Ns = np.append(Ns,dtest[0][-1])
    P1s = np.append(P1s,dtest[1][-1])
    P2s = np.append(P2s,dtest[2][-1])
    Z1s = np.append(Z1s,dtest[3][-1])
    Z2s = np.append(Z2s,dtest[4][-1])

# show convergene to equilibrium for special case
ax1[0].plot(times[:int(ndayss/delt)],dtestshow[0][:int(ndayss/delt)],ncp,label=mnl)
ax1[1].plot(times[:int(ndayss/delt)],dtestshow[1][:int(ndayss/delt)],p1cp,label=mp1l)
ax1[1].plot(times[:int(ndayss/delt)],dtestshow[2][:int(ndayss/delt)],p2cp,label=mp2l)
ax1[2].plot(times[:int(ndayss/delt)],dtestshow[3][:int(ndayss/delt)],zcp,label=mz1l)
ax1[2].plot(times[:int(ndayss/delt)],dtestshow[4][:int(ndayss/delt)],z2cp,label=mz2l)

# plot equilibrium
ax4[1,0].plot(xdim,P1s,p1cp,label=mp1l,lw=3)
ax4[1,0].plot(xdim,P2s,p2cp,label=mp2l,lw=3)
ax4[1,0].plot(xdim,Ns,ncp,label=mnl,lw=3)
ax4[1,1].plot(xdim,Z1s,zcp,label=mz1l,lw=3)
ax4[1,1].plot(xdim,Z2s,z2cp,label=mz2l,lw=3)

ax9[1,0].plot(times[:int(ndayss/delt)],dtestshow[1][:int(ndayss/delt)],p1cp,label=mp1l,lw=2)
ax9[1,0].plot(times[:int(ndayss/delt)],dtestshow[2][:int(ndayss/delt)],p2cp,label=mp2l,lw=2)
ax9[1,0].plot(times[:int(ndayss/delt)],dtestshow[0][:int(ndayss/delt)],ncp,label=mnl,lw=2)
ax9[1,1].plot(times[:int(ndayss/delt)],dtestshow[3][:int(ndayss/delt)],zcp,label=mz1l,lw=2)
ax9[1,1].plot(times[:int(ndayss/delt)],dtestshow[4][:int(ndayss/delt)],z2cp,label=mz2l,lw=2)

# plot dynamic - weak cost
u = np.r_[[Ns[-5],P1s[-5],P2s[-5],Z1s[-5],Z2s[-5]]]
Ns,P1s,P2s,Z1s,Z2s = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
for t in timesl:
    tloc = dyn_temp(t,tmin,tmax)
    r = get_r(a,b,tloc)
    params = get_params(r,cor=cors[0],bor=bors[0])
    u = u + deltl*parallel(u,t,*params)
    Ns = np.append(Ns,u[0])
    P1s = np.append(P1s,u[1])
    P2s = np.append(P2s,u[2])
    Z1s = np.append(Z1s,u[3])
    Z2s = np.append(Z2s,u[4])
ax7[0].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Ns[int(-365*nyearsl/deltl)::freq]),nca,lw=1.5)
ax7[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P1s[int(-365*nyearsl/deltl)::freq]),p1ca,lw=1.5)
ax7[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P2s[int(-365*nyearsl/deltl)::freq]),p2ca,lw=1.5)
ax7[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Z1s[int(-365*nyearsl/deltl)::freq]),zca,lw=1.5)
ax7[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Z2s[int(-365*nyearsl/deltl)::freq]),z2ca,lw=1.5)
ax11[1,0].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P1s[int(-365*nyears/deltl)::freq]),p1ca,lw=2.5)
ax11[1,0].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P2s[int(-365*nyears/deltl)::freq]),p2ca,lw=2.5)
ax7[0].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Ns)[int(-365*nyearsl/deltl):],ncb,label=mnl,lw=1.5)
ax7[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P1s)[int(-365*nyearsl/deltl):],p1cb,label=mp1l,lw=1.5)
ax7[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P2s)[int(-365*nyearsl/deltl):],p2cb,label=mp2l,lw=1.5)
ax7[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Z1s)[int(-365*nyearsl/deltl):],zcb,label=mz1l,lw=1.5)
ax7[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Z2s)[int(-365*nyearsl/deltl):],z2cb,label=mz2l,lw=1.5)
ax11[1,0].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P1s)[int(-365*nyears/deltl):],p1cb,label=mp1l,lw=2.5)
ax11[1,0].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P2s)[int(-365*nyears/deltl):],p2cb,label=mp2l,lw=2.5)

# plot dynamic - strong cost
u = np.r_[[Ns[-5],P1s[-5],P2s[-5],Z1s[-5],Z2s[-5]]]
Ns,P1s,P2s,Z1s,Z2s = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
for t in timesl:
    tloc = dyn_temp(t,tmin,tmax)
    r = get_r(a,b,tloc)
    params = get_params(r,cor=cors[1],bor=bors[1])
    u = u + deltl*parallel(u,t,*params)
    Ns = np.append(Ns,u[0])
    P1s = np.append(P1s,u[1])
    P2s = np.append(P2s,u[2])
    Z1s = np.append(Z1s,u[3])
    Z2s = np.append(Z2s,u[4])
ax8[0].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Ns[int(-365*nyearsl/deltl)::freq]),nca,lw=1.5)
ax8[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P1s[int(-365*nyearsl/deltl)::freq]),p1ca,lw=1.5)
ax8[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P2s[int(-365*nyearsl/deltl)::freq]),p2ca,lw=1.5)
ax8[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Z1s[int(-365*nyearsl/deltl)::freq]),zca,lw=1.5)
ax8[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Z2s[int(-365*nyearsl/deltl)::freq]),z2ca,lw=1.5)
ax11[1,1].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P1s[int(-365*nyears/deltl)::freq]),p1ca,lw=2.5)
ax11[1,1].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P2s[int(-365*nyears/deltl)::freq]),p2ca,lw=2.5)
ax8[0].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Ns)[int(-365*nyearsl/deltl):],ncb,label=mnl,lw=1.5)
ax8[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P1s)[int(-365*nyearsl/deltl):],p1cb,label=mp1l,lw=1.5)
ax8[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P2s)[int(-365*nyearsl/deltl):],p2cb,label=mp2l,lw=1.5)
ax8[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Z1s)[int(-365*nyearsl/deltl):],zcb,label=mz1l,lw=1.5)
ax8[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Z2s)[int(-365*nyearsl/deltl):],z2cb,label=mz2l,lw=1.5)
ax11[1,1].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P1s)[int(-365*nyears/deltl):],p1cb,label=mp1l,lw=2.5)
ax11[1,1].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P2s)[int(-365*nyears/deltl):],p2cb,label=mp2l,lw=2.5)

#ax1[1].semilogy()

################################################
# final fig formatting
################################################

for axes in [ax1,ax2,ax5,ax6,ax7,ax8]:
    for a in axes.flatten()[1:]:
        l = a.legend(prop={'size':18,'family': 'Times New Roman'})
        l.draw_frame(False)

for axes in [ax4,ax9,ax11]:
    for a in axes.flatten():
        l = a.legend(prop={'size':18,'family': 'Times New Roman'})
        l.draw_frame(False)

################################################
# save figs
################################################

f1.savefig('../figures/tdep_parallel', bbox_inches='tight',dpi=300)
f2.savefig('../figures/tdep_diamond', bbox_inches='tight',dpi=300)
f3.savefig('../figures/temp_response', bbox_inches='tight',dpi=300)
f4.savefig('../figures/equilibrium',  bbox_inches='tight',dpi=300)
f5.savefig('../figures/alldyn1', bbox_inches='tight',dpi=300)
f6.savefig('../figures/alldyn2', bbox_inches='tight',dpi=300)
f7.savefig('../figures/alldyn3', bbox_inches='tight',dpi=300)
f8.savefig('../figures/alldyn4',  bbox_inches='tight',dpi=300)
f9.savefig('../figures/t_dep_equil', bbox_inches='tight',dpi=300)
f10.savefig('../figures/median_dat', bbox_inches='tight',dpi=300)
f11.savefig('../figures/phyto_dyn', bbox_inches='tight',dpi=300)
f12.savefig('../figures/cor_bor',  bbox_inches='tight',dpi=300)

plt.close(f1)
plt.close(f2)
plt.close(f3)
plt.close(f4)
plt.close(f5)
plt.close(f6)
plt.close(f7)
plt.close(f8)
plt.close(f9)
plt.close(f10)
plt.close(f11)
plt.close(f12)
