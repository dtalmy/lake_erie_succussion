import numpy as np
import pandas as pd
import pylab as plt
import seaborn as sns
import models as mds
import sys
import scipy.integrate as si

def get_inits(d,Ntot,NtoP=40):
    if d == 'parallel':
        Ns = np.r_[[20,10,1,0.1,0.1]]
        return Ns
    if d == 'diamond':
        Ns = np.r_[[20,10,1,0.2]]
        return Ns

################################################
# make figures
################################################

# font size
fs=16

# make all font Times New Roman
ex = {'font.family': "Times New Roman",'font.size': 14}
plt.rcParams.update(**ex)
#sns.set_theme(rc=ex)

# colors and line styles
pc,nc,p1c,p2c,zc,z2c = 'bo--','bo-','g^--','g*-','ro--','ro-' # for intermittent model samples
pcp,ncp,p1cp,p2cp,zcp,z2cp = 'b--','b','g--','g','r--','r' # for continious plots
pca,nca,p1ca,p2ca,zca,z2ca = 'bo--','bo','g^','g*','ro','ro' # for intermittent model samples
pcb,ncb,p1cb,p2cb,zcb,z2cb = 'b--','b-','g--','g-','r--','r-' # for intermittent model samples

# labels
p1lab,p2lab = 'Non-cyanobacteria chlorophyll','Phycocyanin'

mnl,mpl,mp1l,mp2l,mz1l,mz2l = r'$R_N$',r'$R_P$',r'$P_{1}$',r'$P_{2}$',r'$Z_{1}$',r'$Z_{2}$'

f1,ax1 = plt.subplots(3,1,figsize=[6,12])
f2,ax2 = plt.subplots(3,1,figsize=[6,12])
f3,ax3 = plt.subplots(1,2,figsize=[9,4])
f4,ax4 = plt.subplots(2,3,figsize=[18,12])
f5,ax5 = plt.subplots(3,1,figsize=[6,12])
f7,ax7 = plt.subplots(3,1,figsize=[6,12])
f9,ax9 = plt.subplots(2,3,figsize=[18,12])
f10,ax10 = plt.subplots(1,3,figsize=[18,6])
f11,ax11 = plt.subplots(4,2,figsize=[15,18])
        
for (a,lab) in zip(ax4.flatten(),'abcdef'):
    a.set_xlabel('Temperature (Celcius)',fontsize=fs)
    a.set_ylabel('Nitrogen (mmol N m$^{-3}$)',fontsize=fs)
    a.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=a.transAxes)

for (a,lab) in zip(ax9.flatten(),'abcdef'):
    a.set_xlabel('Time (days)',fontsize=fs)
    a.set_ylabel('Nitrogen (mmol N m$^{-3}$)',fontsize=fs)
    a.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=a.transAxes)

ax4 = ax4[:,[1,2,0]]
ax9 = ax9[:,[1,2,0]]

ax4[1,0].text(0.2,1.08,'Parallel model',fontsize=20,transform=ax4[1,0].transAxes)
ax4[0,0].text(0.1,1.08,'Shared predation model',fontsize=20,transform=ax4[0,0].transAxes)
f4.subplots_adjust(hspace=0.4,wspace=0.4)

ax9[1,0].text(0.2,1.08,'Parallel model',fontsize=20,transform=ax9[1,0].transAxes)
ax9[0,0].text(0.1,1.08,'Shared predation model',fontsize=20,transform=ax9[0,0].transAxes)
f9.subplots_adjust(hspace=0.32,wspace=0.4)

ax4[0,0].semilogy()
ax4[1,0].semilogy()
ax4[0,1].semilogy()
ax4[1,1].semilogy()

ax9[0,0].semilogy()
ax9[1,0].semilogy()
ax9[0,1].semilogy()
ax9[1,1].semilogy()

for (axes,f) in zip([ax1,ax2],[f1,f2]):
    f.subplots_adjust(hspace=0.3)
    for (lab,ax,tit) in zip('abc',axes,['Nutrient','Phytoplankton','Zooplankton']):
        ax.set_xlabel('Time (days)')
        ax.set_ylabel(r'log mmol N m$^{-3}$')
        #ax.set_title(tit)
        ax.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=ax.transAxes)

for (axes,f) in zip([ax5,ax7],[f5,f7]):
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

# percentiles for data plots
pmin,pmax=5,95

# rsquareds
rsquared_parallel = pd.read_csv('../data/rsquared_parallel.csv')
rsquared_diamond = pd.read_csv('../data/rsquared_diamond.csv')

################################################
# viz data
################################################

import taihu_data
import viz_data as vd
import zoo_viz as zv

lt2 = ax10[0].plot(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.phosname,50),'b^--',linewidth = 3,label='Phosphorus')
ax10[0].fill_between(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.phosname,pmin),\
        vd.get_percentiles(vd.sds,vd.phosname,pmax),color='b',alpha=0.1)
lns = lt2
labs = [l.get_label() for l in lns]
leg = ax10[0].legend(lns, labs, prop={'size':12,'family': 'Times New Roman'},fontsize=16,loc='best')
leg.draw_frame(False)

ax10[1].plot(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.chlname,50),p1c,label=p1lab,linewidth = 3)
ax10[1].plot(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.pycname,50),p2c,label=p2lab,linewidth = 3)
ax10[1].fill_between(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.chlname,pmin),\
        vd.get_percentiles(vd.sds,vd.chlname,pmax),color='g',alpha=0.1)
ax10[1].fill_between(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.pycname,pmin),\
        vd.get_percentiles(vd.sds,vd.pycname,pmax),color='g',alpha=0.3)


g = sns.boxplot(x = zv.lez[zv.lez['YEAR']>2011]['Month'],\
        y = np.log(zv.lez[zv.lez['YEAR']>2011]['total_zooplankton_biomass']),\
        hue = zv.lez[zv.lez['YEAR']>2011]['Month'], ax=ax10[2])
ax10[2].set_ylabel('log Total Zooplankton Biomass\n (mmol C m$^{-3}$)',fontsize=fs)
leg = ax10[2].legend(prop={'size':12,'family': 'Times New Roman'})
new_labels = ['April', 'August']
for t, lab in zip(leg.texts, new_labels):
    t.set_text(lab)
leg = ax10[2].legend(prop={'size':12,'family': 'Times New Roman'},fontsize=16,loc='best')

ax10[0].set_ylabel('log '+vd.phosnamed,fontsize=fs)
ax10[1].set_ylabel('log Pigment (mmol N m$^{-3}$)',fontsize=fs)

################################################
#  dump medians to file for fitting elsewhere
################################################

# phyt data medians
nit_med = vd.get_percentiles(vd.sds,vd.nitname,50)
phos_med = vd.get_percentiles(vd.sds,vd.phosname,50)
phyt_months = vd.get_months(vd.sds).astype(int)
chl_med = vd.get_percentiles(vd.sds,vd.chlname,50)
phyc_med = vd.get_percentiles(vd.sds,vd.pycname,50)

# zoo data
zoo_months = vd.get_months(zv.lez).astype(int)
zoo_med = vd.get_percentiles(zv.lez,'total_zooplankton_biomass',50) 

# state variable names
nit_org = np.r_[['N' for i in nit_med]]
phos_org = np.r_[['P' for i in phos_med]]
chl_org = np.r_[['P1' for i in chl_med]]
phyc_org = np.r_[['P2' for i in phyc_med]]
zoo_org = np.r_[['Z' for i in zoo_med]]

# format for dict
organism = np.concatenate((nit_org,phos_org,chl_org,phyc_org,zoo_org))
log_abundance = np.concatenate((nit_med,phos_med,chl_med,phyc_med,zoo_med))
time = np.concatenate((phyt_months,phyt_months,phyt_months,phyt_months,zoo_months))
log_sigma = np.r_[[0.6 for i in time]]

# dump
for_fitting = pd.DataFrame({'organism':organism,'log_abundance':log_abundance,'time':time,'log_sigma':log_sigma})
for_fitting.to_csv('../data/erie_tseries_all.csv')

################################################
#  show 5th and 95th quantiles for model-data comparison
################################################

# phyt data 
nit_95 = vd.get_percentiles(vd.sds,vd.nitname,95)
phos_95 = vd.get_percentiles(vd.sds,vd.phosname,95)
chl_95 = vd.get_percentiles(vd.sds,vd.chlname,95)
phyc_95 = vd.get_percentiles(vd.sds,vd.pycname,95)
nit_5 = vd.get_percentiles(vd.sds,vd.nitname,5)
phos_5 = vd.get_percentiles(vd.sds,vd.phosname,5)
chl_5 = vd.get_percentiles(vd.sds,vd.chlname,5)
phyc_5 = vd.get_percentiles(vd.sds,vd.pycname,5)

# zoo data
zoo_95 = vd.get_percentiles(zv.lez,'total_zooplankton_biomass',95) 
zoo_5 = vd.get_percentiles(zv.lez,'total_zooplankton_biomass',5) 

ax11[0,1].fill_between(phyt_months,phos_5,phos_95,color='b',alpha=0.1)
ax11[1,1].fill_between(phyt_months,chl_5,chl_95,color='g',alpha=0.1)
ax11[2,1].fill_between(phyt_months,phyc_5,phyc_95,color='g',alpha=0.1)
ax11[3,1].fill_between(zoo_months,zoo_5,zoo_95,color='r',alpha=0.1)

ax11[0,0].fill_between(phyt_months,phos_5,phos_95,color='b',alpha=0.1)
ax11[1,0].fill_between(phyt_months,chl_5,chl_95,color='g',alpha=0.1)
ax11[2,0].fill_between(phyt_months,phyc_5,phyc_95,color='g',alpha=0.1)
ax11[3,0].fill_between(zoo_months,zoo_5,zoo_95,color='r',alpha=0.1)

for (a1,a2,r1,r2) in zip(ax11[:,0],ax11[:,1],rsquared_diamond.values[0][1:],rsquared_parallel.values[0][1:]):
    a1.text(0.8,0.9,r'R$^2$='+str(np.round(r1,2)),ha='center',va='center',color='k',transform=a1.transAxes,fontsize=fs)
    a2.text(0.8,0.9,r'R$^2$='+str(np.round(r2,2)),ha='center',va='center',color='k',transform=a2.transAxes,fontsize=fs)

################################################
# run model
################################################

# time array
ndays,delt = 2000,0.1
times = np.linspace(0,ndays,int(ndays/delt))
ndayss = 2000

# frequency for symbol plotting of model
freq = 1000

# shorter time array for dyn sims
ndaysl,deltl = 365*7,0.05
timesl = np.linspace(0,ndaysl,int(ndaysl/deltl))
fyd1 = ndaysl - 365*6

# indices for single year plots
ayear = 7300
offset = int(0/12*ayear)
imax = timesl.shape[0] - offset
imin = imax - ayear

# nutrient array
Ntotmin,Ntotmax = 0.003,0.6
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
b = 0.06 # eppley exponent

# number years to show
nyears,nyearsl = 1,3
ncut = (ndaysl/365-nyears)*365
ncutl = (ndaysl/365-nyearsl)*365

# iteration to visualize time dependent dynamics
ndyn = 25

# default trade-off valis
ndef = 1
weakcost = 0.5
strongcost = 2.0

#################
# viz temp
#################

dtemp = mds.dyn_temp(timesl[-int(365/deltl):],tmin,tmax)
ax3[0].plot((timesl[-int(365/deltl):]-(365*5)),dtemp,lw=2)

rs = mds.get_r(temps)
ax3[1].plot(temps,rs,lw=2)

btemp = pd.read_csv('../data/burns_temp.csv')
ax3[0].scatter(btemp.month*30,btemp.temp)

#################
# diamond model
#################
Nns,P1ns,P2ns,Zns = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
for temp in temps:
    inits = get_inits('diamond',Ntot)
    params = pd.read_csv('../data/chain_inits_default.csv')
    params['temp'] = temp
    params['Tdyn'] = 0
    dtest = si.odeint(mds.diamond,inits,times,(tuple(params.values[0][1:]),)).T
    if temp == temps[ndyn]:
        dtestshow = dtest        
    Nns = np.append(Nns,dtest[0][-1])
    P1ns = np.append(P1ns,dtest[1][-1])
    P2ns = np.append(P2ns,dtest[2][-1])
    Zns = np.append(Zns,dtest[3][-1])

# show convergence to equilbrium for a special case
ax2[0].plot(times[:int(ndayss/delt)],dtestshow[0][:int(ndayss/delt)],ncp,label=mnl,lw=1.5)
ax2[1].plot(times[:int(ndayss/delt)],dtestshow[1][:int(ndayss/delt)],p1cp,label=mp1l,lw=1.5)
ax2[1].plot(times[:int(ndayss/delt)],dtestshow[2][:int(ndayss/delt)],p2cp,label=mp2l,lw=1.5)
ax2[2].plot(times[:int(ndayss/delt)],dtestshow[3][:int(ndayss/delt)],z2cp,label=mz1l,lw=1.5)

# plot equilibrium
ax4[0,0].plot(xdim,P1ns,p1cp,label=mp1l,lw=3)
ax4[0,0].plot(xdim,P2ns,p2cp,label=mp2l,lw=3)
ax4[0,2].plot(xdim,Nns,ncp,label=mnl,lw=3)
ax4[0,1].plot(xdim,Zns,zcp,label=mz1l,lw=3)

ax9[0,0].plot(times[:int(ndayss/delt)],dtestshow[1][:int(ndayss/delt)],p1cp,label=mp1l,lw=2)
ax9[0,0].plot(times[:int(ndayss/delt)],dtestshow[2][:int(ndayss/delt)],p2cp,label=mp2l,lw=2)
ax9[0,2].plot(times[:int(ndayss/delt)],dtestshow[0][:int(ndayss/delt)],ncp,label=mnl,lw=2)
ax9[0,1].plot(times[:int(ndayss/delt)],dtestshow[3][:int(ndayss/delt)],zcp,label=mz1l,lw=2)

tmin = 4

# plot dynamic
params = pd.read_csv('../data/chain_inits_diamond_test.csv')
params['Tdyn'] = True
inits = get_inits('diamond',Ntot)
u = si.odeint(mds.diamond,inits,timesl,(tuple(params.values[0][1:]),)).T
Nns,P1ns,P2ns,Zns = u[0],u[1],u[2],u[3]

# multi-year plots for SI
ax5[0].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Nns[int(-365*nyearsl/deltl)::freq]),nca,lw=1.5)
ax5[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P1ns[int(-365*nyearsl/deltl)::freq]),p1ca,lw=1.5)
ax5[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P2ns[int(-365*nyearsl/deltl)::freq]),p2ca,lw=1.5)
ax5[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Zns[int(-365*nyearsl/deltl)::freq]),zca,lw=1.5)
ax5[0].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Nns)[int(-365*nyearsl/deltl):],ncb,label=mnl,lw=1.5)
ax5[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P1ns)[int(-365*nyearsl/deltl):],p1cb,label=mp1l,lw=1.5)
ax5[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P2ns)[int(-365*nyearsl/deltl):],p2cb,label=mp2l,lw=1.5)
ax5[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Zns)[int(-365*nyearsl/deltl):],zcb,label=mz1l,lw=1.5)

# single year plots for main text
ax11[0,0].plot((timesl[imin:imax]-365*6)/30.5,np.log(Nns[imin:imax]),'b',label=mnl,lw=2)
ax11[1,0].plot((timesl[imin:imax]-365*6)/30.5,np.log(P1ns[imin:imax]),'g',label=mp1l,lw=2)
ax11[2,0].plot((timesl[imin:imax]-365*6)/30.5,np.log(P2ns[imin:imax]),'g--',label=mp2l,lw=2)
ax11[3,0].plot((timesl[imin:imax]-365*6)/30.5,np.log(Zns[imin:imax]),'r',label=mz1l,lw=2)

#################
# parallel model
#################
Nns,P1ns,P2ns,Z1ns,Z2ns = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
for temp in temps:
    inits = get_inits('parallel',Ntot)
    params = pd.read_csv('../data/chain_inits_default.csv')
    params['temp'] = temp
    params['Tdyn'] = 0
    dtest = si.odeint(mds.parallel,inits,times,(tuple(params.values[0][1:]),)).T
    if temp == temps[ndyn]:
        dtestshow = dtest
    Nns = np.append(Nns,dtest[0][-1])
    P1ns = np.append(P1ns,dtest[1][-1])
    P2ns = np.append(P2ns,dtest[2][-1])
    Z1ns = np.append(Z1ns,dtest[3][-1])
    Z2ns = np.append(Z2ns,dtest[4][-1])

# show convergene to equilibrium for special case
ax1[0].plot(times[:int(ndayss/delt)],dtestshow[0][:int(ndayss/delt)],ncp,label=mnl)
ax1[1].plot(times[:int(ndayss/delt)],dtestshow[1][:int(ndayss/delt)],p1cp,label=mp1l)
ax1[1].plot(times[:int(ndayss/delt)],dtestshow[2][:int(ndayss/delt)],p2cp,label=mp2l)
ax1[2].plot(times[:int(ndayss/delt)],dtestshow[3][:int(ndayss/delt)],zcp,label=mz1l)
ax1[2].plot(times[:int(ndayss/delt)],dtestshow[4][:int(ndayss/delt)],z2cp,label=mz2l)

# plot equilibrium
ax4[1,0].plot(xdim,P1ns,p1cp,label=mp1l,lw=3)
ax4[1,0].plot(xdim,P2ns,p2cp,label=mp2l,lw=3)
ax4[1,2].plot(xdim,Nns,ncp,label=mnl,lw=3)
ax4[1,1].plot(xdim,Z1ns,zcp,label=mz1l,lw=3)
ax4[1,1].plot(xdim,Z2ns,z2cp,label=mz2l,lw=3)

ax9[1,0].plot(times[:int(ndayss/delt)],dtestshow[1][:int(ndayss/delt)],p1cp,label=mp1l,lw=2)
ax9[1,0].plot(times[:int(ndayss/delt)],dtestshow[2][:int(ndayss/delt)],p2cp,label=mp2l,lw=2)
ax9[1,2].plot(times[:int(ndayss/delt)],dtestshow[0][:int(ndayss/delt)],ncp,label=mnl,lw=2)
ax9[1,1].plot(times[:int(ndayss/delt)],dtestshow[3][:int(ndayss/delt)],zcp,label=mz1l,lw=2)
ax9[1,1].plot(times[:int(ndayss/delt)],dtestshow[4][:int(ndayss/delt)],z2cp,label=mz2l,lw=2)

# plot dynamic
params = pd.read_csv('../data/chain_inits_parallel_test.csv')
params['Tdyn'] = True
inits = get_inits('parallel',Ntot)
u = si.odeint(mds.parallel,inits,timesl,(tuple(params.values[0][1:]),)).T
Nns,P1ns,P2ns,Z1ns,Z2ns = u[0],u[1],u[2],u[3],u[4]
Zns = Z1ns+Z2ns

ax7[0].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Nns[int(-365*nyearsl/deltl)::freq]),nca,lw=1.5)
ax7[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P1ns[int(-365*nyearsl/deltl)::freq]),p1ca,lw=1.5)
ax7[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P2ns[int(-365*nyearsl/deltl)::freq]),p2ca,lw=1.5)
ax7[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Z1ns[int(-365*nyearsl/deltl)::freq]),zca,lw=1.5)
ax7[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Z2ns[int(-365*nyearsl/deltl)::freq]),z2ca,lw=1.5)
ax7[0].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Nns)[int(-365*nyearsl/deltl):],ncb,label=mnl,lw=1.5)
ax7[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P1ns)[int(-365*nyearsl/deltl):],p1cb,label=mp1l,lw=1.5)
ax7[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P2ns)[int(-365*nyearsl/deltl):],p2cb,label=mp2l,lw=1.5)
ax7[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Z1ns)[int(-365*nyearsl/deltl):],zcb,label=mz1l,lw=1.5)
ax7[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Z2ns)[int(-365*nyearsl/deltl):],z2cb,label=mz2l,lw=1.5)

# single year plots for main text
ax11[0,1].plot((timesl[imin:imax]-365*6)/30.5,np.log(Nns[imin:imax]),'b',label=mnl,lw=2)
ax11[1,1].plot((timesl[imin:imax]-365*6)/30.5,np.log(P1ns[imin:imax]),'g',label=mp1l,lw=2)
ax11[2,1].plot((timesl[imin:imax]-365*6)/30.5,np.log(P2ns[imin:imax]),'g--',label=mp2l,lw=2)
ax11[3,1].plot((timesl[imin:imax]-365*6)/30.5,np.log(Zns[imin:imax]),'r',label=mz1l,lw=2)

################################################
# final fig formatting
################################################

for axes in [ax1,ax2,ax5,ax7]:
    for a in axes.flatten()[1:]:
        l = a.legend(prop={'size':18,'family': 'Times New Roman'})
        l.draw_frame(False)

for axes in [ax4[:,:2],ax9[:,:2]]:
    for a in axes.flatten():
        l = a.legend(prop={'size':18,'family': 'Times New Roman'})
        l.draw_frame(False)

for a in [ax4[0,2],ax4[1,2],ax9[0,2],ax9[1,2]]:
    a.set_ylabel('Phosphorus (mmol P m$^{-3}$)',fontsize=fs)

for (a,lab) in zip(ax11.flatten(),'abcdefgh'):
    a.set_xlabel('Month',fontsize=fs)
    a.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=a.transAxes,fontsize=fs)
ax11[0,0].set_ylabel('log Dissolved P (mmol P m$^{-3}$)',fontsize=fs)
ax11[1,0].set_ylabel('log Phytoplankton N (mmol N m$^{-3}$)',fontsize=fs)
ax11[2,0].set_ylabel('log Phytoplankton N (mmol N m$^{-3}$)',fontsize=fs)
ax11[3,0].set_ylabel('log Zooplankton N (mmol N m$^{-3}$)',fontsize=fs)

################################################
# save figs
################################################

f1.savefig('../figures/tdep_parallel', bbox_inches='tight',dpi=300)
f2.savefig('../figures/tdep_diamond', bbox_inches='tight',dpi=300)
f3.savefig('../figures/ltemp_response', bbox_inches='tight',dpi=300)
f4.savefig('../figures/equilibrium',  bbox_inches='tight',dpi=300)
f5.savefig('../figures/alldyn1', bbox_inches='tight',dpi=300)
f7.savefig('../figures/alldyn3', bbox_inches='tight',dpi=300)
f9.savefig('../figures/t_dep_equil', bbox_inches='tight',dpi=300)
f10.savefig('../figures/median_dat', bbox_inches='tight',dpi=300)
f11.savefig('../figures/model_data_comp',bbox_inches='tight',dpi=300)

plt.close(f1)
plt.close(f2)
plt.close(f3)
plt.close(f4)
plt.close(f5)
plt.close(f7)
plt.close(f9)
plt.close(f10)
plt.close(f11)
