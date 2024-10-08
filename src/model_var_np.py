import numpy as np
import pandas as pd
import pylab as plt
import seaborn as sns
import models as mds
import sys
import scipy.integrate as si

def get_inits(d,Ntot,NtoP=40):
    fnut = 0.1
    felse = 1 - fnut
    if d == 'parallel':
        Ns = np.r_[[fnut,0.25,0.25,0.25,0.25]]*Ntot
        Ps = Ns / NtoP
        return np.append(np.append(Ns,Ps),np.r_[[0,0]])
    if d == 'diamond':
        Ns = np.r_[[fnut,felse/4.0,felse/4.0,felse/2.0]]*Ntot
        Ps = Ns / NtoP
        return np.append(np.append(Ns,Ps),np.r_[[0,0]])

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

mnl,mpl,mp1l,mp2l,mz1l,mz2l = r'$R_N$',r'$R_P$',r'$P_{1,N}$',r'$P_{2,N}$',r'$Z_{1,N}$',r'$Z_{2,N}$'

f1,ax1 = plt.subplots(3,1,figsize=[6,12])
f2,ax2 = plt.subplots(3,1,figsize=[6,12])
f3,ax3 = plt.subplots(1,2,figsize=[9,4])
f4,ax4 = plt.subplots(2,3,figsize=[18,12])
f5,ax5 = plt.subplots(3,1,figsize=[6,12])
f6,ax6 = plt.subplots(3,1,figsize=[6,12])
f7,ax7 = plt.subplots(3,1,figsize=[6,12])
f8,ax8 = plt.subplots(3,1,figsize=[6,12])
f9,ax9 = plt.subplots(2,3,figsize=[18,12])
f10,ax10 = plt.subplots(2,2,figsize=[15,12])
f11,ax11 = plt.subplots(2,2,figsize=[12,12])
f12,ax12 = plt.subplots()
f13,ax13 = plt.subplots(2,2,figsize=[12,12])
f14,ax14 = plt.subplots(2,2,figsize=[12,12])

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

ax3 = ax3[::-1]

ax4[1,0].text(0.2,1.08,'Parallel model',fontsize=20,transform=ax4[1,0].transAxes)
ax4[0,0].text(0.1,1.08,'Shared predation model',fontsize=20,transform=ax4[0,0].transAxes)
f4.subplots_adjust(hspace=0.4,wspace=0.4)

ax9[1,0].text(0.2,1.08,'Parallel model',fontsize=20,transform=ax9[1,0].transAxes)
ax9[0,0].text(0.1,1.08,'Shared predation model',fontsize=20,transform=ax9[0,0].transAxes)
f9.subplots_adjust(hspace=0.32,wspace=0.4)

ax11[0,0].text(0.2,1.1,'High cost of resistance',fontsize=20,transform=ax11[0,0].transAxes)
ax11[0,1].text(0.2,1.1,'Low cost of resistance',fontsize=20,transform=ax11[0,1].transAxes)
ax11[0,1].text(1.1,0.3,'Shared predation model',fontsize=20,transform=ax11[0,1].transAxes,rotation=-90)
ax11[1,1].text(1.1,0.3,'Parallel model',fontsize=20,transform=ax11[1,1].transAxes,rotation=-90)

ax13[0,0].text(0.2,1.1,'High cost of resistance',fontsize=20,transform=ax13[0,0].transAxes)
ax13[0,1].text(0.2,1.1,'Low cost of resistance',fontsize=20,transform=ax13[0,1].transAxes)
ax13[0,1].text(1.3,0.3,'Shared predation model',fontsize=20,transform=ax13[0,1].transAxes,rotation=-90)
ax13[1,1].text(1.3,0.3,'Parallel model',fontsize=20,transform=ax13[1,1].transAxes,rotation=-90)

ax14[0,0].text(0.2,1.1,'High cost of resistance',fontsize=20,transform=ax14[0,0].transAxes)
ax14[0,1].text(0.2,1.1,'Low cost of resistance',fontsize=20,transform=ax14[0,1].transAxes)
ax14[0,1].text(1.1,0.3,'Shared predation model',fontsize=20,transform=ax14[0,1].transAxes,rotation=-90)
ax14[1,1].text(1.1,0.3,'Parallel model',fontsize=20,transform=ax14[1,1].transAxes,rotation=-90)

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
f10.subplots_adjust(wspace=0.5)
f13.subplots_adjust(wspace=0.4)

ax3[0].set_xlabel('Time (days)')
ax3[0].set_ylabel('Temperature (Celcius)')
ax3[1].set_xlabel('Temperature (Celcius)')
ax3[1].set_ylabel('Temperature response')
for (lab,ax) in zip('ba',ax3):
    ax.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=ax.transAxes)

for (a,lab,b) in zip(ax10.flatten(),'abcd',ax14.flatten()):
    a.set_xlabel('Month',fontsize=fs)
    a.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=a.transAxes)
    b.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=b.transAxes)
for (a,lab) in zip(ax11.flatten(),'abcd'):
    a.set_xlabel('Time (days)',fontsize=fs)
    a.set_ylabel('log Nitrogen (mmol N m$^{-3}$)',fontsize=fs)
    a.text(0.07,0.9,lab,ha='center',va='center',color='k',transform=a.transAxes)
for (a,lab) in zip(ax13.flatten(),'abcd'):
    a.set_xlabel('Month',fontsize=fs)
for (a,lab) in zip(ax14.flatten(),'abcd'):
    a.set_xlabel('Month',fontsize=fs)

ax10[0,1].set_ylabel('Molar N:P ratio')
ax12.set_ylabel('Defense ability')
ax12.set_xlabel('Growth ability')

# twin axes
ax10twin = ax10[0,0].twinx()
ax13t1 = ax13[0,0].twinx()
ax13t2 = ax13[0,1].twinx()
ax13t3 = ax13[1,0].twinx()
ax13t4 = ax13[1,1].twinx()

ax4t1 = ax4[0,2].twinx()
ax4t2 = ax4[1,2].twinx()
ax9t1 = ax9[0,2].twinx()
ax9t2 = ax9[1,2].twinx()

# percentiles for data plots
pmin,pmax=25,75

################################################
# viz data
################################################

import taihu_data
import viz_data as vd
import zoo_viz as zv

lt1 = ax10[0,0].plot(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.nitname,50),nc,linewidth = 3,label='Nitrogen')
ax10[0,0].fill_between(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.nitname,pmin),\
        vd.get_percentiles(vd.sds,vd.nitname,pmax),color='b',alpha=0.3)
lt2 = ax10twin.plot(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.phosname,50),'b^--',linewidth = 3,label='Phosphorus')
ax10twin.fill_between(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.phosname,pmin),\
        vd.get_percentiles(vd.sds,vd.phosname,pmax),color='b',alpha=0.1)
lns = lt1+lt2
labs = [l.get_label() for l in lns]
leg = ax10twin.legend(lns, labs, prop={'size':12,'family': 'Times New Roman'},fontsize=16,loc='best')
leg.draw_frame(False)

ax10[0,1].plot(vd.get_months(vd.sds),np.exp(vd.get_percentiles(vd.sds,vd.ntopname,50)),label=p1lab,linewidth = 3)
ax10[0,1].fill_between(vd.get_months(vd.sds),np.exp(vd.get_percentiles(vd.sds,vd.ntopname,pmin)),\
        np.exp(vd.get_percentiles(vd.sds,vd.ntopname,pmax)),color='b',alpha=0.2)

ax10[1,0].plot(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.chlname,50),p1c,label=p1lab,linewidth = 3)
ax10[1,0].plot(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.pycname,50),p2c,label=p2lab,linewidth = 3)
ax10[1,0].fill_between(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.chlname,pmin),\
        vd.get_percentiles(vd.sds,vd.chlname,pmax),color='g',alpha=0.1)
ax10[1,0].fill_between(vd.get_months(vd.sds),vd.get_percentiles(vd.sds,vd.pycname,pmin),\
        vd.get_percentiles(vd.sds,vd.pycname,pmax),color='g',alpha=0.3)


g = sns.boxplot(x = zv.lez[zv.lez['YEAR']>2011]['Month'],\
        y = np.log(zv.lez[zv.lez['YEAR']>2011]['total_zooplankton_biomass']),\
        hue = zv.lez[zv.lez['YEAR']>2011]['Month'], ax=ax10[1,1])
ax10[1,1].set_ylabel('log Total Zooplankton Biomass\n (mg m$^{-3}$)',fontsize=fs)
leg = ax10[1,1].legend(prop={'size':12,'family': 'Times New Roman'})
new_labels = ['April', 'August']
for t, lab in zip(leg.texts, new_labels):
    t.set_text(lab)
leg = ax10[1,0].legend(prop={'size':12,'family': 'Times New Roman'},fontsize=16,loc='best')

ax10[0,0].set_ylabel('log '+vd.nitnamed,fontsize=fs)
ax10twin.set_ylabel('log '+vd.phosnamed,fontsize=fs)
ax10[1,0].set_ylabel('log Pigment (mg m$^{-3}$)',fontsize=fs)

################################################
#  dump medians to file for fitting elsewhere
################################################

# phyt data medians
nit_med = vd.get_percentiles(vd.sds,vd.nitname,50)
phos_med = vd.get_percentiles(vd.sds,vd.phosname,50)
phyt_months = vd.get_months(vd.sds)
chl_med = vd.get_percentiles(vd.sds,vd.chlname,50)
phyc_med = vd.get_percentiles(vd.sds,vd.pycname,50)

# zoo data
zoo_months = vd.get_months(zv.lez)
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
# run model
################################################

# time array
ndays,delt = 1000,0.1
times = np.linspace(0,ndays,int(ndays/delt))
ndayss = 2000

# frequency for symbol plotting of model
freq = 1000

# shorter time array for dyn sims
ndaysl,deltl = 365*5,0.05
timesl = np.linspace(0,ndaysl,int(ndaysl/deltl))
fyd1 = ndaysl - 365*4

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
ndyn = -1

# default trade-off valis
ndef = 1
weakcost = 0.5
strongcost = 2.0

#################
# viz temp
#################

dtemp = mds.dyn_temp(times[-int(365/delt):],tmin,tmax)
ax3[0].plot(times[-int(365/delt):],dtemp,lw=2)

rs = mds.get_r(a,b,temps)
ax3[1].plot(temps,rs,lw=2)

#################
# viz trade-off function
#################

cors = np.linspace(0,1,100)
for (n,c,lab) in zip((weakcost,ndef,strongcost),('--','-.','-'),('High cost of defense','Medium cost of defense','Low cost of defense')):
    bors = mds.corbor(cors,n)
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
Nns,P1ns,P2ns,Zns = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
Nps,P1ps,P2ps,Zps = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
ONs,OPs = np.r_[[]],np.r_[[]]
for temp in temps:
    params = mds.get_params(temp,cor=cors[0],bor=bors[0])
    params = tuple([[np.array(p) for p in params]])
    inits = get_inits('diamond',Ntot)
    dtest = si.odeint(mds.diamond,inits,times,args=params).T
    if temp == temps[ndyn]:
        dtestshow = dtest        
    Nns = np.append(Nns,dtest[0][-1])
    P1ns = np.append(P1ns,dtest[1][-1])
    P2ns = np.append(P2ns,dtest[2][-1])
    Zns = np.append(Zns,dtest[3][-1])
    Nps = np.append(Nps,dtest[4][-1])
    P1ps = np.append(P1ps,dtest[5][-1])
    P2ps = np.append(P2ps,dtest[6][-1])
    Zps = np.append(Zps,dtest[7][-1])
    ONs = np.append(ONs,dtest[8][-1])
    OPs = np.append(OPs,dtest[9][-1])

# show convergence to equilbrium for a special case
ax2[0].plot(times[:int(ndayss/delt)],dtestshow[0][:int(ndayss/delt)],ncp,label=mnl,lw=1.5)
ax2[1].plot(times[:int(ndayss/delt)],dtestshow[1][:int(ndayss/delt)],p1cp,label=mp1l,lw=1.5)
ax2[1].plot(times[:int(ndayss/delt)],dtestshow[2][:int(ndayss/delt)],p2cp,label=mp2l,lw=1.5)
ax2[2].plot(times[:int(ndayss/delt)],dtestshow[3][:int(ndayss/delt)],z2cp,label=mz1l,lw=1.5)

# plot equilibrium
ax4[0,0].plot(xdim,P1ns,p1cp,label=mp1l,lw=3)
ax4[0,0].plot(xdim,P2ns,p2cp,label=mp2l,lw=3)
lt1 = ax4[0,2].plot(xdim,Nns,ncp,label=mnl,lw=3)
lt2 = ax4t1.plot(xdim,Nps,pcp,label=mpl,lw=3,ls='--')
ax4[0,1].plot(xdim,Zns,zcp,label=mz1l,lw=3)
lns = lt1+lt2
labs = [l.get_label() for l in lns]
leg = ax4t1.legend(lns, labs, prop={'size':18,'family': 'Times New Roman'},fontsize=18,loc='best')
leg.draw_frame(False)

ax9[0,0].plot(times[:int(ndayss/delt)],dtestshow[1][:int(ndayss/delt)],p1cp,label=mp1l,lw=2)
ax9[0,0].plot(times[:int(ndayss/delt)],dtestshow[2][:int(ndayss/delt)],p2cp,label=mp2l,lw=2)
lt1 = ax9[0,2].plot(times[:int(ndayss/delt)],dtestshow[0][:int(ndayss/delt)],ncp,label=mnl,lw=2)
lt2 = ax9t1.plot(times[:int(ndayss/delt)],dtestshow[4][:int(ndayss/delt)],pcp,label=mpl,lw=2,ls='--')
ax9[0,1].plot(times[:int(ndayss/delt)],dtestshow[3][:int(ndayss/delt)],zcp,label=mz1l,lw=2)
lns = lt1+lt2
labs = [l.get_label() for l in lns]
leg = ax9t1.legend(lns, labs, prop={'size':18,'family': 'Times New Roman'},fontsize=18,loc='best')
leg.draw_frame(False)

tmin = 4

# plot dynamic - weak cost
u = np.r_[[Nns[-5],P1ns[-5],P2ns[-5],Zns[-5],Nps[-5],P1ps[-5],P2ps[-5],Zps[-5],ONs[-5],OPs[-5]]]
Nns,P1ns,P2ns,Zns = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
Nps,P1ps,P2ps,Zps = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
ONs,OPs = np.r_[[]],np.r_[[]]
NtoPst, NtoPs = np.r_[[]],np.r_[[]]
for t in timesl:
    params = mds.get_params(mds.dyn_temp(t),cor=cors[0],bor=bors[0],Tdyn=True)
    params = tuple([[np.array(p) for p in params]])
    u = u + deltl*mds.diamond(u,t,params[0])
    Nns = np.append(Nns,u[0])
    P1ns = np.append(P1ns,u[1])
    P2ns = np.append(P2ns,u[2])
    Zns = np.append(Zns,u[3])
    Nps = np.append(Nps,u[4])
    P1ps = np.append(P1ps,u[5])
    P2ps = np.append(P2ps,u[6])
    Zps = np.append(Zps,u[7])
    ONs = np.append(ONs,u[8])
    OPs = np.append(OPs,u[9])
    if t > ndaysl - 365:
        NtoPs = np.append(NtoPs,Nns[-1]/Nps[-1])
        NtoPst = np.append(NtoPst,t - (ndaysl - 365))

lt1 = ax13[0,0].plot(timesl[-7300:]-fyd1,Nns[-7300:],label='Nitrogen')
lt2 = ax13t1.plot(timesl[-7300:]-fyd1,Nps[-7300:],c='b',label='Phosphorus',ls='--')
lns = lt1+lt2
labs = [l.get_label() for l in lns]
leg = ax13t1.legend(lns, labs, prop={'size':18,'family': 'Times New Roman'},fontsize=18,loc='best')
leg.draw_frame(False)
ax14[0,0].plot(NtoPst,NtoPs,c='b')

ax5[0].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Nns[int(-365*nyearsl/deltl)::freq]),nca,lw=1.5)
ax5[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P1ns[int(-365*nyearsl/deltl)::freq]),p1ca,lw=1.5)
ax5[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P2ns[int(-365*nyearsl/deltl)::freq]),p2ca,lw=1.5)
ax5[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Zns[int(-365*nyearsl/deltl)::freq]),zca,lw=1.5)
ax11[0,0].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P1ns[int(-365*nyears/deltl)::freq]),p1ca,lw=2.5)
ax11[0,0].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P2ns[int(-365*nyears/deltl)::freq]),p2ca,lw=2.5)
ax5[0].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Nns)[int(-365*nyearsl/deltl):],ncb,label=mnl,lw=1.5)
ax5[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P1ns)[int(-365*nyearsl/deltl):],p1cb,label=mp1l,lw=1.5)
ax5[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P2ns)[int(-365*nyearsl/deltl):],p2cb,label=mp2l,lw=1.5)
ax5[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Zns)[int(-365*nyearsl/deltl):],zcb,label=mz1l,lw=1.5)
ax11[0,0].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P1ns)[int(-365*nyears/deltl):],p1cb,label=mp1l,lw=2.5)
ax11[0,0].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P2ns)[int(-365*nyears/deltl):],p2cb,label=mp2l,lw=2.5)

# plot dynamic - strong cost
u = np.r_[[Nns[-5],P1ns[-5],P2ns[-5],Zns[-5],Nps[-5],P1ps[-5],P2ps[-5],Zps[-5]],ONs[-5],OPs[-5]]
Nns,P1ns,P2ns,Zns = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
Nps,P1ps,P2ps,Zps = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
ONs,OPs = np.r_[[]],np.r_[[]]
NtoPst, NtoPs = np.r_[[]],np.r_[[]]
for t in timesl:
    params = mds.get_params(mds.dyn_temp(t),cor=cors[1],bor=bors[1],Tdyn=True)
    params = tuple([[np.array(p) for p in params]])
    u = u + deltl*mds.diamond(u,t,params[0])
    Nns = np.append(Nns,u[0])
    P1ns = np.append(P1ns,u[1])
    P2ns = np.append(P2ns,u[2])
    Zns = np.append(Zns,u[3])
    Nps = np.append(Nps,u[4])
    P1ps = np.append(P1ps,u[5])
    P2ps = np.append(P2ps,u[6])
    Zps = np.append(Zps,u[7])
    ONs = np.append(ONs,u[8])
    OPs = np.append(OPs,u[9])
    if t > ndaysl - 365:
        NtoPs = np.append(NtoPs,Nns[-1]/Nps[-1])
        NtoPst = np.append(NtoPst,t - (ndaysl - 365))

ax13[0,1].plot(timesl[-7300:]-fyd1,Nns[-7300:])
ax13t2.plot(timesl[-7300:]-fyd1,Nps[-7300:],c='b',ls='--')
ax14[0,1].plot(NtoPst,NtoPs,c='b')

ax6[0].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Nns[int(-365*nyearsl/deltl)::freq]),nca,lw=1.5)
ax6[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P1ns[int(-365*nyearsl/deltl)::freq]),p1ca,lw=1.5)
ax6[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P2ns[int(-365*nyearsl/deltl)::freq]),p2ca,lw=1.5)
ax6[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Zns[int(-365*nyearsl/deltl)::freq]),zca,lw=1.5)
ax11[0,1].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P1ns[int(-365*nyears/deltl)::freq]),p1ca,lw=2.5)
ax11[0,1].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P2ns[int(-365*nyears/deltl)::freq]),p2ca,lw=2.5)
ax6[0].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Nns)[int(-365*nyearsl/deltl):],ncb,label=mnl,lw=1.5)
ax6[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P1ns)[int(-365*nyearsl/deltl):],p1cb,label=mp1l,lw=1.5)
ax6[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P2ns)[int(-365*nyearsl/deltl):],p2cb,label=mp2l,lw=1.5)
ax6[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Zns)[int(-365*nyearsl/deltl):],zcb,label=mz1l,lw=1.5)
ax11[0,1].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P1ns)[int(-365*nyears/deltl):],p1cb,label=mp1l,lw=2.5)
ax11[0,1].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P2ns)[int(-365*nyears/deltl):],p2cb,label=mp2l,lw=2.5)

#################
# parallel model
#################
Nns,P1ns,P2ns,Z1ns,Z2ns = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
Nps,P1ps,P2ps,Z1ps,Z2ps = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
ONs,OPs = np.r_[[]],np.r_[[]]
for temp in temps:
    params = mds.get_params(temp,cor=cors[0],bor=bors[0])
    params = tuple([[np.array(p) for p in params]])
    inits = get_inits('parallel',Ntot)
    dtest = si.odeint(mds.parallel,inits,times,params).T
    if temp == temps[ndyn]:
        dtestshow = dtest
    Nns = np.append(Nns,dtest[0][-1])
    P1ns = np.append(P1ns,dtest[1][-1])
    P2ns = np.append(P2ns,dtest[2][-1])
    Z1ns = np.append(Z1ns,dtest[3][-1])
    Z2ns = np.append(Z2ns,dtest[4][-1])
    Nps = np.append(Nps,dtest[5][-1])
    P1ps = np.append(P1ps,dtest[6][-1])
    P2ps = np.append(P2ps,dtest[7][-1])
    Z1ps = np.append(Z1ps,dtest[8][-1])
    Z2ps = np.append(Z2ps,dtest[9][-1])
    ONs = np.append(ONs,dtest[10][-1])
    OPs = np.append(OPs,dtest[11][-1])

# show convergene to equilibrium for special case
ax1[0].plot(times[:int(ndayss/delt)],dtestshow[0][:int(ndayss/delt)],ncp,label=mnl)
ax1[1].plot(times[:int(ndayss/delt)],dtestshow[1][:int(ndayss/delt)],p1cp,label=mp1l)
ax1[1].plot(times[:int(ndayss/delt)],dtestshow[2][:int(ndayss/delt)],p2cp,label=mp2l)
ax1[2].plot(times[:int(ndayss/delt)],dtestshow[3][:int(ndayss/delt)],zcp,label=mz1l)
ax1[2].plot(times[:int(ndayss/delt)],dtestshow[4][:int(ndayss/delt)],z2cp,label=mz2l)

# plot equilibrium
ax4[1,0].plot(xdim,P1ns,p1cp,label=mp1l,lw=3)
ax4[1,0].plot(xdim,P2ns,p2cp,label=mp2l,lw=3)
lt1 = ax4[1,2].plot(xdim,Nns,ncp,label=mnl,lw=3)
lt2 = ax4t2.plot(xdim,Nps,pcp,label=mpl,lw=3,ls='--')
ax4[1,1].plot(xdim,Z1ns,zcp,label=mz1l,lw=3)
ax4[1,1].plot(xdim,Z2ns,z2cp,label=mz2l,lw=3)
lns = lt1+lt2
labs = [l.get_label() for l in lns]
leg = ax4t2.legend(lns, labs, prop={'size':18,'family': 'Times New Roman'},fontsize=18,loc='best')
leg.draw_frame(False)

ax9[1,0].plot(times[:int(ndayss/delt)],dtestshow[1][:int(ndayss/delt)],p1cp,label=mp1l,lw=2)
ax9[1,0].plot(times[:int(ndayss/delt)],dtestshow[2][:int(ndayss/delt)],p2cp,label=mp2l,lw=2)
lt1 = ax9[1,2].plot(times[:int(ndayss/delt)],dtestshow[0][:int(ndayss/delt)],ncp,label=mnl,lw=2)
lt2 = ax9t2.plot(times[:int(ndayss/delt)],dtestshow[5][:int(ndayss/delt)],pcp,label=mpl,lw=2,ls='--')
ax9[1,1].plot(times[:int(ndayss/delt)],dtestshow[3][:int(ndayss/delt)],zcp,label=mz1l,lw=2)
ax9[1,1].plot(times[:int(ndayss/delt)],dtestshow[4][:int(ndayss/delt)],z2cp,label=mz2l,lw=2)
lns = lt1+lt2
labs = [l.get_label() for l in lns]
leg = ax9t2.legend(lns, labs, prop={'size':18,'family': 'Times New Roman'},fontsize=18,loc='best')
leg.draw_frame(False)

# plot dynamic - weak cost
u = np.r_[[Nns[-5],P1ns[-5],P2ns[-5],Z1ns[-5],Z2ns[-5],Nps[-5],P1ps[-5],P2ps[-5],Z1ps[-5],Z2ps[-5]],ONs[-5],OPs[-5]]
Nns,P1ns,P2ns,Z1ns,Z2ns = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
Nps,P1ps,P2ps,Z1ps,Z2ps = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
ONs,OPs = np.r_[[]],np.r_[[]]
NtoPst, NtoPs = np.r_[[]],np.r_[[]]
for t in timesl:
    params = mds.get_params(mds.dyn_temp(t),cor=cors[0],bor=bors[0])
    params = tuple([[np.array(p) for p in params]])
    u = u + deltl*mds.parallel(u,t,params[0])
    Nns = np.append(Nns,u[0])
    P1ns = np.append(P1ns,u[1])
    P2ns = np.append(P2ns,u[2])
    Z1ns = np.append(Z1ns,u[3])
    Z2ns = np.append(Z2ns,u[4])
    Nps = np.append(Nps,u[5])
    P1ps = np.append(P1ps,u[6])
    P2ps = np.append(P2ps,u[7])
    Z1ps = np.append(Z1ps,u[8])
    Z2ps = np.append(Z2ps,u[9])
    ONs = np.append(ONs,u[10])
    OPs = np.append(OPs,u[11])
    if t > ndaysl - 365:
        NtoPs = np.append(NtoPs,Nns[-1]/Nps[-1])
        NtoPst = np.append(NtoPst,t - (ndaysl - 365))

ax13[1,0].plot(timesl[-7300:]-fyd1,Nns[-7300:])
ax13t3.plot(timesl[-7300:]-fyd1,Nps[-7300:],c='b',ls='--')
ax14[1,0].plot(NtoPst,NtoPs,c='b')

ax7[0].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Nns[int(-365*nyearsl/deltl)::freq]),nca,lw=1.5)
ax7[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P1ns[int(-365*nyearsl/deltl)::freq]),p1ca,lw=1.5)
ax7[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P2ns[int(-365*nyearsl/deltl)::freq]),p2ca,lw=1.5)
ax7[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Z1ns[int(-365*nyearsl/deltl)::freq]),zca,lw=1.5)
ax7[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Z2ns[int(-365*nyearsl/deltl)::freq]),z2ca,lw=1.5)
ax11[1,0].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P1ns[int(-365*nyears/deltl)::freq]),p1ca,lw=2.5)
ax11[1,0].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P2ns[int(-365*nyears/deltl)::freq]),p2ca,lw=2.5)
ax7[0].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Nns)[int(-365*nyearsl/deltl):],ncb,label=mnl,lw=1.5)
ax7[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P1ns)[int(-365*nyearsl/deltl):],p1cb,label=mp1l,lw=1.5)
ax7[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P2ns)[int(-365*nyearsl/deltl):],p2cb,label=mp2l,lw=1.5)
ax7[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Z1ns)[int(-365*nyearsl/deltl):],zcb,label=mz1l,lw=1.5)
ax7[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Z2ns)[int(-365*nyearsl/deltl):],z2cb,label=mz2l,lw=1.5)
ax11[1,0].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P1ns)[int(-365*nyears/deltl):],p1cb,label=mp1l,lw=2.5)
ax11[1,0].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P2ns)[int(-365*nyears/deltl):],p2cb,label=mp2l,lw=2.5)

# plot dynamic - strong cost
u = np.r_[[Nns[-5],P1ns[-5],P2ns[-5],Z1ns[-5],Z2ns[-5],Nps[-5],P1ps[-5],P2ps[-5],Z1ps[-5],Z2ps[-5]],ONs[-5],OPs[-5]]
Nns,P1ns,P2ns,Z1ns,Z2ns = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
Nps,P1ps,P2ps,Z1ps,Z2ps = np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]],np.r_[[]]
ONs,OPs = np.r_[[]],np.r_[[]]
NtoPst, NtoPs = np.r_[[]],np.r_[[]]
for t in timesl:
    params = mds.get_params(mds.dyn_temp(t),cor=cors[1],bor=bors[1])
    params = tuple([[np.array(p) for p in params]])
    u = u + deltl*mds.parallel(u,t,params[0])
    Nns = np.append(Nns,u[0])
    P1ns = np.append(P1ns,u[1])
    P2ns = np.append(P2ns,u[2])
    Z1ns = np.append(Z1ns,u[3])
    Z2ns = np.append(Z2ns,u[4])
    Nps = np.append(Nps,u[5])
    P1ps = np.append(P1ps,u[6])
    P2ps = np.append(P2ps,u[7])
    Z1ps = np.append(Z1ps,u[8])
    Z2ps = np.append(Z2ps,u[9])
    ONs = np.append(ONs,u[10])
    OPs = np.append(OPs,u[11])
    if t > ndaysl - 365:
        NtoPs = np.append(NtoPs,Nns[-1]/Nps[-1])
        NtoPst = np.append(NtoPst,t - (ndaysl - 365))

ax13[1,1].plot(timesl[-7300:]-fyd1,Nns[-7300:])
ax13t4.plot(timesl[-7300:]-fyd1,Nps[-7300:],c='b',ls='--')
ax14[1,1].plot(NtoPst,NtoPs,c='b')

ax8[0].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Nns[int(-365*nyearsl/deltl)::freq]),nca,lw=1.5)
ax8[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P1ns[int(-365*nyearsl/deltl)::freq]),p1ca,lw=1.5)
ax8[1].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(P2ns[int(-365*nyearsl/deltl)::freq]),p2ca,lw=1.5)
ax8[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Z1ns[int(-365*nyearsl/deltl)::freq]),zca,lw=1.5)
ax8[2].plot(timesl[int(-365*nyearsl/deltl)::freq]-ncutl,np.log(Z2ns[int(-365*nyearsl/deltl)::freq]),z2ca,lw=1.5)
ax11[1,1].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P1ns[int(-365*nyears/deltl)::freq]),p1ca,lw=2.5)
ax11[1,1].plot(timesl[int(-365*nyears/deltl)::freq]-ncut,np.log(P2ns[int(-365*nyears/deltl)::freq]),p2ca,lw=2.5)
ax8[0].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Nns)[int(-365*nyearsl/deltl):],ncb,label=mnl,lw=1.5)
ax8[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P1ns)[int(-365*nyearsl/deltl):],p1cb,label=mp1l,lw=1.5)
ax8[1].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(P2ns)[int(-365*nyearsl/deltl):],p2cb,label=mp2l,lw=1.5)
ax8[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Z1ns)[int(-365*nyearsl/deltl):],zcb,label=mz1l,lw=1.5)
ax8[2].plot(timesl[int(-365*nyearsl/deltl):]-ncutl,np.log(Z2ns)[int(-365*nyearsl/deltl):],z2cb,label=mz2l,lw=1.5)
ax11[1,1].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P1ns)[int(-365*nyears/deltl):],p1cb,label=mp1l,lw=2.5)
ax11[1,1].plot(timesl[int(-365*nyears/deltl):]-ncut,np.log(P2ns)[int(-365*nyears/deltl):],p2cb,label=mp2l,lw=2.5)

#ax1[1].semilogy()

################################################
# final fig formatting
################################################

for axes in [ax1,ax2,ax5,ax6,ax7,ax8]:
    for a in axes.flatten()[1:]:
        l = a.legend(prop={'size':18,'family': 'Times New Roman'})
        l.draw_frame(False)

for axes in [ax4[:,:2],ax9[:,:2],ax11]:
    for a in axes.flatten():
        l = a.legend(prop={'size':18,'family': 'Times New Roman'})
        l.draw_frame(False)

for a in [ax4t1,ax4t2,ax9t1,ax9t2]:
    a.set_ylim([0,0.1])
    a.set_ylabel('Phosphorus (mmol P m$^{-3}$)',fontsize=fs)

for a in [ax4[0,2],ax4[1,2],ax9[0,2],ax9[1,2]]:
    a.set_ylim([0,5])

for a in ax13.flatten():
    a.set_ylabel('Nitrogen (mmol N m$^{-3}$)',fontsize=fs)
    a.set_xlabel('Time (days)',fontsize=fs)

for a in [ax13t1,ax13t2,ax13t3,ax13t4]:
    a.set_ylabel('Phosphorus (mmol P m$^{-3}$)',fontsize=fs)
    a.set_xlabel('Time (days)',fontsize=fs)

for a in ax14.flatten():
    a.set_ylabel('molar N:P',fontsize=fs)
    a.set_xlabel('Time (days)',fontsize=fs)

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
f13.savefig('../figures/nandp_vtime',  bbox_inches='tight',dpi=300)
f14.savefig('../figures/ntops',  bbox_inches='tight',dpi=300)

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
plt.close(f13)
plt.close(f14)

