import numpy as np
import pandas as pd
import models as mds
import scipy as sc
import pylab as py
from matplotlib.backends.backend_pdf import PdfPages
import ODElib

#########################################################
# load_data
#########################################################

df = pd.read_csv('../data/erie_tseries_all_repeat.csv')
df['abundance'] = np.exp(df['log_abundance'])
df['log_abundance'] = np.log(df['abundance'])
df['time'] = 365*4 + 30*(df['time']+0.5)
df['log_sigma'] = 0.2

dN = df[df.organism == 'P']
dP1 = df[df.organism == 'P1']
dP2 = df[df.organism == 'P2']
dZ = df[df.organism == 'Z']

#########################################################
# load time array and other params
#########################################################

delt = 1.0/240.0
ndays = 365*5
times = np.linspace(0,ndays,int(ndays/delt))

cits = pd.read_csv('../data/chain_inits_diamond_test.csv')

#########################################################
# priors
#########################################################

pw = 1.0

mu1max_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                      hyperparameters={'s':pw,'scale':0.5})
alpha1_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                      hyperparameters={'s':pw,'scale':0.12})
mp_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                      hyperparameters={'s':pw,'scale':0.15})
mz_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                      hyperparameters={'s':pw,'scale':0.15})
F_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                      hyperparameters={'s':pw,'scale':0.01})
w1_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                      hyperparameters={'s':pw,'scale':1.0})
tau1_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                      hyperparameters={'s':pw,'scale':10.0})
w2_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                      hyperparameters={'s':pw,'scale':0.1})
Sin_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                      hyperparameters={'s':pw,'scale':26.04382314})
temp_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                      hyperparameters={'s':pw,'scale':1.0})
Tdyn_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                            hyperparameters={'s':pw,'scale':1})
mu2max_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                      hyperparameters={'s':pw,'scale':0.25})
alpha2_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                      hyperparameters={'s':pw,'scale':1.0})
rho_prior=ODElib.parameter(stats_gen=sc.stats.lognorm,
                      hyperparameters={'s':pw,'scale':0.056732246})

snames = ['P','P1','P2','Z']

d1=ODElib.ModelFramework(ODE=mds.diamond,
                          parameter_names=['mu1max','alpha1','mp','mz','F','w1','tau1','w2','Sin',\
                                  'temp','Tdyn','mu2max','alpha2','rho'],
                          state_names = snames,
                          dataframe=df,
                          t_steps=1000000,
                          mu1max = mu1max_prior.copy(),
                          alpha1 = alpha1_prior.copy(),
                          mp = mp_prior.copy(),
                          mz = mz_prior.copy(),
                          F = F_prior.copy(),
                          w1 = w1_prior.copy(),
                          tau1 = tau1_prior.copy(),
                          w2 = w2_prior.copy(),
                          Sin = Sin_prior.copy(),
                          temp = temp_prior.copy(),
                          Tdyn = Tdyn_prior.copy(),
                          mu2max = mu2max_prior.copy(),
                          alpha2 = alpha2_prior.copy(),
                          rho = rho_prior.copy(),
                          P = 20,
                          P1 = 10,
                          P2 = 1,
                          Z = 0.2,
                         )

#########################################################
# setup figs
#########################################################

f1,ax1 = py.subplots(4,1,figsize=[10,8])

#########################################################
# run model
#########################################################

fp = ['temp','Tdyn','mu1max','mu2max','Sin','rho','F']
posteriors = d1.MCMC(chain_inits=cits,static_parameters=fp,iterations_per_chain=10000,cpu_cores=1)
mdat = d1.integrate()

pd.DataFrame(d1.get_Rsqrd(d1.integrate(predict_obs=True,as_dataframe=False),print_each=True),index=[0]).to_csv('../data/rsquared_diamond.csv')

pframe = pd.DataFrame(d1.get_parameters(),columns=d1.get_pnames())
pframe.to_csv('../data/chain_inits_diamond_test.csv')
#posteriors.to_csv('../data/posteriors_diamond.csv')

ax1[0].plot(dN.time/365.0,np.log(dN.abundance))
ax1[1].plot(dP1.time/365.0,np.log(dP1.abundance))
ax1[2].plot(dP2.time/365.0,np.log(dP2.abundance))
ax1[3].plot(dZ.time/365.0,np.log(dZ.abundance))

ax1[0].plot(mdat.time/365.0,np.log(mdat.P))
ax1[1].plot(mdat.time/365.0,np.log(mdat.P1))
ax1[2].plot(mdat.time/365.0,np.log(mdat.P2))
ax1[3].plot(mdat.time/365.0,np.log(mdat.Z))

ax1[0].set_ylabel('log Dissolved P \n (mmol P m$^{-3}$)')
ax1[1].set_ylabel('log Phytoplankton N \n (mmol N m$^{-3}$)')
ax1[2].set_ylabel('log Phytoplankton N \n (mmol N m$^{-3}$)')
ax1[3].set_ylabel('log Zooplankton N \n (mmol N m$^{-3}$)')
ax1[3].set_xlabel('Time (years)')

py.show()
f1.savefig('../figures/diamond_spinup')
py.close(f1)


