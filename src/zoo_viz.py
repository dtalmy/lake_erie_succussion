import numpy as np
import pandas as pd
import pylab as plt
import seaborn as sns

################################################
# Data wrangling
################################################

# read data
df = pd.read_csv('../data/Zooplankton_Data_1997-2018_v4.2.csv', encoding = 'latin1')

# clean data
df.dropna(axis=0, inplace=True)
df['YEAR']=df["ï»¿YEAR"]
df.drop(columns=["ï»¿YEAR"], inplace=True)

# format time
df.loc[:,'SAMPLE DATE'] = pd.to_datetime(df['SAMPLE DATE'])
df.loc[:,'Month'] = pd.DatetimeIndex(df['SAMPLE DATE']).month

# column headings for zooplankton groups
kstrings = ['CAL_ugDW/m3','CALIM_ugDW/m3','CYC_ugDW/m3', \
            'CYCIM_ugDW/m3','DAP_ugDW/m3','NDAP_ugDW/m3','PRED_ugDW/m3']

# display names for zooplankton groups
dnames = ['Calanoid Copepod Adults','Calanoid Copepod Copepodites',\
            'Cyclopoid Copepod Adults','Cyclopoid Copepod Copepodites',\
            'Daphnia', 'Non-Daphnia Herbivorous Cladocerans','Predatory Cladocerans']

# calculate total zooplankton biomass
df['total_zooplankton_biomass'] = (df['CAL_ugDW/m3'] + df['CALIM_ugDW/m3'] + \
                                  df['CYC_ugDW/m3'] + df['CYCIM_ugDW/m3'] + \
                                  df['DAP_ugDW/m3'] + df['NDAP_ugDW/m3'] + \
                                  df['PRED_ugDW/m3']) / (1000.0*12.0*6.6)

# lake erie zoo - lez
lez = df[df['LAKE']=='ER']

################################################
# plot zoo
################################################

# setup figs
f1a, axsa = plt.subplots(4,1,figsize=(30,36))
f1b, axsb = plt.subplots(3,1,figsize=(30,36))
axs = np.append(axsa,axsb)
f2, axq = plt.subplots(figsize=(14,6))

fs = 40

# totals
g = sns.boxplot(x = lez[lez['YEAR']>2011]['YEAR'], y = np.log(lez[lez['YEAR']>2011]['total_zooplankton_biomass']),hue = lez[lez['YEAR']>2011]['Month'], ax=axq)
axq.set_title('Total Zooplankton Biomass')
axq.set_ylabel('log mmol N m$^{-3}$')
leg = g.axes.get_legend()
new_title = 'Month'
leg.set_title(new_title)
new_labels = ['April', 'August']
for t, l in zip(leg.texts, new_labels):
    t.set_text(l)


# fontsize
sns.set_theme(font_scale=3)

# by group
for (a,ks,dn) in zip(axs,kstrings,dnames):
    sns.boxplot(x = lez['YEAR'], y = np.log(lez[ks]),hue = lez['Month'], ax = a)
    l = a.legend(labels=['March','April','August'])
    l.draw_frame(False)
    a.set_xlabel('Year',fontsize=fs)
    a.set_ylabel('μg/m$^{3}$',fontsize=fs)
    a.tick_params(axis='both', which='major', labelsize=22)
    a.set_title(dn,fontsize=fs)

################################################
# format and save
################################################

f1a.subplots_adjust(hspace=0.3)

f1a.savefig('../figures/zooplankton_groups_a', bbox_inches='tight')
f1b.savefig('../figures/zooplankton_groups_b', bbox_inches='tight')
f2.savefig('../figures/zooplankton_totals', bbox_inches='tight', dpi=300)

plt.close(f1a)
plt.close(f1b)
plt.close(f2)

