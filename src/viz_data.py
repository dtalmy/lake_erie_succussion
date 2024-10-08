import numpy as np
import pandas as pd
import pylab as plt

################################################
# local functions for dataframe queries
################################################

def get_months(df):
    df = df.sort_values(by='Month')
    return df['Month'].unique()

def get_percentiles(df,variable,percentile):
    plist = []
    for m in get_months(df):
        plist.append(np.nanpercentile(df[df['Month']== m][variable],percentile))
    return np.log(plist)

def plot_by_year(a,df,year,varname):
    sdsyear = df[df['Year']==year].sort_values(by='Date')
    a.scatter(sdsyear['Month'],np.log(df[df['Year']==year][varname]) , color = 'gray')
    for (p,c,l) in zip([25,50,75],[':r','black',':g'], ('25th percentile', 'Median', '75th percentile')):
        a.plot(get_months(sdsyear),get_percentiles(sdsyear,varname,p),c,label=l)
    a.legend(loc='best')

def plot_by_quantile(a,df,quantile,varname):
    a.plot(get_months(df),get_percentiles(df,varname,quantile),'b--',linewidth = 3)
    for (year,c,l) in zip(years,['gray','k','b','r','g','orange'], ('2012', '2013', '2014', '2015', '2016','2017')):
        sdsyear = df[df['Year']==year].sort_values(by='Date')
        a.plot(get_months(sdsyear),get_percentiles(sdsyear,varname,quantile),c,label=l)
    a.legend(loc='best')

################################################
# Data wrangling
################################################

#import data from Lake Erie from bouys between 2012 to 2018
#HABs occur in Western Basin
df = pd.read_csv('../data/lake_erie_habs_field_sampling_results_2012_2018_v2.csv', encoding = 'latin1')

#change date to year first for sorting purposes
df['Date'] = pd.to_datetime(df['Date'])

#sort data by site and date - sds; 'site, data, sorted'
sds = df.sort_values(by =['Site', 'Date'])

# Add year and day of year columns
sds['Year'] = pd.DatetimeIndex(sds['Date']).year
sds['DayofYear'] = pd.DatetimeIndex(sds['Date']).dayofyear
sds['Month'] = sds['Date'].dt.strftime('%m')

# strings as variables for shorthand
chlname = 'Extracted Chlorophyll a (mmol N m$^{-3}$)'
pycname = 'Extracted Phycocyanin (mmol N m$^{-3}$)'
phosname = 'Total Dissolved Phosphorus (mmol P m$^{-3}$)'
nitname = 'Nitrate + Nitrite + Ammonia (mmol N m$^{-3}$)'

# strings for display
chlnamed = 'Extracted Chlorophyll a\n'+ '(mmol N m$^{-3}$)'
pycnamed = 'Extracted Phycocyanin \n'+ '(mmol N m$^{-3}$)'
phosnamed = 'Total Dissolved Phosphorus\n'+ '(mmol P m$^{-3}$)'
nitnamed = 'Nitrate + Nitrite + Ammonia\n' + '(mmol N m$^{-3}$)'
ntopname = 'ntop'

# to numeric
vnames = ['Extracted Phycocyanin (µg/L)','Nitrate + Nitrite (mg N/L)','Ammonia (µg N/L)']
for vname in vnames:
    sds[vname] = sds[vname].str.replace(r'<', '')
    sds[vname] = pd.to_numeric(sds[vname])

########################################################
# units conversions - convert all to mmol N m-3
########################################################

# nutrient unit conversion
sds['Nitrate + Nitrite + Ammonia (µg N/L)']=sds['Nitrate + Nitrite (mg N/L)']*1000 + sds['Ammonia (µg N/L)']
sds[nitname]=sds['Nitrate + Nitrite + Ammonia (µg N/L)'] * 14.0 # grams to mols
sds[phosname] = sds['Total Dissolved Phosphorus (µg P/L)'] * 30.97

# pigment unit conversion
# macintyre 2002 has chl:C g g-1 of microcystis 0.015
# same ref, C:N 3
chltoN = 0.015*3 # g/g
chltoN = 0.1 # placeholder
phyctochl = 4 # g/g
sds[chlname] = sds['Extracted Chlorophyll a (µg/L)']*(1/chltoN)*(1/14.0)
sds[pycname] = sds['Extracted Phycocyanin (µg/L)']/4

# non-cyanonbacteria chlorophyll
sds[chlname] = sds[chlname] - sds[pycname]

# ntop
sds[ntopname] = sds[nitname]/sds[phosname]

# years
years = np.sort(sds.Year.unique())[:-1]

###############################################
# plot all data
###############################################

# axes font size
fs = 16

# create figs
f1, axs = plt.subplots(4,6,figsize=(30,18))
f2, axq = plt.subplots(3,4,figsize=(30,18))

# plot all data by year and variable
for (i,vname) in zip(range(4),[chlname,pycname,phosname,nitname]):
    for (j,year) in zip(range(6),years):
        plot_by_year(axs[i,j],sds,year,vname)

# plot all data by variable and quantile
for (j,vname) in zip(range(4),[phosname,nitname,chlname,pycname]):
    for (i,quantile) in zip(range(3),[25,50,75]):
        plot_by_quantile(axq[i,j],sds,quantile,vname)
        axq[i,j].plot(get_months(sds),get_percentiles(sds,vname,quantile),'b--',linewidth = 3)

for (l,a) in zip([chlnamed,pycnamed,phosnamed,nitnamed],axs[:,0]):
    a.set_ylabel('log '+l,fontsize=fs)
    a.set_ylabel('log '+l,fontsize=fs)
    a.set_ylabel('log '+l,fontsize=fs)
    a.set_ylabel('log '+l,fontsize=fs)

for (a,year) in zip(axs.flatten(),years):
    a.set_title(year)

for a in axs.flatten():
    #a.set_ylim([-7,11])
    a.set_xlabel('Month',fontsize=fs)

for a in axq.flatten():
    #a.set_ylim([-4,9])
    a.set_xlabel('Month',fontsize=fs)

for (j,vname) in zip(range(4),[phosnamed,nitnamed,chlnamed,pycnamed]):
    for (i,quantile) in zip(range(3),[25,50,75]):
        axq[i,j].set_ylabel('log '+vname,fontsize=fs)
        axq[i,j].set_title(str(quantile)+'th quantile',fontsize=fs)

f1.subplots_adjust(hspace=0.3)
f2.subplots_adjust(hspace=0.3)

f1.savefig('../figures/noaa_by_year', bbox_inches='tight')
f2.savefig('../figures/noaa_by_percentiles', bbox_inches='tight')

plt.close(f1)
plt.close(f2)
