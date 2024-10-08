from scipy.integrate import *
from pandas import *
from numpy import *
from pylab import *
from scipy import *

# make all font Times New Roman
plt.rcParams["font.family"] = "Times New Roman"

# import data
d2004 = read_excel("../data/ke_2008_data.xlsx", sheet_name = "t2004")
d2005 = read_excel("../data/ke_2008_data.xlsx", sheet_name = "t2005")
T4,M4,D4,Z4 =array(d2004['T']),array(d2004['M']),array(d2004['E']),array(d2004['Z'])
T5,M5,D5,Z5 =array(d2005['T']),array(d2005['M']),array(d2005['E']),array(d2005['Z'])

f2,ax2 = subplots(1,2,figsize=[11,5])

labs = ['Nutrient','Eukaryotic algae','Microcystis','Viruses of Eukaryotes','Cyanophage','Microzooplankton']
cols = ['C0','C1','C2','C3','C4','C5']

ax2[0].plot(T4,D4,'-o',label=labs[1],c=cols[1])
ax2[0].plot(T4,M4,'-^',label=labs[2],c=cols[2])
ax2[0].plot(T4,Z4,'-*',label=labs[5],c=cols[5])
ax2[1].plot(T5,D5,'-o',label=labs[1],c=cols[1])
ax2[1].plot(T5,M5,'-^',label=labs[2],c=cols[2])
ax2[1].plot(T5,Z5,'-*',label=labs[5],c=cols[5])

ax2[0].set_xlim([amin(T4)-10,amax(T4)+10])
ax2[1].set_xlim([amin(T5)-10,amax(T5)+10])

ax2[0].set_ylim([0,16])
ax2[1].set_ylim([0,16])

l2 = ax2[0].legend()
l2.draw_frame(False)

ax2[0].set_xlabel('Time (julian day)')
ax2[1].set_xlabel('Time (julian day)')

ax2[0].set_ylabel(r'Biomass (mg L$^{-1}$)')
ax2[1].set_ylabel(r'Biomass (mg L$^{-1}$)')

ax2[0].set_title('Lake Taihu 2004')
ax2[1].set_title('Lake Taihu 2005')

ax2[0].text(0.07,0.9,'a',ha='center',va='center',color='k',transform=ax2[0].transAxes)
ax2[1].text(0.07,0.9,'b',ha='center',va='center',color='k',transform=ax2[1].transAxes)

f2.savefig('../figures/taihu_data', bbox_inches='tight')
plt.close(f2)
