import numpy as np
import pandas as pd
import models as mds
import scipy as sc
import pylab as py
from matplotlib.backends.backend_pdf import PdfPages

#########################################################
# setup figures
#########################################################

ppdf = PdfPages('../figures/chain_view_combined')
fps,axps = py.subplots(4,3,figsize=[16,16])
ft1,axt1 = py.subplots(4,1,figsize=[10,16])
ft2,axt2 = py.subplots(4,1,figsize=[10,16])
ft3,axt3 = py.subplots(4,1,figsize=[10,16])

fc1,act1 = py.subplots(4,4,figsize=[16,14])
fc2,act2 = py.subplots(4,4,figsize=[16,14])
fc3,act3 = py.subplots(4,4,figsize=[16,14])
fc4,act4 = py.subplots(4,4,figsize=[16,14])
fc5,act5 = py.subplots(4,4,figsize=[16,14])
fc6,act6 = py.subplots(4,4,figsize=[16,14])
fc7,act7 = py.subplots(4,4,figsize=[16,14])
fc8,act8 = py.subplots(4,4,figsize=[16,14])
fc9,act9 = py.subplots(4,4,figsize=[16,14])

fa1,axa1 = py.subplots(2,3,figsize=[16,12])
axa1 = axa1.flatten()

a = np.concatenate((act1,act2,act3))
b = np.concatenate((act4,act5,act6))
c = np.concatenate((act7,act8,act9))
acts = np.concatenate((a,b,c),axis=1)
cfigs = [fc1,fc2,fc3,fc4,fc5,fc6,fc7,fc8,fc9]

fd1,adt1 = py.subplots(4,4,figsize=[16,14])
fd2,adt2 = py.subplots(4,4,figsize=[16,14])
fd3,adt3 = py.subplots(4,4,figsize=[16,14])
fd4,adt4 = py.subplots(4,4,figsize=[16,14])
fd5,adt5 = py.subplots(4,4,figsize=[16,14])
fd6,adt6 = py.subplots(4,4,figsize=[16,14])
fd7,adt7 = py.subplots(4,4,figsize=[16,14])
fd8,adt8 = py.subplots(4,4,figsize=[16,14])
fd9,adt9 = py.subplots(4,4,figsize=[16,14])

a = np.concatenate((adt1,adt2,adt3))
b = np.concatenate((adt4,adt5,adt6))
c = np.concatenate((adt7,adt8,adt9))
adts = np.concatenate((a,b,c),axis=1)
dfigs = [fd1,fd2,fd3,fd4,fd5,fd6,fd7,fd8,fd9]

#########################################################
# load data
#########################################################

print('reading')
posteriors_parallel = pd.read_csv('../data/posteriors_parallel.csv', skiprows=lambda i: i % 100)
posteriors_diamond = pd.read_csv('../data/posteriors_diamond.csv', skiprows=lambda i: i % 100)
print('read')

#########################################################
# plot output
#########################################################

# for output
output = np.r_[[posteriors_diamond.keys()[i] for i in (1, 2, 3, 4, 5, 6, 7, 8, 9, 12, 13, 14, 15)]]

# comp
axa1[0].hist(posteriors_parallel.alpha1,label='P1')
axa1[0].hist(posteriors_parallel.alpha2,label='P2')
axa1[1].hist(posteriors_diamond.alpha1,label='P1')
axa1[1].hist(posteriors_diamond.alpha2,label='P2')
axa1[2].hist(posteriors_parallel.mu1max)
axa1[2].hist(posteriors_parallel.mu2max)
axa1[3].hist(posteriors_diamond.mu1max)
axa1[3].hist(posteriors_diamond.mu2max)
axa1[4].hist(posteriors_parallel.w1)
axa1[4].hist(posteriors_parallel.w2)
axa1[5].hist(posteriors_diamond.w1)
axa1[5].hist(posteriors_diamond.w2)
axa1[0].set_xlabel('alphas parallel')
axa1[1].set_xlabel('alphas diamond')
axa1[2].set_xlabel('mumax parallel')
axa1[3].set_xlabel('mumax diamond')
axa1[4].set_xlabel('ws parallel')
axa1[5].set_xlabel('ws diamond')
l = axa1[0].legend()
l.draw_frame(False)

print('hists')
# posteriors
for (a,b,v) in zip(axps.flatten(),np.concatenate((axt1,axt2,axt3)),output):
    a.hist(posteriors_parallel[v],rasterized=True,alpha=0.5)
    a.hist(posteriors_diamond[v],rasterized=True,alpha=0.5)
    b.plot(posteriors_parallel[v],rasterized=True)
    b.plot(posteriors_diamond[v],rasterized=True)
    a.set_xlabel(v)
    b.set_ylabel(v)

print('correlations')
# correlation plots
for i in range(12):
    for j in range(12):
        acts[i,j].scatter(posteriors_parallel[output[i]],posteriors_parallel[output[j]],rasterized=True)
        acts[i,j].set_xlabel(output[i])
        acts[i,j].set_ylabel(output[j])
        adts[i,j].scatter(posteriors_diamond[output[i]],posteriors_diamond[output[j]],rasterized=True,color='orange')
        adts[i,j].set_xlabel(output[i])
        adts[i,j].set_ylabel(output[j])

fa1.savefig(ppdf,format='pdf')
fps.savefig(ppdf,format='pdf')
ft1.savefig(ppdf,format='pdf')
ft2.savefig(ppdf,format='pdf')
ft3.savefig(ppdf,format='pdf')
py.close(fa1)
py.close(fps)
py.close(ft1)
py.close(ft2)
py.close(ft3)

for f in cfigs:
    f.subplots_adjust(hspace=0.35,wspace=0.35)
    f.savefig(ppdf,format='pdf')
    py.close(f)

for f in dfigs:
    f.subplots_adjust(hspace=0.35,wspace=0.35)
    f.savefig(ppdf,format='pdf')
    py.close(f)

ppdf.close()
