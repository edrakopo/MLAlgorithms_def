import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='Times', size=20)
import pylab
pylab.rcParams['figure.figsize'] = 6, 3

from sklearn import linear_model, ensemble
import sys
sys.argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
from root_numpy import root2array, tree2array, fill_hist

import random
from sklearn.metrics import mean_squared_error
from sklearn.utils import shuffle
import seaborn as sns
from root_pandas import read_root
from matplotlib import rc, font_manager
from matplotlib.pyplot import figure, axes, plot, xlabel, ylabel, title, \
grid, savefig, show

#df = read_root('/Users/edrakopo/work/work_energy_reco/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root','nu_eneNEW',columns=['trueKE','total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2'])
#df=sns.load_dataset("iris")
df = read_root('/Users/edrakopo/work/work_energy_reco/ener_reco_scikit/nu_nue_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB.root','nu_eneNEW',columns=['trueKE','total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2'])

#df.loc[:,'quantity'] *= -1
#df['quantity'] = df['quantity'].apply(lambda x: x*-1)
df['trueKE'] = df['trueKE'].apply(lambda x: x*0.001)
df['total_hits2'] = df['total_hits2'].apply(lambda x: x*40)
df['total_ring_PEs2'] = df['total_ring_PEs2'].apply(lambda x: x*20)
df['recoDWallR2'] = df['recoDWallR2'].apply(lambda x: x*5.5)
df['recoDWallZ2'] = df['recoDWallZ2'].apply(lambda x: x*12)
df['lambda_max_2'] = df['lambda_max_2'].apply(lambda x: x*25)

def hexbin(x, y, color, **kwargs):
    cmap = sns.light_palette(color, as_cmap=True)
    plt.hexbin(x, y, gridsize=200, cmap=cmap, **kwargs)

#def pairgrid_heatmap(x, y, color, **kws):
def pairgrid_heatmap(x, y, **kws):
    cmap = sns.light_palette(kws.pop("color"), as_cmap=True)
    #color="navy"
    #cmap = sns.light_palette(color, as_cmap=True)
    plt.hist2d(x, y, cmap=cmap, cmin=0.01, **kws)


#sns.palplot(sns.color_palette("coolwarm", 7))
#sns.palplot(sns.color_palette("RdBu_r", 7))
sns.set(style="ticks")
sns. set_context(font_scale=5)
#sns.set_style("white")
#sns.set(font="Times New Roman")
g = sns.PairGrid(df, diag_sharey=False)
#g = sns.PairGrid(df, diag_sharey=False, palette="Reds")
#g.map_lower(sns.kdeplot, cmap="Blues_d")
#g.map_upper(pairgrid_heatmap, bins=100)
g.map_offdiag(pairgrid_heatmap, bins=100)
#g.map_upper(hexbin);
#g.map_upper(plt.scatter,s=30, edgecolor="white")
g.map_diag(sns.kdeplot, lw=3)

plt.rc('text', usetex=True)
plt.rc('font', family='times')
xlabels,ylabels = ['$E_{MC, muon}\ [GeV]$','$N_{C}^{Hit} \ [x10^{3}]$','$N_{Ring}^{Hit}\ [x10^{3}]$','$D_{R}^{Wall}\ [m]$','$D_{Z}^{Wall}\ [m]$','$L\ [m]$'],['$E_{MC, muon}\ [GeV]$','$N_{C}^{Hit}\ [x10^{3}]$','$N_{Ring}^{Hit}\ [x10^{3}]$','$D_{R}^{Wall}\ [m]$','$D_{Z}^{Wall}\ [m]$','$L\ [m]$']


for i in range(len(xlabels)):
    for j in range(len(ylabels)):
        g.axes[j,i].xaxis.set_label_text(xlabels[i],fontdict={'fontsize' : 22, 'family' :'Times New Roman'})
        g.axes[j,i].yaxis.set_label_text(ylabels[j],fontdict={'fontsize' : 22, 'family' :'Times New Roman'})
        #g.axes[j,i].xaxis.set_label_text(xlabels[i],fontdict={'family' :'Times New Roman'})
        #g.axes[j,i].yaxis.set_label_text(ylabels[j],fontdict={'family' :'Times New Roman'})

for u in range(len(ylabels)):
    for v in range(len(ylabels)):
        if u>0 and v>0:
            if u+1<(len(ylabels)) and v+1<(len(ylabels)):
               g.axes[u, v].xaxis.get_label().set_visible(False)
               g.axes[u, v].yaxis.get_label().set_visible(False)
        if u<1 and v+1<(len(ylabels)):
           g.axes[v, u].xaxis.get_label().set_visible(False)
        if u+1==(len(ylabels)) and v>0:
           #print('u: ' ,u)
           g.axes[u, v].yaxis.get_label().set_visible(False)

g.axes[0,0].set_xlim(0,2)
labels = [item.get_text() for item in g.axes[5,0].get_xticklabels()]
#for l in labels:
#    print('l: ' ,l)
new_labels=[0, 0.5, 1, 1.5, 2]
#new_labels=[0, 1, 2]
    #new_labels = ["%d" % int(float(l)) if '.0' in l else '']
g.axes[5,0].set_xticklabels(new_labels)
g.axes[0,0].set_yticklabels(new_labels)

new_labelsHits=[0, 5, 10, 15, 20]
g.axes[1,0].set_yticks(new_labelsHits)
g.axes[1,0].set_yticklabels(new_labelsHits)
g.axes[2,0].set_yticks(new_labelsHits)
g.axes[2,0].set_yticklabels(new_labelsHits)
g.axes[5,0].set_yticks(new_labelsHits)
g.axes[5,0].set_yticklabels(new_labelsHits)

g.axes[5,1].set_xticks(new_labelsHits)
g.axes[5,1].set_xticklabels(new_labelsHits)
g.axes[5,2].set_xticks(new_labelsHits)
g.axes[5,2].set_xticklabels(new_labelsHits)
g.axes[5,5].set_xticks(new_labelsHits)
g.axes[5,5].set_xticklabels(new_labelsHits)

new_labelsDR=[0, 2, 4, 6]
g.axes[3,0].set_yticks(new_labelsDR)
g.axes[3,0].set_yticklabels(new_labelsDR)
g.axes[5,3].set_xticks(new_labelsDR)
g.axes[5,3].set_xticklabels(new_labelsDR)
new_labelsDZ=[0, 5, 10, 15, 20]
g.axes[4,0].set_yticks(new_labelsDZ)
g.axes[4,0].set_yticklabels(new_labelsDZ)
g.axes[5,4].set_xticks(new_labelsDZ)
g.axes[5,4].set_xticklabels(new_labelsDZ)

#ticks_font = font_manager.FontProperties(family='Helvetica', style='normal',
#    size=sizeOfFont, weight='normal', stretch='normal')
ticks_font = font_manager.FontProperties(family='Times New Roman',size=22, weight='normal', stretch='normal',style='normal')
for k in range(len(ylabels)):
    g.axes[0,k].set_ylim(0,2)
    g.axes[1,k].set_ylim(0,20)
    g.axes[2,k].set_ylim(0,20)
    g.axes[3,k].set_ylim(0,6)
    g.axes[4,k].set_ylim(0,12)
    g.axes[5,k].set_ylim(0,22)
    g.axes[k,1].set_xlim(0,20)
    g.axes[k,2].set_xlim(0,20)
    g.axes[k,3].set_xlim(0,6)
    g.axes[k,4].set_xlim(0,12)
    g.axes[k,5].set_xlim(0,22)

    for ticksety in [g.axes[k,0].yaxis.get_major_ticks()]:
        #[(tick.label.set_fontsize(20)) for tick in ticksety]
        [tick.label.set_fontproperties(ticks_font) for tick in ticksety]
    
    for ticksetx in [g.axes[5,k].xaxis.get_major_ticks()]:
#        #[(tick.label.set_fontsize(20)) for tick in tickset]
        [tick.label.set_fontproperties(ticks_font) for tick in ticksetx]
    #for label in g.axes[5,k].get_xticklabels():
    #    label.set_fontproperties(ticks_font)

#    for l in range(len(ylabels)):
#        if l+1<len(ylabels):
i#           g.axes[0,l+1].set_xlim(0,1)
#           print('l+1: ',l+1)
#draw only lower plots
for i, j in zip(*np.triu_indices_from(g.axes, 1)):
    g.axes[i, j].set_visible(False)

g.savefig("diag_plot_nue.png")
#g.savefig("diag_plotC.png")
#savefig('diag_plotC.pdf')
#g.savefig("thescat_plot.png")
#plt.show()



