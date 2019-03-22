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

#rfile = ROOT.TFile('/Users/edrakopo/work/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
#intree = rfile.Get('nu_eneNEW')

#E_threshold = 0 #100
#arr_lo_E = tree2array(intree,selection='trueKE<'+str(E_threshold))
#arr_hi_E = tree2array(intree,selection='trueKE>0')#+str(E_threshold))

#arr2_hi_E = arr_hi_E[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]#,'hits_pot_length2']]
#arr2_hi_E_n = arr2_hi_E.view(arr2_hi_E.dtype[0]).reshape(arr2_hi_E.shape + (-1,))
#arr3_hi_E = arr_hi_E['trueKE']

sns.set(style="ticks")

df = read_root('/Users/edrakopo/work/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root','nu_eneNEW',columns=['trueKE','total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2'])
#df=sns.load_dataset("iris")

def hexbin(x, y, color, **kwargs):
    cmap = sns.light_palette(color, as_cmap=True)
    plt.hexbin(x, y, gridsize=1000, cmap=cmap, **kwargs)

#k = kde.gaussian_kde(data.T)
#xi, yi = np.mgrid[x.min():x.max():nbins*1j, y.min():y.max():nbins*1j]
#zi = k(np.vstack([xi.flatten(), yi.flatten()]))

g = sns.PairGrid(df, diag_sharey=False, palette="GnBu_d")
#g.map_lower(sns.kdeplot, cmap="Blues_d")
#g = g.map(plt.scatter)
#g.map(hexbin);
#g.map(sns.regplot, color=".3")
#g.set(ylim=(-1, 11), yticks=[0, 5, 10]);
#g.map_upper(plt.scatter)
#g.map_upper(plt.pcolormesh)
g.map_upper(hexbin);
#g.map_upper(plt.scatter,s=30, edgecolor="white")
g.map_diag(sns.kdeplot, lw=3)
#g.map_lower(plt.scatter)
g.savefig("diag_plot3.png")
#g.savefig("thescat_plot.png")
#plt.show()



