import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='Times New Roman', size=20)
import pylab
pylab.rcParams['figure.figsize'] = 10, 6

import sys
sys.argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
from root_numpy import root2array, tree2array, fill_hist

import seaborn as sns

################ DE/E vs E scatter plots
#rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/FINALTMVARegApp_BDTG_ener_lilia_CCQE_1040_1100_NEWlookupsTEST_MYselectedVars_large_statTEST10000treesB2.root') #TMVA BDTG/TITUS
#rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/FINALTMVARegApp_BDTG_ener_lilia_CCQE_1040_1100_NEWlookupsTEST_MYselectedVars_nue.root')#TMVA BDTG/TITUS
#regTree=rfile.Get("outTree")
#arr_hi_E = tree2array(regTree,selection='trueKE1>0')#and 'trueKE1<3000' and 'recoDwall*550.>100.' and 'recoToWall*2500.>100.')
#trueKE1 = arr_hi_E['trueKE1']
#recoE_lookup1 = arr_hi_E['recoE_lookup1']
##recoKE = arr_hi_E['recoKE']
#recoDwall = arr_hi_E['recoDWall']
#recoToWall = arr_hi_E['recoToWall']

rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/scikit_recoLAST.root')
#rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/scikit_recoLASTnue.root')
regTree=rfile.Get("tuple")
arr_hi_E = tree2array(regTree,selection='trueKE>0')# and 'trueKE<3000' and 'recoDwall*550.>100.' and 'recoToWall*2500.>100.')
trueKE1 = arr_hi_E['trueKE']
recoE_lookup1 = arr_hi_E['recoKE_BDTG']
recoDwall = arr_hi_E['recoDwall']
recoToWall = arr_hi_E['recoToWall']

v=0
E_l=200
E_h=1000
for i in range(0,len(trueKE1)):
    if trueKE1[i]>E_l and trueKE1[i]<E_h and recoDwall[i]*550.>100. and recoToWall[i]*2500.>100.:
       v=v+1
print('v: ', v)
x = [0 for j in range (0,v)]
y = [0 for j in range (0,v)]

k=0
for i in range(0,len(trueKE1)):
#for i in range(0,10):
    if trueKE1[i]>E_l and trueKE1[i]<E_h and recoDwall[i]*550.>100. and recoToWall[i]*2500.>100.:
       x[k] = trueKE1[i]
       y[k] = 100.*(trueKE1[i]-recoE_lookup1[i])/(1.*trueKE1[i])
       k=k+1


#y = 100.*(trueKE1-recoE_lookup1)/(1.*trueKE1)
#print('entries: ', len(y), '/',len(x))
print('len(recoDwall): ', len(recoDwall))
v2=0
for i in range(0,len(trueKE1)):
    if trueKE1[i]>E_l and trueKE1[i]<E_h:
       v2=v2+1
recoDwallB = [0 for j in range (0,v2)]
recoTowallB = [0 for j in range (0,v2)]
yB = [0 for j in range (0,v2)]
k2=0
for i in range(0,v2): 
    if trueKE1[i]>E_l and trueKE1[i]<E_h:
       recoDwallB[k2]=recoDwall[i]*550.
       recoTowallB[k2]=recoToWall[i]*2500.
       yB[k2]=100.*(trueKE1[i]-recoE_lookup1[i])/(1.*trueKE1[i])
       k2=k2+1
print(recoDwallB[1], yB[1])

##fig, axs = plt.subplots(ncols=1, sharey=True, figsize=(8, 4))
##fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)
##ax = axs[0]
fig,ax=plt.subplots(ncols=1, sharey=True)#, figsize=(8, 6))
#cmap = sns.light_palette('xkcd:navy',as_cmap=True)
cmap = sns.light_palette('darkblue',as_cmap=True)

#hb = ax.hexbin(x, y, gridsize=522, mincnt=1, cmap=cmap)#,bins='log')#TITUS
#hb = ax.hexbin(x, y, gridsize=900, mincnt=1, cmap=cmap)    #scikit
#hb = ax.hexbin(x, y, gridsize=300, mincnt=1, cmap=cmap)#scikit nue
#hb = ax.hexbin(x, y, gridsize=150, mincnt=1, cmap=cmap)#TITUS nue
hb = ax.hexbin(recoDwallB, yB, gridsize=150, mincnt=1, cmap=cmap)
#hb = ax.hexbin(recoTowallB, yB, gridsize=150, mincnt=1, cmap=cmap)
#ax.axis([x.min(), x.max(), y.min(), y.max()])
##ax.set_title("Hexagon binning")
ax.set_ylim(-100.,100.)
ax.set_xlim(0.,600.)
ax.set_xlabel('$D{R}^{Wall}$ [cm]')
#ax.set_xlabel('$E_{MC, muon}$ [MeV]')
#ax.set_xlabel('$E_{MC, electron}$ [MeV]')
ax.set_ylabel('$\Delta E/E$ [%]')
ax.yaxis.set_label_coords(-0.1, 0.9) 
ax.xaxis.set_label_coords(0.85, -0.08) 
#ax.text(left, bottom, 'left bottom', horizontalalignment='left', verticalalignment='bottom',
cb = fig.colorbar(hb, ax=ax)
##cb.set_label('counts')

#plt.subplot(1, 2, 1)
##plt.scatter(trueKE1,100.*(trueKE1-recoE_lookup1)/trueKE1)
##plt.plot(y)
#n, bins, patches = plt.hist(y, 50, normed=1, facecolor='green', alpha=0.75)
##plt.ylim(-100.,100.)
##plt.xlim(0.,3000.)
#plt.savefig("enerResolTITUSBDTGInfid.png")
plt.savefig("enerResolScikitBDTGInfid_DWall.png")
#plt.savefig("enerResolScikitBDTGInfid_ToWall.png")
#plt.savefig("enerResolScikitBDTGInfid.png")
#plt.savefig("enerResolScikitBDTGInfidnue.png")
#plt.savefig("enerResolTITUSnueInfid.png")


##########################################
