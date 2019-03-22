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
#rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/EnergyRecoLUT.root')
#rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/mysimpleLUT3.root')

#regTree=rfile.Get("outTree")
#arr_hi_E = tree2array(regTree,selection='trueKE1>0')# and 'trueKE1<3000' and 'recoDwall*550.>100.' and 'recoToWall*2500.>100.')
#trueKE1 = arr_hi_E['trueKE1']
#recoE_lookup1 = arr_hi_E['recoE_lookup1']
##recoE_lookup1 = arr_hi_E['recoKE']
#recoDwall = arr_hi_E['recoDWall']
#recoToWall = arr_hi_E['recoToWall']

#rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/scikit_recoLAST.root')
rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/scikit_recoLASTnue.root')
#rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/WChscikit_reco4vars.root')#BDTG with only DWAllR,Z Nhits, Nringpes
regTree=rfile.Get("tuple")
arr_hi_E = tree2array(regTree,selection='trueKE>0')# and 'trueKE<3000' and 'recoDwall*550.>100.' and 'recoToWall*2500.>100.')
trueKE1 = arr_hi_E['trueKE']
recoE_lookup1 = arr_hi_E['recoKE_BDTG']
recoDwall = arr_hi_E['recoDwall']
recoToWall = arr_hi_E['recoToWall']

v=0
E_l=0
E_h=3000
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
       y[k] = 100.*(recoE_lookup1[i]-trueKE1[i])/(1.*trueKE1[i])
       k=k+1

#x = trueKE1
#y = 100.*(trueKE1-recoE_lookup1)/(1.*trueKE1)
print('entries: ', len(y), '/',len(x))
#for i in range(0,50): 
#    print(y[i])

##fig, axs = plt.subplots(ncols=1, sharey=True, figsize=(8, 4))
##fig.subplots_adjust(hspace=0.5, left=0.07, right=0.93)
##ax = axs[0]
fig,ax=plt.subplots(ncols=1, sharey=True)#, figsize=(8, 6))
#cmap = sns.light_palette('xkcd:navy',as_cmap=True)
cmap = sns.light_palette('darkblue',as_cmap=True)

#hb = ax.hexbin(x, y, gridsize=522, mincnt=1, cmap=cmap)#,bins='log')#TITUS
#hb = ax.hexbin(x, y, gridsize=1260, mincnt=1, cmap=cmap)#,bins='log')# NEW lookups
#hb = ax.hexbin(x, y, gridsize=900, mincnt=1, cmap=cmap)    #scikit
#hb = ax.hexbin(x, y, gridsize=300, mincnt=1, cmap=cmap)#scikit nue
#hb = ax.hexbin(x, y, gridsize=150, mincnt=1, cmap=cmap)#TITUS nue
#ax.axis([x.min(), x.max(), y.min(), y.max()])
##ax.set_title("Hexagon binning")
ax.set_ylim(-100.,100.)
ax.set_xlim(0.,3000.)
#ax.set_xlabel('$E_{MC, muon}$ [MeV]')
ax.set_xlabel('$E_{MC, electron}$ [MeV]')
ax.set_ylabel('$\Delta E/E$ [%]')
ax.yaxis.set_label_coords(-0.1, 0.9) 
ax.xaxis.set_label_coords(0.85, -0.08) 
ax.text(2250, 85, 'Scikit BDT', bbox={'facecolor':'white','alpha': 0.2, 'pad': 5})
#ax.text(2050, 85, 'Lookup Tables', bbox={'facecolor':'white','alpha': 0.2, 'pad': 5})
#ax.text(left, bottom, 'left bottom', horizontalalignment='left', verticalalignment='bottom',
cb = fig.colorbar(hb, ax=ax)
##cb.set_label('counts')

#plt.subplot(1, 2, 1)
##plt.scatter(trueKE1,100.*(trueKE1-recoE_lookup1)/trueKE1)
##plt.plot(y)
#n, bins, patches = plt.hist(y, 50, normed=1, facecolor='green', alpha=0.75)
##plt.ylim(-100.,100.)
##plt.xlim(0.,3000.)
#plt.savefig("TESTScikitBDTGInfid4vars.png")
#plt.savefig("enerResolTITUSBDTGInfid.png")
#plt.savefig("enerResolScikitBDTGInfid.png")
#plt.savefig("enerResolScikitBDTGInfidnue.png")
#plt.savefig("enerResolTITUSnueInfid.png")
#plt.savefig("enerResolNEWlookupsInfid.png")

twod_GBR_abs_hi_E = np.dstack((x,y))
hist_GBR_abs_hi_E = ROOT.TH2D('name_GBR_abs_hi_E', 'title', 30, 0, 3000, 200, -100, 100)
fill_hist(hist_GBR_abs_hi_E, twod_GBR_abs_hi_E[0])
profile_GBR_abs_hi_E = hist_GBR_abs_hi_E.ProfileX()
canvas_prof = ROOT.TCanvas()
canvas_prof.cd(1)
style = ROOT.TStyle()
style.cd()
profile_GBR_abs_hi_E.Draw("ColZ")
profile_GBR_abs_hi_E.GetYaxis().SetRangeUser(-100.,100.)
canvas_prof.Draw()
canvas_prof.SaveAs("DE_E_profX.png")
##########################################
