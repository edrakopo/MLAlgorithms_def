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
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

#rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/FINALTMVARegApp_BDTG_ener_lilia_CCQE_1040_1100_NEWlookupsTEST_MYselectedVars_large_statTEST10000treesB.root')

#t=rfile.Get("EnerResoBDTG_withcutsTitus")
#c1 = ROOT.TCanvas()
#c1.cd()
#t.Draw()
#c1.SaveAs("enerresol_TITUSINFID.png")

################ DE/E vs E scatter plots
rfile1 = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/FINALTMVARegApp_BDTG_ener_lilia_CCQE_1040_1100_NEWlookupsTEST_MYselectedVars_large_statTEST10000treesB2.root') #TMVA BDTG/TITUS

regTree1=rfile1.Get("outTree")
#arr_hi_E = tree2array(regTree1,selection=('trueKE1>200' and 'trueKE1<600' and 'recoDwall*550.>100.' and 'recoToWall*2500.>100.'))
arr_hi_E1 = tree2array(regTree1,selection="trueKE1>0")#200 && trueKE1<600")# && (recoDwall*550.>100. && recoToWall*2500.>100.)")
#a = rnp.tree2array(chain,
#        branches=['d_double'],
#        selection="f_float < 100 && n_int%2 == 1")
trueKE2 = arr_hi_E1['trueKE1']
recoE_lookup1 = arr_hi_E1['recoE_lookup1']
#recoKE_TMVABDTG = arr_hi_E1['recoKE']

#rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/EnergyRecoLUT.root')
rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/mysimpleLUT3.root')
regTree=rfile.Get("outTree")
arr_hi_E = tree2array(regTree,selection="trueKE1>0")
trueKE = arr_hi_E['trueKE1']
#recoKE_TMVABDTG = arr_hi_E['recoE_lookup1']
recoKE_TMVABDTG = arr_hi_E['recoKE']
recoDwall = arr_hi_E['recoDWall']
recoToWall = arr_hi_E['recoToWall']

#rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/scikit_recoLAST.root')
#regTree=rfile.Get("tuple")
#arr_hi_E = tree2array(regTree,selection="trueKE>0")#200 && trueKE<600")# and "recoDwall>0.1818 && recoToWall>0.04")
##arr_hi_E = tree2array(regTree,selection='trueKE>200' and 'trueKE<600' and 'recoDwall*550.>100.' and 'recoToWall*2500.>100.')
#trueKE = arr_hi_E['trueKE']
#recoKE_BDTG = arr_hi_E['recoKE_BDTG']
#recoDwall = arr_hi_E['recoDwall']
#recoToWall = arr_hi_E['recoToWall']

####low energy
#E_l=200
#E_h=600
####high energy 
E_l=600
E_h=1400

#dwall_thres=80.
#towall_thres=90.
#def:
dwall_thres=100.
towall_thres=100.

v=0
for i in range(0,len(trueKE2)):
    if trueKE2[i]>E_l and trueKE2[i]<E_h and recoDwall[i]*550.>dwall_thres and recoToWall[i]*2500.>towall_thres:
       v=v+1
print('v: ', v) 
x = [0 for j in range (0,v)]
y = [0 for j in range (0,v)] 
y1 = [0 for j in range (0,v)] 
yTitus = [0 for j in range (0,v)]

k=0
for i in range(0,len(trueKE2)):
#for i in range(0,10):
    if trueKE2[i]>E_l and trueKE2[i]<E_h and recoDwall[i]*550.>dwall_thres and recoToWall[i]*2500.>towall_thres:
       #print('trueKE: ',trueKE[i], ' recoDwall*550.: ',recoDwall[i]*550.)
       #print(' DE y: ',100.*(trueKE[i]-recoKE_BDTG[i])/(1.*trueKE[i]), ' ytitus: ', 100.*(trueKE2[i]-recoE_lookup1[i])/(1.*trueKE2[i])) 
       x[k] = trueKE[i]
       #y[k] = 100.*(trueKE[i]-recoKE_BDTG[i])/(1.*trueKE[i]) 
       y1[k] = 100.*(trueKE2[i]-recoKE_TMVABDTG[i])/(1.*trueKE2[i]) 
       yTitus[k] = 100.*(trueKE2[i]-recoE_lookup1[i])/(1.*trueKE2[i])
       k=k+1
       #print('x: ',x[i] , ' y: ',y[i], ' y1: ' ,yTitus[i])


#nbins=[-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100]
nbins=np.arange(-100,100,1)
print('entries: ', len(y), '/',len(y1),':',len(yTitus) )

fig,ax=plt.subplots(ncols=1, sharey=True)#, figsize=(8, 6))
cmap = sns.light_palette('b',as_cmap=True)
#n, bins, patches = ax.hist(y, 20000, facecolor='b')    #scikit
#f0=ax.hist(y, nbins, histtype='step', fill=True, color='gold',alpha=0.75)
f1=ax.hist(y1, nbins, histtype='step', fill=False, edgecolor='red', hatch="/",alpha=0.75)#normed=1
f2=ax.hist(yTitus, nbins, histtype='step', fill=False, edgecolor='blue', hatch="\\\\",alpha=0.75)
#ax.bar(range(-100,100), range(0,13000),edgecolor='black', hatch="/")
ax.set_xlim(-100.,100.)
ax.set_xlabel('$\Delta E/E$ [%]')
ax.set_ylabel('Number of Entries')
ax.xaxis.set_label_coords(0.95, -0.08)
ax.yaxis.set_label_coords(-0.1, 0.75)
#ax.text(-95, 3500, 'TITUS lookup tables', bbox={'facecolor':'white','alpha': 0.2, 'pad': 5})

plt.yscale('log', nonposy='clip')
### box filled
#b_line = mpatches.Patch(fill=True, color='gold', label='Scikit BDT')
#r_line = mpatches.Patch(fill=False, edgecolor='red', hatch="/", label='TMVA BDT')
r_line = mpatches.Patch(fill=False, edgecolor='red', hatch="/", label='New lookup tables')
bl_line = mpatches.Patch(fill=False, edgecolor='blue', hatch="\\\\",label='TITUS lookup tables')
### solid line:
#b_line = mlines.Line2D([], [], color='black', label='Scikit BDT')#marker='*',markersize=15,
#r_line = mlines.Line2D([], [], color='red', label='TMVA BDT')
#bl_line = mlines.Line2D([], [], color='blue', label='TITUS lookup tables')
#plt.legend(handles=[b_line, r_line, bl_line])
#plt.legend(handles=[r_line, bl_line])
plt.legend(prop={'size': 19},handles=[r_line, bl_line])

plt.savefig("AINfidener_reso_red_newlookups_blueTITUS.png")
#plt.savefig("AINfidener_reso_ScikitBDTG_redTMVABDTG_blueTITUS.png")
#plt.savefig("AINfidener_reso_ScikitBDTG_redTMVABDTG_blueTITUSHE.png")
#plt.savefig("AINfidener_reso_ScikitBDTG_redTMVABDTG_blueTITUS_newFid.png")
#plt.savefig("AINfidener_reso_ScikitBDTG_redTMVABDTG_blueTITUSHE_newFid.png")

hist_trueKE_hi_E = ROOT.TH1D('trueKE', 'title', 200, -100, 100)
hist_recoKE_hi_E = ROOT.TH1D('recoKE', 'title', 200, -100, 100)
hist_trueKE_hi_E.SetLineColor(ROOT.kBlue)
hist_recoKE_hi_E.SetLineColor(ROOT.kRed)
fill_hist(hist_trueKE_hi_E, yTitus)
fill_hist(hist_recoKE_hi_E, y1)
c2 = ROOT.TCanvas()
style = ROOT.TStyle()
style.cd()
style.SetOptStat(0)
 #style.SetLabelSize(1.2)
hist_recoKE_hi_E.SetTitle("")
#hist_recoKE_hi_E.GetYaxis().SetRangeUser(-100.,100.)
hist_recoKE_hi_E.GetYaxis().SetTitle('Number of Entries')
hist_recoKE_hi_E.GetXaxis().SetTitle('#Delta E/E [%]')
hist_recoKE_hi_E.Draw()
hist_trueKE_hi_E.Draw("same")
c2.SetLogy()
c2.Draw()
c2.SaveAs("NEWlookupsred_TITUSblue.png")

