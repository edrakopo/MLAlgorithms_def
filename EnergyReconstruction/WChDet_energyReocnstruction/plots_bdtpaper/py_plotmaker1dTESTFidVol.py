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
recoKE_TMVABDTG = arr_hi_E1['recoKE']

rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/scikit_recoLAST.root')
regTree=rfile.Get("tuple")
arr_hi_E = tree2array(regTree,selection="trueKE>0")#200 && trueKE<600")# and "recoDwall>0.1818 && recoToWall>0.04")
#arr_hi_E = tree2array(regTree,selection='trueKE>200' and 'trueKE<600' and 'recoDwall*550.>100.' and 'recoToWall*2500.>100.')
trueKE = arr_hi_E['trueKE']
recoKE_BDTG = arr_hi_E['recoKE_BDTG']
recoDwall = arr_hi_E['recoDwall']
recoToWall = arr_hi_E['recoToWall']

####low energy
E_l=200
E_h=600
####high energy 
#E_l=600
#E_h=1400

dwall_thres=80.
towall_thres=90.

defdwall_thres=100.
deftowall_thres=100.

v=0
u=0
for i in range(0,len(trueKE)):
    if trueKE[i]>E_l and trueKE[i]<E_h and recoDwall[i]*550.>dwall_thres and recoToWall[i]*2500.>towall_thres:
       v=v+1
       if (recoDwall[i]*550.<defdwall_thres or recoToWall[i]*2500.<deftowall_thres):
          u=u+1
#print('v: ', v) 
x = [0 for j in range (0,v)]
y = [0 for j in range (0,v)] 
y1 = [0 for j in range (0,u)] 
y1Titus = [0 for j in range (0,u)]
yTitus = [0 for j in range (0,v)]

v0=0
for i in range(0,len(trueKE)):
    if trueKE[i]>E_l and trueKE[i]<E_h and recoDwall[i]*550.>defdwall_thres and recoToWall[i]*2500.>deftowall_thres:
       v0=v0+1
#print('v0-def cut: ', v0)
x0 = [0 for j in range (0,v0)]
y0 = [0 for j in range (0,v0)]
y0Titus = [0 for j in range (0,v0)]
#-----------#
good=0  
bad=0 
k=0
l=0
for i in range(0,len(trueKE)):
#for i in range(0,10):
    if trueKE[i]>E_l and trueKE[i]<E_h and recoDwall[i]*550.>dwall_thres and recoToWall[i]*2500.>towall_thres:
       #print('trueKE: ',trueKE[i], ' recoDwall*550.: ',recoDwall[i]*550.)
       #print(' DE y: ',100.*(trueKE[i]-recoKE_BDTG[i])/(1.*trueKE[i]), ' ytitus: ', 100.*(trueKE2[i]-recoE_lookup1[i])/(1.*trueKE2[i])) 
       x[k] = trueKE[i]
       y[k] = 100.*(recoKE_BDTG[i]-trueKE[i])/(1.*trueKE[i]) 
       yTitus[k] = 100.*(recoE_lookup1[i]-trueKE2[i])/(1.*trueKE2[i])
       if (recoDwall[i]*550.<defdwall_thres or recoToWall[i]*2500.<deftowall_thres):
           y1[l] = 100.*(recoKE_BDTG[i]-trueKE[i])/(1.*trueKE[i]) 
           y1Titus[l] = 100.*(recoE_lookup1[i]-trueKE2[i])/(1.*trueKE2[i])
           l=l+1
#       if y[k]>10. and (recoDwall[i]*550.<defdwall_thres or recoToWall[i]*2500.<deftowall_thres): 
#          print('bad evet: ', y[k], '|', recoDwall[i]*550., ',' , recoToWall[i]*2500.)
       if y[k]>10.:
          bad=bad+1
       if y[k]<=10.:
          good=good+1
       k=k+1
       #print('x: ',x[i] , ' y: ',y[i], ' y1: ' ,yTitus[i])
print('good evts: ', good, ' bad evnts: ', bad)

good0=0
bad0=0       
k0=0
good0TITUS=0
bad0TITUS=0
for i in range(0,len(trueKE)):
    if trueKE[i]>E_l and trueKE[i]<E_h and recoDwall[i]*550.>defdwall_thres and recoToWall[i]*2500.>deftowall_thres:
       #print('trueKE: ',trueKE[i], ' recoDwall*550.: ',recoDwall[i]*550.)
       x0[k0] = trueKE[i]
       y0[k0] = 100.*(recoKE_BDTG[i]-trueKE[i])/(1.*trueKE[i])
       y0Titus[k0] = 100.*(recoE_lookup1[i]-trueKE2[i])/(1.*trueKE2[i])
       if y0[k0]>10.:
          bad0=bad0+1
       if y0[k0]<=10.:
          good0=good0+1
       if yTitus[k0]>10.:
          bad0TITUS=bad0TITUS+1
       if yTitus[k0]<=10.:
          good0TITUS=good0TITUS+1
       k0=k0+1
print('Def: good evts: ', good0, ' bad evnts: ', bad0)
print('new/def in good: ', 100.*(good-good0)/good0, ' in bad: ', 100.*(bad-bad0)/bad0)
#------------#
print('in def entries: ', len(y0), ' in new cut entries: ', len(y) , ' TITUS: ', len(yTitus))
print('Efficiency increase in evts %: ', 100*(len(y)-len(y0))/len(y0))
print(' Eff %: good/all def: ', 100*good0/len(y0) , ' and bad/all: ', 100*bad0/len(y0) )
print(' Eff %: good/all TITUS def: ', 100*good0TITUS/len(yTitus) , ' and bad/all: ', 100*bad0TITUS/len(yTitus) )
print(' Eff %: good/all new cut: ', 100*good/len(y) , ' and bad/all: ', 100*bad/len(y) )
print('std dev in def: ', np.std(y0), ' and with new cut: ',np.std(y), ' and in TITUS def: ', np.std(yTitus))
print('mean in def: ', np.mean(y0), ' and with new cut: ',np.mean(y), ' and in TITUS def: ', np.mean(yTitus))
print('median in def: ', np.median(y0), ' and with new cut: ',np.median(y), ' and in TITUS def: ', np.median(yTitus))
V1=( (550.-defdwall_thres)*(550.-defdwall_thres))
V2=( (550.-dwall_thres)*(550.-dwall_thres))
print('v1: ',V1, ' Increase in fid vol: % ', 100.*(V2-V1)/V1 )

#nbins=[-100,-90,-80,-70,-60,-50,-40,-30,-20,-10,0,10,20,30,40,50,60,70,80,90,100]
nbins=np.arange(-100,100,1)

fig,ax=plt.subplots(ncols=1, sharey=True)#, figsize=(8, 6))
cmap = sns.light_palette('b',as_cmap=True)
#n, bins, patches = ax.hist(y, 20000, facecolor='b')    #scikit
#f0=ax.hist(y0, nbins, histtype='step', fill=False, edgecolor='black')#, hatch="|",alpha=0.75)
f0=ax.hist(y0, nbins, histtype='step', fill=True, color='gold',alpha=0.75)
f1=ax.hist(y, nbins, histtype='step', fill=False, edgecolor='red')#, hatch="/",alpha=0.75)#normed=1
f3=ax.hist(y1, nbins, histtype='step', fill=False, edgecolor='green', hatch="\\",alpha=0.75)#normed=1
#f2=ax.hist(yTitus, nbins, histtype='step', fill=False, edgecolor='blue', hatch="\\",alpha=0.75)
#ax.bar(range(-100,100), range(0,13000),edgecolor='black', hatch="/")
ax.set_xlim(-100.,100.)
ax.set_ylim(1.,7000.)
ax.set_xlabel('$\Delta E/E$ [%]')
ax.set_ylabel('Number of Entries')
ax.xaxis.set_label_coords(0.95, -0.08)
ax.yaxis.set_label_coords(-0.1, 0.75)
ax.text(-95, 3500, 'Scikit BDT', bbox={'facecolor':'white','alpha': 0.2, 'pad': 5})
plt.savefig("testFidVol.png")

plt.yscale('log', nonposy='clip')
### box filled
#b_line = mpatches.Patch(fill=False, edgecolor='black', label='Fiducial volume')
b_line = mpatches.Patch(fill=True, color='gold', label='Fiducial volume')
r_line = mpatches.Patch(fill=False, edgecolor='red', label='New fiducial volume')
bl_line = mpatches.Patch(fill=False, edgecolor='green', hatch="\\", label='Events added')
plt.legend(prop={'size': 19},handles=[b_line, r_line, bl_line])
plt.savefig("testFidVol_log.png")
#plt.savefig("testFidVol_logHE.png")
################################################
fig,ax=plt.subplots(ncols=1, sharey=True)#, figsize=(8, 6))
cmap = sns.light_palette('b',as_cmap=True)
f0=ax.hist(y0Titus, nbins, histtype='step', fill=False, edgecolor='blue',alpha=0.75)
f1=ax.hist(yTitus, nbins, histtype='step', fill=False, edgecolor='red')#, hatch="/",alpha=0.75)#normed=1
f3=ax.hist(y1Titus, nbins, histtype='step', fill=False, edgecolor='green', hatch="\\",alpha=0.75)#normed=1
ax.set_xlim(-100.,100.)
ax.set_ylim(1.,7000.)
ax.set_xlabel('$\Delta E/E$ [%]')
ax.set_ylabel('Number of Entries')
ax.xaxis.set_label_coords(0.95, -0.08)
ax.yaxis.set_label_coords(-0.1, 0.75)
ax.text(-95, 3500, 'Lookup Tables', bbox={'facecolor':'white','alpha': 0.2, 'pad': 5})

plt.yscale('log', nonposy='clip')
b_lineT = mpatches.Patch(fill=False, edgecolor='blue', label='Fiducial volume')
r_lineT = mpatches.Patch(fill=False, edgecolor='red', label='New fiducial volume')
bl_lineT = mpatches.Patch(fill=False, edgecolor='green', hatch="\\", label='Events added')
plt.legend(prop={'size': 19},handles=[b_lineT, r_lineT, bl_line])
plt.savefig("testFidVolTITUS_log.png")

#b_line = mpatches.Patch(fill=False, edgecolor='black', hatch="|", label='Scikit BDT')
#r_line = mpatches.Patch(fill=False, edgecolor='red', hatch="/", label='TMVA BDT')
#bl_line = mpatches.Patch(fill=False, edgecolor='blue', hatch="\\",label='TITUS lookup tables')
### solid line:
#b_line = mlines.Line2D([], [], color='black', label='Scikit BDT')#marker='*',markersize=15,
#r_line = mlines.Line2D([], [], color='red', label='TMVA BDT')
#bl_line = mlines.Line2D([], [], color='blue', label='TITUS lookup tables')
#plt.legend(handles=[b_line, r_line, bl_line])

#plt.savefig("AINfidener_reso_ScikitBDTG_redTMVABDTG_blueTITUS.png")
#plt.savefig("AINfidener_reso_ScikitBDTG_redTMVABDTG_blueTITUSHE.png")

