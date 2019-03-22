
# coding: utf-8

# In[1]:

#from IPython import get_ipython
#ipython_shell = get_ipython()
from IPython import get_ipython
ipython = get_ipython()
import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif', size=80)
import matplotlib.cm as cm
import pylab
#get_ipython().magic(u'matplotlib inline')
pylab.rcParams['figure.figsize'] = 8, 8


# In[ ]:

import ROOT
from root_numpy import root2array, tree2array, fill_hist
from sklearn import linear_model, ensemble
from ROOT import TFile, TTree, TBranch


# In[1]:

#rfile = ROOT.TFile('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree.root')
rfile = ROOT.TFile('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
intree = rfile.Get('nu_eneNEW')


# In[ ]:

f = ROOT.TFile( "enerReco0.root", "recreate" )
tr = ROOT.TTree( "tr", "tree with histos" )
arr=tree2array(intree)


# In[ ]:

arr2=arr[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]#,'hits_pot_length2']]
arr2_n=arr2.view(arr2.dtype[0]).reshape(arr2.shape + (-1,))
arr3=arr['trueKE']


# In[ ]:

clf = linear_model.SGDRegressor()
clf.fit(arr2_n,arr3)
clf


# In[ ]:

plt.scatter(arr3,clf.predict(arr2_n)-arr3)
plt.xlabel("trueKE [MeV]")
plt.ylabel("recoKE - trueKE [MeV]")


# In[ ]:

plt.scatter(arr3,clf.predict(arr2_n)-arr3)
plt.xlabel("trueKE [MeV]")
plt.ylabel("recoKE - trueKE [MeV]")


# In[ ]:

chain = ROOT.TChain('nu_eneNEW')
for i in range(1040,1099):
    if i==1051: i=i+1
    if i==1075: i=i+1    
    chain.Add('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_'+str(i)+'_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
              #CCQE_12in_energy_studies_recoquant_tree_NEWlookups.root')
test_data = tree2array(chain)


# In[ ]:

test_data_reduced = test_data[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]#,'hits_pot_length2']]
test_data_reduced_n = test_data_reduced.view(test_data_reduced.dtype[0]).reshape(test_data_reduced.shape + (-1,))
test_data_trueKE = test_data['trueKE']


# In[ ]:

plt.scatter(test_data_trueKE,clf.predict(test_data_reduced_n)-test_data_trueKE)
plt.ylim((-5000,3000))
plt.xlabel("trueKE [MeV]")
plt.ylabel("recoKE - trueKE [MeV]")


# In[ ]:

plt.scatter(test_data_trueKE,(clf.predict(test_data_reduced_n)-test_data_trueKE)/test_data_trueKE)
plt.ylim((0,1))
plt.xlabel("trueKE [MeV]")
plt.ylabel("DeltaE/E")
res_twod_SGD = np.dstack((test_data_trueKE, (clf.predict(test_data_reduced_n)-test_data_trueKE)/test_data_trueKE))


# In[ ]:

hist_SGD = ROOT.TH2D('name', 'title', 100, 0, 5000, 100, 0, 100)
fill_hist(hist_SGD, res_twod_SGD[0])
hist_SGD.Draw()
ROOT.gPad.Draw()


# In[ ]:

profile_SGD = hist_SGD.ProfileX()
profile_SGD.SetLineColor(ROOT.kBlue)
profile_SGD.Draw()
ROOT.gPad.Draw()


# In[ ]:

params = {'n_estimators': 1000, 'max_depth': 10, 'min_samples_split': 1,
          'learning_rate': 0.01, 'loss': 'lad'}
net = ensemble.GradientBoostingRegressor(**params)
net.fit(arr2_n,arr3)
net


# In[ ]:

plt.scatter(arr3,net.predict(arr2_n)-arr3,c='r')
plt.xlabel("trueKE [MeV]")
plt.ylabel("recoKE - trueKE [MeV]")
plt.ylim(-50,50)


# In[ ]:

plt.scatter(test_data_trueKE,(net.predict(test_data_reduced_n)-test_data_trueKE), c='r')
plt.xlabel("trueKE [MeV]")
plt.ylabel("recoKE - trueKE [MeV]")
matrix = np.dstack((test_data_trueKE, (net.predict(test_data_reduced_n)-test_data_trueKE)))


# In[ ]:

plt.scatter(test_data_trueKE,((net.predict(test_data_reduced_n)-test_data_trueKE)/test_data_trueKE),c='r')
plt.xlabel("trueKE [MeV]")
plt.ylabel("DeltaE/E")
plt.ylim(-2,2)
twod_GBR_abs = np.dstack((test_data_trueKE, np.abs(net.predict(test_data_reduced_n)-test_data_trueKE)/test_data_trueKE))


# In[ ]:

hist_GBR_abs = ROOT.TH2D('name_GBR_abs', 'title', 100, 0, 2000, 200, -100, 100)
fill_hist(hist_GBR_abs, twod_GBR_abs[0])
canvas = ROOT.TCanvas()
hist_GBR_abs.Draw()
hist_GBR_abs.GetXaxis().SetTitle('true KE [MeV]')
hist_GBR_abs.GetYaxis().SetTitle('abs(#Delta E)/E')
canvas.Draw()


# In[ ]:

profile_GBR_abs = hist_GBR_abs.ProfileX()
profile_GBR_abs.SetLineColor(ROOT.kBlue+2)
profile_GBR_abs.SetMarkerColor(ROOT.kBlue+2)
profile_GBR_abs.SetLineWidth(1)
canvas_prof = ROOT.TCanvas()
profile_GBR_abs.Draw()
profile_GBR_abs.SetMinimum(0)
profile_GBR_abs.SetMaximum(1)
profile_GBR_abs.GetXaxis().SetTitle('true KE [MeV]')
profile_GBR_abs.GetYaxis().SetTitle('abs(#Delta E)/E')
canvas_prof.Draw()


# In[ ]:

twodemu_GBR = np.dstack((test_data_trueKE, 100*(test_data_trueKE-net.predict(test_data_reduced_n))/test_data_trueKE))
hist_GBR_emu = ROOT.TH2D('name_GBR_emu', 'title', 100, 0, 3000, 200, -100, 100)
fill_hist(hist_GBR_emu, twodemu_GBR[0])
canvas = ROOT.TCanvas()
hist_GBR_emu.Draw()
hist_GBR_emu.GetXaxis().SetTitle('E_{MC,muon} [MeV]')
hist_GBR_emu.GetYaxis().SetTitle('#Delta E/E')
canvas.Draw()


# In[ ]:

profile_GBR_emu = hist_GBR_emu.ProfileX("pf",0, 3000,"s")
profile_GBR_emu.SetLineColor(ROOT.kBlue+2)
profile_GBR_emu.SetMarkerColor(ROOT.kBlue+2)
profile_GBR_emu.SetLineWidth(1)
canvas_prof2 = ROOT.TCanvas()
profile_GBR_emu.Draw()
#profile_GBR_emu.SetMinimum(0)
#profile_GBR_emu.SetMaximum(1)
profile_GBR_emu.SetTitle('Energy Resolution using Scikit BDTG')
profile_GBR_emu.SetStats(0)
#profile_GBR_emu.GetYaxis().SetRangeUser(-100,100)
profile_GBR_emu.GetYaxis().SetRangeUser(-50,50)
profile_GBR_emu.GetXaxis().SetRangeUser(0,3000)
profile_GBR_emu.GetXaxis().SetTitle('E_{MC,muon} [MeV]')
profile_GBR_emu.GetYaxis().SetTitle('#Delta E/E [%]')
canvas_prof2.Draw()


# In[ ]:

#plt.scatter(test_data_trueKE, net.predict(test_data_reduced_n),c='r',cmap='viridis')
#plt.xlabel('E_{MC,muon} [MeV]')
#plt.ylabel('E_{reco,muon} [MeV]')
#plt.xlim(0,3000)
#plt.ylim(0,3000)
emu_GBR = np.dstack((test_data_trueKE, net.predict(test_data_reduced_n)))
hist_GBR_emusca = ROOT.TH2D('name_GBR_emuscat','', 150, 0, 3000, 150, 0, 3000)
fill_hist(hist_GBR_emusca, emu_GBR[0])
canvas = ROOT.TCanvas()
ROOT.gStyle.SetPalette(1)
hist_GBR_emusca.SetStats(0)
hist_GBR_emusca.Draw("ColZ")
hist_GBR_emusca.GetXaxis().SetTitle('E_{MC,muon} [MeV]')
hist_GBR_emusca.GetYaxis().SetTitleOffset(1.2)
hist_GBR_emusca.GetYaxis().SetTitle('E_{reco,muon} [MeV]')
canvas.Draw()


# In[ ]:

hist_trueKE = ROOT.TH1D('trueKE', 'title', 100, 0, 5000)
hist_recoKE = ROOT.TH1D('recoKE', 'title', 100, 0, 5000)
hist_recoKE_GBR = ROOT.TH1D('recoKE_GBR', 'title', 100, 0, 5000)
hist_trueKE.SetLineColor(ROOT.kBlack)
hist_recoKE.SetLineColor(ROOT.kRed)
hist_recoKE_GBR.SetLineColor(ROOT.kBlue+2)
hist_trueKE.SetLineWidth(2)
hist_recoKE_GBR.SetLineWidth(2)
fill_hist(hist_trueKE, test_data_trueKE)
fill_hist(hist_recoKE, clf.predict(test_data_reduced_n))
fill_hist(hist_recoKE_GBR, net.predict(test_data_reduced_n))
hist_trueKE.Draw()
#hist_recoKE.Draw("same")
hist_recoKE_GBR.Draw("same")
hist_trueKE.GetXaxis().SetTitle('true or reco KE [MeV]')
hist_trueKE.GetYaxis().SetTitle('Events')
ROOT.gPad.SetLogy()
ROOT.gPad.Draw()


# In[ ]:

#tr.Branch('truKE', test_data_trueKE, 'trueKE/D')
#tr.Branch('recoKEL',clf.predict(test_data_reduced_n),'recoKEL/F')
#tr.Branch('recoKEGB',net.predict(test_data_reduced_n), 'recoKEGB/F')
#tr.Fill()
#tr.Write()


# In[ ]:

hist_trueKE_zoom = ROOT.TH1D('trueKE_zoom', 'title', 100, 0, 2000)
hist_recoKE_zoom = ROOT.TH1D('recoKE_zoom', 'title', 100, 0, 2000)
hist_recoKE_GBR_zoom = ROOT.TH1D('recoKE_GBR_zoom', 'title', 100, 0, 2000)
hist_trueKE_zoom.SetLineColor(ROOT.kBlack)
hist_recoKE_zoom.SetLineColor(ROOT.kRed)
hist_recoKE_GBR_zoom.SetLineColor(ROOT.kBlue+2)
hist_trueKE_zoom.SetLineWidth(2)
hist_recoKE_GBR_zoom.SetLineWidth(2)
fill_hist(hist_trueKE_zoom, test_data_trueKE)
fill_hist(hist_recoKE_zoom, clf.predict(test_data_reduced_n))
fill_hist(hist_recoKE_GBR_zoom, net.predict(test_data_reduced_n))
hist_trueKE_zoom.Draw()
#hist_recoKE_zoom.Draw("same")
hist_recoKE_GBR_zoom.Draw("same")
hist_trueKE_zoom.GetXaxis().SetTitle('true or reco KE [MeV]')
hist_trueKE_zoom.GetYaxis().SetTitle('Events')
ROOT.gPad.SetLogy()
ROOT.gPad.Draw()


# In[1]:

#tr.Write()
#hist_trueKE_zoom.Write()
f.Write()
f.Close()


# In[ ]:



