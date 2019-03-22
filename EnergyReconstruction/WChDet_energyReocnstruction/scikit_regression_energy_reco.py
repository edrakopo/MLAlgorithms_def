
# coding: utf-8

# In[2]:

import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif', size=80)
import pylab
get_ipython().magic(u'matplotlib inline')
pylab.rcParams['figure.figsize'] = 8, 8


# In[3]:

import ROOT
from root_numpy import root2array, tree2array, fill_hist
from sklearn import linear_model, ensemble


# In[4]:

rfile = ROOT.TFile('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree.root')
intree = rfile.Get('nu_eneNEW')


# In[5]:

arr=tree2array(intree)


# In[6]:

arr2=arr[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','hits_pot_length2','lambda_max_2']]#,'hits_pot_length2']]
arr2_n=arr2.view(arr2.dtype[0]).reshape(arr2.shape + (-1,))
arr3=arr['trueKE']


# In[7]:

clf = linear_model.SGDRegressor()
clf.fit(arr2_n,arr3)
clf


# In[8]:

plt.scatter(arr3,clf.predict(arr2_n)-arr3)
plt.xlabel("trueKE [MeV]")
plt.ylabel("recoKE - trueKE [MeV]")


# In[9]:

chain = ROOT.TChain('nu_eneNEW')
for i in range(1040,1099):
    chain.Add('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_'+str(i)+'_CCQE_12in_energy_studies_recoquant_tree_NEWlookups.root')
test_data = tree2array(chain)


# In[10]:

test_data_reduced = test_data[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','hits_pot_length2','lambda_max_2']]#,'hits_pot_length2']]
test_data_reduced_n = test_data_reduced.view(test_data_reduced.dtype[0]).reshape(test_data_reduced.shape + (-1,))
test_data_trueKE = test_data['trueKE']


# In[11]:

plt.scatter(test_data_trueKE,clf.predict(test_data_reduced_n)-test_data_trueKE)
plt.ylim((-5000,3000))
plt.xlabel("trueKE [MeV]")
plt.ylabel("recoKE - trueKE [MeV]")


# In[12]:

plt.scatter(test_data_trueKE,(clf.predict(test_data_reduced_n)-test_data_trueKE)/test_data_trueKE)
plt.ylim((0,1))
plt.xlabel("trueKE [MeV]")
plt.ylabel("DeltaE/E")
res_twod_SGD = np.dstack((test_data_trueKE, (clf.predict(test_data_reduced_n)-test_data_trueKE)/test_data_trueKE))


# In[13]:

hist_SGD = ROOT.TH2D('name', 'title', 100, 0, 5000, 100, -1, 10)
fill_hist(hist_SGD, res_twod_SGD[0])
hist_SGD.Draw()
ROOT.gPad.Draw()


# In[14]:

profile_SGD = hist_SGD.ProfileX()
profile_SGD.SetLineColor(ROOT.kBlue)
profile_SGD.Draw()
ROOT.gPad.Draw()


# In[15]:

params = {'n_estimators': 1000, 'max_depth': 10, 'min_samples_split': 1,
          'learning_rate': 0.01, 'loss': 'lad'}
net = ensemble.GradientBoostingRegressor(**params)
net.fit(arr2_n,arr3)
net


# In[16]:

plt.scatter(arr3,net.predict(arr2_n)-arr3,c='r')
plt.xlabel("trueKE [MeV]")
plt.ylabel("recoKE - trueKE [MeV]")
plt.ylim(-50,50)


# In[17]:

plt.scatter(test_data_trueKE,(net.predict(test_data_reduced_n)-test_data_trueKE), c='r')
plt.xlabel("trueKE [MeV]")
plt.ylabel("recoKE - trueKE [MeV]")
matrix = np.dstack((test_data_trueKE, (net.predict(test_data_reduced_n)-test_data_trueKE)))


# In[25]:

plt.scatter(test_data_trueKE,((net.predict(test_data_reduced_n)-test_data_trueKE)/test_data_trueKE),c='r')
plt.xlabel("trueKE [MeV]")
plt.ylabel("DeltaE/E")
plt.ylim(-2,2)
twod_GBR_abs = np.dstack((test_data_trueKE, np.abs(net.predict(test_data_reduced_n)-test_data_trueKE)/test_data_trueKE))


# In[42]:

hist_GBR_abs = ROOT.TH2D('name_GBR_abs', 'title', 100, 0, 2000, 100, -1, 10)
fill_hist(hist_GBR_abs, twod_GBR_abs[0])
canvas = ROOT.TCanvas()
hist_GBR_abs.Draw()
hist_GBR_abs.GetXaxis().SetTitle('true KE [MeV]')
hist_GBR_abs.GetYaxis().SetTitle('abs(#Delta E)/E')
canvas.Draw()


# In[43]:

profile_GBR_abs = hist_GBR_abs.ProfileX()
profile_GBR_abs.SetLineColor(ROOT.kBlue)
canvas_prof = ROOT.TCanvas()
profile_GBR_abs.Draw()
profile_GBR_abs.SetMinimum(0)
profile_GBR_abs.SetMaximum(1)
profile_GBR_abs.GetXaxis().SetTitle('true KE [MeV]')
profile_GBR_abs.GetYaxis().SetTitle('abs(#Delta E)/E')
canvas_prof.Draw()


# In[47]:

hist_trueKE = ROOT.TH1D('trueKE', 'title', 100, 0, 5000)
hist_recoKE = ROOT.TH1D('recoKE', 'title', 100, 0, 5000)
hist_recoKE_GBR = ROOT.TH1D('recoKE_GBR', 'title', 100, 0, 5000)
hist_trueKE.SetLineColor(ROOT.kBlack)
hist_recoKE.SetLineColor(ROOT.kBlue)
hist_recoKE_GBR.SetLineColor(ROOT.kRed)
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


# In[48]:

hist_trueKE_zoom = ROOT.TH1D('trueKE_zoom', 'title', 100, 0, 2000)
hist_recoKE_zoom = ROOT.TH1D('recoKE_zoom', 'title', 100, 0, 2000)
hist_recoKE_GBR_zoom = ROOT.TH1D('recoKE_GBR_zoom', 'title', 100, 0, 2000)
hist_trueKE_zoom.SetLineColor(ROOT.kBlack)
hist_recoKE_zoom.SetLineColor(ROOT.kBlue)
hist_recoKE_GBR_zoom.SetLineColor(ROOT.kRed)
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


# In[ ]:



