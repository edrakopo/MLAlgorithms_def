import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif', size=80)
import pylab
pylab.rcParams['figure.figsize'] = 16, 8


# In[3]:

from sklearn import linear_model, ensemble
import sys
sys.argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
from root_numpy import root2array, tree2array, fill_hist

import random
from sklearn.metrics import mean_squared_error
from sklearn.utils import shuffle

# In[4]:

#rfile = ROOT.TFile('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree.root')
rfile = ROOT.TFile('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
intree = rfile.Get('nu_eneNEW')

E_threshold = 0 #100
#arr_lo_E = tree2array(intree,selection='neutrinoE<'+str(E_threshold))
arr_hi_E = tree2array(intree,selection='neutrinoE>0')#+str(E_threshold))

# In[7]:

#arr2_lo_E = arr_lo_E[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]#,'hits_pot_length2']]
#arr2_lo_E_n = arr2_lo_E.view(arr2_lo_E.dtype[0]).reshape(arr2_lo_E.shape + (-1,))
#arr3_lo_E = arr_lo_E['neutrinoE']
arr2_hi_E = arr_hi_E[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]#,'hits_pot_length2']]
arr2_hi_E_n = arr2_hi_E.view(arr2_hi_E.dtype[0]).reshape(arr2_hi_E.shape + (-1,))
arr3_hi_E = arr_hi_E['neutrinoE']


# In[8]:

#clf_lo_E = linear_model.SGDRegressor()
#clf_lo_E.fit(arr2_lo_E_n,arr3_lo_E)
#clf_lo_E
clf_hi_E = linear_model.SGDRegressor()
clf_hi_E.fit(arr2_hi_E_n,arr3_hi_E)
clf_hi_E


# In[9]:
#f, (ax1, ax2) = plt.subplots(1, 2)
f, (ax1, ax2) = plt.subplots(1, 2)
#ax1.scatter(arr3_lo_E,clf_lo_E.predict(arr2_lo_E_n)-arr3_lo_E)
ax1.set_xlabel("neutrinoE [MeV]")
ax1.set_ylabel("recoKE - neutrinoE [MeV]")
ax2.scatter(arr3_hi_E,clf_hi_E.predict(arr2_hi_E_n)-arr3_hi_E)
ax2.set_xlabel("neutrinoE [MeV]")
ax2.set_ylabel("recoKE - neutrinoE [MeV]")


# In[10]:

chain = ROOT.TChain('nu_eneNEW')
for i in range(1040,1099):
    chain.Add('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_'+str(i)+'_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
#test_data_lo_E = tree2array(chain, selection='neutrinoE<'+str(E_threshold))
test_data_hi_E = tree2array(chain, selection='neutrinoE')#>'+str(E_threshold))


# In[11]:

#test_data_reduced_lo_E = test_data_lo_E[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]#,'hits_pot_length2']]
#test_data_reduced_lo_E_n = test_data_reduced_lo_E.view(test_data_reduced_lo_E.dtype[0]).reshape(test_data_reduced_lo_E.shape + (-1,))
#test_data_neutrinoE_lo_E = test_data_lo_E['neutrinoE']
test_data_reduced_hi_E = test_data_hi_E[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]#,'hits_pot_length2']]
test_data_reduced_hi_E_n = test_data_reduced_hi_E.view(test_data_reduced_hi_E.dtype[0]).reshape(test_data_reduced_hi_E.shape + (-1,))
test_data_neutrinoE_hi_E = test_data_hi_E['neutrinoE']

# In[12]:

f, (ax1, ax2) = plt.subplots(1, 2)
#ax1.scatter(test_data_neutrinoE_lo_E,clf_lo_E.predict(test_data_reduced_lo_E_n)-test_data_neutrinoE_lo_E)
ax2.scatter(test_data_neutrinoE_hi_E,clf_hi_E.predict(test_data_reduced_hi_E_n)-test_data_neutrinoE_hi_E)
ax1.set_ylim((-100,100))
ax1.set_xlabel("neutrinoE [MeV]")
ax1.set_ylabel("recoKE - neutrinoE [MeV]")
ax2.set_ylim((-5000,1000))
ax2.set_xlabel("neutrinoE [MeV]")
ax2.set_ylabel("recoKE - neutrinoE [MeV]")

test_data_recoKE_hi_E = clf_hi_E.predict(test_data_reduced_hi_E_n)
outputFile = ROOT.TFile.Open("scikit_recoLAST_Eneu.root", "RECREATE")
#outputTuple = ROOT.TNtuple("tuple", "tuple", "neutrinoE:recoKE")
#for i in range(len(test_data_recoKE_hi_E)):
#    outputTuple.Fill(test_data_neutrinoE_hi_E[i], test_data_recoKE_hi_E[i])
#outputTuple.Write()
#outputFile.Close()


#res_twod_SGD_lo_E = np.dstack((test_data_neutrinoE_lo_E, np.abs(clf_lo_E.predict(test_data_reduced_lo_E_n)-test_data_neutrinoE_lo_E)/test_data_neutrinoE_lo_E))
#res_twod_SGD_hi_E = np.dstack((test_data_neutrinoE_hi_E, np.abs(clf_hi_E.predict(test_data_reduced_hi_E_n)-test_data_neutrinoE_hi_E)/test_data_neutrinoE_hi_E))
res_twod_SGD_hi_E = np.dstack((test_data_neutrinoE_hi_E, 100.*(clf_hi_E.predict(test_data_reduced_hi_E_n)-test_data_neutrinoE_hi_E)/test_data_neutrinoE_hi_E))

# In[14]:

#hist_SGD_lo_E = ROOT.TH2D('name_lo', 'title', 100, 0, E_threshold, 100, -1, 10)
#fill_hist(hist_SGD_lo_E, res_twod_SGD_lo_E[0])
hist_SGD_hi_E = ROOT.TH2D('name_hi', 'title', 100, E_threshold, 5000, 200, -100, 100)
fill_hist(hist_SGD_hi_E, res_twod_SGD_hi_E[0])
c = ROOT.TCanvas()
c.Divide(2,1)
c.SetLogy(0)
c.cd(1)
#hist_SGD_lo_E.Draw()
c.cd(2)
hist_SGD_hi_E.Draw("ColZ")
c.Draw()
c.SaveAs("plots/SGDRegDE_E.png")

# In[15]:

#profile_SGD_lo_E = hist_SGD_lo_E.ProfileX()
#profile_SGD_lo_E.SetLineColor(ROOT.kBlue)
profile_SGD_hi_E = hist_SGD_hi_E.ProfileX()
profile_SGD_hi_E.SetLineColor(ROOT.kBlue)
c1 = ROOT.TCanvas()
c1.Divide(2,1)
c1.SetLogy(0)
c1.cd(1)
#profile_SGD_lo_E.Draw()
c1.cd(2)
profile_SGD_hi_E.Draw()
c1.Draw()
c1.SetLogy(0)
c1.Draw()
c1.SaveAs("plots/SGDRegDE_E_profX.png")


# In[16]: ########## BDTG ###########

params = {'n_estimators': 1000, 'max_depth': 10, 'min_samples_split': 1,
          'learning_rate': 0.01, 'loss': 'lad'} 
#net_lo_E = ensemble.GradientBoostingRegressor(**params)
#net_lo_E.fit(arr2_lo_E_n,arr3_lo_E)
#net_lo_E
net_hi_E = ensemble.GradientBoostingRegressor(**params)
net_hi_E.fit(arr2_hi_E_n,arr3_hi_E)
net_hi_E

#adding train/test sample with random selection and check deviations:
#arr_hi_E_Rn = random.shuffle(arr_hi_E)  
arr_hi_E_Rn = shuffle(arr_hi_E, random_state=13) #split randomly the elements
arr2_hi_E_Rn = arr_hi_E_Rn[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]
arr2_hi_E_Rn_n = arr2_hi_E_Rn.view(arr2_hi_E_Rn.dtype[0]).reshape(arr2_hi_E_Rn.shape + (-1,))
arr3_hi_E_Rn = arr_hi_E_Rn['neutrinoE']

offset = int(arr2_hi_E_Rn_n.shape[0] * 0.7) #Y.shape[0] returns the number of rows in Y 
arr2_hi_E_train, arr3_hi_E_train = arr2_hi_E_Rn_n[:offset], arr3_hi_E_Rn[:offset]  # train sample
arr2_hi_E_test, arr3_hi_E_test   = arr2_hi_E_Rn_n[offset:], arr3_hi_E_Rn[offset:]  # test sample

net_hi_E2 = ensemble.GradientBoostingRegressor(**params)
net_hi_E2.fit(arr2_hi_E_train, arr3_hi_E_train)

mse = mean_squared_error(arr3_hi_E_test, net_hi_E2.predict(arr2_hi_E_test))
print("MSE: %.4f" % mse)

test_score = np.zeros((params['n_estimators'],), dtype=np.float64)

for i, y_pred in enumerate(net_hi_E2.staged_predict(arr2_hi_E_test)):
    test_score[i] = net_hi_E2.loss_(arr3_hi_E_test, y_pred)

plt.figure(figsize=(10, 6))
#plt.subplot(1, 2, 1)
#plt.title('Deviance')
plt.plot(np.arange(params['n_estimators']) + 1, net_hi_E2.train_score_, 'b-',
         label='Training Set Deviance')
plt.plot(np.arange(params['n_estimators']) + 1, test_score, 'r-',
         label='Test Set Deviance')
#plt.legend(loc='upper right')
#plt.xlabel('Boosting Iterations')
#plt.ylabel('Deviance')
plt.savefig("plots/deviation_train_test.png")

params01 = {'n_estimators': 2000, 'max_depth': 6, 'min_samples_split': 1,
           'learning_rate': 0.01, 'loss': 'lad'}
params02 = {'n_estimators': 2000, 'max_depth': 7, 'min_samples_split': 1,
           'learning_rate': 0.01, 'loss': 'lad'}

net_hi_E_par01 = ensemble.GradientBoostingRegressor(**params01)
net_hi_E_par01.fit(arr2_hi_E_n,arr3_hi_E)
net_hi_E_par01
net_hi_E_par02 = ensemble.GradientBoostingRegressor(**params02)
net_hi_E_par02.fit(arr2_hi_E_n,arr3_hi_E)
net_hi_E_par02
#params = {'n_estimators': 1000, 'max_depth': 10, 'min_samples_split': 1,'learning_rate': 0.01, 'loss': 'lad'}
# params0 = {'n_estimators': 1000, 'max_depth': 5, 'min_samples_split': 1,
#           'learning_rate': 0.01, 'loss': 'lad'}
# print 'params0'
# params1 = {'n_estimators': 1000, 'max_depth': 7, 'min_samples_split': 1,
#           'learning_rate': 0.01, 'loss': 'lad'}
# print 'params1'
# params2 = {'n_estimators': 2000, 'max_depth': 10, 'min_samples_split': 1,
#           'learning_rate': 0.01, 'loss': 'lad'}
# print 'params2'
# params3 = {'n_estimators': 5000, 'max_depth': 10, 'min_samples_split': 1,
#           'learning_rate': 0.01, 'loss': 'lad'}
# print 'params3'
# params4 = {'n_estimators': 10000, 'max_depth': 10, 'min_samples_split': 1,
#           'learning_rate': 0.01, 'loss': 'lad'}
# print 'params4'
# params5 = {'n_estimators': 1000, 'max_depth': 9, 'min_samples_split': 1,
#           'learning_rate': 0.01, 'loss': 'lad'}
# print 'params5'
# 
# net_hi_E_par0 = ensemble.GradientBoostingRegressor(**params0)
# net_hi_E_par0.fit(arr2_hi_E_n,arr3_hi_E)
# net_hi_E_par0
# net_hi_E_par1 = ensemble.GradientBoostingRegressor(**params1)
# net_hi_E_par1.fit(arr2_hi_E_n,arr3_hi_E)
# net_hi_E_par1
# net_hi_E_par2 = ensemble.GradientBoostingRegressor(**params2)
# net_hi_E_par2.fit(arr2_hi_E_n,arr3_hi_E)
# net_hi_E_par2
# net_hi_E_par3 = ensemble.GradientBoostingRegressor(**params3)
# net_hi_E_par3.fit(arr2_hi_E_n,arr3_hi_E)
# net_hi_E_par3
# net_hi_E_par4 = ensemble.GradientBoostingRegressor(**params4)
# net_hi_E_par4.fit(arr2_hi_E_n,arr3_hi_E)
# net_hi_E_par4
# net_hi_E_par5 = ensemble.GradientBoostingRegressor(**params5)
# net_hi_E_par5.fit(arr2_hi_E_n,arr3_hi_E)
# net_hi_E_par5
# 
test_data_recoKE_hi_E_BDTG = net_hi_E.predict(test_data_reduced_hi_E_n)
test_data_recoKE_hi_E_BDTG_par01 = net_hi_E_par01.predict(test_data_reduced_hi_E_n)
test_data_recoKE_hi_E_BDTG_par02 = net_hi_E_par02.predict(test_data_reduced_hi_E_n)
#test_data_recoKE_hi_E_BDTG2 = net_hi_E2.predict(test_data_reduced_hi_E_n)
#test_data_recoKE_hi_E_BDTG_par0 = net_hi_E_par0.predict(test_data_reduced_hi_E_n)
#test_data_recoKE_hi_E_BDTG_par1 = net_hi_E_par1.predict(test_data_reduced_hi_E_n)
#test_data_recoKE_hi_E_BDTG_par2 = net_hi_E_par2.predict(test_data_reduced_hi_E_n)
#test_data_recoKE_hi_E_BDTG_par3 = net_hi_E_par3.predict(test_data_reduced_hi_E_n)
#test_data_recoKE_hi_E_BDTG_par4 = net_hi_E_par4.predict(test_data_reduced_hi_E_n)
#test_data_recoKE_hi_E_BDTG_par5 = net_hi_E_par5.predict(test_data_reduced_hi_E_n)

#write everything in tree:
#outputTuple = ROOT.TNtuple("tuple", "tuple", "neutrinoE:recoKE:recoKE_BDTG:recoKE_BDTG2:recoKE_BDTG_par0:recoKE_BDTG_par1:recoKE_BDTG_par2:recoKE_BDTG_par3:recoKE_BDTG_par4:recoKE_BDTG_par5")
outputTuple = ROOT.TNtuple("tuple", "tuple", "neutrinoE:recoKE:recoKE_BDTG:recoKE_BDTG01:recoKE_BDTG_par02")

for i in range(len(test_data_recoKE_hi_E)):
     outputTuple.Fill(test_data_neutrinoE_hi_E[i], test_data_recoKE_hi_E[i], test_data_recoKE_hi_E_BDTG[i], test_data_recoKE_hi_E_BDTG_par01[i], test_data_recoKE_hi_E_BDTG_par02[i])
 #   outputTuple.Fill(test_data_neutrinoE_hi_E[i], test_data_recoKE_hi_E[i], test_data_recoKE_hi_E_BDTG[i], test_data_recoKE_hi_E_BDTG2[i], test_data_recoKE_hi_E_BDTG_par0[i], test_data_recoKE_hi_E_BDTG_par1[i], test_data_recoKE_hi_E_BDTG_par2[i], test_data_recoKE_hi_E_BDTG_par3[i], test_data_recoKE_hi_E_BDTG_par4[i], test_data_recoKE_hi_E_BDTG_par5[i],)

outputTuple.Write()
outputFile.Close()
#######################
# In[17]:

f, (ax1, ax2) = plt.subplots(1, 2)
#ax1.scatter(arr3_lo_E,net_lo_E.predict(arr2_lo_E_n)-arr3_lo_E,c='r')
ax2.scatter(arr3_hi_E,net_hi_E.predict(arr2_hi_E_n)-arr3_hi_E,c='r')
ax1.set_xlabel("neutrinoE [MeV]")
ax1.set_ylabel("recoKE - neutrinoE [MeV]")
ax2.set_xlabel("neutrinoE [MeV]")
ax2.set_ylabel("recoKE - neutrinoE [MeV]")


# In[18]:

f, (ax1, ax2) = plt.subplots(1, 2)
#ax1.scatter(test_data_neutrinoE_lo_E,(net_lo_E.predict(test_data_reduced_lo_E_n)-test_data_neutrinoE_lo_E), c='r')
ax2.scatter(test_data_neutrinoE_hi_E,(net_hi_E.predict(test_data_reduced_hi_E_n)-test_data_neutrinoE_hi_E), c='r')
ax1.set_xlabel("neutrinoE [MeV]")
ax1.set_ylabel("recoKE - neutrinoE [MeV]")
ax2.set_xlabel("neutrinoE [MeV]")
ax2.set_ylabel("recoKE - neutrinoE [MeV]")
#matrix_lo_E = np.dstack((test_data_neutrinoE_lo_E , (net_lo_E.predict(test_data_reduced_lo_E_n)-test_data_neutrinoE_lo_E )))
matrix_hi_E = np.dstack((test_data_neutrinoE_hi_E , (net_hi_E.predict(test_data_reduced_hi_E_n)-test_data_neutrinoE_hi_E )))


# In[29]:

f, (ax1, ax2) = plt.subplots(1, 2)
#ax1.scatter(test_data_neutrinoE_lo_E,(np.abs(net_lo_E.predict(test_data_reduced_lo_E_n)-test_data_neutrinoE_lo_E)/test_data_neutrinoE_lo_E),c='r')
ax2.scatter(test_data_neutrinoE_hi_E,(np.abs(net_hi_E.predict(test_data_reduced_hi_E_n)-test_data_neutrinoE_hi_E)/test_data_neutrinoE_hi_E),c='r')
ax1.set_xlabel("neutrinoE [MeV]")
ax1.set_ylabel("DeltaE/E")
ax1.set_ylim(0,2)
ax2.set_xlabel("neutrinoE [MeV]")
ax2.set_ylabel("DeltaE/E")
ax2.set_ylim(0,2)
#twod_GBR_abs_lo_E = np.dstack((test_data_neutrinoE_lo_E, np.abs(net_lo_E.predict(test_data_reduced_lo_E_n)-test_data_neutrinoE_lo_E)/test_data_neutrinoE_lo_E))
#twod_GBR_abs_hi_E = np.dstack((test_data_neutrinoE_hi_E, np.abs(net_hi_E.predict(test_data_reduced_hi_E_n)-test_data_neutrinoE_hi_E)/test_data_neutrinoE_hi_E))
twod_GBR_abs_hi_E = np.dstack((test_data_neutrinoE_hi_E, 100.*(net_hi_E.predict(test_data_reduced_hi_E_n)-test_data_neutrinoE_hi_E)/test_data_neutrinoE_hi_E))
twod_GBR_abs_hi_E2 = np.dstack((test_data_neutrinoE_hi_E, 100.*(net_hi_E2.predict(test_data_reduced_hi_E_n)-test_data_neutrinoE_hi_E)/test_data_neutrinoE_hi_E))

# In[30]:

#hist_GBR_abs_lo_E = ROOT.TH2D('name_GBR_abs_lo_E', 'title', 50, 0, E_threshold, 100, -1, 10)
hist_GBR_abs_hi_E2 = ROOT.TH2D('name_GBR_abs_hi_E2', 'title', 100, E_threshold, 5000, 200, -100, 100)
hist_GBR_abs_hi_E = ROOT.TH2D('name_GBR_abs_hi_E', 'title', 100, E_threshold, 5000, 200, -100, 100)
#fill_hist(hist_GBR_abs_lo_E, twod_GBR_abs_lo_E[0])
fill_hist(hist_GBR_abs_hi_E2, twod_GBR_abs_hi_E2[0])
fill_hist(hist_GBR_abs_hi_E, twod_GBR_abs_hi_E[0])
canvas = ROOT.TCanvas()
canvas.Divide(2,1)
canvas.cd(1)
hist_GBR_abs_hi_E2.Draw()
hist_GBR_abs_hi_E2.GetXaxis().SetTitle('true KE [MeV]')
hist_GBR_abs_hi_E2.GetYaxis().SetTitle('abs(#Delta E)/E')
#hist_GBR_abs_lo_E.Draw()
#hist_GBR_abs_lo_E.GetXaxis().SetTitle('true KE [MeV]')
#hist_GBR_abs_lo_E.GetYaxis().SetTitle('abs(#Delta E)/E')
canvas.cd(2)
hist_GBR_abs_hi_E.Draw()
hist_GBR_abs_hi_E.GetXaxis().SetTitle('true KE [MeV]')
hist_GBR_abs_hi_E.GetYaxis().SetTitle('abs(#Delta E)/E')
canvas.Draw()
canvas.SaveAs("plots/BDTGscikitDE_E.png")

# In[31]:

#profile_GBR_abs_lo_E = hist_GBR_abs_lo_E.ProfileX()
#profile_GBR_abs_lo_E.SetLineColor(ROOT.kBlue+2)
#profile_GBR_abs_lo_E.SetMarkerColor(ROOT.kBlue+2)
#profile_GBR_abs_lo_E.SetLineWidth(1)
profile_GBR_abs_hi_E2 = hist_GBR_abs_hi_E2.ProfileX()
profile_GBR_abs_hi_E2.SetLineColor(ROOT.kBlue+2)
profile_GBR_abs_hi_E2.SetMarkerColor(ROOT.kBlue+2)
profile_GBR_abs_hi_E2.SetLineWidth(1)
profile_GBR_abs_hi_E = hist_GBR_abs_hi_E.ProfileX()
profile_GBR_abs_hi_E.SetLineColor(ROOT.kBlue+2)
profile_GBR_abs_hi_E.SetMarkerColor(ROOT.kBlue+2)
profile_GBR_abs_hi_E.SetLineWidth(1)
canvas_prof = ROOT.TCanvas()
canvas_prof.Divide(2,1)
canvas_prof.cd(1)
profile_GBR_abs_hi_E2.Draw("ColZ")
#profile_GBR_abs_lo_E.Draw()
#profile_GBR_abs_lo_E.SetMinimum(0)
#profile_GBR_abs_lo_E.SetMaximum(1)
#profile_GBR_abs_lo_E.GetXaxis().SetTitle('true KE [MeV]')
#profile_GBR_abs_lo_E.GetYaxis().SetTitle('abs(#Delta E)/E')
canvas_prof.cd(2)
profile_GBR_abs_hi_E.Draw("ColZ")
#profile_GBR_abs_hi_E.SetMinimum(0)
#profile_GBR_abs_hi_E.SetMaximum(1)
profile_GBR_abs_hi_E.GetXaxis().SetTitle('true KE [MeV]')
profile_GBR_abs_hi_E.GetYaxis().SetTitle('abs(#Delta E)/E')
canvas_prof.Draw()
canvas_prof.SaveAs("plots/BDTGscikitDE_E_profX.png")


# In[24]:

hist_neutrinoE = ROOT.TH1D('neutrinoE', 'title', 100, 0, 5000)
#hist_recoKE_lo_E = ROOT.TH1D('recoKE', 'title', 4, 0, E_threshold)
hist_recoKE_hi_E = ROOT.TH1D('recoKE_GBR', 'title', 96, E_threshold, 5000)
hist_neutrinoE.SetLineColor(ROOT.kBlack)
#hist_recoKE_lo_E.SetLineColor(ROOT.kRed+2)
hist_recoKE_hi_E.SetLineColor(ROOT.kBlue+2)
hist_neutrinoE.SetLineWidth(2)
#hist_recoKE_lo_E.SetLineWidth(2)
hist_recoKE_hi_E.SetLineWidth(2)
fill_hist(hist_neutrinoE, test_data_neutrinoE_hi_E)
#fill_hist(hist_recoKE_lo_E, net_lo_E.predict(test_data_reduced_lo_E_n))
fill_hist(hist_recoKE_hi_E, net_hi_E.predict(test_data_reduced_hi_E_n))
c2 = ROOT.TCanvas()
hist_neutrinoE.Draw()
#hist_recoKE_lo_E.Draw("same")
hist_recoKE_hi_E.Draw("same")
hist_neutrinoE.GetXaxis().SetTitle('true or reco KE [MeV]')
hist_neutrinoE.GetYaxis().SetTitle('Events')
c2.SetLogy()
c2.Draw()
c2.SaveAs("plots/BDTG_MCEmu.png")


# In[ ]:

hist_neutrinoE_zoom = ROOT.TH1D('neutrinoE_zoom', 'title', 100, 0, 2000)
hist_recoKE_zoom = ROOT.TH1D('recoKE_zoom', 'title', 100, 0, 2000)
hist_recoKE_GBR_zoom = ROOT.TH1D('recoKE_GBR_zoom', 'title', 100, 0, 2000)
hist_neutrinoE_zoom.SetLineColor(ROOT.kBlack)
hist_recoKE_zoom.SetLineColor(ROOT.kRed)
hist_recoKE_GBR_zoom.SetLineColor(ROOT.kBlue+2)
hist_neutrinoE_zoom.SetLineWidth(2)
hist_recoKE_GBR_zoom.SetLineWidth(2)
#fill_hist(hist_neutrinoE_zoom, test_data_neutrinoE)
#fill_hist(hist_recoKE_zoom, clf.predict(test_data_reduced_n))
#fill_hist(hist_recoKE_GBR_zoom, net.predict(test_data_reduced_n))
hist_neutrinoE_zoom.Draw()
#hist_recoKE_zoom.Draw("same")
hist_recoKE_GBR_zoom.Draw("same")
hist_neutrinoE_zoom.GetXaxis().SetTitle('true or reco KE [MeV]')
hist_neutrinoE_zoom.GetYaxis().SetTitle('Events')
ROOT.gPad.SetLogy()
ROOT.gPad.Draw()



# In[ ]:

#net.feature_importances_








