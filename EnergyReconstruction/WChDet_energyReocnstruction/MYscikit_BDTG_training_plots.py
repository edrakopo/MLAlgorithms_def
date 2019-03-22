import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='Times New Roman', size=20)
import pylab
pylab.rcParams['figure.figsize'] = 10, 6


from sklearn import linear_model, ensemble
import sys
sys.argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
from root_numpy import root2array, tree2array, fill_hist

import random
from sklearn.metrics import mean_squared_error
from sklearn.utils import shuffle


#rfile = ROOT.TFile('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
intree = rfile.Get('nu_eneNEW')

E_threshold = 0 #100
#arr_lo_E = tree2array(intree,selection='trueKE<'+str(E_threshold))
arr_hi_E = tree2array(intree,selection='trueKE>0')#+str(E_threshold))

arr2_hi_E = arr_hi_E[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]#,'hits_pot_length2']]
arr2_hi_E_n = arr2_hi_E.view(arr2_hi_E.dtype[0]).reshape(arr2_hi_E.shape + (-1,))
arr3_hi_E = arr_hi_E['trueKE']


# In[16]: ########## BDTG ###########

params = {'n_estimators': 600, 'max_depth': 10, 'random_state': 100, #'min_samples_split': 1,
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
arr3_hi_E_Rn = arr_hi_E_Rn['trueKE']

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

#----------------
#myhist = ROOT.TH2D('myhist', 'title', 100, E_threshold, 5000, 200, -100, 100)
#fill_hist=(myhist, )
##fill_hist=np.dstack()
#a_train=np.dstack((np.arange(params['n_estimators']) + 1, net_hi_E2.train_score_)
#a_test=np.dstack((np.arange(params['n_estimators']) + 1, test_score)

#c = ROOT.TCanvas()
#c.Divide(1,1)
##c.SetLogy(0)
#c.cd(1)
##c.cd(2)
#myhist.Draw("ColZ")
#c.Draw()
#c.SaveAs("plots/deviation_train_testA.png")
#----------------

plt.plot(np.arange(params['n_estimators']) + 1, net_hi_E2.train_score_, 'b-',
         label='Training Set Deviation')
plt.plot(np.arange(params['n_estimators']) + 1, test_score, 'r-',
         label='Test Set Deviation')
plt.ylim(0,300)
plt.legend(loc='upper right')
plt.xlabel('Number of Estimators')
plt.ylabel('Least Absolute Deviations [MeV]')#,'right')
#plt.ylabel('Mean Squared Error [MeV]','horizontalalignment' : 'right')
plt.savefig("plots/deviation_train_test.png")

fig,ax=plt.subplots(ncols=1, sharey=True)
ax.plot(np.arange(params['n_estimators']) + 1, net_hi_E2.train_score_, 'b-',
       label='Training Set Deviation')
ax.plot(np.arange(params['n_estimators']) + 1, test_score, 'r-',
         label='Test Set Deviation')
ax.set_ylim(0.,300.)
ax.set_xlim(0.,600.)
ax.legend(loc='upper right')
ax.set_ylabel('Least Absolute Deviations [MeV]')
ax.set_xlabel('Number of Estimators')
ax.yaxis.set_label_coords(-0.1, 0.6)
ax.xaxis.set_label_coords(0.85, -0.08)
plt.savefig("plots/deviation_train_test2.png")
# 
#######################








