import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model, ensemble
from sklearn.model_selection import cross_val_score, cross_val_predict
import xgboost as xgb
#import sys
#sys.argv.append( '-b-' )
import ROOT
ROOT.gROOT.SetBatch(True)
from root_numpy import root2array, tree2array, fill_hist
import os

from sklearn import metrics

chain = ROOT.TChain('nu_eneNEW')
#chain.Add('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree.root')
chain.Add('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
for i in range(1040,1099):
    input_file = '/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_'+str(i)+'_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root'
    if os.path.exists(input_file):
        chain.Add(input_file)
#data = tree2array(chain)
#data_reduced   = data[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]#,'hits_pot_length2']]
#data_reduced_n = data_reduced.view(data_reduced.dtype[0]).reshape(data_reduced.shape + (-1,))
#data_target    = data['trueKE']/1e3

E_threshold_lo = 200
E_threshold_hi = 600
select='trueKE>'+str(E_threshold_lo)+'&&'+'trueKE<'+str(E_threshold_hi)
data = tree2array(chain,selection=select )
data_reduced   = data[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]#,'hits_pot_length2']]
data_reduced_n = data_reduced.view(data_reduced.dtype[0]).reshape(data_reduced.shape + (-1,))
data_target    = data['trueKE']/1e3

##########
#params = {'n_estimators': 1000, 'max_depth': 10, 'min_samples_split': 2,
#          'learning_rate': 0.01, 'loss': 'lad'}
params = {'n_estimators': 1000, 'max_depth': 10, 'learning_rate': 0.01, 'loss': 'lad'}
net_hi_E = ensemble.GradientBoostingRegressor(**params)
#net_hi_E.fit(data_reduced_n,data_target)
#net_hi_E
#test_data_recoKE_hi_E_BDTG = net_hi_E.predict(data_reduced_n)

params01 = {'n_estimators': 1000, 'max_depth': 6, 'learning_rate': 0.01, 'loss': 'lad'}
params02 = {'n_estimators': 1000, 'max_depth': 7, 'learning_rate': 0.01, 'loss': 'lad'}
net_hi_E_par01 = ensemble.GradientBoostingRegressor(**params01)
#net_hi_E_par01.fit(arr2_hi_E_n,arr3_hi_E)
#net_hi_E_par01
net_hi_E_par02 = ensemble.GradientBoostingRegressor(**params02)
#net_hi_E_par02.fit(arr2_hi_E_n,arr3_hi_E)
#net_hi_E_par02

scoresBDTG = cross_val_score(net_hi_E, data_reduced_n, data_target, cv=5, scoring='neg_mean_squared_error')
print scoresBDTG
print("BDTG MSE: %0.3f (+/- %0.3f)" % (scoresBDTG.mean(), scoresBDTG.std() * 2))

predictedBDTG = cross_val_predict(net_hi_E, data_reduced_n, data_target, cv=5)

scoresBDTG2 = metrics.mean_squared_error(data_target, predictedBDTG)
print scoresBDTG2
print("BDTG MSE metrics: %0.3f (+/- %0.3f)" % (scoresBDTG2.mean(), scoresBDTG2.std() * 2))

scoresBDTG3 = metrics.median_absolute_error(data_target, predictedBDTG)
print scoresBDTG3
print("BDTG MSE metrics-median: %0.3f (+/- %0.3f)" % (scoresBDTG3.mean(), scoresBDTG3.std() * 2))
#------------------------#
scoresBDTG_01 = cross_val_score(net_hi_E_par01, data_reduced_n, data_target, cv=5, scoring='neg_mean_squared_error')
print scoresBDTG_01
print("BDTG MSE: %0.3f (+/- %0.3f)" % (scoresBDTG_01.mean(), scoresBDTG_01.std() * 2))

predictedBDTG01 = cross_val_predict(net_hi_E_par01, data_reduced_n, data_target, cv=5)
scoresBDTG2_01 = metrics.mean_squared_error(data_target, predictedBDTG01)
print scoresBDTG2_01
print("BDTG MSE metrics: %0.3f (+/- %0.3f)" % (scoresBDTG2_01.mean(), scoresBDTG2_01.std() * 2))

scoresBDTG3_01 = metrics.median_absolute_error(data_target, predictedBDTG01)
print scoresBDTG3_01
print("BDTG MSE metrics-median: %0.3f (+/- %0.3f)" % (scoresBDTG3_01.mean(), scoresBDTG3_01.std() * 2))
#------------------------#
scoresBDTG_02 = cross_val_score(net_hi_E_par02, data_reduced_n, data_target, cv=5, scoring='neg_mean_squared_error')
print scoresBDTG_02
print("BDTG MSE: %0.3f (+/- %0.3f)" % (scoresBDTG_02.mean(), scoresBDTG_02.std() * 2))

predictedBDTG02 = cross_val_predict(net_hi_E_par02, data_reduced_n, data_target, cv=5)
scoresBDTG2_02 = metrics.mean_squared_error(data_target, predictedBDTG02)
print scoresBDTG2_02
print("BDTG MSE metrics: %0.3f (+/- %0.3f)" % (scoresBDTG2_02.mean(), scoresBDTG2_02.std() * 2))

scoresBDTG3_02 = metrics.median_absolute_error(data_target, predictedBDTG02)
print scoresBDTG3_02
print("BDTG MSE metrics-median: %0.3f (+/- %0.3f)" % (scoresBDTG3_02.mean(), scoresBDTG3_02.std() * 2))

#-----------------------------------------------------------------------
ROOT.gStyle.SetOptStat(1)
hist1 = ROOT.TH1D("hist1","hist1", 100, -1, 1)
diff1 = (data_target - predictedBDTG)/data_target
fill_hist(hist1, diff1)
hist1.GetXaxis().SetTitle("#DeltaE/E")
canvas = ROOT.TCanvas()
hist1.Draw()
canvas.SaveAs("xgb_cross_val_DeltaE_200_600MeV_BDTG.pdf")

ROOT.gStyle.SetOptStat(1)
hist01 = ROOT.TH1D("hist01","hist01", 100, -1, 1)
diff01 = (data_target - predictedBDTG01)/data_target
fill_hist(hist01, diff01)
hist01.GetXaxis().SetTitle("#DeltaE/E")
canvas = ROOT.TCanvas()
hist01.Draw()
canvas.SaveAs("xgb_cross_val_DeltaE_200_600MeV_BDTG01.pdf")

ROOT.gStyle.SetOptStat(1)
hist02 = ROOT.TH1D("hist02","hist02", 100, -1, 1)
diff02 = (data_target - predictedBDTG02)/data_target
fill_hist(hist02, diff02)
hist02.GetXaxis().SetTitle("#DeltaE/E")
canvas = ROOT.TCanvas()
hist02.Draw()
canvas.SaveAs("xgb_cross_val_DeltaE_200_600MeV_BDTG02.pdf")
