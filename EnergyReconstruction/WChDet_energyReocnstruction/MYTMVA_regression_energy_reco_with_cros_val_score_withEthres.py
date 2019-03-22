import ROOT
ROOT.gROOT.SetBatch(True)
import random
#from root_numpy.tmva import add_regression_events, evaluate_reader
#from root_numpy import ROOT_VERSION
from ROOT import TMVA, TFile, TCut, TString

import sys
sys.path.insert(0,'/usr/local/lib/python2.7/site-packages/')
import sklearn
import numpy as np
import matplotlib.pyplot as plt
from sklearn import linear_model, ensemble, metrics
from sklearn.model_selection import cross_val_score, cross_val_predict
from root_numpy import root2array, tree2array, fill_hist
import os


chain = ROOT.TChain('nu_eneNEW')
#chain.Add('/Users/edrakopo/work/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
chain.Add('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
for i in range(1040,1099):
    #input_file = '/Users/edrakopo/work/ener_reco_scikit/nu_numu_'+str(i)+'_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root'
    input_file = '/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_'+str(i)+'_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root'
    if os.path.exists(input_file):
        chain.Add(input_file)

E_threshold_lo = 200
E_threshold_hi = 600
select='trueKE>'+str(E_threshold_lo)+'&&'+'trueKE<'+str(E_threshold_hi)
data = tree2array(chain,selection=select )
data_reduced   = data[['total_hits2','total_ring_PEs2','recoDWallR2','recoDWallZ2','lambda_max_2']]#,'hits_pot_length2']]
print 'data_reduced.shape: ', data_reduced.shape
data_reduced_n = data_reduced.copy(data_reduced.dtype[0]).reshape(data_reduced.shape + (-1,))
print 'data_reduced_n.shape: ', data_reduced_n.shape
data_target    = data['trueKE']/1e3

#outputFile = ROOT.TFile("TMVAReginpy_test.root","RECREATE")
#factory = TMVA.Factory("TMVARegression", outputFile, "!V:!Silent:Color:DrawProgressBar")
#tmva.evaluate_method("BDTG",data_reduced)
#tmva.evaluate_method("TMVA::MethodBDT",data_reduced)
outputFile = ROOT.TFile("TMVAReginpytest.root","RECREATE")
factory = TMVA.Factory("TMVARegression", outputFile, "!V:!Silent:Color:DrawProgressBar")

dataloader = TMVA.DataLoader("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training")
#dataloader = TMVA.DataLoader("dl")

#factory.SetVerbose( True )
dataloader.AddVariable("total_hits2", 'F')
dataloader.AddVariable("total_ring_PEs2", 'F')
dataloader.AddVariable("recoDWallR2", 'F')
dataloader.AddVariable("recoDWallZ2", 'F')
dataloader.AddVariable("lambda_max_2", 'F')
dataloader.AddTarget( "trueKE" )

#factory.AddVariable("total_hits2", 'F')
#factory.AddVariable("total_ring_PEs2", 'F')
#factory.AddVariable("recoDWallR2", 'F')
#factory.AddVariable("recoDWallZ2", 'F')
#factory.AddVariable("lambda_max_2", 'F')
#factory.AddTarget( "trueKE" )

#root_data = ROOT.TFile.Open('/Users/edrakopo/work/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
root_data =ROOT.TFile.Open('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')

regTree = root_data.Get('nu_eneNEW')
regWeight=1.
factory->AddRegressionTree( regTree, regWeight )
mycut=ROOT.TCut("")
factory.PrepareTrainingAndTestTree( mycut,"nTrain_Regression=44420:nTest_Regression=19038:SplitMode=Random:NormMode=NumEvents:!V" )

tmvaBDTG=factory.BookMethod( TMVA.Types.kBDT, "BDTG", "!H:!V:NegWeightTreatment=NoNegWeightsInTraining:NTrees=10000:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3")

scoresBDTG = cross_val_score(tmvaBDTG, data_reduced_n, data_target, cv=5, scoring='neg_mean_squared_error')
print scoresBDTG
print("BDTG MSE: %0.3f (+/- %0.3f)" % (scoresBDTG.mean(), scoresBDTG.std() * 2))

predictedBDTG = cross_val_predict(tmvaBDTG, data_reduced_n, data_target, cv=5)

scoresBDTG2 = metrics.mean_squared_error(data_target, predictedBDTG)
print scoresBDTG2
print("BDTG MSE metrics: %0.3f (+/- %0.3f)" % (scoresBDTG2.mean(), scoresBDTG2.std() * 2))

scoresBDTG3 = metrics.median_absolute_error(data_target, predictedBDTG)
print scoresBDTG3
print("BDTG MSE metrics-median: %0.3f (+/- %0.3f)" % (scoresBDTG3.mean(), scoresBDTG3.std() * 2))


#factory.EvaluateRegression()

################################
#ROOT.TMVA.Tools.Instance()
#ROOT.TMVA.Types.Instance()
#root_data = ROOT.TFile.Open('/Users/edrakopo/work/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root').Get('nu_eneNEW')

########### BDTG-TMVA ###########
#outputFile = ROOT.TFile("TMVAReginpy.root","RECREATE")

#factory = TMVA.Factory("TMVARegression", outputFile, "!V:!Silent:Color:DrawProgressBar")

#factory.AddVariable("total_hits2", 'F')
#factory.AddVariable("total_ring_PEs2", 'F')
#factory.AddVariable("recoDWallR2", 'F')
#factory.AddVariable("recoDWallZ2", 'F')
#factory.AddVariable("lambda_max_2", 'F')
#factory.AddTarget( "trueKE" ) 

#root_data = ROOT.TFile.Open('/Users/edrakopo/work/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')

#regTree = root_data.Get('nu_eneNEW')
#regWeight=1.
#factory.AddRegressionTree( regTree, regWeight )
#mycut=ROOT.TCut("")
#factory.PrepareTrainingAndTestTree( mycut,"nTrain_Regression=44420:nTest_Regression=19038:SplitMode=Random:NormMode=NumEvents:!V" )

#factory.BookMethod( TMVA.Types.kBDT, "BDTG", "!H:!V:NegWeightTreatment=NoNegWeightsInTraining:NTrees=10000:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3")


#factory.TrainAllMethods()
#factory.TestAllMethods()

#outputFile.Close()
#######################

