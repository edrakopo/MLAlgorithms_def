import ROOT
ROOT.gROOT.SetBatch(True)
#from root_numpy import root2array, tree2array, fill_hist

import random

#from root_numpy.tmva import add_regression_events, evaluate_reader
#from root_numpy import ROOT_VERSION
from ROOT import TMVA, TFile, TCut, TString

#ROOT.TMVA.Tools.Instance()
#ROOT.TMVA.Types.Instance()
#root_data = ROOT.TFile.Open('/Users/edrakopo/work/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root').Get('nu_eneNEW')
root_data = ROOT.TFile.Open('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root').Get('nu_eneNEW')

########### BDTG-TMVA ###########
outputFile = ROOT.TFile("TMVAReginpy.root","RECREATE")

#factory = ROOT.TMVA.Factory("TMVARegression", fout,
#                            ":".join([
#                                "!V",
#                                "!Silent",
#                                "Color",
#                                "DrawProgressBar",
#                                ] ))

factory = TMVA.Factory("TMVARegression", outputFile, "!V:!Silent:Color:DrawProgressBar")

#dataloader = TMVA.DataLoader("nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training")
#dataloader = TMVA.DataLoader("dl")

#factory.SetVerbose( True )
#dataloader.AddVariable("total_hits2", 'F')
#dataloader.AddVariable("total_ring_PEs2", 'F')
#dataloader.AddVariable("recoDWallR2", 'F')
#dataloader.AddVariable("recoDWallZ2", 'F')
#dataloader.AddVariable("lambda_max_2", 'F')
#dataloader.AddTarget( "trueKE" )

factory.AddVariable("total_hits2", 'F')
factory.AddVariable("total_ring_PEs2", 'F')
factory.AddVariable("recoDWallR2", 'F')
factory.AddVariable("recoDWallZ2", 'F')
factory.AddVariable("lambda_max_2", 'F')
factory.AddTarget( "trueKE" ) 

#root_data = ROOT.TFile.Open('/Users/edrakopo/work/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')
root_data = ROOT.TFile.Open('/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root')

regTree = root_data.Get('nu_eneNEW')
regWeight=1.
factory.AddRegressionTree( regTree, regWeight )
mycut=ROOT.TCut("")
factory.PrepareTrainingAndTestTree( mycut,"nTrain_Regression=44420:nTest_Regression=19038:SplitMode=Random:NormMode=NumEvents:!V" )

factory.BookMethod( TMVA.Types.kBDT, "BDTG", "!H:!V:NegWeightTreatment=NoNegWeightsInTraining:NTrees=10000:BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3")

#factory.BookMethod( dataloader, ROOT.TString('kBDT'), ROOT.TString("BDTG") )
#factory.BookMethod( dataloader, 'BDT', 'BDTG', '' )
#factory.BookMethod( dataloader, TMVA.Types.kBDT, "BDTG", "!H,!V,NegWeightTreatment=NoNegWeightsInTraining,NTrees=10000,BoostType=Grad,Shrinkage=0.1,UseBaggedBoost,BaggedSampleFraction=0.5,nCuts=20,MaxDepth=3" )
#factory.BookMethod( dataloader, TMVA.Types.kBDT, "BDTG",":".join(["!H,!V,NegWeightTreatment=NoNegWeightsInTraining,NTrees=10000,BoostType=Grad,Shrinkage=0.1,UseBaggedBoost,BaggedSampleFraction=0.5,nCuts=20,MaxDepth=3"]) )
#factory.BookMethod( dataloader, TMVA.Types.kBDT, "BDTG", "!H,!V,NegWeightTreatment=NoNegWeightsInTraining,NTrees=10000,BoostType=Grad,Shrinkage=0.1,UseBaggedBoost,BaggedSampleFraction=0.5,nCuts=20,MaxDepth=3" )
#factory.BookMethod( dataloader, BDT,"BDTG","!H,!V,NegWeightTreatment=NoNegWeightsInTraining,NTrees=10000,BoostType=Grad,Shrinkage=0.1,UseBaggedBoost,BaggedSampleFraction=0.5,nCuts=20,MaxDepth=3" )

factory.TrainAllMethods()
factory.TestAllMethods()

outputFile.Close()

#######################








