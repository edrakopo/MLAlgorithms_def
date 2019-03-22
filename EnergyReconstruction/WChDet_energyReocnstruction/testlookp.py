import sys
import ROOT
ROOT.gROOT.SetBatch(True)
from root_numpy import root2array, tree2array, fill_hist
import numpy as np
import matplotlib.pyplot as plt

##function to return value from lookup table:
def lookup(value, dict):
  nearest = sys.maxsize
  result = ""

  for k,v in dict.items():
    if abs(value - k) < nearest:
      nearest = abs(value - k)
      result = v

  return result

rfile = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquan    t_tree_NEWlookupsB_for_training.root");//nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree.root')
intree = rfile.Get('nu_eneNEW');



rfile2 = ROOT.TFile('/Users/edrakopo/work/work_energy_reco/ener_reco_scikit/nu_numu_%d_CCQE_12in_energy_studies_recoquant_tree_NEWl    ookupsB_for_training.root')
intree2 = rfile2.Get('nu_eneNEW');

#print(lookup(60, {1:'a', 100:'b', 1024:'c'}) )
