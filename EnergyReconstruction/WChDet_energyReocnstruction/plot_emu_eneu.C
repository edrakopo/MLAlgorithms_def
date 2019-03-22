#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TF1.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif


void plot_emu_eneu()
{
 
 TH1D *MuEn =new TH1D("MuEn","MuEn; E_{MC, muon} [MeV]", 40, 0., 2000.);
 TH1D *NeuEn =new TH1D("NeuEn","NeuEn; E_{MC, neutrino} [MeV]", 40, 0., 2000.); 

 TFile *input(0);
   char *fname = Form("/Users/edrakopo/work/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root");

  if (!gSystem->AccessPathName( fname )){
    cout<<"Open file "<<endl;
    input = TFile::Open( fname );
  }

  TTree *regTree = (TTree*)input->Get("nu_eneNEW");
  Double_t trueKE, neutrinoE;
  regTree->SetBranchAddress("trueKE", &trueKE);
  regTree->SetBranchAddress("neutrinoE", &neutrinoE);

  for (Long64_t ievt=0; ievt<regTree->GetEntries();ievt++) {
  //for (Long64_t ievt=0; ievt<100;ievt++) {
     if (ievt%1000 == 0) {
        std::cout << "--- ... Processing event: " << ievt << std::endl;
     }
     regTree->GetEntry(ievt);
     MuEn->Fill(trueKE);
     NeuEn->Fill(neutrinoE);
  }

  new TCanvas();
  MuEn->Draw();
  new TCanvas();
  NeuEn->Draw();
}
