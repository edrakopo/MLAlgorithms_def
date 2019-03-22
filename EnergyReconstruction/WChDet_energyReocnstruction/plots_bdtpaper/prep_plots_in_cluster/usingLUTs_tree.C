#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TF1.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

void usingLUTs_tree()
{
  cout<<"starting.... "<<endl;

  TH2D *EnerResoLUTs =new TH2D("EnerResoLUTs","Energy Resolution using LUTs;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
 TH2D* emu_EnerLUTs = new TH2D("emu_EnerLUTs", "emu_EnerLUTs;E_{MC, muon} [MeV];E_{reco, muon} [MeV]", 250, 0., 5000., 250, 0., 5000.);;

   TFile *input(0);
   for (int i = 1040; i < 1100; i++) {
      char *fname = Form("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_%i_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root",i);
  
  if (!gSystem->AccessPathName( fname )){ 
    cout<<"Open file "<<i<<endl;
    input = TFile::Open( fname ); 
  }else{ 
  cout<<"file"<<i<<" does not exist. Next one is: "<<i+1<<endl;
  i=i+1;
  char *fname = Form("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_%i_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root",i);
  input = TFile::Open( fname );
 }

  TTree *regTree = (TTree*)input->Get("nu_eneNEW");

  Double_t trueKE, neutrinoE, recoE_lookup;
  Float_t total_PMTs_hits2, total_hits2, total_ring_PEs2, pot_length2, hits_pot_length2, recoDWallR2, recoDWallZ2, lambda_max_2, myweight;
  regTree->SetBranchAddress("trueKE", &trueKE);
  regTree->SetBranchAddress("neutrinoE", &neutrinoE);
  regTree->SetBranchAddress("total_hits2", &total_hits2);
  regTree->SetBranchAddress("total_ring_PEs2", &total_ring_PEs2);
  regTree->SetBranchAddress("pot_length2", &pot_length2);
  regTree->SetBranchAddress("recoDWallR2", &recoDWallR2);
  regTree->SetBranchAddress("recoDWallZ2", &recoDWallZ2);
  regTree->SetBranchAddress("lambda_max_2", &lambda_max_2);
   
  for (Long64_t ievt=0; ievt<regTree->GetEntries();ievt++) {
     if (ievt%1000 == 0) {
        std::cout << "--- ... Processing event: " << ievt << std::endl;
     }
     regTree->GetEntry(ievt);

     TF1* f = new TF1("a","pol3");
     f->SetParameters( 1.10176e+02,  3.09111e+03 , -2.13215e+03, 1.08625e+04);
     double elut = f->Eval(total_ring_PEs2,0,0);
     //cout<<"total_ring_PEs2= "<<total_ring_PEs2<<" elut= "<<elut<<" / trueKE= "<<trueKE<<endl;
     EnerResoLUTs->Fill(trueKE, (100*(trueKE-elut)/trueKE));
     emu_EnerLUTs->Fill(trueKE, elut);
  }
 }//end of files
 regTree->Delete();
 input->Close();
 
 new TCanvas();
 TProfile *EnerResoLUTs_pfx= EnerResoLUTs->ProfileX("EnerResoLUTs_pfx", 0., 3000., "s");
 EnerResoLUTs_pfx->SetLineColor(1);
 EnerResoLUTs_pfx->Draw(); 

 TFile* outputFile = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoLUT.root", "RECREATE" );

 EnerResoLUTs->Write();
 EnerResoLUTs_pfx->Write();
 emu_EnerLUTs->Write();

 //regTree->Delete();
 //input->Close();
 outputFile->Close();
}
