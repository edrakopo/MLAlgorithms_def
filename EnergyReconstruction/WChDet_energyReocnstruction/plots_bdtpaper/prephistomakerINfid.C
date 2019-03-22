#include "TFile.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPolyLine.h"

void prephistomakerINfid(){

  gROOT -> ProcessLine( ".x ./mattStyle.C" );

  //TFile* file = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/FINALTMVARegApp_MLP_ener_lilia_CCQE_1040_1100_NEWlookupsTEST_MYselectedVars_large_statTEST10000treesB.root");
  TFile* file = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/FINALTMVARegApp_BDTG_ener_lilia_CCQE_1040_1100_NEWlookupsTEST_MYselectedVars_large_statTEST10000treesB2.root");
  
  TH1F * ener_resolution = new TH1F("ener_resolution","ener_resolution; (E_{#mu}-E_{reco})/E_{#mu}", 200, -100., 100.);
  TH1F * ener_resolutionTITUS= new TH1F("ener_resolutionTITUS","ener_resolutionTITUS; (E_{#mu}-E_{reco})/E_{#mu}", 200, -100., 100.);
  
  TTree *regTree = (TTree*)file->Get("outTree");

  Float_t trueKE1, neutrinoE, recoE_lookup1, recoKE, recoKE_BDTG, recoKE_BDTG2, recoDWall, recoToWall;

  regTree->SetBranchAddress("trueKE1", &trueKE1);
  regTree->SetBranchAddress("recoKE", &recoKE);
  regTree->SetBranchAddress("recoE_lookup1", &recoE_lookup1);
  regTree->SetBranchAddress("recoDWall", &recoDWall);
  regTree->SetBranchAddress("recoToWall", &recoToWall);

  for (Long64_t ievt=0; ievt<regTree->GetEntries();ievt++) {
  //for (Long64_t ievt=0; ievt<100;ievt++) {
     if (ievt%1000 == 0) {
        std::cout << "--- ... Processing event: " << ievt <<" recoKE: "<<recoKE<< std::endl;
     }
     regTree->GetEntry(ievt);
     float trueKE=trueKE1;

    if( (recoDWall*550.)>100. && (recoToWall*2500.)>100.){
     if(trueKE>200. && trueKE<600.){ 
       ener_resolution -> Fill(100*(trueKE-recoKE)/trueKE);
       ener_resolutionTITUS -> Fill(100*(trueKE-recoE_lookup1)/trueKE);
    }}
   
  }

  TFile *target  = new TFile( "/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/TMVABDTG_INFID_enerRsol.root","RECREATE" );
  //TFile *target  = new TFile( "/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/TMVAMLP_INFID_enerRsol.root","RECREATE" );

  ener_resolution->Write();
  ener_resolutionTITUS->Write();
  
  target->Close();
}
