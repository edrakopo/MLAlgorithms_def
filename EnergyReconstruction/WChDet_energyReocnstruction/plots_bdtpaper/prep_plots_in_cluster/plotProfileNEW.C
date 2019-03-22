
#include "TFile.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPolyLine.h"

typedef struct{

  TPolyLine line1;
  TPolyLine line2;
  
}Envelope;

Envelope* makeEnvelope(TProfile* prof ){

  Envelope* env = new Envelope();
  int nbin = prof->GetNbinsX();
  for (auto i = 1; i <= nbin; ++i){
    env->line1.SetPoint(i, prof->GetBinCenter(i), prof->GetBinContent(i));
  }
  
  return env;
}

void plotProfileNEW(){

  gROOT -> ProcessLine( ".x ./mattStyleB.C" );
 
  //TFile* file = new TFile("TMVARegApp_BDTG_ener_lilia_CCQE_1040_1100_NEWlookupsTEST_MYselectedVars.root");
  //TFile* file = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoScikitNEW.root");
  //TFile* file = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoScikitNEW2.root");
  //TFile* file = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoLUT.root");
 //TFile* file = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/outputs/FINALTMVARegApp_BDTG_ener_lilia_CCQE_1040_1100_NEWlookupsTEST_MYselectedVars_large_statTEST10000treesB.root");
  //TFile* file = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoKeras_200_600MeV.root");
  //TFile* file = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoKeras_200_600MeV_2_100_1_10.root");
  //TFile* file = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoKeras_200_600MeV_2_100_1_10_1_1.root");
  TFile* file = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoKeras_200_600MeV_3_100_1_10_1_1.root");

  //TProfile* profBDT = (TProfile*)file->Get("EnerResoLUTs_pfx");
  TProfile* profBDT = (TProfile*)file->Get("EnerResoBDTG_pfx");
  //TProfile* profBDT = (TProfile*)file->Get("EnerResoBDTG_pfx_withcuts");
  //TProfile* profBDT = (TProfile*)file->Get("EnerResoBDTG_Eneu_withcuts_pfx");
  //TProfile* profBDT = (TProfile*)file->Get("pfx_BDTG2");
  //TProfile* profBDT = (TProfile*)file->Get("pfx_TITUS2");
  //TProfile* profBDT = (TProfile*)file->Get("EnerResoBDTG_pfx_withcutsTitus");
  //profBDT->Rebin(5);
  profBDT->GetXaxis()->SetRangeUser(0,3000);
  profBDT->GetXaxis()->SetTitle("E_{MC, muon} [MeV]");
  profBDT->GetYaxis()->SetTitle("#Delta E/E [%]");
  profBDT->SetLineColor(1);
  //profBDT->SetFillColor(434); //17,2,4,29,31,36,38,46 
  //profBDT->SetFillColor(424);
  //profBDT->Draw("E3L");
  //profBDT->Draw("E0");
  //profBDT->Draw();

  gStyle->SetPalette(87);
  TH2D *EnerResoBDTG2 = (TH2D*)file->Get("EnerResoBDTG");
  EnerResoBDTG2->Draw("ColZ");
  //EnerResoBDTG2->Draw("E0 same");

 /* Envelope* env = makeEnvelope(profBDT);
  env->line1.Draw();*/
}
