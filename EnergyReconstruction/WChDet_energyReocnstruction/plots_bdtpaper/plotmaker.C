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
    //cout<<"i: "<<i<<" prof->GetXaxis()->GetBinCenter(i): "<<prof->GetXaxis()->GetBinCenter(i)<<" prof->GetBinContent(i): "<<prof->GetBinContent(i)<<endl;
    env->line1.SetPoint(i, prof->GetBinCenter(i), prof->GetBinContent(i));
  }
  
  return env;
}
/*Envelope* makeEnvelope2(TH2D* prof2 ){

  Envelope* env2 = new Envelope();
  int nbin = prof2->GetNbinsX();
  for (auto i = 1; i <= nbin; ++i){
 //   cout<<"i: "<<i<<" prof2->GetXaxis()->GetBinCenter(i): "<<prof2->GetXaxis()->GetBinCenter(i)<<" prof2->GetBinContent(i): "<<prof2->GetBinContent(i)<<endl; 
    //env2->line2.SetPoint(i, prof2->GetXaxis()->GetBinCenter(i), prof2->GetBinContent(i));
  }

  return env2;
}
*/

//_______________

void plotmaker(){

  gROOT -> ProcessLine( ".x ./mattStyle.C" );
  //if (isColor == false){
  //  gStyle->SetPalette(kGreyScale,0);
 /* }
  else {*/
  // gStyle->SetPalette(kLightTemperature,0);
  //} 
  //gStyle->SetPalette(kInvertedDarkBodyRadiator,0);
  gStyle->SetPalette(kCMYK,0);
  //gStyle->SetPalette(kSolar);
  //gStyle->SetPalette(kPearl,0);

  //TFile* file = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/EnergyRecoKeras_TrainTestSampleG.root");
  //TFile* file = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/EnergyRecoScikitNEW2_200_600MeV.root");
  TFile* file = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/FINALTMVARegApp_BDTG_ener_lilia_CCQE_1040_1100_NEWlookupsTEST_MYselectedVars_large_statTEST10000treesB.root");
  //TFile* file = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/FINALTMVARegApp_MLP_ener_lilia_CCQE_1040_1100_NEWlookupsTEST_MYselectedVars_large_statTEST10000treesB.root");

  //TFile* file = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/EnergyRecoScikitNEW2_200_600MeVnue.root");
  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorScikit_nue.root");
  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorScikit_nueINfid.root");
  //TFile* file = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/FINALTMVARegApp_BDTG_ener_lilia_CCQE_1040_1100_NEWlookupsTEST_MYselectedVars_nue.root");
  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorTITUSnueINfid.root"); 
  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorTITUSnue.root");

  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorTITUS.root");
  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorTMVAINfid.root"); 
  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorTMVA.root");
  TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorTITUSINfid.root");
  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorScikitINfid.root");
  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorScikit.root");
  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorKeras_TrainTestSampleGINfid.root");
  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorKeras_TrainTestSampleG.root");
  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorTMVAINfid_MLP.root");
  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/SigmaErrorTMVA_MLP.root");

  TGraphErrors* ge = (TGraphErrors*)file2->Get("Graph");

  new TCanvas();
  //TH2D *EnerResoBDTG =(TH2D*)file->Get("EnerResoBDTG");
  //TH2D *EnerResoBDTG =(TH2D*)file->Get("hDiffEnemlpmuTitus2");
  //TH2D *EnerResoBDTG =(TH2D*)file->Get("EnerResoBDTG_withcuts");
  TH2D *EnerResoBDTG =(TH2D*)file->Get("EnerResoBDTG_withcutsTitus");
  //EnerResoBDTG->Draw("ColZ");
  //EnerResoBDTG->DrawClone("col");
  EnerResoBDTG->SetMinimum(0);
  EnerResoBDTG->SetMaximum(1300);  //for numu
  //EnerResoBDTG->SetMaximum(70); //for nue 
  EnerResoBDTG->GetYaxis()->SetRangeUser(-100, 100);
  EnerResoBDTG->GetXaxis()->SetTitle("E_{MC, muon} [MeV]");
  //EnerResoBDTG->GetXaxis()->SetTitle("E_{MC, electron} [MeV]");
  EnerResoBDTG->GetYaxis()->SetTitle("#Delta E/E [%]");
  EnerResoBDTG->Draw("colz1");


  ge->SetLineWidth(5);
  ge->SetMarkerStyle(1);
  ge->Draw("L || same");
  //EnerResoBDTG->Draw("colz1 same ");


  //TProfile* profBDT = (TProfile*)file->Get("EnerResoLUTs_pfx");
  //TProfile* profBDT = (TProfile*)file->Get("EnerResoBDTG_pfx");
  TProfile* profBDT = (TProfile*)file->Get("EnerResoBDTG_pfx_withcuts");
  //TProfile* profBDT = (TProfile*)file->Get("pfx_BDTG2");
  //TProfile* profBDT = (TProfile*)file->Get("pfx_TITUS2");
  //TProfile* profBDT = (TProfile*)file->Get("EnerResoBDTG_pfx_withcutsTitus");
  ////profBDT->Rebin(5);
  ////profBDT->GetXaxis()->SetRangeUser(0,3000);
  //profBDT->GetXaxis()->SetTitle("E_{MC, muon} [MeV]");
  //profBDT->GetYaxis()->SetTitle("#Delta E/E [%]");
  profBDT->SetMarkerStyle(1);
  //profBDT->SetLineColor(1);
  profBDT->SetLineWidth(5);
  ////profBDT->SetFillColor(1); //17,2,4,29,31,36,38,46 
  ////profBDT->SetFillColor(424);
  //profBDT->SetFillStyle(3354);
  
  profBDT->SetFillColorAlpha(1, 0.); //add transparency
  //profBDT->Draw("APL || same");
  profBDT->Draw("E3L same");

  Envelope* env = makeEnvelope(profBDT);
  env->line1.Draw();


}
