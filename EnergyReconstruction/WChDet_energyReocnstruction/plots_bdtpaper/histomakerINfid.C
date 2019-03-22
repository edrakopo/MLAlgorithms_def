#include "TFile.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPolyLine.h"

void histomakerINfid(){

  gROOT -> ProcessLine( ".x ./mattStyle.C" );

  //TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/EnergyRecoKeras_TrainTestSampleG.root");
  TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/WChEnergyRecoScikit600.root");//--TEST
  TFile* file1 = new TFile("/Users/edrakopo/work/work_energy_reco/WChEnergyRecoScikit.root");//--TEST
  //TFile* file3 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/TMVAMLP_INFID_enerRsol.root");
  TFile* file3 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/TMVAMLP_INFID_enerRsolTEST.root");
  //TFile* file1 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/EnergyRecoScikitNEW2_200_600MeV.root");
  TFile* file = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/TMVABDTG_INFID_enerRsol.root");

  new TCanvas();
  TH1D *plot_emu_erecoBDTG =(TH1D*)file1->Get("plot_emu_erecoBDTG_Infid");
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file1->Get("plot_emu_erecoBDTG_Infid_he");  //600-1400
 
  TH1D *plot_emu_erecoMLP =(TH1D*)file3->Get("ener_resolution");
  TH1D *plot_emu_erecoTMVABDTG =(TH1D*)file->Get("ener_resolution");
  TH1D *plot_emu_erecoKERAS =(TH1D*)file2->Get("plot_emu_erecoBDTG_Infid");
  TH1D *ener_resolutionTITUS=(TH1D*)file->Get("ener_resolutionTITUS");
  /*
  TH1D *plot_emu_erecoMLP =(TH1D*)file3->Get("plot_emu_erecoBDTG_Infid_he");// MLP //600-1400
  TH1D *plot_emu_erecoTMVABDTG =(TH1D*)file->Get("plot_emu_erecoBDTG_Infid_he");//TMVABDTG
  TH1D *plot_emu_erecoKERAS =(TH1D*)file2->Get("plot_emu_erecoBDTG_Infid_he"); //keras 600-1400
  TH1D *ener_resolutionTITUS=(TH1D*)file->Get("plot_emu_erecoBDTG_Infid_heTITUS"); //TITUS
  */

  //profBDT->GetXaxis()->SetRangeUser(0,3000);
  plot_emu_erecoBDTG->GetYaxis()->SetTitle("Number of Entries");
  plot_emu_erecoBDTG->GetXaxis()->SetTitle("#Delta E/E [%]");
  plot_emu_erecoBDTG->SetLineColor(1);
  plot_emu_erecoBDTG->SetLineWidth(1);
  plot_emu_erecoBDTG->SetFillColor(18);
  plot_emu_erecoBDTG->SetMarkerStyle(1);
  plot_emu_erecoBDTG->Draw();

  plot_emu_erecoKERAS->SetMarkerStyle(1);
  plot_emu_erecoKERAS->SetFillColorAlpha(kCyan-6, 0.4); //0.6);
  plot_emu_erecoKERAS->Draw("same");
/*
  plot_emu_erecoTMVABDTG->SetFillColorAlpha(kRed-4, 0.4);
  plot_emu_erecoTMVABDTG->SetMarkerColor(kBlue);
  plot_emu_erecoTMVABDTG->SetMarkerSize(0.5);
  plot_emu_erecoTMVABDTG->SetMarkerStyle(1);
  plot_emu_erecoTMVABDTG->Draw("same");
*/
  /*plot_emu_erecoMLP->SetMarkerStyle(1);
  plot_emu_erecoMLP->SetFillColorAlpha(kCyan-6, 0.8);
  plot_emu_erecoMLP->Draw("same");
*/
 
  TLegend* legend = new TLegend(0.1,0.7,0.48,0.6);
  legend->SetBorderSize(1);
  legend->SetTextFont(22);
  legend->SetTextSize(0.035);
  legend->AddEntry(plot_emu_erecoBDTG,"Scikit Boosted Decision Tree");
  legend->AddEntry(plot_emu_erecoTMVABDTG,"TMVA Boosted Decision Tree");
  legend->AddEntry(plot_emu_erecoMLP,"TMVA Neural Network");
  legend->Draw("same");
 
/*
  ener_resolutionTITUS->SetMarkerStyle(1);
  ener_resolutionTITUS->SetFillColorAlpha(kAzure+2, 0.8);
  ener_resolutionTITUS->Draw("same");

  plot_emu_erecoKERAS->SetMarkerStyle(1);
  plot_emu_erecoKERAS->SetFillColorAlpha(kOrange+8, 0.6);
  plot_emu_erecoKERAS->Draw("same");

  TLegend* legend2 = new TLegend(0.1,0.7,0.48,0.6);
  legend2->SetBorderSize(1);
  legend2->SetTextFont(22);
  legend2->SetTextSize(0.035);
  legend2->AddEntry(plot_emu_erecoBDTG,"Scikit Boosted Decision Tree");
  legend2->AddEntry(ener_resolutionTITUS,"TITUS lookup tables");
  legend2->AddEntry(plot_emu_erecoKERAS,"Tensorflow Neural Network");
  legend2->Draw();
 
*/
}
