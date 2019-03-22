#include "TFile.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPolyLine.h"

void histomaker2(){

  gROOT -> ProcessLine( ".x ./mattStyle.C" );
  gStyle->SetPalette(kCMYK,0);

  TFile* file = new TFile("/Users/edrakopo/work/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root");
  
  TH2D *plot_emu_hits =(TH2D*)file->Get("hMuEn2_HITS");
  TH2D *plot_emu_pes =(TH2D*)file->Get("hMuEn2__ring_PEs");

  new TCanvas();
  plot_emu_hits->GetXaxis()->SetTitle("E_{MC, muon} [MeV]");
  plot_emu_hits->GetYaxis()->SetTitle("$N_{C}^{Hit}$");
  plot_emu_hits->GetYaxis()->SetRange(0,20000);
  //plot_emu_hits->SetMinimum(0);
  //plot_emu_hits->SetMaximum(20000);
  plot_emu_hits->Draw("colz1");

  //profBDT->SetMarkerStyle(1);

  new TCanvas();
  plot_emu_pes->GetXaxis()->SetTitle("E_{MC, muon} [MeV]");
  plot_emu_pes->GetYaxis()->SetTitle("$N_{Ring}^{Hit}$");
  plot_emu_pes->GetYaxis()->SetRange(0,20000);
  //plot_emu_pes->SetMinimum(0);
  //plot_emu_pes->SetMaximum(20000);
  plot_emu_pes->Draw("colz1");

/*  TLegend* legend2 = new TLegend(0.1,0.7,0.48,0.6);
  legend2->SetBorderSize(1);
  legend2->AddEntry("plot_emu_erecoBDTG","Scikit Boosted Decision Tree");
  //legend2->AddEntry("ener_resolution","ROOT-TMVA BDTG");
  legend2->AddEntry("ener_resolutionTITUS","TITUS lookup tables");
  legend2->AddEntry("plot_emu_erecoKERAS","Deep Learning");
  legend2->Draw();
 */

}
