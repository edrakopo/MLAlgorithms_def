#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TF1.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

void usingReco_scikitKeras()
{
  cout<<"starting.... "<<endl;


  //TH2D *EnerResoBDTG =new TH2D("EnerResoBDTG","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  TH2D *EnerResoBDTG =new TH2D("EnerResoBDTG","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",60, 0., 3000., 200, -100., 100.);
  TH2D *EnerResoBDTG2 =new TH2D("EnerResoBDTG2","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  TH2D *emu_mlp = new TH2D("emu_mlp", "emu_MLP;E_{MC, muon} [MeV];E_{reco, muon} [MeV]", 250, 0., 5000., 250, 0., 5000.);
  TH2D *plot_emuEneuMC= new TH2D("plot_emuEneuMC", "plot_emuEneuMC;E_{MC, muon} [MeV];E_{MC, neutrino} [MeV]", 250, 0., 5000., 250, 0., 5000.);
  
  TH2D *EnerResoBDTG_Eneu=new TH2D("EnerResoBDTG_Eneu","Energy Resolution using Scikit BDTG;E_{MC, neutrino} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  TH2D *EnerResoBDTG_Eneu_withcuts=new TH2D("EnerResoBDTG_Eneu_withcuts","Energy Resolution using Scikit BDTG;E_{MC, neutrino} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);

  //TH2D *EnerResoBDTG_withcuts=new TH2D("EnerResoBDTG_withcuts","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  TH2D *EnerResoBDTG_withcuts=new TH2D("EnerResoBDTG_withcuts","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",60, 0., 3000., 200, -100., 100.);

  TH1D *plot_emu_erecoBDTG= new TH1D("plot_emu_erecoKERAS", "Energy Resolution;#Delta E/E [%]", 200, -100., 100.);
  TH1D *plot_emu_erecoBDTG_Infid= new TH1D("plot_emu_erecoBDTG_Infid", "Energy Resolution;#Delta E/E [%]", 200, -100., 100.);
  TH1D *plot_emu_erecoBDTG_Infid_he= new TH1D("plot_emu_erecoBDTG_Infid_he", "Energy Resolution;#Delta E/E [%]", 200, -100., 100.);
  TH1D *plot_emu_erecoBDTG_he= new TH1D("plot_emu_erecoBDTG_he", "Energy Resolution;#Delta E/E [%]", 200, -100., 100.);
  TH1D *plot_emu_erecoBDTG2= new TH1D("plot_emu_erecoKERAS2", "plot_emu_erecoBDTG;#Delta E/E", 200, -1., 1.);
  TH2D * plot_y= new TH2D("plot_y", "plot_y;E_{MC, neutrino} [MeV]; Bjorken y", 250, 0., 5000., 100, 0., 1.);


   TFile *input(0);  
    //for (int i = 1040; i < 1100; i++) {
   //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco.root");
   //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco.root");
   //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco_3_10.root");
   //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco_1_50_2_10.root");
   // char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco_2_50_1_10.root");  
   //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco_2_100_2_10.root"); 
   //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco_1_100_2_10.root");  
   //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco_2_100_1_10_1_1.root"); 
   //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco_3_100_1_10_1_1.root"); //GOOD   
   //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco_3_100_2_10_1_1.root");//overtrained
   //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco_4_100_1_10_1_1.root");//overtained
   //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco_3_100_3_10.root");
   //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco_TrainTestSample.root");
   //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco_TrainTestSample100ep.root");
   char *fname = Form("/home/edrakopo/work/ener_reco_scikit/keras_reco_TrainTestSamplePAPER.root"); //keras_reco_TrainTestSampleG.root");

  if (!gSystem->AccessPathName( fname )){ 
    cout<<"Open file "<<endl;
    input = TFile::Open( fname ); 
  }/*else{ 
  cout<<"file"<<i<<" does not exist. Next one is: "<<i+1<<endl;
  i=i+1;
  char *fname = Form("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/editing_ene/outputs/nu_numu_%i_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root",i);
  input = TFile::Open( fname );
 }*/

  TTree *regTree = (TTree*)input->Get("tuple");

  Float_t trueKE, neutrinoE, recoE_lookup, recoKE, recoKE_BDTG, recoKE_BDTG2, recoKE_BDTG_par0, recoKE_BDTG_par1, recoKE_BDTG_par2, recoKE_BDTG_par3, recoKE_BDTG_par4, recoKE_BDTG_par5, recoKE_BDTG01, recoKE_BDTG_par02, recoDwall, recoToWall;
  Float_t total_PMTs_hits2, total_hits2, total_ring_PEs2, pot_length2, hits_pot_length2, recoDWallR2, recoDWallZ2, lambda_max_2, myweight;
  regTree->SetBranchAddress("trueKE", &trueKE);
  regTree->SetBranchAddress("neutrinoE", &neutrinoE);
  regTree->SetBranchAddress("recoKE", &recoKE);
  regTree->SetBranchAddress("recoDwall", &recoDwall);
  regTree->SetBranchAddress("recoToWall", &recoToWall);

  for (Long64_t ievt=0; ievt<regTree->GetEntries();ievt++) {
  //for (Long64_t ievt=0; ievt<100;ievt++) {
     if (ievt%1000 == 0) {
        std::cout << "--- ... Processing event: " << ievt << std::endl;
     }
     regTree->GetEntry(ievt);

     //cout<<"neutrinoE: "<<neutrinoE<<" trueKE= "<<trueKE<<" recoKE_BDTG: "<<recoKE_BDTG<<endl;
     EnerResoBDTG->Fill(trueKE, (100*(trueKE-recoKE)/trueKE));
     EnerResoBDTG_Eneu->Fill(neutrinoE, (100*(neutrinoE-recoKE)/neutrinoE));
 
    //cout<<"recoDwall*550.: "<<recoDwall*550.<<" recoToWall*2500.: "<<recoToWall*2500.<<endl; 
    if( (recoDwall*550.)>100. && (recoToWall*2500.)>100.){
       EnerResoBDTG_withcuts->Fill(trueKE, (100*(trueKE-recoKE)/trueKE));
       EnerResoBDTG_Eneu_withcuts->Fill(neutrinoE, (100*(neutrinoE-recoKE)/neutrinoE));
       if(trueKE>200. && trueKE<600.){
         plot_emu_erecoBDTG_Infid->Fill(100.*(trueKE-recoKE)/trueKE);
      }
      if(trueKE>=600. && trueKE<1400.){
        plot_emu_erecoBDTG_Infid_he->Fill(100.*(trueKE-recoKE)/trueKE);
      } 
    }
     
     emu_mlp->Fill(trueKE,recoKE);
     plot_emuEneuMC->Fill(trueKE,neutrinoE);
     plot_y->Fill( neutrinoE, (1.-trueKE/neutrinoE) );

     if(trueKE>200. && trueKE<600.){ //cout<<"in"<<" DE/E: "<<trueKE-recoKE_BDTG/trueKE<<endl;
     //if(trueKE>500. && trueKE<550.){
       //if((100*(trueKE-recoKE_BDTG)/trueKE)<80. && (100*(trueKE-recoKE_BDTG)/trueKE)>-80.){
         plot_emu_erecoBDTG->Fill(100.*(trueKE-recoKE)/trueKE);
       //}
     }
     if(trueKE>=600. && trueKE<1400.){
        plot_emu_erecoBDTG_he->Fill(100.*(trueKE-recoKE)/trueKE);
      }
     
  }
 //}//end of files
 //regTree->Delete();
 //input->Close();
 
 new TCanvas();
 //plot_y->Draw("ColZ");
 EnerResoBDTG->Draw("ColZ");

 TProfile * EnerResoBDTG_Eneu_pfx= EnerResoBDTG_Eneu->ProfileX("EnerResoBDTG_Eneu_pfx", 0., 3000., "s");
 TProfile * EnerResoBDTG_Eneu_withcuts_pfx= EnerResoBDTG_Eneu_withcuts->ProfileX("EnerResoBDTG_Eneu_withcuts_pfx", 0., 3000., "s");
 //new TCanvas();
 TProfile *EnerResoBDTG_pfx= EnerResoBDTG->ProfileX("EnerResoBDTG_pfx", 0., 3000., "s");
 //EnerResoBDTG_pfx->GetXaxis()->SetRangeUser(0,3000);
/* EnerResoBDTG_pfx->GetXaxis()->SetTitle("E_{MC, muon} [MeV]");
 EnerResoBDTG_pfx->GetYaxis()->SetTitle("#Delta E/E [%]");
 EnerResoBDTG_pfx->SetLineColor(1);
 EnerResoBDTG_pfx->SetFillColor(4); //17,2,4,29
 */
 //EnerResoBDTG_pfx->Draw("E3L");
 //EnerResoBDTG_pfx->SetLineStyle(2);
 //EnerResoBDTG_pfx->Draw("same");
 
 new TCanvas();
 TProfile *EnerResoBDTG_pfx_withcuts= EnerResoBDTG_withcuts->ProfileX("EnerResoBDTG_pfx_withcuts", 0., 3000., "s");
 EnerResoBDTG_pfx->GetXaxis()->SetTitle("E_{MC, muon} [MeV]");
 EnerResoBDTG_pfx->GetYaxis()->SetTitle("#Delta E/E [%]");
 EnerResoBDTG_pfx->SetLineColor(1);
 EnerResoBDTG_pfx->Draw();
 EnerResoBDTG_pfx_withcuts->SetLineColor(2); 
 EnerResoBDTG_pfx_withcuts->SetLineStyle(2);
 EnerResoBDTG_pfx_withcuts->Draw("same");

 new TCanvas();
 plot_emu_erecoBDTG->Draw();

  //TFile* outputFile = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoKeras_200_600MeV_4vars.root", "RECREATE" );
  //TFile* outputFile = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoKeras_200_600MeV_4_100_1_10_1_1.root", "RECREATE" );
  //TFile* outputFile = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoKeras_3_100_1_10_1_1B.root", "RECREATE" );
 TFile* outputFile = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoKeras_TrainTestSampleG.root", "RECREATE" );
  //TFile* outputFile = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoKeras_TrainTestSampleG_100epochs.root", "RECREATE" );

 plot_y->Write();
 EnerResoBDTG->Write();
 emu_mlp->Write();
 plot_emuEneuMC->Write();

 EnerResoBDTG_Eneu_pfx->Write();
 EnerResoBDTG_Eneu_withcuts_pfx->Write();
 EnerResoBDTG_pfx->Write();
 EnerResoBDTG_pfx_withcuts->Write();
 EnerResoBDTG_withcuts->Write();

 plot_emu_erecoBDTG->Write();
 plot_emu_erecoBDTG_Infid->Write();
 plot_emu_erecoBDTG_Infid_he->Write();
 plot_emu_erecoBDTG_he->Write();
 //regTree->Delete();
 //input->Close();
 outputFile->Close();
}
