#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TF1.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

void usingReco_scikit()
{
  cout<<"starting.... "<<endl;


  TH2D *EnerResoLUTs =new TH2D("EnerResoLUTs","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  //TH2D *EnerResoBDTG =new TH2D("EnerResoBDTG","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  TH2D *EnerResoBDTG =new TH2D("EnerResoBDTG","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",60, 0., 3000., 200, -100., 100.);
  TH2D *EnerResoBDTG2 =new TH2D("EnerResoBDTG2","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  TH2D *emu_mlp = new TH2D("emu_mlp", "emu_MLP;E_{MC, muon} [MeV];E_{reco, muon} [MeV]", 250, 0., 5000., 250, 0., 5000.);
  TH2D *plot_emuEneuMC= new TH2D("plot_emuEneuMC", "plot_emuEneuMC;E_{MC, muon} [MeV];E_{MC, neutrino} [MeV]", 250, 0., 5000., 250, 0., 5000.);
 
  TH2D *EnerResoBDTG_par0=new TH2D("EnerResoBDTG_par0","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  TH2D *EnerResoBDTG_par1=new TH2D("EnerResoBDTG_par1","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  TH2D *EnerResoBDTG_par2=new TH2D("EnerResoBDTG_par2","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.); 
  TH2D *EnerResoBDTG_par3=new TH2D("EnerResoBDTG_par3","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  TH2D *EnerResoBDTG_par4=new TH2D("EnerResoBDTG_par4","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  TH2D *EnerResoBDTG_par5=new TH2D("EnerResoBDTG_par5","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);

  TH2D *EnerResoBDTG_par01=new TH2D("EnerResoBDTG_par01","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  TH2D *EnerResoBDTG_par02=new TH2D("EnerResoBDTG_par02","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
 
  TH2D *EnerResoBDTG_Eneu=new TH2D("EnerResoBDTG_Eneu","Energy Resolution using Scikit BDTG;E_{MC, neutrino} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  TH2D *EnerResoBDTG_Eneu_withcuts=new TH2D("EnerResoBDTG_Eneu_withcuts","Energy Resolution using Scikit BDTG;E_{MC, neutrino} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);

  //TH2D *EnerResoBDTG_withcuts=new TH2D("EnerResoBDTG_withcuts","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",100, 0., 5000., 200, -100., 100.);
  TH2D *EnerResoBDTG_withcuts=new TH2D("EnerResoBDTG_withcuts","Energy Resolution using Scikit BDTG;E_{MC, muon} [MeV]; #Delta E/E [%]",60, 0., 3000., 200, -100., 100.);

  TH1D *plot_emu_erecoBDTG= new TH1D("plot_emu_erecoBDTG", "Energy Resolution;#Delta E/E [%]", 200, -100., 100.);
  TH1D *plot_emu_erecoBDTG_Infid= new TH1D("plot_emu_erecoBDTG_Infid", "Energy Resolution;#Delta E/E [%]", 200, -100., 100.);
  TH1D *plot_emu_erecoBDTG_Infid_he= new TH1D("plot_emu_erecoBDTG_Infid_he", "Energy Resolution;#Delta E/E [%]", 200, -100., 100.);
  TH1D *plot_emu_erecoBDTG_he= new TH1D("plot_emu_erecoBDTG_he", "Energy Resolution;#Delta E/E [%]", 200, -100., 100.);

  TH1D *plot_emu_erecoBDTG2= new TH1D("plot_emu_erecoBDTG2", "plot_emu_erecoBDTG;#Delta E/E", 200, -1., 1.);
  TH1D *plot_emu_erecoBDTG_par0= new TH1D("plot_emu_erecoBDTG_par0", "plot_emu_erecoBDTG_par0;#Delta E/E", 200, -1., 1.);
  TH1D *plot_emu_erecoBDTG_par1= new TH1D("plot_emu_erecoBDTG_par1", "plot_emu_erecoBDTG_par1;#Delta E/E", 200, -1., 1.);
  TH1D *plot_emu_erecoBDTG_par2= new TH1D("plot_emu_erecoBDTG_par2", "plot_emu_erecoBDTG_par0;#Delta E/E", 200, -1., 1.);
  TH1D *plot_emu_erecoBDTG_par3= new TH1D("plot_emu_erecoBDTG_par3", "plot_emu_erecoBDTG_par0;#Delta E/E", 200, -1., 1.);
  TH1D *plot_emu_erecoBDTG_par4= new TH1D("plot_emu_erecoBDTG_par4", "plot_emu_erecoBDTG_par0;#Delta E/E", 200, -1., 1.);
  TH1D *plot_emu_erecoBDTG_par5= new TH1D("plot_emu_erecoBDTG_par5", "plot_emu_erecoBDTG_par0;#Delta E/E", 200, -1., 1.);

  TH1D *plot_emu_erecoBDTG_par01= new TH1D("plot_emu_erecoBDTG_par01", "plot_emu_erecoBDTG_par0;#Delta E/E", 200, -1., 1.);
  TH1D *plot_emu_erecoBDTG_par02= new TH1D("plot_emu_erecoBDTG_par02", "plot_emu_erecoBDTG_par0;#Delta E/E", 200, -1., 1.);

  TH2D * plot_y= new TH2D("plot_y", "plot_y;E_{MC, neutrino} [MeV]; Bjorken y", 250, 0., 5000., 100, 0., 1.);


   TFile *input(0);  
    //for (int i = 1040; i < 1100; i++) {
     //char *fname = Form("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/scikit_reco-2.root");
     //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/scikit_reco.root");  
    //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/scikit_recoLAST.root"); //-paper
    //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/scikit_recoLASTnue.root"); //-paper
    //char *fname = Form("/home/edrakopo/work/ener_reco_scikit/scikit_recoLASTtest.root"); //100 estimators
    
    char *fname = Form("/home/edrakopo/WChscikit_reco.root");

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
  //regTree->SetBranchAddress("recoKE", &recoKE);
  regTree->SetBranchAddress("recoKE_BDTG", &recoKE_BDTG);
  //regTree->SetBranchAddress("recoKE_BDTG01", &recoKE_BDTG01);
  //regTree->SetBranchAddress("recoKE_BDTG_par02", &recoKE_BDTG_par02);
  regTree->SetBranchAddress("recoDwall", &recoDwall);
  regTree->SetBranchAddress("recoToWall", &recoToWall);
  /*regTree->SetBranchAddress("recoKE_BDTG2", &recoKE_BDTG2);
  regTree->SetBranchAddress("recoKE_BDTG_par0", &recoKE_BDTG_par0);
  regTree->SetBranchAddress("recoKE_BDTG_par1", &recoKE_BDTG_par1);
  regTree->SetBranchAddress("recoKE_BDTG_par2", &recoKE_BDTG_par2);
  regTree->SetBranchAddress("recoKE_BDTG_par3", &recoKE_BDTG_par3);
  regTree->SetBranchAddress("recoKE_BDTG_par4", &recoKE_BDTG_par4); 
  regTree->SetBranchAddress("recoKE_BDTG_par5", &recoKE_BDTG_par5);
  */

  for (Long64_t ievt=0; ievt<regTree->GetEntries();ievt++) {
  //for (Long64_t ievt=0; ievt<100;ievt++) {
     if (ievt%1000 == 0) {
        std::cout << "--- ... Processing event: " << ievt << std::endl;
     }
     regTree->GetEntry(ievt);

     //cout<<"neutrinoE: "<<neutrinoE<<" trueKE= "<<trueKE<<" recoKE_BDTG: "<<recoKE_BDTG<<endl;
     //EnerResoLUTs->Fill(trueKE, (100*(trueKE-recoKE)/trueKE));
     EnerResoBDTG->Fill(trueKE, (100*(trueKE-recoKE_BDTG)/trueKE));
     EnerResoBDTG_Eneu->Fill(neutrinoE, (100*(neutrinoE-recoKE_BDTG)/neutrinoE));
     //EnerResoBDTG_par01->Fill(trueKE, (100*(trueKE-recoKE_BDTG01)/trueKE));
     //EnerResoBDTG_par02->Fill(trueKE, (100*(trueKE-recoKE_BDTG_par02)/trueKE));
 
    //cout<<"recoDwall*550.: "<<recoDwall*550.<<" recoToWall*2500.: "<<recoToWall*2500.<<endl; 
    if( (recoDwall*550.)>100. && (recoToWall*2500.)>100.){
       EnerResoBDTG_withcuts->Fill(trueKE, (100*(trueKE-recoKE_BDTG)/trueKE));
       EnerResoBDTG_Eneu_withcuts->Fill(neutrinoE, (100*(neutrinoE-recoKE_BDTG)/neutrinoE));
      if(trueKE>200. && trueKE<600.){
        plot_emu_erecoBDTG_Infid->Fill(100.*(trueKE-recoKE_BDTG)/trueKE);
      }
      if(trueKE>=600. && trueKE<1400.){
        plot_emu_erecoBDTG_Infid_he->Fill(100.*(trueKE-recoKE_BDTG)/trueKE);
      } 
    }
  /*   EnerResoBDTG2->Fill(trueKE, (100*(trueKE-recoKE_BDTG2)/trueKE));
     EnerResoBDTG_par0->Fill(trueKE, (100*(trueKE-recoKE_BDTG_par0)/trueKE));
     EnerResoBDTG_par1->Fill(trueKE, (100*(trueKE-recoKE_BDTG_par1)/trueKE));
     EnerResoBDTG_par2->Fill(trueKE, (100*(trueKE-recoKE_BDTG_par2)/trueKE));
     EnerResoBDTG_par3->Fill(trueKE, (100*(trueKE-recoKE_BDTG_par3)/trueKE));
     EnerResoBDTG_par4->Fill(trueKE, (100*(trueKE-recoKE_BDTG_par4)/trueKE));
     EnerResoBDTG_par5->Fill(trueKE, (100*(trueKE-recoKE_BDTG_par5)/trueKE));
  */
     emu_mlp->Fill(trueKE,recoKE);
     plot_emuEneuMC->Fill(trueKE,neutrinoE);
     plot_y->Fill( neutrinoE, (1.-trueKE/neutrinoE) );

     if(trueKE>200. && trueKE<600.){ //cout<<"in"<<" DE/E: "<<trueKE-recoKE_BDTG/trueKE<<endl;
       //if((100*(trueKE-recoKE_BDTG)/trueKE)<80. && (100*(trueKE-recoKE_BDTG)/trueKE)>-80.){
         plot_emu_erecoBDTG->Fill(100.*(trueKE-recoKE_BDTG)/trueKE);
       //}
       //plot_emu_erecoBDTG_par01->Fill((trueKE-recoKE_BDTG01)/trueKE);
       //plot_emu_erecoBDTG_par02->Fill((trueKE-recoKE_BDTG_par02)/trueKE);
      /* plot_emu_erecoBDTG2->Fill((trueKE-recoKE_BDTG2)/trueKE);
       plot_emu_erecoBDTG_par0->Fill((trueKE-recoKE_BDTG_par0)/trueKE);
       plot_emu_erecoBDTG_par1->Fill((trueKE-recoKE_BDTG_par1)/trueKE);
       plot_emu_erecoBDTG_par2->Fill((trueKE-recoKE_BDTG_par2)/trueKE);
       plot_emu_erecoBDTG_par3->Fill((trueKE-recoKE_BDTG_par3)/trueKE);
       plot_emu_erecoBDTG_par4->Fill((trueKE-recoKE_BDTG_par4)/trueKE);
       plot_emu_erecoBDTG_par5->Fill((trueKE-recoKE_BDTG_par5)/trueKE);  */
     }
     if(trueKE>=600. && trueKE<1400.){
        plot_emu_erecoBDTG_he->Fill(100.*(trueKE-recoKE_BDTG)/trueKE);
      }

  }
 //}//end of files
 //regTree->Delete();
 //input->Close();
 
 new TCanvas();
 plot_y->Draw("ColZ");

 TProfile *EnerResoLUTs_pfx= EnerResoLUTs->ProfileX("EnerResoLUTs_pfx", 0., 3000., "s");
/* EnerResoLUTs_pfx->GetXaxis()->SetRangeUser(0,3000);
 EnerResoLUTs_pfx->GetXaxis()->SetTitle("E_{MC, muon} [MeV]");
 EnerResoLUTs_pfx->GetYaxis()->SetTitle("#Delta E/E [%]");
 EnerResoLUTs_pfx->SetLineColor(1);
 EnerResoLUTs_pfx->SetFillColor(29); //17,2,4,29
 EnerResoLUTs_pfx->Draw("E3L");
 //EnerResoLUTs_pfx->Draw(); 
*/
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

 //
 TProfile *EnerResoBDTG_par02_pfx= EnerResoBDTG_par02->ProfileX("EnerResoBDTG_par02_pfx", 0., 3000., "s");
 /*//EnerResoBDTG2_pfx->GetXaxis()->SetRangeUser(0,3000);
 EnerResoBDTG_par02_pfx->GetXaxis()->SetTitle("E_{MC, muon} [MeV]");
 EnerResoBDTG_par02_pfx->GetYaxis()->SetTitle("#Delta E/E [%]");
 EnerResoBDTG_par02_pfx->SetLineColor(1);
 EnerResoBDTG_par02_pfx->SetFillColor(2); //17,2,4,29
 EnerResoBDTG_par02_pfx->Draw("E3L");
 //EnerResoBDTG2_pfx->SetLineStyle(2);
 //EnerResoBDTG2_pfx->Draw("same");
*/
 TProfile *EnerResoBDTG_par01_pfx= EnerResoBDTG_par01->ProfileX("EnerResoBDTG_par01_pfx", 0., 3000., "s");
 TProfile *EnerResoBDTG_par02_pfx= EnerResoBDTG_par02->ProfileX("EnerResoBDTG_par02_pfx", 0., 3000., "s");
 /*TProfile *EnerResoBDTG_par0_pfx= EnerResoBDTG_par0->ProfileX("EnerResoBDTG_par0_pfx", 0., 3000., "s");
 TProfile *EnerResoBDTG_par1_pfx= EnerResoBDTG_par1->ProfileX("EnerResoBDTG_par1_pfx", 0., 3000., "s");
 TProfile *EnerResoBDTG_par2_pfx= EnerResoBDTG_par2->ProfileX("EnerResoBDTG_par2_pfx", 0., 3000., "s");
 TProfile *EnerResoBDTG_par3_pfx= EnerResoBDTG_par3->ProfileX("EnerResoBDTG_par3_pfx", 0., 3000., "s");
 TProfile *EnerResoBDTG_par4_pfx= EnerResoBDTG_par4->ProfileX("EnerResoBDTG_par4_pfx", 0., 3000., "s");
 TProfile *EnerResoBDTG_par5_pfx= EnerResoBDTG_par5->ProfileX("EnerResoBDTG_par5_pfx", 0., 3000., "s");
*/


  TFile* outputFile = new TFile("/home/edrakopo/WChEnergyRecoScikit.root", "RECREATE" );
 //TFile* outputFile = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoScikitNEW2.root", "RECREATE" );
 //TFile* outputFile = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoScikitNEW2_200_600MeV.root", "RECREATE" ); //-paper
 //TFile* outputFile = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoScikitNEW2_200_600MeVnue.root", "RECREATE" ); //-paper
 //TFile* outputFile = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/EnergyRecoScikitNEW2_200_600MeVtest.root", "RECREATE" );

 plot_y->Write();
 EnerResoLUTs->Write();
 EnerResoLUTs_pfx->Write();
 EnerResoBDTG->Write();
 EnerResoBDTG2->Write();
 emu_mlp->Write();
 plot_emuEneuMC->Write();
 plot_emu_erecoBDTG_Infid->Write();
 plot_emu_erecoBDTG_Infid_he->Write();
 plot_emu_erecoBDTG_he->Write();
 EnerResoBDTG_withcuts->Write();

 EnerResoBDTG_Eneu_pfx->Write();
 EnerResoBDTG_Eneu_withcuts_pfx->Write();
 EnerResoBDTG_pfx->Write();
 EnerResoBDTG_par01_pfx->Write();
 EnerResoBDTG_par02_pfx->Write();
 EnerResoBDTG_pfx_withcuts->Write();

/* EnerResoBDTG2_pfx->Write();
 EnerResoBDTG_par0_pfx->Write();
 EnerResoBDTG_par1_pfx->Write();
 EnerResoBDTG_par2_pfx->Write();
 EnerResoBDTG_par3_pfx->Write();
 EnerResoBDTG_par4_pfx->Write();
 EnerResoBDTG_par5_pfx->Write();
*/
 plot_emu_erecoBDTG->Write();
 plot_emu_erecoBDTG_par01->Write();
 plot_emu_erecoBDTG_par02->Write();
/* plot_emu_erecoBDTG2->Write();
 plot_emu_erecoBDTG_par0->Write();
 plot_emu_erecoBDTG_par1->Write();
 plot_emu_erecoBDTG_par2->Write();
 plot_emu_erecoBDTG_par3->Write();
 plot_emu_erecoBDTG_par4->Write();
 plot_emu_erecoBDTG_par5->Write();
*/
 //regTree->Delete();
 //input->Close();
 outputFile->Close();
}
