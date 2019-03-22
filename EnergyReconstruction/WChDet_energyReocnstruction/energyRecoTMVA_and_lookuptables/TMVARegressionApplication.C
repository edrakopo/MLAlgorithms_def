/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVARegressionApplication                                          *
 *                                                                                *
 * This macro provides a simple example on how to use the trained regression MVAs *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif

using namespace TMVA;

void TMVARegressionApplication( TString myMethodList = "" ) 
{
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDEFoam"]         = 0; 
   Use["KNN"]             = 0;
   // 
   // --- Linear Discriminant Analysis
   Use["LD"]		  = 0;
   // 
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   // 
   // --- Neural Network
   Use["MLP"]             = 0; 
   // 
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 0;
   Use["BDTG"]            = 1;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVARegressionApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   /*Float_t var1, var2;
   reader->AddVariable( "var1", &var1 );
   reader->AddVariable( "var2", &var2 );*/
   Float_t pot_length2, hits_pot_length2, recoDWallR2, recoDWallZ2, lambda_max_2, total_hits2;
   Float_t recoDWall_2, recoToWall_2, vtxTrackBias_2;
   Float_t total_ring_PEs2=0.;
   reader->AddVariable( "total_hits2", &total_hits2);
   reader->AddVariable( "total_ring_PEs2", &total_ring_PEs2);
   reader->AddVariable( "recoDWallR2", &recoDWallR2);
   reader->AddVariable( "recoDWallZ2", &recoDWallZ2);
   reader->AddVariable( "lambda_max_2", &lambda_max_2);
   /*Float_t spec1,spec2;
   reader->AddSpectator( "spec1:=var1*2",  &spec1 );
   reader->AddSpectator( "spec2:=var1*3",  &spec2 );
   */
   // --- Book the MVA methods

   TString dir    = "weights/";
   TString prefix = "TMVARegression";
   TString THEmethod;

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = it->first + " method";
	 THEmethod=it->first;
         TString weightfile = dir + prefix + "_" + TString(it->first) + ".weights.xml";
         reader->BookMVA( methodName, weightfile ); 
      }
   }
   // Book output histograms
   TH1 *eneu = new TH1F("eneu", "eneu;E_{#nu}", 250, 0., 5000.);
   TH1 *emu = new TH1F("emu", "emu;E_{#mu}", 250, 0., 5000.);
   TH1 *recoTitus= new TH1F("recoTitus", "recoTitus;E_{#mu,reco}", 250, 0., 5000.);
   TH2D *emu_recoTitus  = new TH2D("emu_recoTitus", "emu_recoTitus;E_{#mu};E_{reco}", 250, 0., 5000., 250, 0., 5000.);
   TH2D * hDiffEnemlpmuTitus=new TH2D("hDiffEnemlpmuTitus","hDiffEnemlpmuTitus;E_{#mu};|E_{#mu}-E_{reco}|/E_{#mu}",250, 0., 5000., 100, 0., 100.);
   TH2D * hDiffEnemlpneuTitus=new TH2D("hDiffEnemlpneuTitus","hDiffEnemlpneuTitus;E_{#nu};|E_{#nu}-E_{reco}|/E_{#nu}",250, 0., 5000., 100, 0., 100.);
   TH1 *recoMLP= new TH1F("recoMLP", "recoMLP;E_{#mu}", 250, 0., 5000.);
   TH2D *eneu_mlp = new TH2D("eneu_mlp", "eneu_MLP;E_{#nu};E_{reco}", 250, 0., 5000., 250, 0., 5000.);
   TH2D *emu_mlp = new TH2D("emu_mlp", "emu_MLP;E_{#mu};E_{MLP}", 250, 0., 5000., 250, 0., 5000.);
   TH2D * hDiffEnemlpmu=new TH2D("hDiffEnemlpmu","hDiffEnemlpmu;E_{#mu};|E_{#mu}-E_{MLP}|/E_{#mu}",250, 0., 5000., 100, 0., 100.);
   TH2D * hDiffEnemlpneu=new TH2D("hDiffEnemlpneu","hDiffEnemlpneu;E_{#nu};|E_{#nu}-E_{MLP}|E_{#nu}",250, 0., 5000., 100, 0., 100.);

   TH1* hists[100];  
   Int_t nhists = -1;
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
     TH1* h = new TH1F( it->first.c_str(), TString(it->first) + " method", 500, 0., 5000. );
     //cout<<"it->first.c_str()= "<<it->first.c_str()<<" TString(it->first)= "<<TString(it->first)<<endl;  
     if (it->second){ 
       hists[++nhists] = h;	
     }
   }
   nhists++;
   
   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   
   TFile *input(0);
   //TString fname = "./tmva_reg_example.root";
   for (int i = 1040; i < 1100; i++) {
      char *fname = Form("nu_numu_%i_CCQE_12in_energy_studies_recoquant_tree.root",i);


   if (!gSystem->AccessPathName( fname )) {
      input = TFile::Open( fname ); // check if file in local directory exists
   } 
   else {
      i=i+1;
      char *fname = Form("nu_numu_%i_CCQE_12in_energy_studies_recoquant_tree.root",i);
       input = TFile::Open( fname );  
      //input = TFile::Open( "http://root.cern.ch/files/tmva_reg_example.root" ); // if not: download from ROOT server
   }
   
   if (!input) {
      i=i+1;
      //std::cout << "ERROR: could not open data file" << std::endl;
      //exit(1);
   }
   std::cout << "--- TMVARegressionApp        : Using input file: " << input->GetName() << std::endl;

   // --- Event loop

   // Prepare the tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   /*TTree* theTree = (TTree*)input->Get("TreeR");
   std::cout << "--- Select signal sample" << std::endl;
   theTree->SetBranchAddress( "var1", &var1 );
   theTree->SetBranchAddress( "var2", &var2 );*/
   Double_t trueKE, neutrinoE, recoE_lookup;
   TTree *theTree = (TTree*)input->Get("nu_eneNEW");
   theTree->SetBranchAddress("total_hits2", &total_hits2);
   theTree->SetBranchAddress("total_ring_PEs2", &total_ring_PEs2);
   theTree->SetBranchAddress("recoDWallR2", &recoDWallR2);
   theTree->SetBranchAddress("recoDWallZ2", &recoDWallZ2);
   theTree->SetBranchAddress("lambda_max_2", &lambda_max_2);
   theTree->SetBranchAddress("trueKE", &trueKE);
   theTree->SetBranchAddress("neutrinoE", &neutrinoE);
   theTree->SetBranchAddress("recoE_lookup",&recoE_lookup);

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<theTree->GetEntries();ievt++) {

      if (ievt%1000 == 0) {
         std::cout << "--- ... Processing event: " << ievt << std::endl;
      }

      theTree->GetEntry(ievt);

      // Retrieve the MVA target values (regression outputs) and fill into histograms
      // NOTE: EvaluateRegression(..) returns a vector for multi-target regression

      for (Int_t ih=0; ih<nhists; ih++) {
	TString title = hists[ih]->GetTitle();
	Float_t val = (reader->EvaluateRegression( title ))[0];
	hists[ih]->Fill( val ); 
      }
      	eneu_mlp ->Fill(neutrinoE, val);
	emu_mlp->Fill(trueKE, val);
	hDiffEnemlpmu->Fill(trueKE, (100*abs(trueKE-val)/trueKE));
	hDiffEnemlpneu->Fill(neutrinoE, (100*abs(neutrinoE-val)/neutrinoE));
	eneu ->Fill(neutrinoE);
	emu->Fill(trueKE);
	recoTitus->Fill(recoE_lookup);
	emu_recoTitus->Fill(trueKE, recoE_lookup);
	hDiffEnemlpmuTitus->Fill(trueKE, (100*abs(trueKE-recoE_lookup)/trueKE));
	hDiffEnemlpneuTitus->Fill(neutrinoE, (100*abs(neutrinoE-recoE_lookup)/neutrinoE));
   }//end of entries
   
    }//files
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();
    
    new TCanvas();
    TProfile *mlp= hDiffEnemlpmu->ProfileX("pfx_BDTG", 0., 100., "s");
    //mlp->Reset("s");
    mlp->SetLineColor(2);
    mlp->Draw();
    new TCanvas();
    TProfile *titus= hDiffEnemlpmuTitus->ProfileX("pfx_TITUS", 0., 100., "s");
    //titus->Reset("s");
    titus->SetLineColor(4);
    titus->Draw();
    
    // --- Write histograms
  TFile *target  = new TFile("TMVARegApp_"+THEmethod+"_ener_CCQE_1040_1100_selectedVars.root","RECREATE" );

    for (Int_t ih=0; ih<nhists; ih++){ 
      hists[ih]->Write();
    }
    mlp->Write();
    titus->Write();
    eneu->Write();
    emu->Write();
    recoTitus->Write();
    emu_recoTitus->Write();
    hDiffEnemlpmuTitus->Write();
    hDiffEnemlpneuTitus->Write();
    eneu_mlp ->Write();
    emu_mlp->Write();
    hDiffEnemlpmu->Write();
    hDiffEnemlpneu->Write();
    
    target->Close();
    
    
    std::cout << "--- Created root file: \"" << target->GetName() 
	      << "\" containing the MVA output histograms" << std::endl;
    
    delete reader;
    
    std::cout << "==> TMVARegressionApplication is done!" << std::endl << std::endl;
}
