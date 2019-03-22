#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TF1.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#endif


void MYcalcSigma_scikit_nue()
{
  cout<<"starting.... "<<endl;

  /*TFile * outRootFile = new TFile("test.root", "RECREATE");

  vector<double> r;
  double dr;
  r.push_back(50);
  r.push_back(100);


  vector<TH1F*> hAllPe;
  char nbuff[500];
  char buff[500];
  int nBins = 25;
  TH1F * temp;

  for (int i =0; i < r.size(); i++){
    cout<<"i= "<<i<<" r.at(i): "<<r.at(i)<<endl;
    sprintf(nbuff, "angRespAll_%d", (int) r.at(i));

    sprintf(buff, "Angular response function");
    temp = new TH1F(nbuff, buff, nBins, 0, 1);
    hAllPe.push_back(temp);
  }
*/
  TH2F * plot_emu_sigma_err = new TH2F("plot_emu_sigma_err","MLP;Log(MC Muon Energy (GeV));RMS;", 31, 1.8, 8.,100, 0., 1.);
  TH2F * plot_emu_LOGE_MLP_diff = new TH2F("plot_emu_LOGE_MLP_diff","MLP;Log(MC Muon Energy (GeV));LogEmlp-LogEmc;", 31, 1.8, 8.,100, -1., 1.);

  TH1D *EnerResoBDTG50 =new TH1D("EnerResoBDTG50","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG100 =new TH1D("EnerResoBDTG100","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG150 =new TH1D("EnerResoBDTG150","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG200 =new TH1D("EnerResoBDTG200","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG250 =new TH1D("EnerResoBDTG250","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG300 =new TH1D("EnerResoBDTG300","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG350 =new TH1D("EnerResoBDTG350","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG400 =new TH1D("EnerResoBDTG400","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG450 =new TH1D("EnerResoBDTG450","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG500 =new TH1D("EnerResoBDTG500","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG550 =new TH1D("EnerResoBDTG550","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG600 =new TH1D("EnerResoBDTG600","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG650 =new TH1D("EnerResoBDTG650","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG700 =new TH1D("EnerResoBDTG700","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG750 =new TH1D("EnerResoBDTG750","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG800 =new TH1D("EnerResoBDTG800","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG850 =new TH1D("EnerResoBDTG850","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG900 =new TH1D("EnerResoBDTG900","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG950 =new TH1D("EnerResoBDTG950","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1000 =new TH1D("EnerResoBDTG1000","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);   

  TH1D *EnerResoBDTG1050 =new TH1D("EnerResoBDTG1050","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1100 =new TH1D("EnerResoBDTG1100","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1150 =new TH1D("EnerResoBDTG1150","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1200 =new TH1D("EnerResoBDTG1200","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1250 =new TH1D("EnerResoBDTG1250","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1300 =new TH1D("EnerResoBDTG1300","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1350 =new TH1D("EnerResoBDTG1350","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1400 =new TH1D("EnerResoBDTG1400","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1450 =new TH1D("EnerResoBDTG1450","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1500 =new TH1D("EnerResoBDTG1500","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1550 =new TH1D("EnerResoBDTG1550","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1600 =new TH1D("EnerResoBDTG1600","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1650 =new TH1D("EnerResoBDTG1650","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1700 =new TH1D("EnerResoBDTG1700","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1750 =new TH1D("EnerResoBDTG1750","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1800 =new TH1D("EnerResoBDTG1800","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1850 =new TH1D("EnerResoBDTG1850","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1900 =new TH1D("EnerResoBDTG1900","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG1950 =new TH1D("EnerResoBDTG1950","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2000 =new TH1D("EnerResoBDTG2000","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
 
  TH1D *EnerResoBDTG2050 =new TH1D("EnerResoBDTG2050","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2100 =new TH1D("EnerResoBDTG2100","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2150 =new TH1D("EnerResoBDTG2150","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2200 =new TH1D("EnerResoBDTG2200","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2250 =new TH1D("EnerResoBDTG2250","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2300 =new TH1D("EnerResoBDTG2300","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2350 =new TH1D("EnerResoBDTG2350","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2400 =new TH1D("EnerResoBDTG2400","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2450 =new TH1D("EnerResoBDTG2450","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2500 =new TH1D("EnerResoBDTG2500","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2550 =new TH1D("EnerResoBDTG2550","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2600 =new TH1D("EnerResoBDTG2600","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2650 =new TH1D("EnerResoBDTG2650","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2700 =new TH1D("EnerResoBDTG2700","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2750 =new TH1D("EnerResoBDTG2750","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2800 =new TH1D("EnerResoBDTG2800","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2850 =new TH1D("EnerResoBDTG2850","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2900 =new TH1D("EnerResoBDTG2900","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG2950 =new TH1D("EnerResoBDTG2950","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);
  TH1D *EnerResoBDTG3000 =new TH1D("EnerResoBDTG3000","Energy Resolution using Scikit BDTG; #Delta E/E [%]", 200, -100., 100.);



   TFile *input(0);  
   char *fname = Form("/home/edrakopo/work/ener_reco_scikit/scikit_recoLASTnue.root");

  if (!gSystem->AccessPathName( fname )){ 
    cout<<"Open file "<<endl;
    input = TFile::Open( fname ); 
  }

  TTree *regTree = (TTree*)input->Get("tuple");

  Float_t trueKE, neutrinoE, recoE_lookup, recoKE, recoKE_BDTG, recoKE_BDTG2, recoKE_BDTG_par0, recoKE_BDTG_par1, recoKE_BDTG_par2, recoKE_BDTG_par3, recoKE_BDTG_par4, recoKE_BDTG_par5, recoKE_BDTG01, recoKE_BDTG_par02, recoDwall, recoToWall;
  Float_t total_PMTs_hits2, total_hits2, total_ring_PEs2, pot_length2, hits_pot_length2, recoDWallR2, recoDWallZ2, lambda_max_2, myweight;
  regTree->SetBranchAddress("trueKE", &trueKE);
  regTree->SetBranchAddress("neutrinoE", &neutrinoE);
  //regTree->SetBranchAddress("recoKE", &recoKE);
  regTree->SetBranchAddress("recoKE_BDTG", &recoKE_BDTG);
  regTree->SetBranchAddress("recoDwall", &recoDwall);
  regTree->SetBranchAddress("recoToWall", &recoToWall);


   Float_t sigma_err[60]; Float_t meanDlogE[60]; Float_t n[60]; 

  for (Long64_t ievt=0; ievt<regTree->GetEntries();ievt++) {
  //for (Long64_t ievt=0; ievt<100;ievt++) {
     if (ievt%1000 == 0) {
        std::cout << "--- ... Processing event: " << ievt << std::endl;
     }
     regTree->GetEntry(ievt);
      recoKE=recoKE_BDTG;

    if( (recoDwall*550.)>100. && (recoToWall*2500.)>100.){

     if(trueKE<=50.){
       EnerResoBDTG50->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>50. && trueKE<=100.){
       EnerResoBDTG100->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>100. && trueKE<=150.){
       EnerResoBDTG150->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>150. && trueKE<=200.){
       EnerResoBDTG200->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>200. && trueKE<=250.){
       EnerResoBDTG250->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>250. && trueKE<=300.){
       EnerResoBDTG300->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>300. && trueKE<=350.){
       EnerResoBDTG350->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>350. && trueKE<=400.){
       EnerResoBDTG400->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>400. && trueKE<=450.){
       EnerResoBDTG450->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>450. && trueKE<=500.){
       EnerResoBDTG500->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>500. && trueKE<=550.){
       EnerResoBDTG550->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>550. && trueKE<=600.){
       EnerResoBDTG600->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>600. && trueKE<=650.){
       EnerResoBDTG650->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>650. && trueKE<=700.){
       EnerResoBDTG700->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>700. && trueKE<=750.){
       EnerResoBDTG750->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>750. && trueKE<=800.){
       EnerResoBDTG800->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>800. && trueKE<=850.){
       EnerResoBDTG850->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>850. && trueKE<=900.){
       EnerResoBDTG900->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>900. && trueKE<=950.){
       EnerResoBDTG950->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>950. && trueKE<=1000.){
       EnerResoBDTG1000->Fill(100*(trueKE-recoKE)/trueKE);
     }
   ///----
     if(trueKE>1000. && trueKE<=1050.){
       EnerResoBDTG1050->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1050. && trueKE<=1100.){
       EnerResoBDTG1100->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1100. && trueKE<=1150.){
       EnerResoBDTG1150->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1150. && trueKE<=1200.){
       EnerResoBDTG1200->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1200. && trueKE<=1250.){
       EnerResoBDTG1250->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1250. && trueKE<=1300.){
       EnerResoBDTG1300->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1300. && trueKE<=1350.){
       EnerResoBDTG1350->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1350. && trueKE<=1400.){
       EnerResoBDTG1400->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1400. && trueKE<=1450.){
       EnerResoBDTG1450->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1450. && trueKE<=1500.){
       EnerResoBDTG1500->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1500. && trueKE<=1550.){
       EnerResoBDTG550->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1550. && trueKE<=1600.){ 
       EnerResoBDTG1600->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1600. && trueKE<=1650.){
       EnerResoBDTG1650->Fill(100*(trueKE-recoKE)/trueKE);
     } 
     if(trueKE>1650. && trueKE<=1700.){
       EnerResoBDTG1700->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1700. && trueKE<=1750.){
       EnerResoBDTG1750->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1750. && trueKE<=1800.){
       EnerResoBDTG1800->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1800. && trueKE<=1850.){
       EnerResoBDTG1850->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1850. && trueKE<=1900.){
       EnerResoBDTG1900->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1900. && trueKE<=1950.){
       EnerResoBDTG1950->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>1950. && trueKE<=2000.){
       EnerResoBDTG2000->Fill(100*(trueKE-recoKE)/trueKE);
     }
    //------
    if(trueKE>2000. && trueKE<=2050.){
       EnerResoBDTG2050->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2050. && trueKE<=2100.){
       EnerResoBDTG2100->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2100. && trueKE<=2150.){
       EnerResoBDTG2150->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2150. && trueKE<=2200.){
       EnerResoBDTG2200->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2200. && trueKE<=2250.){
       EnerResoBDTG2250->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2250. && trueKE<=2300.){
       EnerResoBDTG2300->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2300. && trueKE<=2350.){
       EnerResoBDTG2350->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2350. && trueKE<=2400.){
       EnerResoBDTG2400->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2400. && trueKE<=2450.){
       EnerResoBDTG2450->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2450. && trueKE<=2500.){
       EnerResoBDTG2500->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2500. && trueKE<=2550.){
       EnerResoBDTG2550->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2550. && trueKE<=2600.){
       EnerResoBDTG2600->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2600. && trueKE<=2650.){
       EnerResoBDTG2650->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2650. && trueKE<=2700.){
       EnerResoBDTG2700->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2700. && trueKE<=2750.){
       EnerResoBDTG2750->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2750. && trueKE<=2800.){
       EnerResoBDTG2800->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2800. && trueKE<=2850.){
       EnerResoBDTG2850->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2850. && trueKE<=2900.){
       EnerResoBDTG2900->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2900. && trueKE<=2950.){
       EnerResoBDTG2950->Fill(100*(trueKE-recoKE)/trueKE);
     }
     if(trueKE>2950. && trueKE<=3000.){
       EnerResoBDTG3000->Fill(100*(trueKE-recoKE)/trueKE);
     }

    //
    /*for (int shell = 0; shell < r.size(); shell++){
      
      hAllPe[shell]->Fill(100*(trueKE-recoKE)/trueKE);
    }*/
     /*
     if(trueKE>0. && trueKE<=50.){ 
        n_50++;
        sigma_50 +=((recoKE-trueKE)*(recoKE-trueKE)); 
        DlogE_50 += (recoKE-trueKE);
     }*/

    if( (recoDwall*550.)>100. && (recoToWall*2500.)>100.){
     }
     

     /*if(trueKE>200. && trueKE<600.){ //cout<<"in"<<" DE/E: "<<trueKE-recoKE_BDTG/trueKE<<endl;
         //plot_emu_erecoBDTG->Fill(100.*(trueKE-recoKE)/trueKE);
     }*/

   }//Infid
  }//entries

  sigma_err[0]=EnerResoBDTG50->GetStdDev();  //0.5*sqrt(sigma_50/n_50);
  meanDlogE[0]=EnerResoBDTG50->GetMean();  //
  sigma_err[1]=EnerResoBDTG100->GetStdDev();  //0.5*sqrt(sigma_50/n_50);
  meanDlogE[1]=EnerResoBDTG100->GetMean(); 
  sigma_err[2]=EnerResoBDTG150->GetStdDev();  //0.5*sqrt(sigma_50/n_50);
  meanDlogE[2]=EnerResoBDTG150->GetMean();
  sigma_err[3]=EnerResoBDTG200->GetStdDev();  
  meanDlogE[3]=EnerResoBDTG200->GetMean();
  sigma_err[4]=EnerResoBDTG250->GetStdDev();  
  meanDlogE[4]=EnerResoBDTG250->GetMean();
  sigma_err[5]=EnerResoBDTG300->GetStdDev();  
  meanDlogE[5]=EnerResoBDTG300->GetMean();
  sigma_err[6]=EnerResoBDTG350->GetStdDev();  
  meanDlogE[6]=EnerResoBDTG350->GetMean();
  sigma_err[7]=EnerResoBDTG400->GetStdDev();  
  meanDlogE[7]=EnerResoBDTG400->GetMean();
  sigma_err[8]=EnerResoBDTG450->GetStdDev();  
  meanDlogE[8]=EnerResoBDTG450->GetMean();
  sigma_err[9]=EnerResoBDTG500->GetStdDev();  
  meanDlogE[9]=EnerResoBDTG500->GetMean();
  sigma_err[10]=EnerResoBDTG550->GetStdDev();
  meanDlogE[10]=EnerResoBDTG550->GetMean();
  sigma_err[11]=EnerResoBDTG600->GetStdDev();
  meanDlogE[11]=EnerResoBDTG600->GetMean(); 
  sigma_err[12]=EnerResoBDTG650->GetStdDev();
  meanDlogE[12]=EnerResoBDTG650->GetMean();
  sigma_err[13]=EnerResoBDTG700->GetStdDev();
  meanDlogE[13]=EnerResoBDTG700->GetMean();
  sigma_err[14]=EnerResoBDTG750->GetStdDev();
  meanDlogE[14]=EnerResoBDTG750->GetMean();
  sigma_err[15]=EnerResoBDTG800->GetStdDev();
  meanDlogE[15]=EnerResoBDTG800->GetMean();
  sigma_err[16]=EnerResoBDTG850->GetStdDev();
  meanDlogE[16]=EnerResoBDTG850->GetMean();
  sigma_err[17]=EnerResoBDTG900->GetStdDev();
  meanDlogE[17]=EnerResoBDTG900->GetMean();
  sigma_err[18]=EnerResoBDTG950->GetStdDev();
  meanDlogE[18]=EnerResoBDTG950->GetMean();
  sigma_err[19]=EnerResoBDTG1000->GetStdDev();
  meanDlogE[19]=EnerResoBDTG1000->GetMean();
  //
  sigma_err[20]=EnerResoBDTG1050->GetStdDev();  
  meanDlogE[20]=EnerResoBDTG1050->GetMean();  
  sigma_err[21]=EnerResoBDTG1100->GetStdDev();  
  meanDlogE[21]=EnerResoBDTG1100->GetMean();
  sigma_err[22]=EnerResoBDTG1150->GetStdDev();  
  meanDlogE[22]=EnerResoBDTG1150->GetMean();
  sigma_err[23]=EnerResoBDTG1200->GetStdDev();
  meanDlogE[23]=EnerResoBDTG1200->GetMean();
  sigma_err[24]=EnerResoBDTG1250->GetStdDev();
  meanDlogE[24]=EnerResoBDTG1250->GetMean();
  sigma_err[25]=EnerResoBDTG1300->GetStdDev();
  meanDlogE[25]=EnerResoBDTG1300->GetMean();
  sigma_err[26]=EnerResoBDTG1350->GetStdDev();
  meanDlogE[26]=EnerResoBDTG1350->GetMean();
  sigma_err[27]=EnerResoBDTG1400->GetStdDev();
  meanDlogE[27]=EnerResoBDTG1400->GetMean();
  sigma_err[28]=EnerResoBDTG1450->GetStdDev();
  meanDlogE[28]=EnerResoBDTG1450->GetMean();
  sigma_err[29]=EnerResoBDTG1500->GetStdDev();
  meanDlogE[29]=EnerResoBDTG1500->GetMean();
  sigma_err[30]=EnerResoBDTG1550->GetStdDev();
  meanDlogE[30]=EnerResoBDTG1550->GetMean();
  sigma_err[31]=EnerResoBDTG1600->GetStdDev();
  meanDlogE[31]=EnerResoBDTG1600->GetMean();
  sigma_err[32]=EnerResoBDTG1650->GetStdDev();
  meanDlogE[32]=EnerResoBDTG1650->GetMean();
  sigma_err[33]=EnerResoBDTG1700->GetStdDev();
  meanDlogE[33]=EnerResoBDTG1700->GetMean();
  sigma_err[34]=EnerResoBDTG1750->GetStdDev();
  meanDlogE[34]=EnerResoBDTG1750->GetMean();
  sigma_err[35]=EnerResoBDTG1800->GetStdDev();
  meanDlogE[35]=EnerResoBDTG1800->GetMean();
  sigma_err[36]=EnerResoBDTG1850->GetStdDev();
  meanDlogE[36]=EnerResoBDTG1850->GetMean();
  sigma_err[37]=EnerResoBDTG1900->GetStdDev();
  meanDlogE[37]=EnerResoBDTG1900->GetMean();
  sigma_err[38]=EnerResoBDTG1950->GetStdDev();
  meanDlogE[38]=EnerResoBDTG1950->GetMean();
  sigma_err[39]=EnerResoBDTG2000->GetStdDev();
  meanDlogE[39]=EnerResoBDTG2000->GetMean();
  //
  sigma_err[40]=EnerResoBDTG2050->GetStdDev();  
  meanDlogE[40]=EnerResoBDTG2050->GetMean();  
  sigma_err[41]=EnerResoBDTG2100->GetStdDev();  
  meanDlogE[41]=EnerResoBDTG2100->GetMean();
  sigma_err[42]=EnerResoBDTG2150->GetStdDev();  
  meanDlogE[42]=EnerResoBDTG2150->GetMean();
  sigma_err[43]=EnerResoBDTG2200->GetStdDev();
  meanDlogE[43]=EnerResoBDTG2200->GetMean();
  sigma_err[44]=EnerResoBDTG2250->GetStdDev();
  meanDlogE[44]=EnerResoBDTG2250->GetMean();
  sigma_err[45]=EnerResoBDTG2300->GetStdDev();
  meanDlogE[45]=EnerResoBDTG2300->GetMean();
  sigma_err[46]=EnerResoBDTG2350->GetStdDev();
  meanDlogE[46]=EnerResoBDTG2350->GetMean();
  sigma_err[47]=EnerResoBDTG2400->GetStdDev();
  meanDlogE[47]=EnerResoBDTG2400->GetMean();
  sigma_err[48]=EnerResoBDTG2450->GetStdDev();
  meanDlogE[48]=EnerResoBDTG2450->GetMean();
  sigma_err[49]=EnerResoBDTG2500->GetStdDev();
  meanDlogE[49]=EnerResoBDTG2500->GetMean();
  sigma_err[50]=EnerResoBDTG2550->GetStdDev();
  meanDlogE[50]=EnerResoBDTG2550->GetMean();
  sigma_err[51]=EnerResoBDTG2600->GetStdDev();
  meanDlogE[51]=EnerResoBDTG2600->GetMean();
  sigma_err[52]=EnerResoBDTG2650->GetStdDev();
  meanDlogE[52]=EnerResoBDTG2650->GetMean();
  sigma_err[53]=EnerResoBDTG2700->GetStdDev();
  meanDlogE[53]=EnerResoBDTG2700->GetMean();
  sigma_err[54]=EnerResoBDTG2750->GetStdDev();
  meanDlogE[54]=EnerResoBDTG2750->GetMean();
  sigma_err[55]=EnerResoBDTG2800->GetStdDev();
  meanDlogE[55]=EnerResoBDTG2800->GetMean();
  sigma_err[56]=EnerResoBDTG2850->GetStdDev();
  meanDlogE[56]=EnerResoBDTG2850->GetMean();
  sigma_err[57]=EnerResoBDTG2900->GetStdDev();
  meanDlogE[57]=EnerResoBDTG2900->GetMean();
  sigma_err[58]=EnerResoBDTG2950->GetStdDev();
  meanDlogE[58]=EnerResoBDTG2950->GetMean();
  sigma_err[59]=EnerResoBDTG3000->GetStdDev();
  meanDlogE[59]=EnerResoBDTG3000->GetMean();

  //------------------------------------
   Float_t emu_=0.; Float_t muonE[60];

   for(int em=0; em<60;){ 
    n[em]=0.;
    muonE[em]=emu_;
    cout<<"n[em]: "<<n[em]<<" em: "<<em<<" emu_ "<<muonE[em]<<" sigma_err[em]: "<<sigma_err[em]<<" meanDlogE[em]: "<<meanDlogE[em]<<endl;
    //plot_emu_sigma_err -> Fill(emu_,2*sigma_err[em]);
    //plot_emu_LOGE_MLP_diff -> Fill(emu_, meanDlogE[em]);
    TGraphErrors* ge = new TGraphErrors(60, muonE, meanDlogE, n, sigma_err);
 
    em=em+1;
    emu_+=50.;
   }
  //- - - - - - - -
  /* const Int_t nq = 100;
   const Int_t nshots = 100;
   Double_t xq[nq];  // position where to compute the quantiles in [0,1]
   Double_t yq[nq];  // array to contain the quantiles
   for (Int_t i=0;i<nq;){ 
     xq[i] = Float_t(i+1)/nq;
     //cout<<"i: "<<i<<" xq[i]: "<<xq[i]<<endl;
     i=i+100;
   }
   TGraph *gr68 = new TGraph(nshots);
   //TH1F *h = new TH1F("h","demo quantiles",50,-3,3);

   for (Int_t shot=0;shot<nshots;) {
      EnerResoBDTG50->GetQuantiles(nq,yq,xq);
      gr68->SetPoint(shot,shot+1,yq[68]);
      cout<<"shot: "<<shot<<" yq[68]: "<<yq[68]<<endl;
      shot=shot+100;
  }
  */
//_______________________________________
 new TCanvas();
 ge->SetFillColor(4);
 ge->SetFillStyle(3010);
 ge->SetMarkerStyle(1);
 ge->Draw("ALP ||");
 //ge->Draw();  
 //EnerResoBDTG200->Draw(); 
 //EnerResoBDTG2000->Draw();
 
 //new TCanvas();
 //EnerResoBDTG1500->Draw();
 //EnerResoBDTG300->Draw();
 
 //gr68->Draw();
 //plot_emu_LOGE_MLP_diff->Draw(); 

  //TFile* outputFile = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/SigmaErrorScikit_nue.root", "RECREATE" );
  TFile* outputFile = new TFile("/Disk/ds-sopa-group/PPE/titus/ts-WChRecoSandBox/scripts/tmva_ene/SigmaErrorScikit_nueINfid.root", "RECREATE" );

  ge->Write();

  outputFile->Close();
}
