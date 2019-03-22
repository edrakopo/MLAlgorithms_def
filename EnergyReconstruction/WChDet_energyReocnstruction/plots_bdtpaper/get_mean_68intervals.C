#include "TFile.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPolyLine.h"

void get_mean_68intervals(){

  gROOT -> ProcessLine( ".x ./mattStyle.C" );

  /*TFile* file = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/FINALTMVARegApp_BDTG_ener_lilia_CCQE_1040_1100_NEWlookupsTEST_MYselectedVars_large_statTEST10000treesB2.root"); //TMVABDTG
  TFile* file1 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/EnergyRecoScikitNEW2_200_600MeV.root");
  TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/EnergyRecoKeras_TrainTestSampleG.root");
  TFile* file3 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/FINALTMVARegApp_MLP_ener_lilia_CCQE_1040_1100_NEWlookupsTEST_MYselectedVars_large_statTEST10000treesB.root"); //MLP
*/
  TFile* file = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/TMVABDTG_INFID_enerRsol.root");  //TMVABDTG
  TFile* file1 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/EnergyRecoScikitNEW2_200_600MeV.root");
  TFile* file2 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/EnergyRecoKeras_TrainTestSampleG.root");
  TFile* file3 = new TFile("/Users/edrakopo/work/work_energy_reco/plots_bdtpaper/TMVAMLP_INFID_enerRsol.root"); 

  //TH1D *plot_emu_erecoBDTG =(TH1D*)file->Get("ener_resolution");  //TMVABDTG
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file1->Get("plot_emu_erecoBDTG"); //scikit
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file2->Get("plot_emu_erecoKERAS"); //keras
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file3->Get("ener_resolution"); //MLP
  //TH1D *plot_emu_erecoBDTG=(TH1D*)file->Get("ener_resolutionTITUS");
 //----INFID
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file->Get("ener_resolution"); //TMVABDTG
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file1->Get("plot_emu_erecoBDTG_Infid"); //scikit
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file2->Get("plot_emu_erecoBDTG_Infid"); //keras
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file3->Get("ener_resolution"); //MLP
  //TH1D *plot_emu_erecoBDTG=(TH1D*)file->Get("ener_resolutionTITUS");
  // cout<<"-- " <<plot_emu_erecoBDTG->GetNbinsX()<<endl;	
  //------- 600-1400 MeV  
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file->Get("plot_emu_erecoBDTG_he"); //TMVABDTG  
  TH1D *plot_emu_erecoBDTG =(TH1D*)file1->Get("plot_emu_erecoBDTG_he"); //scikit
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file2->Get("plot_emu_erecoBDTG_he"); //keras
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file3->Get("plot_emu_erecoBDTG_he");
  //TH1D *plot_emu_erecoBDTG=(TH1D*)file->Get("plot_emu_erecoBDTG_heTITUS"); //TITUS

  //----INFID 600-1400 MeV
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file->Get("plot_emu_erecoBDTG_Infid_he"); //TMVABDTG 
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file1->Get("plot_emu_erecoBDTG_Infid_he"); //scikit
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file2->Get("plot_emu_erecoBDTG_Infid_he"); //keras
  //TH1D *plot_emu_erecoBDTG =(TH1D*)file3->Get("plot_emu_erecoBDTG_Infid_he"); //MLP
  //TH1D *plot_emu_erecoBDTG=(TH1D*)file->Get("plot_emu_erecoBDTG_Infid_heTITUS"); //TITUS

  int max_Nentries=0; int max_bin=0; double bincenter=0.; int allentries=0; 
  int sumevts1=0; bool first1=1; int sumevts2=0; bool first2=1; int binL=0; int binU=0;
  int sumevts_outL=0; int sumevts_outU=0;

  for(int i=1; i<201; i++){
    double binContent = plot_emu_erecoBDTG->GetBinContent(i);
    allentries+=binContent;
    if(binContent>max_Nentries){
      max_bin=i;
      max_Nentries=binContent;
      bincenter=plot_emu_erecoBDTG->GetBinCenter(i);
    }
    //cout<<"binContent1: "<<binContent1<<endl;     
  } 
  cout<<"allentries: "<<allentries<<" max_Nentries: "<<max_Nentries<<" in: "<<max_bin<<" at MEAN: "<<bincenter<<endl;
  int max_Nentries2=int((1.*max_Nentries)/2.);
  cout<<"max_Nentries2-> "<<max_Nentries2<<" | "<<0.34*allentries<<endl; 
 
  for(int i=(max_bin+1); i<201; i++){
      double binContent1 = plot_emu_erecoBDTG->GetBinContent(i);  
      sumevts1+=binContent1;
     if(sumevts1<0.34*allentries){
      cout<<"i: "<<i<<" binContent1: "<<binContent1<<" sumevts1: "<<sumevts1<<" ALL: "<<(max_Nentries2+sumevts1)<<endl;
     }
     if((sumevts1+max_Nentries2)>=0.34*allentries){ 
      if(first1==1){ binU=i;
         cout<<"Upper Limit at: "<<plot_emu_erecoBDTG->GetBinCenter(i)<<" bin: "<<binU<<" sumevts1: "<<sumevts1<<endl; }
      first1=0;   }
  }

  for(int i=(max_bin-1); i>1; i--){
      double binContent2 = plot_emu_erecoBDTG->GetBinContent(i);
      sumevts2+=binContent2;
      if(sumevts2<(0.34*allentries)){
       cout<<"i: "<<i<<" binContent2: "<<binContent2<<" sumevts2: "<<sumevts2<<" ALL: "<<(max_Nentries2+sumevts2)<<endl;
      }
      if((max_Nentries2+sumevts2)>=0.34*allentries){
       if(first2==1){ binL=i;
          cout<<"Lower Limit at: "<<plot_emu_erecoBDTG->GetBinCenter(i)<<" bin: "<<binL<<" sumevts2: "<<sumevts2<<endl; }
      first2=0;   }
  }

  for(int i=1; i<201; i++){ 
    if(binL>0){
     if(i<binL){ //cout<<"iL: "<<i<<endl;
       sumevts_outL+=plot_emu_erecoBDTG->GetBinContent(i);  } }
  }
  for(int i=1; i<201; i++){
    if(binU>0.){
     if(i>binU){ //cout<<"iU: "<<i<<" binU: "<<binU<<endl;
        sumevts_outU+=plot_emu_erecoBDTG->GetBinContent(i);  } }
  }
  cout<<"#Events below L: "<<sumevts_outL<<" #Events above U: "<<sumevts_outU<<endl;
  cout<<"Percentage of events outside 68%: "<<1.*(sumevts_outL+sumevts_outU)/allentries<<endl;

}
