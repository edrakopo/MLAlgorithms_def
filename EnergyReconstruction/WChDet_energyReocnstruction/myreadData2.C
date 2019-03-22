#include "TFile.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TPolyLine.h"
#include "TTree.h"
#include <iostream>
#include <string>
#include <vector>
#include "MakeBins.C"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TH2D.h"
#include "TCanvas.h"
/*
class StandardTable{

public:  
  
  StandardTable(std::string lookupfilename);
  
  ~StandardTable(){}

  double eval(double rwall, double dwall, double npes, bool isMuon);


  
private:

  TH2D* muonEnergyLookup;
  TH2D *electronEnergyLookup;
  
};

StandardTable::StandardTable(std::string lookupfilename){

  TFile *energyTables = new TFile(lookupfilename.c_str(), "READ");
  electronEnergyLookup = (TH2D *) energyTables->Get("nHitDWallKELookup_e");
  muonEnergyLookup = (TH2D *) energyTables->Get("nHitDWallKELookup_mu");
  
}

double StandardTable::eval(double dWallR, double dWallZ , double lookupPEs ,bool isMuon   ){

 double dWall = dWallR < dWallZ ? dWallR : dWallZ;
 
 double corr_factor_phot_coverage = 0.9; // 36%/40%=0.9
 
 double dWallMin = muonEnergyLookup->GetYaxis()->GetBinCenter(muonEnergyLookup->GetYaxis()->GetFirst());
 double dWallMax = muonEnergyLookup->GetYaxis()->GetBinCenter(muonEnergyLookup->GetYaxis()->GetLast());
 double PEmin = muonEnergyLookup->GetXaxis()->GetBinCenter(muonEnergyLookup->GetXaxis()->GetFirst());
 double PEmax = muonEnergyLookup->GetXaxis()->GetBinCenter(muonEnergyLookup->GetXaxis()->GetLast());
 if (dWall <= dWallMin) dWall = dWallMin + 0.01;
 else if (dWall >= dWallMax) dWall = dWallMax - 0.01;
 double correctionFactor = 1.;
 if (lookupPEs <= PEmin) {
    correctionFactor = lookupPEs / PEmin;
    lookupPEs = PEmin + 1;
 }
 else if (lookupPEs >= PEmax) {
   correctionFactor = lookupPEs / PEmax;
   lookupPEs = PEmax - 1;
 } 

 return isMuon ? muonEnergyLookup->Interpolate(lookupPEs, dWall)*correctionFactor :
                 electronEnergyLookup->Interpolate(lookupPEs, dWall)*correctionFactor;   
}
*/

void plotGraph(TGraph* graph, std::string xTitle ="#PE^{ring}/20000",
               std::string yTitle ="E_{#mu} [MeV]", int marker = 20, int color = 1, int style = 1 ){

  graph->SetMarkerStyle(marker);
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetLineStyle(style);
  graph->SetLineWidth(1);
  //  graph->SetMinimum(0.0);
  graph->SetTitle("");
  graph->GetXaxis()->SetTitleFont(132);
  graph->GetYaxis()->SetTitleFont(132);
  graph->GetXaxis()->SetTitleOffset(1.2);
  graph->GetYaxis()->SetTitleOffset(1.1);
  graph->GetXaxis()->SetLabelFont(132);
  graph->GetYaxis()->SetLabelFont(132);
  graph->GetYaxis()->SetTitle(yTitle.c_str());
  graph->GetXaxis()->SetTitle(xTitle.c_str());
  //graph->GetYaxis()->SetTitle("");
}


MakeBins* divideDatabyE(TTree* tree, int nbin = 15){

  //split up the data into bins
  double trueE;
  tree->SetBranchAddress("trueKE",&trueE);
  
  std::vector<double> energies; energies.resize(tree->GetEntries());
  for (auto i = 0; i < tree->GetEntries(); ++i){
    tree->GetEntry(i);
    energies[i] = trueE;
  } //i

  MakeBins* bins = new MakeBins(energies,nbin);
  return bins;
}

int toRZbin(double recoDWallR2, double recoDWallZ2){

  int bin = 0;
  if (recoDWallZ2 < 0.4 && recoDWallR2 < 0.4){
  //if (recoDWallZ2 < 0.3 && recoDWallR2 < 0.3){
  //if (recoToWall_2 < 0.4 && recoDWall_2 < 0.4){
     bin = 0;
  }
  else if ((recoDWallZ2 > 0.4 && recoDWallZ2 < 0.6 && recoDWallR2 < 0.4) || (recoDWallR2 <0.6 && recoDWallR2 > 0.4 && recoDWallZ2 <0.6)){
  //else if ((recoDWallZ2 > 0.3 && recoDWallZ2 < 0.6 && recoDWallR2 < 0.3) || (recoDWallR2 <0.6 && recoDWallR2 > 0.3 && recoDWallZ2 <0.6)){ 
  //else if ((recoToWall_2 > 0.4 && recoToWall_2 < 0.6 && recoDWall_2 < 0.4) || (recoDWall_2 <0.6 && recoDWall_2 > 0.4 && recoToWa    ll_2 <0.6)){
      bin = 1;
  }
  else {
     bin = 2;
  } 
  
  return bin;
}

MakeBins* divideDataByERZ(int requiredbin, TTree* tree, int nbin = 10){

  //split up the data into bins
  double trueE;
  tree->SetBranchAddress("trueKE",&trueE);

  float recoDWallR2;
  tree->SetBranchAddress("recoDWallR2",&recoDWallR2);

  float recoDWallZ2;
  tree->SetBranchAddress("recoDWallZ2",&recoDWallZ2);

  float recoDWall_2;
  tree->SetBranchAddress("recoDWall_2",&recoDWall_2);
   
  float recoToWall_2;
  tree->SetBranchAddress("recoToWall_2",&recoToWall_2);

  std::vector<double> energies; energies.reserve(tree->GetEntries());
 
 
  
  for (auto i = 0; i < tree->GetEntries(); ++i){
    tree->GetEntry(i);
    //int bin = toRZbin(recoDWallR2, recoDWallZ2);
    int bin = toRZbin(recoDWall_2, recoToWall_2);
    if (bin == requiredbin){
      energies.push_back(trueE);
    }
  } //i
  
  MakeBins* bins = new MakeBins(energies,nbin);
 
  return bins;
}


typedef std::pair<TF1*,MakeBins*> LUT;
LUT makeLUT(TTree* tree,int nbin = 15, int requiredbin = 3){

  
  MakeBins* bins;
  if (requiredbin == 3) {
    bins = divideDatabyE(tree,nbin);
  }
  else {
    bins = divideDataByERZ(requiredbin,tree,nbin);
  }

  float total_ring_PEs2;
  tree->SetBranchAddress("total_ring_PEs2",&total_ring_PEs2);

  double trueE;
  tree->SetBranchAddress("trueKE",&trueE); 
  
  // histos of pes in each bin;
  std::vector<TH1F*> histos;   std::vector<TH1F*> ehistos; 
  for (auto ih= 0; ih < nbin ;++ih){
    std::string name1 = "h_bin" + to_string(requiredbin) + "_" + to_string(ih);
    std::string name2 = "eh_bin" + to_string(requiredbin) + "_"+ to_string(ih); 
    TH1F* h = new TH1F(name1.c_str(),name1.c_str(), 100, 0, 1.5); h->Sumw2();
    TH1F* eh = new TH1F(name2.c_str(),name2.c_str(), 100, 0, 5000); eh->Sumw2();
    histos.push_back(h); ehistos.push_back(eh);
  }
  
  // now look and for each bin get the mean hits or whatever
  for (auto i = 0; i < tree->GetEntries(); ++i){
    tree->GetEntry(i);
    // if (total_ring_PEs2 > 2e-3){
      int binNum = TMath::Min((int)bins->toBin(trueE),nbin);
      //  std::cout << binNum <<" " << trueE << std::endl;
      histos[binNum]->Fill(total_ring_PEs2);
      ehistos[binNum]->Fill(trueE);
      // }
  }

  // now make a tgraph
  TGraphErrors* graph = new TGraphErrors(nbin);
  for (auto i = 0; i < nbin; ++i){
    graph->SetPoint(i,histos[i]->GetMean(), bins->binCenter(i));
    graph->SetPointError(i,histos[i]->GetMeanError(), ehistos[i]->GetMeanError());
  }

  TCanvas* can = new TCanvas("can","can", 600, 400);
  graph->SetMarkerStyle(20);
  graph->Fit("pol3");
  TF1* f = graph->GetFunction("pol3");
  f->SetLineColor(2);
  f->SetLineStyle(2);
  plotGraph(graph);
  graph->Draw("AP");
  f->Draw("SAME");
  can->Print("mylut.png");
  
  std::string funname = "fun" + to_string(requiredbin) ;
  return std::make_pair((TF1*)f->Clone(funname.c_str()), bins); 
}




/*
void readData(int nbin = 15){

  gROOT -> ProcessLine( ".x ./mattStyle.C" );
  
  TFile* file = new TFile("nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree.root");
  TTree* tree = (TTree*)file->Get("nu_eneNEW");

  // opening the file
  LUT theLUT = makeLUT(tree,nbin);
  TF1* f = theLUT.first;
  f->Print();
  TFormula* form = f->GetFormula();
  form->Print();
  
  // setting the branches
  float total_ring_PEs2;
  tree->SetBranchAddress("total_ring_PEs2",&total_ring_PEs2);
  double trueE;
  tree->SetBranchAddress("trueKE",&trueE); 

  //set up histogram
  std::vector<TH1F*> histos;    std::vector<TH1F*> dhistos;  
  for (auto ih= 0; ih < nbin ;++ih){
    std::string name1 = "absdEhisto" + to_string(ih);
    std::string name2 = "dEhisto" + to_string(ih);
    TH1F* h = new TH1F(name1.c_str(),name1.c_str(), 100, 0, 1.); h->Sumw2();
    TH1F* h2 = new TH1F(name2.c_str(),name2.c_str(), 100, -1, 1.); h2->Sumw2();
    histos.push_back(h);
    dhistos.push_back(h2);
  }

  //const double constant = -2.6913;//used maximum bin value
  //const double gradient = 7.66944;

  
  TH2D* histo = new TH2D("histo", "histo", 25, 0, 3000, 50, 0, 1);
  for (auto i = 0; i < tree->GetEntries(); ++i){
    tree->GetEntry(i);
    if (total_ring_PEs2 > 0.001){
      double recoE = f->Eval(total_ring_PEs2);
      //   double recoE = ((total_ring_PEs2*1.25*20000) - constant)/gradient; 
      histo->Fill(trueE,TMath::Abs(trueE-recoE)/trueE);
      int binNum = TMath::Min((int)(theLUT.second)->toBin(trueE),nbin);
      dhistos[binNum]->Fill((trueE-f->Eval(total_ring_PEs2))/trueE);
      histos[binNum]->Fill(TMath::Abs(trueE-recoE)/trueE);
    }
  }


  
  TProfile* prof = histo->ProfileX("_pfx", 1, -1 , "S");
  prof->Draw();

  
  TFile* output = new TFile("simpleLUT.root", "RECREATE");
  for (auto ih= 0; ih < nbin ;++ih){
    histos[ih]->Write();
    dhistos[ih]->Write();
  }
  prof->Write();
  output->Close();
}*/



void myreadData2(int nbin = 10, int nbinE = 15){

  gROOT -> ProcessLine( ".x ./mattStyle.C" );


  TFile* file = new TFile("/Users/edrakopo/work/work_energy_reco/ener_reco_scikit/nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root");//nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree.root");
  TTree* tree = (TTree*)file->Get("nu_eneNEW");

  // opening the file
  std::vector<LUT> luts;
  for (auto i =0; i < 3; ++i){
    LUT theLUT = makeLUT(tree,nbin,i);
    (theLUT.first)->SetLineColor(2);
    (theLUT.first)->SetLineStyle(i+1);
    luts.push_back(theLUT);
  }


  TCanvas* can = new TCanvas("can","can", 600, 400);
  luts[0].first->GetXaxis()->SetTitle("#PE^{ring}/20000");
  luts[0].first->GetYaxis()->SetTitle("E_{#mu} [MeV]");
  luts[0].first->GetYaxis()->SetTitleFont(132);
  luts[0].first->Draw();
  luts[1].first->Draw("SAME");
  luts[2].first->Draw("SAME");
  can->Print("mylut3.png");

  // LUT in energy versus photoelectrons to get the bins
  std::cout << "make my LUT " <<std::endl;
  MakeBins* bins;
  bins = divideDatabyE(tree,nbinE); 
  
/*  // setting the branches
  float total_ring_PEs2;
  tree1->SetBranchAddress("total_ring_PEs2",&total_ring_PEs2);

  float recoDWallR2;
  tree1->SetBranchAddress("recoDWallR2",&recoDWallR2);

  float recoDWallZ2;
  tree1->SetBranchAddress("recoDWallZ2",&recoDWallZ2);
  
  double trueE;
  tree1->SetBranchAddress("trueKE",&trueE); 
*/
  //set up histogram
  std::vector<TH1F*> histos;    std::vector<TH1F*> dhistos;  
  for (auto ih= 0; ih < nbinE ;++ih){
    std::string name1 = "absdEhisto" + to_string(ih);
    std::string name2 = "dEhisto" + to_string(ih);
    TH1F* h = new TH1F(name1.c_str(),name1.c_str(), 100, 0, 1.); h->Sumw2();
    TH1F* h2 = new TH1F(name2.c_str(),name2.c_str(), 100, -1, 1.); h2->Sumw2();
    histos.push_back(h);
    dhistos.push_back(h2);
  }

  //const double constant = -2.6913;//used maximum bin value
  //const double gradient = 7.66944;

  std::cout << histos.size() << " " << dhistos.size() << std::endl;  
  
  TH2D* histo = new TH2D("histo", "histo", 25, 0, 3000, 50, 0, 1);
 
  TFile* output = new TFile("mysimpleLUT3.root", "RECREATE");  

  double recoKE=0.; double trueKE1=0.; float recoDWall=0.; float recoToWall=0.;
  TTree *outTree= new TTree("outTree","output tree");
  outTree->Branch("trueKE1",&trueKE1,"trueKE1/D");
  outTree->Branch("recoKE",&recoKE,"recoKE/D");
  outTree->Branch("recoDWall", &recoDWall, "recoDWall/F");
  outTree->Branch("recoToWall", &recoToWall, "recoToWall/F");

  char filename1[1000];
  for (int i = 1040; i < 1100; i++) {
   if(i!=1051 && i!=1075){
   //TFile* file1 = new TFile("/Users/edrakopo/work/work_energy_reco/ener_reco_scikit/nu_numu_%i_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root");
    sprintf(filename1,"/Users/edrakopo/work/work_energy_reco/ener_reco_scikit/nu_numu_%d_CCQE_12in_energy_studies_recoquant_tree_NEWlookupsB_for_training.root",i);
  TFile *file1=new TFile(filename1,"READONLY");
  cout<<"File: "<<i<<endl;
  TTree* tree1 = (TTree*)file1->Get("nu_eneNEW");
   
   // setting the branches
   float total_ring_PEs2;
   tree1->SetBranchAddress("total_ring_PEs2",&total_ring_PEs2);
 
   float recoDWallR2;
   tree1->SetBranchAddress("recoDWallR2",&recoDWallR2);
 
   float recoDWallZ2;
   tree1->SetBranchAddress("recoDWallZ2",&recoDWallZ2);

   float recoDWall_2;
   tree1->SetBranchAddress("recoDWall_2",&recoDWall_2);
 
   float recoToWall_2;
   tree1->SetBranchAddress("recoToWall_2",&recoToWall_2);
 
   double trueE;
   tree1->SetBranchAddress("trueKE",&trueE);

  cout<<"Entries in evaluation tree: "<<tree1->GetEntries()<<endl;
  for (auto i = 0; i < tree1->GetEntries(); ++i){
  //for (auto i = 0; i < 10; ++i){
    tree1->GetEntry(i);
   // if (total_ring_PEs2 > 0.001){
      //int bin = toRZbin(recoDWallR2,recoDWallZ2);
      int bin = toRZbin(recoDWall_2,recoToWall_2);
      TF1* f = luts[bin].first;
      double recoE = f->Eval(total_ring_PEs2);
      //   double recoE = ((total_ring_PEs2*1.25*20000) - constant)/gradient; 
      histo->Fill(trueE,TMath::Abs(trueE-recoE)/trueE);
      int binNum = TMath::Min((int)bins->toBin(trueE),nbinE);
      dhistos[binNum]->Fill((trueE-recoE)/trueE);
      histos[binNum]->Fill(TMath::Abs(trueE-recoE)/trueE);
      //cout<<"----- recoKE: "<<recoKE<<","<<trueKE1<<"|"<<recoDWall<<","<<recoToWall<<endl;
  /*    if(i<5){
        cout<<"recoE: "<<recoE<<" trueE: "<<trueE<<" recoDWallR2: "<<recoDWallR2<<" recoDWallZ2: "<<recoDWallZ2<<endl;
        cout<<"recoDWall_2: "<<recoDWall_2<<" recoToWall_2: "<<recoToWall_2<<endl;
      }*/
      recoKE=recoE;
      trueKE1=trueE;
      recoDWall=recoDWall_2;
      recoToWall=recoToWall_2;

      outTree->Fill();
   // }
   }
   }//if
  }

  TProfile* prof = histo->ProfileX("_pfx", 1, -1 , "S");
  prof->Draw();

  //TFile* output = new TFile("simpleLUT3.root", "RECREATE");
  for (auto ih= 0; ih < nbinE ;++ih){
    histos[ih]->Write();
    dhistos[ih]->Write();
  }

  output->cd();
  prof->Write();
  outTree->Write();

  output->Close();
}



/*
void readStandard(int nbinE = 15){

  gROOT -> ProcessLine( ".x ./mattStyle.C" );
  
  TFile* file = new TFile("nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree.root");
  TTree* tree = (TTree*)file->Get("nu_eneNEW");

  StandardTable* table = new StandardTable("energyLookups_40_new.root");

  std::cout << "make my LUT " <<std::endl;
  MakeBins* bins;
  bins = divideDatabyE(tree,nbinE); 
  
  
  // setting the branches
  float total_ring_PEs2;
  tree->SetBranchAddress("total_ring_PEs2",&total_ring_PEs2);

  float recoDWallR2;
  tree->SetBranchAddress("recoDWallR2",&recoDWallR2);

  float recoDWallZ2;
  tree->SetBranchAddress("recoDWallZ2",&recoDWallZ2);
  
  double trueE;
  tree->SetBranchAddress("trueKE",&trueE); 

  // LUT in energy versus photoelectrons to get the bins
  
  
  //set up histogram
  std::vector<TH1F*> histos;    std::vector<TH1F*> dhistos;  
  for (auto ih= 0; ih < nbinE ;++ih){
    std::string name1 = "absdEhisto" + to_string(ih);
    std::string name2 = "dEhisto" + to_string(ih);
    TH1F* h = new TH1F(name1.c_str(),name1.c_str(), 100, 0, 1.); h->Sumw2();
    TH1F* h2 = new TH1F(name2.c_str(),name2.c_str(), 100, -1, 1.); h2->Sumw2();
    histos.push_back(h);
    dhistos.push_back(h2);
  }

  //const double constant = -2.6913;//used maximum bin value
  //const double gradient = 7.66944;

  TFile* output = new TFile("standardLUT.root", "RECREATE");
  TTree* otree = new TTree("tree", "RECREATE");
  double rwall; TBranch* branch_r =  otree->Branch("rwall",&rwall, "rwall/D");
  double dwall; TBranch* branch_d =  otree->Branch("dwall",&dwall, "dwall/D");
  double de; TBranch* branch_de =  otree->Branch("de",&de, "e/D");
  double etrue; TBranch* branch_etrue =  otree->Branch("etrue",&etrue, "etrue/D");
  int thebin; TBranch* branch_bin =  otree->Branch("bin",&thebin, "thebin/I");
  
  TH2D* histo = new TH2D("histo", "histo", 25, 0, 3000, 50, 0, 1);
  for (auto i = 0; i < tree->GetEntries(); ++i){
    tree->GetEntry(i);
    //    if (total_ring_PEs2 > 0.001){
      
      double recoE = table->eval(550*recoDWallR2,1100*recoDWallZ2,20000*total_ring_PEs2,true);
      //   double recoE = ((total_ring_PEs2*1.25*20000) - constant)/gradient; 
      histo->Fill(trueE,TMath::Abs(trueE-recoE)/trueE);
      int binNum = TMath::Min((int)bins->toBin(trueE),nbinE);
      dhistos[binNum]->Fill((trueE-recoE)/trueE);
      histos[binNum]->Fill(TMath::Abs(trueE-recoE)/trueE);
      rwall = 550*recoDWallR2; dwall = 1100*recoDWallZ2; etrue = trueE; de = (trueE-recoE)/trueE;
      thebin = binNum;
      otree->Fill();
      //}
  }

  TProfile* prof = histo->ProfileX("_pfx", 1, -1 , "S");
  prof->Draw();

  output->cd();
  otree->Write();
  for (auto ih= 0; ih < nbinE ;++ih){
    histos[ih]->Write();
    dhistos[ih]->Write();
  }
  prof->Write();
 
  output->Close();
}



void addBranch(std::string fileName = "nu_numu_1000_1039_CCQE_12in_energy_studies_recoquant_tree.root", std::string trailer = "_withLUT.root"){

  gROOT -> ProcessLine( ".x ./mattStyle.C" );
  
  TFile* file = new TFile(fileName.c_str());
  TTree* tree = (TTree*)file->Get("nu_eneNEW");

  StandardTable* table = new StandardTable("energyLookups_40_new.root");

  std::cout << "make my LUT " <<std::endl;
  LUT theLUT = makeLUT(tree,15);
  TF1* f = theLUT.first;

   // make the output
  std::string outputName = fileName.substr(0,fileName.size() - 5);
  outputName += trailer;
  TFile* outFile  =new TFile(outputName.c_str(),"RECREATE");

  // clone the tree..
  TTree*  newtree = tree->CloneTree(-1);

  // setting the branches
  float total_ring_PEs2;
  newtree->SetBranchAddress("total_ring_PEs2",&total_ring_PEs2);

  float recoDWallR2;
  newtree->SetBranchAddress("recoDWallR2",&recoDWallR2);

  float recoDWallZ2;
  newtree->SetBranchAddress("recoDWallZ2",&recoDWallZ2);
  
  double trueE;
  newtree->SetBranchAddress("trueKE",&trueE); 

  float standardlutE; TBranch* branch_luts = newtree->Branch("standardHighE",&standardlutE, "standardHighE");
  float lutE; TBranch* branch_lutE = newtree->Branch("lutE",&lutE, "lutE");
  
  
  // const double constant = -2.6913;//used maximum bin value
  //const double gradient = 7.66944;
  
  for (auto i = 0; i < tree->GetEntries(); ++i){
    newtree->GetEntry(i);  
    standardlutE  = table->eval(550*recoDWallR2,1100*recoDWallZ2,20000*total_ring_PEs2,true);
    //lowE = (nphot2*1.25 - constant)/gradient;
    lutE = f->Eval(total_ring_PEs2);
    branch_luts->Fill();
    branch_lutE->Fill();
      //}
  }

  newtree->Write();
  outFile->Close();

}*/

