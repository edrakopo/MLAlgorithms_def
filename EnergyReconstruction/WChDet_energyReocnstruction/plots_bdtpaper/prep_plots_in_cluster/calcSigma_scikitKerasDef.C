#include "TSystem.h"
#include "TProfile.h"
#include "TFile.h"
#include "TROOT.h"
#include "MakeBins.C"
//#include "FitDCBS.C"
#include "TTree.h"
#include "TBranch.h"
#include "TObject.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH1F.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <TString.h>

void loadAtStart(){
gROOT -> ProcessLine( ".x ./mattStyle.C" );
}

typedef struct {

double value;
double error;

} ValueWithError;

typedef struct {


ValueWithError mean;
ValueWithError sigmaEff;

} Result;

Result* fitBin(TH1D* histo, RooPlot* &outRoo, TH1D* &outHist){
Result *res = new Result();

FitDCBS(histo,res->sigmaEff.value, res->sigmaEff.error ,outRoo, outHist, res->mean.value, res->mean.error);


// Create a macro to fit the histograms and output a Result* 
//
// fit.C or something and remember to have a gRoot processline to run MakeBins.C and FitHisto.C

return res;

}

Result* fitBin(TH1D* histo){
Result *res = new Result();


FitDCBS(histo,res->sigmaEff.value, res->sigmaEff.error, res->mean.value, res->mean.error);
//FitSigma(histo,res->sigmaEff.value, res->sigmaEff.error );


// Create a macro to fit the histograms and output a Result* 
//
// fit.C or something and remember to have a gRoot processline to run MakeBins.C and FitHisto.C

return res;

}

void OutResult(std::string fName,TString methName , Result * a){

ofstream myfile;
myfile.open(fName, std::ios::app);

std::cout.setf(ios::fixed, ios::floatfield); // set fixed floating format
std::cout.precision(5);
myfile<<methName<<"\t \t";
myfile<<"Mean: \t";
myfile.width(10);
myfile <<std::internal<<a->mean.value;
myfile<<"\tError: \t";
std::cout.width(10);
myfile<<std::internal<<a->mean.error;
myfile<<"\tSigma: \t";
std::cout.width(10);
myfile<<std::internal<<a->sigmaEff.value;
myfile<<"\tSigmaError: \t";
std::cout.width(10);
myfile<<std::internal<<a->sigmaEff.error<<"\n";
myfile.close();

}



void OutResult(std::string fName,TString methName , Result * a, double percentValD, double percentValU,double percentVal){

ofstream myfile;
myfile.open(fName, std::ios::app);

std::cout.setf(ios::fixed, ios::floatfield); // set fixed floating format
std::cout.precision(5);
myfile<<methName<<"\t \t";
myfile<<"Mean: \t";
myfile.width(10);
myfile <<std::internal<<a->mean.value;
myfile<<"\tSigma: \t";
std::cout.width(10);
myfile<<std::internal<<a->mean.error;
myfile<<"\tLower Percentage: \t";
std::cout.width(10);
myfile<<std::internal<<percentValD;
myfile<<"\tUpper Percentage: \t";
std::cout.width(10);
myfile<<std::internal<<percentValU;
myfile<<"\tPercentage in Sigma: \t";
std::cout.width(10);
myfile<<std::internal<<percentVal<<"\n";
myfile.close();

}




void PrintResult(Result * a){

std::cout.setf(ios::fixed, ios::floatfield); // set fixed floating format
std::cout.precision(5); 

std::cout<<"Mean: \t";
std::cout.width(10);
std::cout <<std::internal<<a->mean.value;
std::cout<<"\tError: \t";
std::cout.width(10);
std::cout<<std::internal<<a->mean.error;
std::cout<<"\tSigma: \t";
std::cout.width(10);
std::cout<<std::internal<<a->sigmaEff.value;
std::cout<<"\tSigmaError: \t";
std::cout.width(10);
std::cout<<std::internal<<a->sigmaEff.error<<std::endl;


}


MakeBins* divideDatabyE(TTree* tree, int nbin){

//split up the data into bins
double trueE;
tree->SetBranchAddress("trueKE",&trueE);
   
std::vector<double> energies; energies.resize(tree->GetEntries());

             for (auto i = 0; i < tree->GetEntries(); ++i){
             tree->GetEntry(i);
             energies[i] = trueE;
             } 
   
MakeBins* bins = new MakeBins(energies,nbin);
return bins;
}
   



//void RunEnergy(int nbin, int minEnergy, int maxEnergy, std::string ineffVal)
void RunEnergy(int nbin, int minEnergy, int maxEnergy)
{

loadAtStart();

std::vector<TString> methods;
/*
methods.push_back("BDTG_10000_Trees_3_Depth");
methods.push_back("BDTG_5000_Trees_3_Depth");
methods.push_back("MLP_6_HiddenLayers_500_cycles");
methods.push_back("BDTG_10000_Trees_4_Depth");
methods.push_back("BDTG_8000_Trees_3_Depth");
methods.push_back("MLP_6_HiddenLayers");
methods.push_back("BDTG_10000_Trees_5_Depth");
methods.push_back("BDTG_9000_Trees_3_Depth");
methods.push_back("MLP_7_HiddenLayers");
methods.push_back("BDTG_1000_Trees_3_Depth");
methods.push_back("MLP_4_HiddenLayers");
methods.push_back("MLP_8_HiddenLayers");
methods.push_back("BDTG_11000_Trees_3_Depth");
methods.push_back("MLP_5_HiddenLayers");
methods.push_back("BDTG_12000_Trees_3_Depth");
methods.push_back("MLP_6_HiddenLayers_1500_cycles");
methods.push_back("BDTG_10000_trees_10_maxdepth");
methods.push_back("BDTG_1000_trees_10_maxdepth");
*/
methods.push_back("BDTG_10000_Trees_3_Depth");
methods.push_back("BDTG_5000_Trees_3_Depth");
methods.push_back("BDTG_10000_Trees_4_Depth");
methods.push_back("BDTG_8000_Trees_3_Depth");
methods.push_back("BDTG_10000_Trees_5_Depth");
methods.push_back("BDTG_9000_Trees_3_Depth");
methods.push_back("BDTG_1000_Trees_3_Depth");
methods.push_back("BDTG_11000_Trees_3_Depth");
methods.push_back("BDTG_12000_Trees_3_Depth");
methods.push_back("LUT");
//methods.push_back("BDTG_10000_trees_10_maxdepth");
//methods.push_back("BDTG_1000_trees_10_maxdepth");
std::string prefix = "/home/edrakopo/work/ener_reco_scikit/keras_reco_TrainTestSample.root";//+ineffVal;
for (int ineff = 6; ineff<7; ineff++){
TString fileName = prefix + methods[ineff] + ".root";
TFile *file = new TFile(fileName);
int totalNum =0;
int numInSig =0;
TTree *tree = (TTree*)file->Get("data");

MakeBins *binnedData = divideDatabyE(tree,nbin);

double trueE;
float val;

tree->SetBranchAddress("trueKE", &trueE);
tree->SetBranchAddress("Val", &val);

TH2D* hist = new TH2D("ReconstructedEnergy","" ,nbin,(float)maxEnergy,(float)minEnergy,200. , -100. , 100.);
	      hist->GetYaxis()->SetTitle("#frac{#Delta E}{E} [%]");
              hist->GetXaxis()->SetTitle("E [MeV]");

TProfile *mlp;
//mlp->GetYaxis()->SetTitle("Events");
//TH2D* hist = new TH2D("ReconstructedEnergy" ,"ReconstructedEnergy" ,nbin,0.,5000.,200. , -100. , 100.);
//std::vector<TH2D*> histos;

double binCenterVal[nbin];

std::vector<TH1D*> projectedhistos;
std::vector<double>vals;

for (auto i = 0; i < nbin; i++){

std::string hname = "fitBin_new" + std::to_string(i);
TH1D* histo = new TH1D(hname.c_str(), hname.c_str(),200,-100 ,100);
binCenterVal[i] = binnedData->binCenter(i);
projectedhistos.push_back(histo);


}


std::vector<double> newVals;
//create 1D histograms to be filled for each bin
//



//fill 2D histogram and fill the 1D histograms for each bin

	      for (auto i =0; i<tree->GetEntries(); i++){
	      tree->GetEntry(i);

//	      hist->Fill(trueE, (100*(trueE-val)/trueE));
if (trueE<=maxEnergy && trueE >= minEnergy){ hist->Fill(trueE, (100*(trueE-val)/trueE));
int binNum = TMath::Min((int)binnedData->toBin(trueE) ,nbin);
totalNum++; vals.push_back((100*(trueE-val)/trueE));
projectedhistos[binNum]->Fill((100*(trueE-val)/trueE));


//Counts the number of data points in the energy range.
}


	     //	      int binNum = TMath::Min((int)binnedData->toBin(trueE),nbin);

//	      projectedhistos[binNum]->Fill((100*(trueE-val)/trueE));


	      }


//	     for(int i =0; i<nbin; i++){

//	     projectedhistos[nbin]=hist->ProjectionY("Projected Y Graph",nbin,nbin);



//	     }



mlp = hist->ProfileX();

	      for (int i = 0; i<nbin; i++){
	      
   	      std::string name1 = "histBin_"+to_string(i);
	      TH1D *h = hist->ProjectionY(name1.c_str(),i,i);
	      h->GetXaxis()->SetTitle("#frac{#Delta E}{E} [%]");
              h->GetYaxis()->SetTitle("Events");
	      h->Sumw2();
//	      projectedhistos.push_back(h);


	      }
 



std::vector<Result*> results;
std::vector<RooPlot*> rooplots ;
std::vector<TH1D*> sigmahists;
TString outnameFile = prefix + "OUT" + methods[ineff] + ".root"; 
TFile* output = new TFile(outnameFile, "RECREATE");

RooPlot* rooplot;
TH1D* tmpHist; 

for (auto ih= 0; ih < nbin ;++ih){
Result *r = new Result;

if (ineff !=9){r = fitBin(projectedhistos[ih], rooplot, tmpHist);}
else {r->mean.value = projectedhistos[ih]->GetMean();
r->mean.error = 0;
r->sigmaEff.value = 0;
r->sigmaEff.error = 0;


}

results.push_back(r);

std::string nameRoo = "FitBin_"+to_string(ih);
if (ineff !=9){
rooplot->SetName(nameRoo.c_str());
rooplot->Write();
std::string nameHist = "SigmaHistBin_"+to_string(ih);
tmpHist->SetName(nameHist.c_str());
tmpHist->GetXaxis()->SetTitle("#sigma (#frac{#Delta E}{E}) [%]");
tmpHist->GetYaxis()->SetTitle("Number of Candidates");
tmpHist->Write();
}
projectedhistos[ih]->Write();

//if(ih == interestBin){ Binterest.push_back(r);}
    
  }
//work out percentage of results in sigma
double sigD =results[0]->mean.value - results[0]->mean.error;
double sigU =results[0]->mean.value + results[0]->mean.error;
int topPC =0;
int botPC =0; 
int totPC = 0;
	      for (auto i =0; i<vals.size(); i++){
	      
if (vals[i]>= sigD && vals[i]<results[1]->mean.value){ botPC ++;totPC++;

}else if (vals[i]<= sigU && vals[i]>=results[1]->mean.value){topPC++;totPC++;}


	      }


double percentValD = (double) botPC/totalNum;
double percentValU = (double) topPC/totalNum;
double percentVal = (double) totPC/totalNum;

//To work out the 68% confidence level

std::sort(vals.begin(),vals.end());

int confidenceLevel = (int) (vals.size()*0.34);
double upperBound;
double lowerBound;
int numberOfMean =0;
std::vector<double> means ; 

bool lessThanMean = true;
int backupMeansVal;
for(int i = 0; i<vals.size(); i++){
if (vals[i] == results[0]->mean.value){means.push_back(i);}
if(lessThanMean && vals[i] > results[0]->mean.value){lessThanMean = false; 

//nested if statement
if(abs(vals[i] - results[0]->mean.value) > abs(vals[i-1] - results[0]->mean.value)){backupMeansVal=i-1;}else{
backupMeansVal = i; }


}
if (vals[i] == results[0]->mean.value){numberOfMean++;}
//std::cout<<"vector member: "<<vals[i]<<std::endl;
}

if (!lessThanMean && means.size()==0){means.push_back(backupMeansVal);}
//std::cout<<"number where mean is: "<<means[(int) means.size()/2] <<"and number of confidence level: "<<confidenceLevel<<std::endl;
upperBound = vals[means[(int) means.size()/2] + confidenceLevel];
lowerBound = vals[means[(int) means.size()/2] - confidenceLevel];




for (auto i =0;i<nbin; i++){

//std::cout<< "Bin "; 
//std::cout.width(4);
//std::cout<<std::left <<to_string(i);
//PrintResult(results[i]);
//OutResult(prefix +"EnergyResPercentInSigma.txt",methods[ineff] ,results[i], percentValD, percentValU, percentVal);
OutResult("EnergyVsIneff.txt",methods[ineff]+ineffVal ,results[i], lowerBound, upperBound, numberOfMean);

}
hist->Write();
mlp->Write();
output->Close();

projectedhistos.clear();
results.clear();
rooplots.clear();
sigmahists.clear();


}
//std::string binofInterest = "Bin"+to_string(interestBin)+"of"+to_string(nbin)+".root";
//TFile* outFile = new TFile(binofInterest.c_str(), "RECREATE");
//double xval[8] =  {0,1,2,3,4,5,6,7};
//double meanval[8], emeanval[8], sigmaval[8], esigmaval[8];
//for (int i= 0; i<Binterest.size(); i++){
//meanval[i] = Binterest[i]->mean.value;
//emeanval[i] = Binterest[i]->mean.error;
//sigmaval[i] = Binterest[i]->sigmaEff.value;
//esigmaval[i] = Binterest[i]->sigmaEff.error;


//}



}

void RunForAll(int nbin, int minEnergy, int maxEnergy){

std::vector<std::string> ineff;
ineff.push_back("QE22/");
ineff.push_back("QE21.78/");
ineff.push_back("QE21.56/");
ineff.push_back("QE21.34/");
ineff.push_back("QE21.12/");
ineff.push_back("QE20.9/");
ineff.push_back("QE20.68/");
ineff.push_back("QE20.46/");

for (int i=0; i<8; i++){

RunEnergy(nbin, minEnergy,maxEnergy, ineff[i]);

}
}

