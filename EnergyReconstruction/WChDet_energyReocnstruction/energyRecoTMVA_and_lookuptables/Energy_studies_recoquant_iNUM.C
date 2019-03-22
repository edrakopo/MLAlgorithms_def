#include "TRandom3.h"
#include <iomanip>
#include <iostream>
#include "TMath.h"
#include "TFile.h"
#include "TH2D.h"
#include "TGraph2D.h"
#include "TH3D.h"
#include "TSpectrum2.h"
#include "../include/WCLTreeReader.hh"
#include "../include/HighEReco.hh"
#include "../include/LowEReco.hh"

using namespace std;

// Performs reconstruction assuming single ring high-E event
//
// filename: root output file from $WCHSANDBOX/analysis/advanced_scripts/SandboxDetSim.C
// genfilename: generatorcardfile.root file from WChSandbox
// outfilename: file to write output
// lnLfilename: file with likelihood tables
//

int* FindClusters(int nHits, double * hitTimes, int * hitCluster, int & nClusters){
    //cout << "Starting cluster search" << endl;
    nClusters = 0;
    int *hitOrder = new int[nHits];
    TMath::Sort(nHits, hitTimes, hitOrder, kFALSE);
    //cout << "Sorted " << nHits << " hits" << endl;
    //Loop over hits
    for(int clusterStart=0; clusterStart+9<nHits; clusterStart++){
        int clusterEnd = clusterStart+9;
        double clusterStartTime = hitTimes[hitOrder[clusterStart]];
        //Check at least 10 hits in 200ns window
        if(hitTimes[hitOrder[clusterEnd]]-clusterStartTime>200)
            continue;
        //Find cluster end time, with hit spacing of less than 10ns, allowing for at most two spacings of 10-20ns
        int nSpacings = 0;
        double clusterEndTime = clusterStartTime;
        for(clusterEnd=clusterStart; clusterEnd<nHits-1; clusterEnd++)
        {
            clusterEndTime = hitTimes[hitOrder[clusterEnd]];
            double nextHitTime = hitTimes[hitOrder[clusterEnd+1]];
            if(nextHitTime - clusterEndTime < 10)
                continue;
            if (nextHitTime - clusterEndTime < 20 && nSpacings < 2) {
                nSpacings++;
                continue;
            }
            break;
        }
        //Check at least 10 hits in cluster
        if(clusterEnd-clusterStart<9)
            continue;
        //Found a cluster
        nClusters++;
        //cout << "Cluster " << nClusters << " found with " << clusterEnd-clusterStart+1 << "hits: "
         //    << clusterStartTime << "ns to " << clusterEndTime << "ns" << endl;
        //Set cluster number for hits in cluster
        for(int iHit = clusterStart; iHit<=clusterEnd; iHit++){
            hitCluster[hitOrder[iHit]] = nClusters;
        }
        //Look for next cluster
        clusterStart = clusterEnd;
    }
    //Count hits in each cluster
    int * clusterHitCounts = new int[nClusters+1](); for(int i=0; i<nClusters+1; i++) clusterHitCounts[i]=0;
    for(int iHit=0; iHit<nHits; iHit++) {
        clusterHitCounts[hitCluster[iHit]]++;
    }
    delete[] hitOrder;
    return clusterHitCounts;
}

//____________________________________
//-------------------------------------------------------------------------------------------------------------------------
std::vector<double> finding_solutions_for_qudratic_eq( double a , double b , double c){
  std::vector<double> resault;resault.clear();
  double discr = b*b -4*a*c;
  if(discr < 0){
    cout<<"ALARM: discr < 0"<<endl;
    return resault;
  }
  if(discr == 0){
    cout<<"ALARM: discr = 0"<<endl;
    double s1 = -b/(2*a);
    resault.push_back(s1);
    return resault;
  }
  if(discr > 0){
    double s1 = (-b + sqrt(discr) )/(2*a);
    resault.push_back(s1);
    double s2 = (-b - sqrt(discr) )/(2*a);
    resault.push_back(s2);
    cout<<"resault.size() = "<<resault.size()<<endl;
    cout<<"sol 1 = "<<resault[0]<<" and sol 2 = "<<resault[1]<<endl;
    return resault;
  }
}

//____________________________________
//-------------------------------------------------------------------------------------------------------------------------
double CalculatingIntersectionPoint( std::vector<double> sol, double lambda_max )
{
  bool solution_found =0;
  double time = 0; 
	 
	 if(sol.size() == 0)
	   cout<<"ALARM if(sol.size() == 0)"<<endl;
  
	 if(sol.size() > 2)
	   cout<<"ALARM sol.size() > 2"<<endl;
  
	 if(sol.size() == 1){
	   if(sol[0] < 0 )
	     cout<<"ALARM if(sol.size() == 1) && if(sol[0] < 0 )"<<endl;
	   else
	     time = sol[0];
	 }
  
	 if(sol.size() == 2){
	   if(sol[0] >0 &&  sol[1] > 0){
	     cout<<"ALARM if(sol.size() == 2) && if(sol[0] >0 &&  sol[1] > 0)"<<endl;
	     solution_found=1;
	       if(sol[0] > lambda_max)
		 time = sol[0];
	       if(sol[1] > lambda_max)
		 time = sol[1];
	   }
	   else{
	     if(sol[0] <0 &&  sol[1] <0){
	       cout<<"ALARM if(sol.size() == 2) && if(sol[0] <0 &&  sol[1] < 0)"<<endl;
	      solution_found=1;
	       if(sol[0] > lambda_max)
		 time = sol[0];
	       if(sol[1] > lambda_max)
		 time = sol[1];
	     }
	     else{ solution_found=1;
	       if(sol[0] > 0)
		 time = sol[0];
	       if(sol[1] > 0)
		 time = sol[1];
	     }}}

	 cout<<"time= "<<time<<endl;
	 return time; 
}

//____________________________________
//-------------------------------------------------------------------------------------------------------------------------
bool FindSolution( std::vector<double> sol,  double lambda_max )
{
  bool solution_found =0;
  double time = 0; 
	 
	 if(sol.size() == 0)
	   cout<<"ALARM if(sol.size() == 0)"<<endl;
  
	 if(sol.size() > 2)
	   cout<<"ALARM sol.size() > 2"<<endl;
  
	 if(sol.size() == 1){
	   if(sol[0] < 0 )
	     cout<<"ALARM if(sol.size() == 1) && if(sol[0] < 0 )"<<endl;
	   else
	     time = sol[0];
	 }
  
	 if(sol.size() == 2){
	   if(sol[0] >0 &&  sol[1] > 0){
	     cout<<"ALARM if(sol.size() == 2) && if(sol[0] >0 &&  sol[1] > 0)"<<endl;
	     solution_found=1;
	       if(sol[0] > lambda_max)
		 time = sol[0];
	       if(sol[1] > lambda_max)
		 time = sol[1];
	   }
	   else{
	     if(sol[0] <0 &&  sol[1] <0){
	       cout<<"ALARM if(sol.size() == 2) && if(sol[0] <0 &&  sol[1] < 0)"<<endl;
	      solution_found=1;
	       if(sol[0] > lambda_max)
		 time = sol[0];
	       if(sol[1] > lambda_max)
		 time = sol[1];
	     }
	     else{ solution_found=1;
	       if(sol[0] > 0)
		 time = sol[0];
	       if(sol[1] > 0)
		 time = sol[1];
	     }}}

	  cout<<"solution_found= "<<solution_found<<endl;
	  return solution_found; 
}

//____________________________________
//-----------------------------------------------------------------------------//
// calculating distance from om (WITH HIT) center to muon track - 2ways
//-----------------------------------------------------------------------------//
double find_vert_dist(double xmu_rec,double ymu_rec,double zmu_rec,double xrecDir,double yrecDir,double zrecDir,double x_pmtpos,double y_pmtpos,double z_pmtpos)
{
  double x_dist = (xmu_rec - x_pmtpos);  //reco vertex - hit position vector
  double y_dist = (ymu_rec - y_pmtpos);
  double z_dist = (zmu_rec - z_pmtpos);
  
  //1st way:
  double cross_x_pulse = y_dist*zrecDir - z_dist*yrecDir;
  double cross_y_pulse = z_dist*xrecDir - x_dist*zrecDir;
  double cross_z_pulse = x_dist*yrecDir - y_dist*xrecDir;
  double magnitude = sqrt(cross_x_pulse*cross_x_pulse + cross_y_pulse*cross_y_pulse + cross_z_pulse*cross_z_pulse);   
  
  //2nd way:
  double dot_pulse = x_dist*xrecDir + y_dist*yrecDir + z_dist*zrecDir;
  double dist_dot = sqrt((x_dist*x_dist + y_dist*y_dist + z_dist*z_dist)-dot_pulse*dot_pulse);

  return magnitude;
}

//____________________________________
//---------------------------------------------------------------------------------------------------------//
// Find Distance between Muon's Vertex and Photon's Production Point (lambda) 
//---------------------------------------------------------------------------------------------------------//
double find_lambda(double xmu_rec,double ymu_rec,double zmu_rec,double xrecDir,double yrecDir,double zrecDir,double x_pmtpos,double y_pmtpos,double z_pmtpos,double theta_cher)
{
     double lambda1 = 0.0;    double lambda2 = 0.0;    double length = 0.0 ;
     double xmupos_t1 = 0.0;  double ymupos_t1 = 0.0;  double zmupos_t1 = 0.0;  
     double xmupos_t2 = 0.0;  double ymupos_t2 = 0.0;  double zmupos_t2 = 0.0; 
     double xmupos_t = 0.0;   double ymupos_t = 0.0;   double zmupos_t = 0.0;  
    
     double theta_muDir_track = 0.0;
     double theta_muDir_track1 = 0.0;  double theta_muDir_track2 = 0.0;
     double cos_thetacher = cos(theta_cher*TMath::DegToRad());
     double xmupos_tmin = 0.0; double ymupos_tmin = 0.0; double zmupos_tmin = 0.0; 
     double xmupos_tmax = 0.0; double ymupos_tmax = 0.0; double zmupos_tmax = 0.0;
     double lambda_min = 10000000;  double lambda_max = -99999999.9;  double lambda = 0.0; 
   
   //---orizw syntelestes ths deyterova8mias alpha, beta, gamma - me agnwsto to lamda(opoy xmupos_t1 = xmu_rec + xrecDir*lambda1;)
     double alpha = (xrecDir*xrecDir + yrecDir*yrecDir + zrecDir*zrecDir) * ( (xrecDir*xrecDir + yrecDir*yrecDir + zrecDir*zrecDir) - (cos_thetacher*cos_thetacher) );
     double beta = ( (-2)*(xrecDir*(x_pmtpos - xmu_rec) + yrecDir*(y_pmtpos - ymu_rec) + zrecDir*(z_pmtpos - zmu_rec) )*((xrecDir*xrecDir + yrecDir*yrecDir + zrecDir*zrecDir) - (cos_thetacher*cos_thetacher)) );
     double gamma = ( ( (xrecDir*(x_pmtpos - xmu_rec) + yrecDir*(y_pmtpos - ymu_rec) + zrecDir*(z_pmtpos - zmu_rec))*(xrecDir*(x_pmtpos - xmu_rec) + yrecDir*(y_pmtpos - ymu_rec) + zrecDir*(z_pmtpos - zmu_rec)) ) - (((x_pmtpos - xmu_rec)*(x_pmtpos - xmu_rec) + (y_pmtpos - ymu_rec)*(y_pmtpos - ymu_rec) + (z_pmtpos - zmu_rec)*(z_pmtpos - zmu_rec))*(cos_thetacher*cos_thetacher)) );
     
     
     double discriminant = ( (beta*beta) - (4*alpha*gamma) );
     
     lambda1 = ( (-beta + sqrt(discriminant))/(2*alpha));
     lambda2 = ( (-beta - sqrt(discriminant))/(2*alpha));
 
     xmupos_t1 = xmu_rec + xrecDir*lambda1;  xmupos_t2 = xmu_rec + xrecDir*lambda2;
     ymupos_t1 = ymu_rec + yrecDir*lambda1;  ymupos_t2 = ymu_rec + yrecDir*lambda2;
     zmupos_t1 = zmu_rec + zrecDir*lambda1;  zmupos_t2 = zmu_rec + zrecDir*lambda2;
     
     double tr1 = sqrt((x_pmtpos - xmupos_t1)*(x_pmtpos - xmupos_t1) + (y_pmtpos - ymupos_t1)*(y_pmtpos - ymupos_t1) + (z_pmtpos - zmupos_t1)*(z_pmtpos - zmupos_t1));      
     double tr2 = sqrt((x_pmtpos - xmupos_t2)*(x_pmtpos - xmupos_t2) + (y_pmtpos - ymupos_t2)*(y_pmtpos - ymupos_t2) + (z_pmtpos - zmupos_t2)*(z_pmtpos - zmupos_t2)); 
     theta_muDir_track1 = (acos( (xrecDir*(x_pmtpos - xmupos_t1) + yrecDir*(y_pmtpos - ymupos_t1) + zrecDir*(z_pmtpos - zmupos_t1))/(tr1) )*TMath::RadToDeg());
     theta_muDir_track2 = (acos( (xrecDir*(x_pmtpos - xmupos_t2) + yrecDir*(y_pmtpos - ymupos_t2) + zrecDir*(z_pmtpos - zmupos_t2))/(tr2) )*TMath::RadToDeg()); 
     
     
     //---------------------------- choose lambda!! ---------------------------------
     if( theta_muDir_track1 < theta_muDir_track2 ){
       lambda = lambda1;
       xmupos_t = xmupos_t1; 
       ymupos_t = ymupos_t1;
       zmupos_t = zmupos_t1;
       theta_muDir_track=theta_muDir_track1;
     }else if( theta_muDir_track2 < theta_muDir_track1 ){
       lambda = lambda2;
       xmupos_t = xmupos_t2; 
       ymupos_t = ymupos_t2;
       zmupos_t = zmupos_t2;
       theta_muDir_track=theta_muDir_track2;
     }
     //cout<<"lambda= "<<lambda<<endl;
     // cout<<"xmupos_t = "<<(xmu_rec + xrecDir*lambda)<<endl; 
     // cout<<"ymupos_t = "<<(ymu_rec + yrecDir*lambda)<<endl; 
     // cout<<"zmupos_t = "<<(zmu_rec + zrecDir*lambda)<<endl;  
    
     return lambda;
}
//_________________________________________________________________________

////////////////////////
void Energy_studies_recoquant_iNUM(TString filename = "nu_numu_iNUM_out_12in_new_fortraining.root",//photonmask output
				   TString genfilename = "generatorcardfile.root",
           		           TString outfilename = "nu_numu_iNUM_CCQE_12in_energy_studies_recoquant_tree.root",
				   TString lnLfilename = "likelihood_tables.root",
				   TString lookupfilename = "energyLookups_40_new.root",
				   int startEvent=0, int maxEvents=0)

{

  TH2F * hMuEn=new TH2F("hMuEn","hMuEn;trueKE;Total pes",60, 0., 3000., 500, 0., 5.);
  TH2F * hMuEn2=new TH2F("hMuEn2","hMuEn2;trueKE;Total pes",40, 0., 2000., 100, 0., 100000.);
  TH2F * hMuEn_HITS=new TH2F("hMuEn_HITS","hMuEn_HITS;trueKE;Total hits",60, 0., 3000., 500, 0., 5.);
  TH2F * hMuEn2_HITS=new TH2F("hMuEn2_HITS","hMuEn2_HITS;trueKE;Total hits",40, 0., 2000., 400, 0., 40000.);
  TH2F * hMuEn2_trueToWall=new TH2F("hMuEn2_trueToWall","hMuEn2_trueToWall;trueKE;Total pes/trueToWall",40, 0., 2000., 100, 0., 500.);
  TH2F * hMuEn2_short_pot_length=new TH2F("hMuEn2_short_pot_length","hMuEn2_short_pot_length;trueKE;Total pes",40, 0., 2000., 100, 0., 100000.);
  TH2F * hMuEn2_large_pot_length=new TH2F("hMuEn2_large_pot_length","hMuEn2_large_pot_length;trueKE;Total pes",40, 0., 2000., 100, 0., 100000.);
  TH2F * hNeuEn2=new TH2F("hNeuEn2","hNeuEn2;E_{#nu};Total pes",40, 0., 2000., 100, 0., 100000.);
  TH2F * hNeuEn2_short_pot_length=new TH2F("hNeuEn2_short_pot_length","hNeuEn2_short_pot_length;E_{#nu};Total pes",40, 0., 2000., 100, 0., 100000.);
  TH2F * hNeuEn2_large_pot_length=new TH2F("hNeuEn2_large_pot_length","hNeuEn2_large_pot_length;E_{#nu};Total pes",40, 0., 2000., 100, 0., 100000.);

  TH2F * hNeuEn2_HITS_short_lambda_max=new TH2F("hNeuEn2_short_lambda_max","hNeuEn2_short_lambda_max;E_{#nu};Total pes",40, 0., 2000., 100, 0., 100000.);
  TH2F * hNeuEn2_HITS_large_lambda_max=new TH2F("hNeuEn2_large_lambda_max","hNeuEn2_large_lambda_max;E_{#nu};Total pes",40, 0., 2000., 100, 0., 100000.);

  TH2F * hNeuEn_HITS=new TH2F("hNeuEn_HITS","hNeuEn_HITS;trueKE;Total hits",60, 0., 3000., 500, 0., 5.);
  TH2F * hNeuEn2_HITS=new TH2F("hNeuEn2_HITS","hNeuEn2_HITS;trueKE;Total hits",40, 0., 2000., 400, 0., 40000.);
  TH2F * hNeuEn_HITS_vs_pot_length=new TH2F("hNeuEn_HITS_vs_pot_length","hNeuEn_HITS_vs_pot_length;E_{#nu};Total hits/pot_length",60, 0., 3000., 200, 0., 2.);
  TH2F * hNeuEn2_HITS_vs_pot_length=new TH2F("hNeuEn2_HITS_vs_pot_length","hNeuEn2_HITS_vs_pot_length;E_{#nu};Total hits/pot_length",40, 0., 2000., 1000, 0., 100.);

  TH2F * hMuEn2_HITS_short_pot_length=new TH2F("hMuEn2_HITS_short_pot_length","hMuEn2_HITS_short_pot_length;trueKE;Total pes",40, 0., 2000., 100, 0., 40000.);
  TH2F * hMuEn2_HITS_large_pot_length=new TH2F("hMuEn2_HITS_large_pot_length","hMuEn2_HITS_large_pot_length;trueKE;Total pes",40, 0., 2000., 100, 0., 40000.);
  TH2F * hNeuEn2_HITS_short_pot_length=new TH2F("hNeuEn2_HITS_short_pot_length","hNeuEn2_HITS_short_pot_length;E_{#nu};Total pes",40, 0., 2000., 100, 0., 40000.);
  TH2F * hNeuEn2_HITS_large_pot_length=new TH2F("hNeuEn2_HITS_large_pot_length","hNeuEn2_HITS_large_pot_length;E_{#nu};Total pes",40, 0., 2000., 100, 0., 40000.);
  // / / 
  TH2F * hMuEn_ring_PEs=new TH2F("hMuEn__ring_PEs","hMuEn_ring_PEs;trueKE;Total hits",60, 0., 3000., 500, 0., 5.);
  TH2F * hMuEn2_ring_PEs=new TH2F("hMuEn2__ring_PEs","hMuEn2_ring_PEs;trueKE;Total hits",40, 0., 2000., 400, 0., 40000.);
  TH2F * hNeuEn_ring_PEs=new TH2F("hNeuEn_ring_PEs","hNeuEn_ring_PEs;trueKE;Total hits",60, 0., 3000., 500, 0., 5.);
  TH2F * hNeuEn2_ring_PEs=new TH2F("hNeuEn2_ring_PEs","hNeuEn2_ring_PEs;trueKE;Total hits",40, 0., 2000., 400, 0., 40000.);
  TH2F * hMuEn2_ring_PEs_short_pot_length=new TH2F("hMuEn2_ring_PEs_short_pot_length","hMuEn2_ring_PEs_short_pot_length;trueKE;Total pes",40, 0., 2000., 100, 0., 40000.);
  TH2F * hMuEn2_ring_PEs_large_pot_length=new TH2F("hMuEn2_ring_PEs_large_pot_length","hMuEn2_ring_PEs_large_pot_length;trueKE;Total pes",40, 0., 2000., 100, 0., 40000.);
  TH2F * hNeuEn2_ring_PEs_short_pot_length=new TH2F("hNeuEn2_ring_PEs_short_pot_length","hNeuEn2_ring_PEs_short_pot_length;E_{#nu};Total pes",40, 0., 2000., 100, 0., 40000.);
  TH2F * hNeuEn2_ring_PEs_large_pot_length=new TH2F("hNeuEn2_ring_PEs_large_pot_length","hNeuEn2_ring_PEs_large_pot_length;E_{#nu};Total pes",40, 0., 2000., 100, 0., 40000.);
  // / /
  TH2F * hMuEn_total_PMTs_hits=new TH2F("hMuEn_total_PMTs_hits","hMuEn_total_PMTs_hits;trueKE;Total PMTs",60, 0., 3000., 400, 0., 4.);
  TH2F * hMuEn2_total_PMTs_hits=new TH2F("hMuEn2_total_PMTs_hits","hMuEn2_total_PMTs_hits;trueKE;Total PMTs",40, 0., 2000., 600, 0., 6000.);
  TH2F * hNeuEn_total_PMTs_hits=new TH2F("hNeuEn_total_PMTs_hits","hNeuEn_total_PMTs_hits;trueKE;Total PMTs",60, 0., 3000., 400, 0., 4.);
  TH2F * hNeuEn2_total_PMTs_hits=new TH2F("hNeuEn2_total_PMTs_hits","hNeuEn2_total_PMTs_hits;trueKE;Total PMTs",40, 0., 2000., 600, 0., 6000.);
  TH2F * hMuEn2_total_PMTs_hits_short_pot_length=new TH2F("hMuEn2_total_PMTs_hits_short_pot_length","hMuEn2_total_PMTs_hits_short_pot_length;trueKE;Total PMT",40, 0., 2000., 600, 0., 6000.);
  TH2F * hMuEn2_total_PMTs_hits_large_pot_length=new TH2F("hMuEn2_total_PMTs_hits_large_pot_length","hMuEn2_total_PMTs_hits_large_pot_length;trueKE;Total PMT",40, 0., 2000., 600, 0., 6000.);
  TH2F * hNeuEn2_total_PMTs_hits_short_pot_length=new TH2F("hNeuEn2_total_PMTs_hits_short_pot_length","hNeuEn2_total_PMTs_hits_short_pot_length;E_{#nu};Total PMT",40, 0., 2000., 600, 0., 6000.);
  TH2F * hNeuEn2_total_PMTs_hits_large_pot_length=new TH2F("hNeuEn2_total_PMTs_hits_large_pot_length","hNeuEn2_total_PMTs_hits_large_pot_length;E_{#nu};Total PMT",40, 0., 2000., 600, 0., 6000.);
  // / / 
  TH1F * plot_trueToWall=new TH1F("plot_trueToWall","plot_trueToWall;trueToWall [cm]",100, 0., 4000.);
  TH1F * plot_trueDWall=new TH1F("plot_trueDWall","plot_trueDWall;trueDWall [cm]",100, 0., 4000.);

  TH1F * plot_trueToWallR=new TH1F("plot_trueToWallR","plot_trueToWallR;trueToWallR [cm]",100, 0., 4000.);
  TH1F * plot_trueDWallR=new TH1F("plot_trueDWallR","plot_trueDWallR;trueDWallR [cm]",100, 0., 4000.);
  TH1F * plot_trueToWallZ=new TH1F("plot_trueToWallZ","plot_trueToWallZ;trueToWallZ [cm]",100, 0., 4000.);
  TH1F * plot_trueDWallZ=new TH1F("plot_trueDWallZ","plot_trueDWallZ;trueDWallZ [cm]",100, 0., 4000.);

  TH1F * plot_pot_length=new TH1F("plot_pot_length","plot_pot_length;potentila length [cm]",100, 0., 4000.);
  TH1F * plot_theta_muDir_track=new TH1F("plot_theta_muDir_track","plot_theta_muDir_track; theta_muDir_track [^{0}]",180, 0., 180.);
  TH1F * plot_vertical_dist=new TH1F("plot_vertical_dist","plot_vertical_dist;vertical distance [cm]",200, 0., 2000.);
  TH1F * plot_vertical_dist_sim=new TH1F("plot_vertical_dist_sim","plot_vertical_dist_sim;vertical distance [cm]",200, 0., 2000.);

// Load Data
    // =========
    // WCLTREEREADER
    //**********************************************************
    // first parameter in the construct: 0=direct output of WCLite
    //                                   1=a "smeared" file made from WCL files
    //                                   2=output from SandboxDetSim
    // second parameter in the constructor: 1 means "include gen level tree"
    //                                      0 means no gen level tree


  TString filename = "nu_numu_iNUM_out_12in_new_fortraining.root";
  TString genfilename = "generatorcardfile.root";
  
  WCLTreeReader *mTR = new WCLTreeReader(2,1);
  mTR->LoadData(filename,genfilename);

  //_____________________________
  TTree * nu_eneNEW = new TTree("nu_eneNEW","nu_eneNEW");
  //double theta_cher = 42.153902; 
  double theta_cher = 45.;
  //_____________________________

  const int maxSubEvts = 20;
  int cluster[maxSubEvts], ring[maxSubEvts];
  int nSubevents, nClusters;
  const int maxPMT = 20000;
  int observedPE[maxSubEvts][maxPMT];
  int allPE[maxSubEvts][maxPMT];
  double houghSpace[2500];

  double trueVtxX;
  double trueVtxY;
  double trueVtxZ;
  double trueTime;
  double trueDirX;
  double trueDirY;
  double trueDirZ;
  double trueKE;
  double diffKE;
  double diffVtxX;
  double diffVtxY;
  double diffVtxZ;
  double diffTime;
  double diffDirX;
  double diffDirY;
  double diffDirZ;
  double diffVtxAbs;
  double diffDirAbs;
  double pointVtxX[maxSubEvts];
  double pointVtxY[maxSubEvts];
  double pointVtxZ[maxSubEvts];
  double pointTime[maxSubEvts];
  double trackVtxX[maxSubEvts];
  double trackVtxY[maxSubEvts];
  double trackVtxZ[maxSubEvts];
  double trackTime[maxSubEvts];
  double trackDirX[maxSubEvts];
  double trackDirY[maxSubEvts];
  double trackDirZ[maxSubEvts];
  double electronTrackCorrection[maxSubEvts];
  double muonTrackCorrection[maxSubEvts];
  int mode;
  double neutrinoE;
  double neutrinoDirX;
  double neutrinoDirY;
  double neutrinoDirZ;
  int neutrinoPID;
  int nNeutrons;
  int neutronCount;
  int nCaptures;
  double trueToWall;
  double trueDWall;
  double recoToWall;
  double recoDWall;
  int recoCaptures;
   
  double recoVtxXLowE[maxSubEvts];
  double recoVtxYLowE[maxSubEvts];
  double recoVtxZLowE[maxSubEvts];
  double recoTimeLowE[maxSubEvts];
  double recoDirXLowE[maxSubEvts];
  double recoDirYLowE[maxSubEvts];
  double recoDirZLowE[maxSubEvts];
  double recoChkvAngleLowE[maxSubEvts];
  double recoEnergyLowE[maxSubEvts];
  double recoVtxXHighEElectron[maxSubEvts];
  double recoVtxYHighEElectron[maxSubEvts];
  double recoVtxZHighEElectron[maxSubEvts];
  double recoTimeHighEElectron[maxSubEvts];
  double recoDirXHighEElectron[maxSubEvts];
  double recoDirYHighEElectron[maxSubEvts];
  double recoDirZHighEElectron[maxSubEvts];
  double recoChkvAngleHighEElectron[maxSubEvts];
  double recoEnergyHighEElectron[maxSubEvts];
  double recoLnLHighEElectron[maxSubEvts];
  double recoVtxXHighEMuon[maxSubEvts];
  double recoVtxYHighEMuon[maxSubEvts];
  double recoVtxZHighEMuon[maxSubEvts];
  double recoTimeHighEMuon[maxSubEvts];
  double recoDirXHighEMuon[maxSubEvts];
  double recoDirYHighEMuon[maxSubEvts];
  double recoDirZHighEMuon[maxSubEvts];
  double recoChkvAngleHighEMuon[maxSubEvts];
  double recoEnergyHighEMuon[maxSubEvts];
  double recoLnLHighEMuon[maxSubEvts];
  int recoNRings[maxSubEvts];
  double recoVtxX[maxSubEvts];
  double recoVtxY[maxSubEvts];
  double recoVtxZ[maxSubEvts];
  double recoTime[maxSubEvts];
  double recoDirX[maxSubEvts];
  double recoDirY[maxSubEvts];
  double recoDirZ[maxSubEvts];
  double recoChkvAngle[maxSubEvts];
  double recoEnergy[maxSubEvts];
  int recoPID[maxSubEvts];
  bool isHighE[maxSubEvts];
  int ringPEs[maxSubEvts];
  double trueExpectedPEHighEMuon[maxSubEvts];
  double recoExpectedPEHighEMuon[maxSubEvts];
  double trueExpectedPEHighEElectron[maxSubEvts];
  double recoExpectedPEHighEElectron[maxSubEvts];
  //---------------------------------
  TFile *likelihoodTables = new TFile(lnLfilename, "READ");
  TH3D *electronPhotons = (TH3D *) likelihoodTables->Get("photons_e");
  TH3D *muonPhotons = (TH3D *) likelihoodTables->Get("photons_mu");
  TH3D *electronTimes = (TH3D *) likelihoodTables->Get("times_e");
  TH3D *muonTimes = (TH3D *) likelihoodTables->Get("times_mu");
 
  TFile *energyTables = new TFile(lookupfilename, "READ");
  TH2D *electronEnergyLookup = (TH2D *) energyTables->Get("nHitDWallKELookup_e");
  TH2D *muonEnergyLookup = (TH2D *) energyTables->Get("nHitDWallKELookup_mu");
  TH2D *electronVtxBiasLookup = (TH2D *) energyTables->Get("nHitDWallVtxTrackBias_e");
  TH2D *muonVtxBiasLookup = (TH2D *) energyTables->Get("nHitDWallVtxTrackBias_mu");
  
  // Get positions, timing res of all PMTs
  //cout << "getting PMT info" << endl;
  int totalPMTs = mTR->GetNpmts();
  //cout << totalPMTs << " PMTs" << endl;
  double * pmtXall = new double[totalPMTs];
  double * pmtYall = new double[totalPMTs];
  double * pmtZall = new double[totalPMTs];
  double * pmtDirXall = new double[totalPMTs];
  double * pmtDirYall = new double[totalPMTs];
  double * pmtDirZall = new double[totalPMTs];
  double * pmtTimeResAll = new double[totalPMTs];
  int * pmtIDall = new int[totalPMTs];
  for(int i=0; i<totalPMTs; i++){
    mTR->LoadPMT(i);
    pmtXall[i] = mTR->get_PMTx()/10.; //pmt positions in mm, convert to cm
    pmtYall[i] = mTR->get_PMTy()/10.;
    pmtZall[i] = mTR->get_PMTz()/10.;
    pmtTimeResAll[i] = mTR->get_PMTtimeRes();
    pmtIDall[i] = mTR->get_PMTid();
    if(TMath::Abs(pmtZall[i]-1100)<0.1){
      pmtDirXall[i] = 0;
      pmtDirYall[i] = 0;
      pmtDirZall[i] = -1;
    }
    else if(TMath::Abs(pmtZall[i]+1100)<0.1){
      pmtDirXall[i] = 0;
      pmtDirYall[i] = 0;
      pmtDirZall[i] = 1;
    }
    else{
      pmtDirXall[i] = -pmtXall[i]/550.;
      pmtDirYall[i] = -pmtYall[i]/550.;
      pmtDirZall[i] = 0;
    }
  }
  //---------------------------------
  int nEntries = mTR->GetEntries();
  cout << "getting nEntries "<<nEntries<<endl;
  
  TFile f_out(outfilename,"recreate");
  int evt;
  for(int i=startEvent; i<(nEntries+1); i++){
  //for(int i=startEvent; i<50; i++){
    cout << "---------------------------------------" << endl;
    cout<<"Loadng event: "<<i<<endl;
    
    mTR->LoadEvent(i);
    evt=mTR->get_genevt();
    int npart = mTR->get_npart();
    Int_t* part_pid = mTR->get_part_pid();
    Double_t* part_KEstart = mTR->get_part_KEstart();
    Double_t* part_KEend = mTR->get_part_KEend();
    Int_t* part_parentid = mTR->get_part_parentid();

    trueVtxX = mTR->get_genvtxx();
    trueVtxY = mTR->get_genvtxy();
    trueVtxZ = mTR->get_genvtxz();
    int trackNumber=0;
    int *pids = mTR->get_genpid();
    while(TMath::Abs(pids[trackNumber])<11 || TMath::Abs(pids[trackNumber])>14)
      trackNumber++;
    trueDirX = mTR->get_genpx()[trackNumber];
    trueDirY = mTR->get_genpy()[trackNumber];
    trueDirZ = mTR->get_genpz()[trackNumber];
    trueKE = mTR->get_genKE()[trackNumber];
    mode = mTR->get_genmode();
    neutrinoE = mTR->get_genE();
    neutrinoDirX = mTR->get_genbeam_px();
    neutrinoDirY = mTR->get_genbeam_py();
    neutrinoDirZ = mTR->get_genbeam_pz();
    neutrinoPID = mTR->get_genbeam_id();
    nNeutrons = mTR->get_gennneutrons();
    neutronCount = mTR->get_neutroncount();
    nCaptures = mTR->get_ncapturecount();
    int totalPEs = mTR->get_nhits();
    int totalPMTs = mTR->GetNpmts();
    int total_hits=0;
    int total_ring_PEs=0;
    int total_PMTs_hits=0;
    int total_reco_Capture_E=0;

    double lambda_min = 10000000;  double lambda_max = -99999999.9; double vertical_dist_max = -99999999.9; 

    if(totalPEs>=10 && (fabs(mode)==1) ){
      double trueVtxR2 = trueVtxX*trueVtxX+trueVtxY*trueVtxY;
      double trueDWallR = 550-TMath::Sqrt(trueVtxR2); //cm
      double trueDWallZ = 1100-TMath::Abs(trueVtxZ);  //cm
      trueDWall = trueDWallR<trueDWallZ ? trueDWallR : trueDWallZ;
      cout<<"trueDWall= "<<trueDWall<<" trueDWallR= "<<trueDWallR<<" trueDWallZ= "<<trueDWallZ<<endl;
      double a = 1-trueDirZ*trueDirZ;
      double b = trueVtxX*trueDirX+trueVtxY*trueDirY;
      double c = trueVtxR2-550*550;
      double trueToWallR = (TMath::Sqrt(b*b-a*c)-b)/a;
      double trueToWallZ = 1100 - trueVtxZ*TMath::Abs(trueDirZ);
      trueToWall = trueToWallR<trueToWallZ ? trueToWallR : trueToWallZ;
      cout<<"trueToWall= "<<trueToWall<<" trueToWallR= "<<trueToWallR<<" trueToWallZ= "<<trueToWallZ<<endl;
      cout<<"trueVtxX= "<<trueVtxX<<" trueVtxY= "<<trueVtxY<<" trueVtxZ= "<<trueVtxZ<<endl;
      cout<<"trueDirX= "<<trueDirX<<" trueDirY= "<<trueDirY<<" truDirZ= "<<trueDirZ<<endl;
      //cout<<"norm= "<<(TMath::Sqrt((trueDirX*trueDirX)+(trueDirY*trueDiY)+(trueDirZ*trueDirZ)))<<endl;
    
      trueTime = 0;
      recoCaptures = 0;
      //    trueDirCosTheta = trueDirZ/TMath::Sqrt(trueDirX*trueDirX+trueDirY*trueDirY+trueDirZ*trueDirZ);
      double trueDirTheta = TMath::ACos(trueDirZ /TMath::Sqrt(trueDirX * trueDirX + trueDirY * trueDirY + trueDirZ * trueDirZ));
      double trueDirPhi = TMath::ATan2(trueDirY, trueDirX);
      
      //cout<<"trueToWall= "<<trueToWall<<" trueToWallR= "<<trueToWallR<<" trueToWallZ= "<<trueToWallZ<<endl;
      //cout<<"trueVtxX= "<<trueVtxX<<" trueVtxY= "<<trueVtxY<<" trueVtxZ= "<<trueVtxZ<<endl;
      //cout<<"trueDirX= "<<trueDirX<<" trueDirY= "<<trueDirY<<" truDirZ= "<<trueDirZ<<endl;
      cout<<"trueKE= "<<trueKE<<" totalPEs= "<<totalPEs<<" trueToWall= "<<trueToWall<<" trueDWall= " <<trueDWall<<endl;
  
      //______________ check hits in rings _______________
      int * hitCluster = new int[totalPEs](); for(int iPE =0; iPE <totalPEs; iPE++) hitCluster[iPE]=0;
      double * PEhitTimes = mTR->get_hitTime();
      int * clusterHitCounts = FindClusters(totalPEs, PEhitTimes, hitCluster, nClusters);

       if(nClusters <1) { //ignore events with no clusters
            //cout << "----- Skipping event with no clusters! -----" << endl << endl;
            nSubevents =0;
            nClusters=0;
            diffVtxX=-9999;
            diffVtxY=-9999;
            diffVtxZ=-9999;
            diffTime=-9999;
            diffDirX=-9999;
            diffDirY=-9999;
            diffDirZ=-9999;
            diffVtxAbs=-9999;
            diffDirAbs=-9999;
            diffKE=-9999;
            continue;
        }
      
      //Loop over sub-events. Cluster 0 is for hits not in any sub-event, so start at 1
      if(nClusters>maxSubEvts) nClusters =maxSubEvts;
      int subevent = 0; 
      nSubevents = nClusters;
      int iHit2=0;
      double *hitx_all = new double[totalPEs];
      double *hity_all = new double[totalPEs];
      double *hitz_all = new double[totalPEs];

      for(int iCluster = 0; subevent< nSubevents; iCluster++, subevent++) {
	cluster[subevent] = iCluster;
	ring[subevent] = 0;
	for(int iPMT=0; iPMT<totalPMTs; iPMT++){
	  observedPE[subevent][iPMT] = 0;
	  allPE[iCluster][iPMT] = 0;
	  //                ring1PE[subevent][iPMT] = 0;
	  //                ring2PE[subevent][iPMT] = 0;
	  //                trueExpectedPEHighEElectron[subevent][iPMT] = 0;
	  //                trueExpectedPEHighEMuon[subevent][iPMT] = 0;
	}
	int nHits = clusterHitCounts[iCluster +1];
	cout << "Hits in sub event " << iCluster +1 << ": " << nHits << endl;
	total_hits+= nHits;

	double *hitx = new double[nHits];
	double *hity = new double[nHits];
	double *hitz = new double[nHits];
	double *hitt = new double[nHits];
	int iClusterHit = 0; 
	double clusterStartTime = 99999;
	for (int iHit = 0; iHit < totalPEs; iHit++) {
	  if (hitCluster[iHit] != iCluster +1) continue;
	  hitx[iClusterHit] = mTR->get_hitX()[iHit];
	  hity[iClusterHit] = mTR->get_hitY()[iHit];
	  hitz[iClusterHit] = mTR->get_hitZ()[iHit];

	  hitx_all[iHit2] = mTR->get_hitX()[iHit];
	  hity_all[iHit2] = mTR->get_hitY()[iHit];
	  hitz_all[iHit2] = mTR->get_hitZ()[iHit];
	  iHit2++;

	  hitt[iClusterHit] = mTR->get_hitTime()[iHit];
	  if(hitt[iClusterHit]<clusterStartTime) clusterStartTime=hitt[iClusterHit];
	  //cout<<"aa_hitx[iClusterHit]= "<<hitx[iClusterHit]<<" hity[iClusterHit]= "<<hity[iClusterHit]<<" hitz[iClusterHit]= "<<hitz[iClusterHit]<<endl;
	  iClusterHit++;
	}
	cout<<"iClusterHit= "<<iClusterHit<<" iHit2= "<<iHit2<<endl;

	// nHits energy estimate.
	LowEReco lowEReco;
	recoEnergyLowE[iCluster] = lowEReco.ReconstructEnergy(nHits, hitx, hity, hitz);
	// Do Low-E reco
	lowEReco.DoLowEReco(i, nHits, hitx, hity, hitz, hitt, recoChkvAngleLowE[iCluster],
			    recoVtxXLowE[iCluster], recoVtxYLowE[iCluster], recoVtxZLowE[iCluster], recoTimeLowE[iCluster],
			    recoDirXLowE[iCluster], recoDirYLowE[iCluster], recoDirZLowE[iCluster]);
	
	delete[] hitx;
	delete[] hity;
	delete[] hitz;
	delete[] hitt;
	
	recoVtxX[subevent] = recoVtxXLowE[iCluster];
	recoVtxY[subevent] = recoVtxYLowE[iCluster];
	recoVtxZ[subevent] = recoVtxZLowE[iCluster];
	recoTime[subevent] = recoTimeLowE[iCluster];
	recoDirX[subevent] = recoDirXLowE[iCluster];
	recoDirY[subevent] = recoDirYLowE[iCluster];
	recoDirZ[subevent] = recoDirZLowE[iCluster];
	recoChkvAngle[subevent] = recoChkvAngleLowE[iCluster];
	recoEnergy[subevent] = recoEnergyLowE[iCluster];
	
	//cout << "--- nhits Energy estimate: " << recoEnergy[subevent];
	//if(subevent == 0) cout << "  (true: " << trueKE << ")";
	//cout << endl;
	cout<<"lowE-- recoEnergy[subevent]= "<<recoEnergy[subevent]<<endl;
	// If < threshold, don't do high E reco
	const int highEThreshold = 60;
	if (recoEnergyLowE[iCluster] < highEThreshold){
	  //Check for neutron capture
	  if(recoEnergy[subevent] > 2 && recoEnergy[subevent] < 10 && recoTime[subevent] > 2000 && recoTime[subevent] < 100000)
	    recoCaptures++;
	  total_reco_Capture_E+= recoEnergy[subevent];
	  recoPID[subevent] = 0;
	  recoVtxXHighEElectron[subevent] = -9999;
	  recoVtxYHighEElectron[subevent] = -9999;
	  recoVtxZHighEElectron[subevent] = -9999;
	  recoTimeHighEElectron[subevent] = -9999;
	  recoDirXHighEElectron[subevent] = -9999;
	  recoDirYHighEElectron[subevent] = -9999;
	  recoDirZHighEElectron[subevent] = -9999;
	  recoChkvAngleHighEElectron[subevent] = -9999;
	  recoEnergyHighEElectron[subevent] = -9999;
	  recoLnLHighEElectron[subevent] = -9999;
	  recoVtxXHighEMuon[subevent] = -9999;
	  recoVtxYHighEMuon[subevent] = -9999;
	  recoVtxZHighEMuon[subevent] = -9999;
	  recoTimeHighEMuon[subevent] = -9999;
	  recoDirXHighEMuon[subevent] = -9999;
	  recoDirYHighEMuon[subevent] = -9999;
	  recoDirZHighEMuon[subevent] = -9999;
	  recoChkvAngleHighEMuon[subevent] = -9999;
	  recoEnergyHighEMuon[subevent] = -9999;
	  recoLnLHighEMuon[subevent] = -9999;
	  recoNRings[iCluster] = 0;
	  isHighE[subevent] = false;
	  continue;
	}
	
	//High-E reco
	HighEReco highEReco;
	highEReco.electronPhotons = electronPhotons;
	highEReco.muonPhotons = muonPhotons;
	highEReco.electronTimes = electronTimes;
	highEReco.muonTimes = muonTimes;
	
	int *hitPMTids = mTR->get_hitPMTid();
	
	//Find number of PEs from current subevent on each PMT
	double maxTime = clusterStartTime+90; //ignore all PEs after 90ns (~22m for photon in water)
	//cout << "Cluster start time: " << clusterStartTime << endl;
	highEReco.nPEs = 0;
	for (int iPE = 0; iPE < totalPEs; iPE++) {
	  if (PEhitTimes[iPE] > maxTime || hitCluster[iPE] != iCluster +1) continue;
	  highEReco.nPEs++;
	  allPE[iCluster][hitPMTids[iPE]]++;
	}
	cout<<"HE: highEReco.nPEs= "<<highEReco.nPEs<<endl;

	//Now create new arrays for PMTs that were actually hit
	highEReco.nHitPMT = 0;
	for (int iPMT = 0; iPMT < totalPMTs; iPMT++) {
	  highEReco.nHitPMT += (allPE[iCluster][iPMT] > 0);
	  total_PMTs_hits+= (allPE[iCluster][iPMT] > 0);
	}
	cout<<"inHE total_PMTs_hits+= "<< total_PMTs_hits<<endl;
	//cout << highEReco.nHitPMT << " PMTs with" << highEReco.nPEs << " hits" << endl;
	highEReco.hitPMTx = new double[highEReco.nHitPMT];
	highEReco.hitPMTy = new double[highEReco.nHitPMT];
	highEReco.hitPMTz = new double[highEReco.nHitPMT];
	highEReco.hitPMTDirX = new double[highEReco.nHitPMT];
	highEReco.hitPMTDirY = new double[highEReco.nHitPMT];
	highEReco.hitPMTDirZ = new double[highEReco.nHitPMT];
	highEReco.hitPMTTimeRes = new double[highEReco.nHitPMT];
	int *pmtID = new int[highEReco.nHitPMT];
	highEReco.hitPMTPEs = new int[highEReco.nHitPMT]();
	for (int iPMT = 0; iPMT < highEReco.nHitPMT; iPMT++) highEReco.hitPMTPEs[iPMT] = 0;
	int *newPMTids = new int[totalPMTs];
	for (int iPMT = 0, iPMT2 = 0; iPMT < totalPMTs; iPMT++) {
	  if (allPE[iCluster][iPMT] == 0) continue;
	  highEReco.hitPMTPEs[iPMT2] = allPE[iCluster][iPMT];
	  highEReco.hitPMTx[iPMT2] = pmtXall[iPMT];
	  highEReco.hitPMTy[iPMT2] = pmtYall[iPMT];
	  highEReco.hitPMTz[iPMT2] = pmtZall[iPMT];
	  highEReco.hitPMTDirX[iPMT2] = pmtDirXall[iPMT];
	  highEReco.hitPMTDirY[iPMT2] = pmtDirYall[iPMT];
	  highEReco.hitPMTDirZ[iPMT2] = pmtDirZall[iPMT];
	  highEReco.hitPMTTimeRes[iPMT2] = pmtTimeResAll[iPMT];
	  pmtID[iPMT2] = pmtIDall[iPMT];
	  newPMTids[iPMT] = iPMT2;
	  iPMT2++;
	}
	////	cout<<"what is iPMT2= "<<iPMT2<<endl; //= total_PMTs_hits
	highEReco.hitPMT = new int[highEReco.nPEs];
	highEReco.hitT = new double[highEReco.nPEs];
	highEReco.hitRing = new int[highEReco.nPEs];
	// Loop over PEs create arrays of hit PMT, hit time, etc, for use later
	for (int iPE = 0, iPE2 = 0; iPE < totalPEs; iPE++) {
	  if (PEhitTimes[iPE] > maxTime || hitCluster[iPE] != iCluster +1) continue;
	  int iPMT = newPMTids[hitPMTids[iPE]];
	  highEReco.hitPMT[iPE2] = iPMT;
	  highEReco.hitT[iPE2] = PEhitTimes[iPE];
	  highEReco.hitRing[iPE2] = 999;
	  iPE2++;
	}
	delete[] newPMTids;
	
	//Initial point fit for starting vertex
	double recoPar[7];
	highEReco.PointFit(recoVtxX[subevent], recoVtxY[subevent], recoVtxZ[subevent], recoTime[subevent], recoChkvAngle[subevent]);
	pointVtxX[iCluster] = recoVtxX[subevent];
	pointVtxY[iCluster] = recoVtxY[subevent];
	pointVtxZ[iCluster] = recoVtxZ[subevent];
	pointTime[iCluster] = recoTime[subevent];
	
	int maxRings = maxSubEvts-subevent;
	double *peakTheta = new double[maxRings];
	double *peakPhi = new double[maxRings];
	int * ringPE = new int[maxRings];
	double * hough = iCluster==0 ? houghSpace : 0;
	recoNRings[iCluster] = highEReco.FindRings(recoVtxX[subevent], recoVtxY[subevent], recoVtxZ[subevent],
						   recoTime[subevent], peakTheta, peakPhi, ringPE, maxRings, false,
						   hough);
	
	nSubevents += recoNRings[iCluster]-1;
	if(nSubevents > maxSubEvts) nSubevents=maxSubEvts;
	
	//cout << "Found " << recoNRings[iCluster] << " rings." << endl;
	for(int iRing = 0; iRing < recoNRings[iCluster]; iRing++)
	  cout << "Ring " << iRing+1 << ": " << ringPE[iRing] << " PEs    theta=" << peakTheta[iRing] << " phi=" << peakPhi[iRing] << endl;

	//Save all cluster hit info
	int * clusterHitPMTPEs = highEReco.hitPMTPEs;
	int * clusterHitPMT = highEReco.hitPMT;
	double * clusterHitT = highEReco.hitT;
	int clusterPEs = highEReco.nPEs;
	double * clusterHitPMTx = highEReco.hitPMTx;
	double * clusterHitPMTy = highEReco.hitPMTy;
	double * clusterHitPMTz = highEReco.hitPMTz;
	double * clusterHitPMTDirX = highEReco.hitPMTDirX;
	double * clusterHitPMTDirY = highEReco.hitPMTDirY;
	double * clusterHitPMTDirZ = highEReco.hitPMTDirZ;
	double * clusterHitPMTTimeRes = highEReco.hitPMTTimeRes;
	int * clusterPMTid = pmtID;
	int clusterHitPMTs = highEReco.nHitPMT;
	
	// Copy point fit seed to each ring
	for(int iRing = 1; iRing < recoNRings[iCluster]; iRing++){
	  recoVtxX[subevent+iRing] = recoVtxX[subevent];
	  recoVtxY[subevent+iRing] = recoVtxY[subevent];
	  recoVtxZ[subevent+iRing] = recoVtxZ[subevent];
	  recoTime[subevent+iRing] = recoTime[subevent];
	  recoChkvAngle[subevent+iRing] = recoChkvAngle[subevent];
	  for(int iPMT=0; iPMT<totalPMTs; iPMT++) observedPE[subevent+iRing][iPMT] = 0;
	}
	//ringPE[0] = clusterPEs;
	for(int iRing = 0; iRing < recoNRings[iCluster]; iRing++, subevent++){
	  ringPEs[subevent] = ringPE[iRing];
	  cluster[subevent] = iCluster;
	  double recoDirTheta = peakTheta[iRing];
	  double recoDirCosTheta = TMath::Cos(recoDirTheta);
	  double recoDirPhi = peakPhi[iRing];

	  //cout << "Ring " << iRing+1 << endl;

	  //Only use PEs from this ring from now on
	  ring[subevent] = iRing+1;
	  double *newHitT = new double[ringPEs[subevent]];
	  int *newHitPMT = new int[ringPEs[subevent]];
	  highEReco.hitPMTPEs = new int[clusterHitPMTs](); for(int iPMT =0; iPMT < clusterHitPMTs; iPMT++) highEReco.hitPMTPEs[iPMT]=0;
	  //                cout << iRing << " " << ringPEs[subevent] << ":" << endl;
	  int nPEnew = 0;
	  for (int iPE = 0; iPE < clusterPEs; iPE++) {
	    if (highEReco.hitRing[iPE] != iRing+1/* && iRing>0*/) continue;
	    //                    cout << iPE << " " << nPEnew << endl;
	    newHitT[nPEnew] = clusterHitT[iPE];
	    newHitPMT[nPEnew] = clusterHitPMT[iPE];
	    highEReco.hitPMTPEs[clusterHitPMT[iPE]]++;
	    nPEnew++;
	  }
	  highEReco.hitT = newHitT;
	  highEReco.hitPMT = newHitPMT;
	  highEReco.nPEs = nPEnew;
	  
	  //Now create new arrays for PMTs that were actually hit
	  int nPMTnew = 0;
	  for (int iPMT = 0; iPMT < clusterHitPMTs; iPMT++) {
	    nPMTnew += (highEReco.hitPMTPEs[iPMT] > 0);
	  }
	  //cout << nPMTnew << " PMTs with hits" << endl;
	  highEReco.hitPMTx = new double[nPMTnew];
	  highEReco.hitPMTy = new double[nPMTnew];
	  highEReco.hitPMTz = new double[nPMTnew];
	  highEReco.hitPMTDirX = new double[nPMTnew];
	  highEReco.hitPMTDirY = new double[nPMTnew];
	  highEReco.hitPMTDirZ = new double[nPMTnew];
	  highEReco.hitPMTTimeRes = new double[nPMTnew];
	  pmtID = new int[nPMTnew];
	  int *pmtPEsNew = new int[nPMTnew](); for(int iPMT =0; iPMT <nPMTnew; iPMT++) pmtPEsNew[iPMT]=0;
	  newPMTids = new int[clusterHitPMTs];
	  for (int iPMT = 0, iPMT2 = 0; iPMT < clusterHitPMTs; iPMT++) {
	    if (highEReco.hitPMTPEs[iPMT] < 1) continue;
	    pmtPEsNew[iPMT2] = highEReco.hitPMTPEs[iPMT];
	    highEReco.hitPMTx[iPMT2] = clusterHitPMTx[iPMT];
	    highEReco.hitPMTy[iPMT2] = clusterHitPMTy[iPMT];
	    highEReco.hitPMTz[iPMT2] = clusterHitPMTz[iPMT];
	    highEReco.hitPMTDirX[iPMT2] = clusterHitPMTDirX[iPMT];
	    highEReco.hitPMTDirY[iPMT2] = clusterHitPMTDirY[iPMT];
	    highEReco.hitPMTDirZ[iPMT2] = clusterHitPMTDirZ[iPMT];
	    highEReco.hitPMTTimeRes[iPMT2] = clusterHitPMTTimeRes[iPMT];
	    pmtID[iPMT2] = clusterPMTid[iPMT];
	    newPMTids[iPMT] = iPMT2;
	    iPMT2++;
	  }
	  delete[] highEReco.hitPMTPEs;
	  highEReco.hitPMTPEs = pmtPEsNew;
	  highEReco.nHitPMT = nPMTnew;
	  for (int iPE = 0; iPE < highEReco.nPEs; iPE++){
	    highEReco.hitPMT[iPE] = newPMTids[highEReco.hitPMT[iPE]];
	  }
	  delete[] newPMTids;
	  cout<<"HE: (ringPEs[subevent])= "<<(ringPEs[subevent])<<endl;
	  total_ring_PEs+=(ringPEs[subevent]);
	  
	  //cout << "Use approx direction of ring as initial track direction" << endl;
	  //cout << "True vtx: (" << trueVtxX << ", " << trueVtxY << ", " << trueVtxZ << ", " << trueTime
	  // << ") dir: (" << TMath::Sin(trueDirTheta) * TMath::Cos(trueDirPhi) << " " <<
	  //  TMath::Sin(trueDirTheta) * TMath::Sin(trueDirPhi) << " " << TMath::Cos(trueDirTheta) << ")" << endl;
	  //cout << "Reco vtx: (" << recoVtxX[subevent] << ", " << recoVtxY[subevent] << ", " << recoVtxZ[subevent] <<
	  // ", " << recoTime[subevent]
	  //<< ") dir: (" << TMath::Sin(recoDirTheta) * TMath::Cos(recoDirPhi) << " " <<
	  // TMath::Sin(recoDirTheta) * TMath::Sin(recoDirPhi) << " " << TMath::Cos(recoDirTheta) << ")" << endl;
	  recoPar[4] = recoDirCosTheta;
	  recoPar[5] = recoDirPhi;
	  
	  //Fit including track:
	  highEReco.TrackFit(recoVtxX[subevent], recoVtxY[subevent], recoVtxZ[subevent], recoTime[subevent],
			     recoDirPhi, recoDirTheta,
			     recoDirCosTheta, recoChkvAngle[subevent]);
	  recoDirX[subevent] = TMath::Sin(recoDirTheta) * TMath::Cos(recoDirPhi);
	  recoDirY[subevent] = TMath::Sin(recoDirTheta) * TMath::Sin(recoDirPhi);
	  recoDirZ[subevent] = TMath::Cos(recoDirTheta);
	  trackVtxX[subevent] = recoVtxX[subevent];
	  trackVtxY[subevent] = recoVtxY[subevent];
	  trackVtxZ[subevent] = recoVtxZ[subevent];
	  trackTime[subevent] = recoTime[subevent];
	  trackDirX[subevent] = recoDirX[subevent];
	  trackDirY[subevent] = recoDirY[subevent];
	  trackDirZ[subevent] = recoDirZ[subevent];

	  //Energy reconstruction and correction to vertex in track direction

	  double dWallZ = 1100 - TMath::Abs(recoVtxZ[subevent]);
	  double dWallR = 550 - TMath::Sqrt(
					    recoVtxX[subevent] * recoVtxX[subevent] + recoVtxY[subevent] * recoVtxY[subevent]);
	  double dWall = dWallR < dWallZ ? dWallR : dWallZ;
	  double corr_factor_phot_coverage = 0.9; // 36%/40%=0.9
	  double lookupPEs = (ringPEs[subevent]);//*corr_factor_phot_coverage;
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
	  recoEnergyHighEMuon[subevent] = muonEnergyLookup->Interpolate(lookupPEs, dWall)*correctionFactor;
	  recoEnergyHighEElectron[subevent] = electronEnergyLookup->Interpolate(lookupPEs, dWall)*correctionFactor;
	  muonTrackCorrection[subevent] = muonVtxBiasLookup->Interpolate(lookupPEs, dWall);
	  electronTrackCorrection[subevent] = electronVtxBiasLookup->Interpolate(lookupPEs, dWall);
	  
	  //cout << "Energy lookup:    Muon: " << recoEnergyHighEMuon[subevent] << "    Electron: " <<
	  //recoEnergyHighEElectron[subevent] << endl;
	  
	  //Need to use unhit PMTs in likelihood
	  for (int iPMT = 0; iPMT < highEReco.nHitPMT; iPMT++) {
	    observedPE[subevent][pmtID[iPMT]] = highEReco.hitPMTPEs[iPMT];
	  }
	  highEReco.nUnhitPMT = totalPMTs - highEReco.nHitPMT;
	  highEReco.unhitPMTx = new double[highEReco.nUnhitPMT];
	  highEReco.unhitPMTy = new double[highEReco.nUnhitPMT];
	  highEReco.unhitPMTz = new double[highEReco.nUnhitPMT];
	  highEReco.unhitPMTDirX = new double[highEReco.nUnhitPMT];
	  highEReco.unhitPMTDirY = new double[highEReco.nUnhitPMT];
	  highEReco.unhitPMTDirZ = new double[highEReco.nUnhitPMT];
	  highEReco.unhitPMTTimeRes = new double[highEReco.nUnhitPMT];
	  //cout << "Total PMTs:" << totalPMTs << " Hit:" << highEReco.nHitPMT << " Unhit:" << highEReco.nUnhitPMT << endl;
	  for (int iPMT = 0, iPMT2 = 0; iPMT < totalPMTs; iPMT++) {
	    if (observedPE[subevent][iPMT] > 0) continue;
	    highEReco.unhitPMTx[iPMT2] = pmtXall[iPMT];
	    highEReco.unhitPMTy[iPMT2] = pmtYall[iPMT];
	    highEReco.unhitPMTz[iPMT2] = pmtZall[iPMT];
	    highEReco.unhitPMTDirX[iPMT2] = pmtDirXall[iPMT];
	    highEReco.unhitPMTDirY[iPMT2] = pmtDirYall[iPMT];
	    highEReco.unhitPMTDirZ[iPMT2] = pmtDirZall[iPMT];
	    highEReco.unhitPMTTimeRes[iPMT2] = pmtTimeResAll[iPMT];
	    iPMT2++;
	  }

	  // Perform maximum likelihood for electron hypothesis
	  //cout << "Fit for electron hypothesis" << endl;
	  recoVtxXHighEElectron[subevent] = recoVtxX[subevent];
	  recoVtxYHighEElectron[subevent] = recoVtxY[subevent];
	  recoVtxZHighEElectron[subevent] = recoVtxZ[subevent];
	  recoTimeHighEElectron[subevent] = recoTime[subevent];
	  recoChkvAngleHighEElectron[subevent] = recoChkvAngleLowE[iCluster];
	  double recoDirPhiHighEElectron = recoDirPhi;
	  double recoDirThetaHighEElectron = recoDirTheta;
	  highEReco.LikelihoodFit(electronTrackCorrection[subevent], recoVtxXHighEElectron[subevent],
				  recoVtxYHighEElectron[subevent], recoVtxZHighEElectron[subevent],
				  recoTimeHighEElectron[subevent],
				  recoDirPhiHighEElectron, recoDirThetaHighEElectron,
				  recoEnergyHighEElectron[subevent], recoLnLHighEElectron[subevent], 11);
	  recoDirXHighEElectron[subevent] =
	    TMath::Sin(recoDirThetaHighEElectron) * TMath::Cos(recoDirPhiHighEElectron);
	  recoDirYHighEElectron[subevent] =
	    TMath::Sin(recoDirThetaHighEElectron) * TMath::Sin(recoDirPhiHighEElectron);
	  recoDirZHighEElectron[subevent] = TMath::Cos(recoDirThetaHighEElectron);
	  trueExpectedPEHighEElectron[subevent] = 0;
	  recoExpectedPEHighEElectron[subevent] = 0;
	  for (int iPMT = 0, iHitPMT = 0, iUnhitPMT = 0; iPMT < totalPMTs; iPMT++) {
	    bool isHit = observedPE[subevent][iPMT] > 0;
	    int thisID = isHit ? iHitPMT++ : iUnhitPMT++;
	    trueExpectedPEHighEElectron[subevent] += highEReco.ExpectedPMTPhotoelectrons(
											 trueVtxX, trueVtxY, trueVtxZ,
											 trueDirX, trueDirY, trueDirZ,
											 thisID, isHit, trueKE, 11);
	    recoExpectedPEHighEElectron[subevent] += highEReco.ExpectedPMTPhotoelectrons(
											 recoVtxXHighEElectron[subevent], recoVtxYHighEElectron[subevent], recoVtxZHighEElectron[subevent],
											 recoDirXHighEElectron[subevent], recoDirYHighEElectron[subevent], recoDirZHighEElectron[subevent],
											 thisID, isHit, recoEnergyHighEElectron[subevent], 11);
	  }
	  //cout << "Expected PE: true: " << trueExpectedPEHighEElectron[subevent] << "  reco: " << recoExpectedPEHighEElectron[subevent] << "  observed: " << ringPEs[subevent] << endl;
	  // Perform maximum likelihood for muon hypothesis
	  //cout << "Fit for muon hypothesis" << endl;
	  recoVtxXHighEMuon[subevent] = recoVtxX[subevent];
	  recoVtxYHighEMuon[subevent] = recoVtxY[subevent];
	  recoVtxZHighEMuon[subevent] = recoVtxZ[subevent];
	  recoTimeHighEMuon[subevent] = recoTime[subevent];
	  recoChkvAngleHighEMuon[subevent] = recoChkvAngleLowE[iCluster];
	  double recoDirPhiHighEMuon = recoDirPhi;
	  double recoDirThetaHighEMuon = recoDirTheta;
	  highEReco.LikelihoodFit(muonTrackCorrection[subevent], recoVtxXHighEMuon[subevent],
				  recoVtxYHighEMuon[subevent], recoVtxZHighEMuon[subevent],
				  recoTimeHighEMuon[subevent],
				  recoDirPhiHighEMuon, recoDirThetaHighEMuon,
				  recoEnergyHighEMuon[subevent], recoLnLHighEMuon[subevent], 13);
	  recoDirXHighEMuon[subevent] = TMath::Sin(recoDirThetaHighEMuon) * TMath::Cos(recoDirPhiHighEMuon);
	  recoDirYHighEMuon[subevent] = TMath::Sin(recoDirThetaHighEMuon) * TMath::Sin(recoDirPhiHighEMuon);
	  recoDirZHighEMuon[subevent] = TMath::Cos(recoDirThetaHighEMuon);
	  trueExpectedPEHighEMuon[subevent] = 0;
	  recoExpectedPEHighEMuon[subevent] = 0;
	  for (int iPMT = 0, iHitPMT = 0, iUnhitPMT = 0; iPMT < totalPMTs; iPMT++) {
	    bool isHit = observedPE[subevent][iPMT] > 0;
	    int thisID = isHit ? iHitPMT++ : iUnhitPMT++;
	    trueExpectedPEHighEMuon[subevent] += highEReco.ExpectedPMTPhotoelectrons(
										     trueVtxX, trueVtxY, trueVtxZ,
										     trueDirX, trueDirY, trueDirZ,
										     thisID, isHit, trueKE, 13);
	    recoExpectedPEHighEMuon[subevent] += highEReco.ExpectedPMTPhotoelectrons(
										     recoVtxXHighEMuon[subevent], recoVtxYHighEMuon[subevent], recoVtxZHighEMuon[subevent],
										     recoDirXHighEMuon[subevent], recoDirYHighEMuon[subevent], recoDirZHighEMuon[subevent],
										     thisID, isHit, recoEnergyHighEMuon[subevent], 13);
	  }
	  //cout << "Expected PE: true: " << trueExpectedPEHighEMuon[subevent] << "  reco: " << recoExpectedPEHighEMuon[subevent] << "  observed: " << ringPEs[subevent] << endl;
	  
	  //cout << "Reco  LnL_mu - LnL_e = " << recoLnLHighEMuon[subevent] - recoLnLHighEElectron[subevent] << endl;
	  //cout << endl << endl;

	  cout<<"recoEnergy[subevent]= "<<recoEnergy[subevent]<<endl;

	  if (recoLnLHighEMuon[subevent] - recoLnLHighEElectron[subevent] > 0) {
	    recoPID[subevent] = 13;
	    recoVtxX[subevent] = recoVtxXHighEMuon[subevent];
	    recoVtxY[subevent] = recoVtxYHighEMuon[subevent];
	    recoVtxZ[subevent] = recoVtxZHighEMuon[subevent];
	    recoTime[subevent] = recoTimeHighEMuon[subevent];
	    recoDirX[subevent] = recoDirXHighEMuon[subevent];
	    recoDirY[subevent] = recoDirYHighEMuon[subevent];
	    recoDirZ[subevent] = recoDirZHighEMuon[subevent];
	    recoChkvAngle[subevent] = recoChkvAngleHighEMuon[subevent];
	    recoEnergy[subevent] = recoEnergyHighEMuon[subevent];
	  }
	  else {
	    recoPID[subevent] = 11;
	    recoVtxX[subevent] = recoVtxXHighEElectron[subevent];
	    recoVtxY[subevent] = recoVtxYHighEElectron[subevent];
	    recoVtxZ[subevent] = recoVtxZHighEElectron[subevent];
	    recoTime[subevent] = recoTimeHighEElectron[subevent];
	    recoDirX[subevent] = recoDirXHighEElectron[subevent];
	    recoDirY[subevent] = recoDirYHighEElectron[subevent];
	    recoDirZ[subevent] = recoDirZHighEElectron[subevent];
	    recoChkvAngle[subevent] = recoChkvAngleHighEElectron[subevent];
	    recoEnergy[subevent] = recoEnergyHighEElectron[subevent];
	  }
	  isHighE[subevent] = true;
	 
	 
 
	  delete[] highEReco.hitPMTx;
	  delete[] highEReco.hitPMTy;
	  delete[] highEReco.hitPMTz;
	  delete[] highEReco.hitPMTDirX;
	  delete[] highEReco.hitPMTDirY;
	  delete[] highEReco.hitPMTDirZ;
	  delete[] highEReco.hitPMTTimeRes;
	  delete[] pmtID;
	  delete[] highEReco.hitT;
	  delete[] highEReco.hitPMT;
	  delete[] highEReco.hitPMTPEs;
	  delete[] highEReco.unhitPMTx;
	  delete[] highEReco.unhitPMTy;
	  delete[] highEReco.unhitPMTz;
	  delete[] highEReco.unhitPMTDirX;
	  delete[] highEReco.unhitPMTDirY;
	  delete[] highEReco.unhitPMTDirZ;
	  delete[] highEReco.unhitPMTTimeRes;
	}
	subevent--;
	delete[] highEReco.hitRing;
	delete[] clusterHitPMTPEs;
	delete[] clusterHitPMT;
	delete[] clusterHitT;
	delete[] clusterHitPMTx;
	delete[] clusterHitPMTy;
	delete[] clusterHitPMTz;
	delete[] clusterHitPMTDirX;
	delete[] clusterHitPMTDirY;
	delete[] clusterHitPMTDirZ;
	delete[] clusterHitPMTTimeRes;
	delete[] clusterPMTid;
	delete[] peakTheta;
	delete[] peakPhi;
	delete[] ringPE;
	//    break;
      
      }//iCluster++, subevent++     

       //_______________________
	  cout<<"i am here: total_hits= "<<total_hits<<" totalPEs= "<< totalPEs<<endl;
	  
	  for (int iClusterHit2=0; iClusterHit2<iHit2; iClusterHit2++){
	   
	    double lambda = find_lambda(recoVtxX[0],recoVtxY[0],recoVtxZ[0],recoDirX[0],recoDirY[0],recoDirZ[0],hitx_all[iClusterHit2]/10.,hity_all[iClusterHit2]/10.,hitz_all[iClusterHit2]/10.,theta_cher);
	    double vertical_dist = find_vert_dist(recoVtxX[0],recoVtxY[0],recoVtxZ[0],recoDirX[0],recoDirY[0],recoDirZ[0],hitx_all[iClusterHit2]/10.,hity_all[iClusterHit2]/10.,hitz_all[iClusterHit2]/10.);
	    double vertical_dist_sim = find_vert_dist(trueVtxX,trueVtxY,trueVtxZ,trueDirX,trueDirY,trueDirZ,hitx_all[iClusterHit2]/10.,hity_all[iClusterHit2]/10.,hitz_all[iClusterHit2]/10.);
	    
	    plot_vertical_dist_sim->Fill(vertical_dist_sim);
	    plot_vertical_dist->Fill(vertical_dist);
	    //cout<<"hitx_all[iClusterHit2]= "<<hitx_all[iClusterHit2]/10.<<" hity_all[iClusterHit2]= "<<hity_all[iClusterHit2]/10.<<" hitz_all[iClusterHit2]= "<<hitz_all[iClusterHit2]/10.<<endl;
	    /* 
	    xmupos_t1 = recoVtxX[0] + recoDirX[0]*lambda; 
	    ymupos_t1 = recoVtxY[0] + recoDirY[0]*lambda; 
	    zmupos_t1 = recoVtxZ[0] + recoDirZ[0]*lambda; 

	    double tr = sqrt((hitx_all[iClusterHit2]/10. - xmupos_t1)*(hitx_all[iClusterHit2]/10. - xmupos_t1) + (hity_all[iClusterHit2]/10. - ymupos_t1)*(hity_all[iClusterHit2]/10. - ymupos_t1) + (hitz_all[iClusterHit2]/10. - zmupos_t1)*(hitz_all[iClusterHit2]/10. - zmupos_t1));      
	    double THE_theta_muDir_track = (acos( (recoDirX[0]*(hitx_all[iClusterHit2]/10. - xmupos_t1) + recoDirY[0]*(hity_all[iClusterHit2]/10. - ymupos_t1) + recoDirZ[0]*(hitz_all[iClusterHit2]/10. - zmupos_t1))/(tr) )*TMath::RadToDeg());
	    plot_theta_muDir_track->Fill(THE_theta_muDir_track); */
	    
	    double photon_length_from_vertex= sqrt( (trueVtxX-hitx_all[iClusterHit2]/10.)*(trueVtxX-hitx_all[iClusterHit2]/10.) + (trueVtxY-hity_all[iClusterHit2]/10.)*(trueVtxY-hity_all[iClusterHit2]/10.) + (trueVtxZ-hitz_all[iClusterHit2]/10.)*(trueVtxZ-hitz_all[iClusterHit2]/10.) );
	  
	    //cout<<" lambda= "<<lambda<<" vertical_dist= "<<vertical_dist<<" photon_length_from_vertex= "<<photon_length_from_vertex<<endl;
	    if( lambda <= lambda_min ){
	      lambda_min = lambda;
	    }
	    if( lambda >= lambda_max ){
	      lambda_max = lambda;
	    }
	    //--------
	    //Find Maximum Magnitude (Vertical Distance between OM & Muon Track):
	    if(lambda>0.){
	      if( vertical_dist_max < vertical_dist){
		vertical_dist_max = vertical_dist;   
		}}
	  }//end of iClusterHit2
	  
	  delete[] hitx_all;
	  delete[] hity_all;
	  delete[] hitz_all;
	  //_____________________
      
      cout<<"------------------------"<<endl;
      cout<<"lambda_min= "<<lambda_min<<" lambda_max= "<<lambda_max<<" vertical_dist_max= "<<vertical_dist_max<<endl;
      cout<<"recoVtxX[0]= "<<recoVtxX[0]<<" recoVtxY[0]= "<<recoVtxY[0]<<" recoVtxZ[0]= "<<recoVtxZ[0]<<endl;
      cout<<"trueVtxX= "<<trueVtxX<<" trueVtxY= "<<trueVtxY<<" trueVtxZ= "<<trueVtxZ<<endl;
      cout<<"Rvtx= "<<sqrt((trueVtxX*trueVtxX)+(trueVtxY*trueVtxY))<<" Rreco_vtx= "<<sqrt((recoVtxX[0]*recoVtxX[0])+(recoVtxY[0]*recoVtxY[0]))<<endl;
      cout<<"trueDirX= "<<trueDirX<<" trueDirY= "<<trueDirY<<" truDirZ= "<<trueDirZ<<endl;
      cout<<"trueToWall= "<<trueToWall<<" trueToWallR= "<<trueToWallR<<" trueToWallZ= "<<trueToWallZ<<endl;
      cout<<"trueKE= "<<trueKE<<" totalPEs= "<<totalPEs<<" trueToWall= "<<trueToWall<<" trueDWall= " <<trueDWall<<endl;


      trueDirX = TMath::Sin(trueDirTheta) * TMath::Cos(trueDirPhi);
      trueDirY = TMath::Sin(trueDirTheta) * TMath::Sin(trueDirPhi);
      trueDirZ = TMath::Cos(trueDirTheta);
      diffVtxX = recoVtxX[0]-trueVtxX;
      diffVtxY = recoVtxY[0]-trueVtxY;
      diffVtxZ = recoVtxZ[0]-trueVtxZ;
      diffVtxAbs = TMath::Sqrt(diffVtxX*diffVtxX + diffVtxY*diffVtxY + diffVtxZ*diffVtxZ);
      diffTime = recoTime[0]-trueTime;
      diffDirX = recoDirX[0]-trueDirX;
      diffDirY = recoDirY[0]-trueDirY;
      diffDirZ = recoDirZ[0]-trueDirZ;
      diffDirAbs = TMath::ACos(recoDirX[0]*trueDirX + recoDirY[0]*trueDirY + recoDirZ[0]*trueDirZ);
      diffKE= recoEnergy[0]-trueKE;
      double recoVtxR2 = recoVtxX[0]*recoVtxX[0]+recoVtxY[0]*recoVtxY[0];
      double recoDWallR = 550-TMath::Sqrt(recoVtxR2);
      double recoDWallZ = 1100-TMath::Abs(recoVtxZ[0]);
      recoDWall = recoDWallR<recoDWallZ ? recoDWallR : recoDWallZ;
      a = 1-recoDirZ[0]*recoDirZ[0];
      b = recoVtxX[0]*recoDirX[0]+recoVtxY[0]*recoDirY[0];
      c = recoVtxR2-550*550;
      double recoToWallR = (TMath::Sqrt(b*b-a*c)-b)/a;
      double recoToWallZ = 1100 - recoVtxZ[0]*TMath::Abs(recoDirZ[0]);
      recoToWall = recoToWallR<recoToWallZ ? recoToWallR : recoToWallZ;
      double vtxTrackBias = diffVtxX*recoDirX[0]+diffVtxY*recoDirY[0]+diffVtxZ*recoDirZ[0];
      
      cout<<" diffVtxAbs= "<< diffVtxAbs<<" diffDirAbs= "<<diffDirAbs<<" diffKE= "<<diffKE<<endl;
      cout<<"recoDWallR= "<<recoDWallR<<" recoDWallZ= "<<recoDWallZ<<endl;

      //---------------------------------------------------------------------------------------
      //-------------- Calcualte Muon Track Intersection Point with Detector -----------------
      
      double det_radious = 550.; double det_z_max = 1100.; double det_z_min = -1100.;
      double lambda_max2=10.;
      std::vector<double> sol = finding_solutions_for_qudratic_eq( (recoDirX[0]*recoDirX[0])+(recoDirY[0]*recoDirY[0]) , 2*(recoDirX[0]*recoVtxX[0]+recoDirY[0]*recoVtxY[0]) , (recoVtxX[0]*recoVtxX[0])+(recoVtxY[0]*recoVtxY[0])-(det_radious*det_radious) );
      double time = CalculatingIntersectionPoint( sol, lambda_max2);
      bool solution_found=FindSolution( sol, lambda_max2 );
      
      if(solution_found==0){cout<<"_______ repeating calc _______"<<endl;
	lambda_max2=lambda_max2+10.;
	time = CalculatingIntersectionPoint( sol, lambda_max2);
	solution_found=FindSolution( sol, lambda_max2 );
      }
      
      double zout = recoVtxZ[0] + time*recoDirZ[0];
      double xout = 0.0; double yout = 0.0;
      
      if(zout < det_z_max && zout > det_z_min){
	xout = recoVtxX[0] + time*recoDirX[0];
	yout = recoVtxY[0] + time*recoDirY[0];
      }
      else{
	if(recoDirZ[0] > 0){
	time = (det_z_max - recoVtxZ[0])/recoDirZ[0];
	zout = det_z_max;
	xout = recoVtxX[0] + time*recoDirX[0];
	yout = recoVtxY[0] + time*recoDirY[0];
	// xout = xmu_rec + time*xrecDir;
	// yout = ymu_rec + time*yrecDir;
	}
	if(recoDirZ[0] < 0){
	  time = (det_z_min - recoVtxZ[0])/recoDirZ[0];
	  zout = det_z_min;
	  xout = recoVtxX[0] + time*recoDirX[0];
	  yout = recoVtxY[0] + time*recoDirY[0];
	}
      }
      cout<<"xout= "<<xout<<" yout= "<<yout<<" zout= "<<zout<<endl;
      cout<<"Rout= "<<sqrt((xout*xout)+(yout*yout))<<endl;
      //---------------------------------------------------------------------------------------
      double pot_length=(TMath::Sqrt( (recoVtxX[0]-xout)*(recoVtxX[0]-xout) + (recoVtxY[0]-yout)*(recoVtxY[0]-yout) + (recoVtxZ[0]-zout)*(recoVtxZ[0]-zout) ));
      cout<<"pot_length= "<<pot_length<<endl;
      
      hMuEn->Fill(trueKE, TMath::Log10(totalPEs) );
      hMuEn2->Fill(trueKE, totalPEs);
      if(pot_length>0.){
	hMuEn2_trueToWall->Fill(trueKE, totalPEs/pot_length);  }
      if(pot_length<300.){
	hMuEn2_short_pot_length->Fill(trueKE, totalPEs);
      }
      if(pot_length>=300.){
	hMuEn2_large_pot_length->Fill(trueKE, totalPEs);
      }
      //----
      hNeuEn2->Fill(neutrinoE, totalPEs);
      if(pot_length<300.){
	hNeuEn2_short_pot_length->Fill(neutrinoE, totalPEs);
      }
      if(pot_length>=300.){
      hNeuEn2_large_pot_length->Fill(neutrinoE, totalPEs);
      }
    
      plot_trueToWall->Fill(trueToWall);
      plot_trueDWall->Fill(trueDWall);
      plot_trueToWallR->Fill(trueToWallR);
      plot_trueDWallR->Fill(trueDWallR);
      plot_trueToWallZ->Fill(trueToWallZ);
      plot_trueDWallZ->Fill(trueDWallZ);
    
      plot_pot_length->Fill(pot_length);
      
      ///_________________________________________________________________________________________
      ///_________________________________________________________________________________________
      cout<<"____final:"<<" mode: " <<mode<<"neutrinoE= "<<neutrinoE<<" trueKE= "<<trueKE<<endl;
      cout<<"total_hits= "<<total_hits<<" total_ring_PEs= "<<total_ring_PEs<<" total_PMTs_hits= "<<total_PMTs_hits<<endl;
      cout<<"total_reco_Capture_E= "<<total_reco_Capture_E<<endl;

      //_______ calculate the total energy of parent-induced particles if you don't have a CCQE interaction _________
      if(fabs(mode)>=31 && fabs(mode)<=60){
	for(int j=0; j<npart; j++){//cout<<"-----------"<<endl;
	  if(part_parentid[j]==0){
	    neutrinoE+=part_KEstart[j];
	    cout<<"part_pid= "<<part_pid[j]<<" part_parentid= "<<part_parentid[j]<<" part_KEstart= "<<part_KEstart[j]<<endl;
	  }}
	cout<<"neutrinoE= "<<neutrinoE<<endl;
      }

      if(pot_length>0.){  cout<<"total_hits/pot_length= "<<total_hits/pot_length<<endl;  }
      //_____________________________________________________________________________________________________________
      hMuEn_HITS->Fill(trueKE, TMath::Log10(total_hits) );
      hMuEn2_HITS->Fill(trueKE, total_hits);
      hNeuEn_HITS->Fill(neutrinoE, TMath::Log10(total_hits) );
      hNeuEn2_HITS->Fill(neutrinoE, total_hits);
      if(pot_length>0.){
	hNeuEn_HITS_vs_pot_length->Fill(neutrinoE, TMath::Log10(total_hits)/pot_length );
	hNeuEn2_HITS_vs_pot_length->Fill(neutrinoE, total_hits/pot_length);    }

      if(pot_length<300.){
	hMuEn2_HITS_short_pot_length->Fill(trueKE, total_hits);
      }
      if(pot_length>=300.){
	hMuEn2_HITS_large_pot_length->Fill(trueKE, total_hits);
      }
      //----
      if(pot_length<300.){
	hNeuEn2_HITS_short_pot_length->Fill(neutrinoE, total_hits);
      }
      if(pot_length>=300.){
	hNeuEn2_HITS_large_pot_length->Fill(neutrinoE, total_hits);
      }
      if(lambda_max<300.){
	hNeuEn2_HITS_short_lambda_max->Fill(neutrinoE, total_hits);
      }
      if(lambda_max>=300.){
	hNeuEn2_HITS_large_lambda_max->Fill(neutrinoE, total_hits);
      }
      //_________ E vs total_ring_PEs ________
      //hMuEn_ring_PEs>Fill(trueKE, TMath::Log10(total_ring_PEs) );
      hMuEn2_ring_PEs->Fill(trueKE, total_ring_PEs);
      hNeuEn_ring_PEs->Fill(neutrinoE, TMath::Log10(total_ring_PEs) );
      hNeuEn2_ring_PEs->Fill(neutrinoE, total_ring_PEs);
      
      if(pot_length<300.){
	hMuEn2_ring_PEs_short_pot_length->Fill(trueKE, total_ring_PEs);
      }
      if(pot_length>=300.){
	hMuEn2_ring_PEs_large_pot_length->Fill(trueKE, total_ring_PEs);
      }
      //----
      if(pot_length<300.){
	hNeuEn2_ring_PEs_short_pot_length->Fill(neutrinoE, total_ring_PEs);
      }
      if(pot_length>=300.){
	hNeuEn2_ring_PEs_large_pot_length->Fill(neutrinoE, total_ring_PEs);
      }
      //___________total_PMTs_hits__________
      hMuEn2_total_PMTs_hits->Fill(trueKE, total_PMTs_hits);
      hNeuEn_total_PMTs_hits->Fill(neutrinoE, TMath::Log10(total_PMTs_hits) );
      hNeuEn2_total_PMTs_hits->Fill(neutrinoE, total_PMTs_hits);
      
      if(pot_length<300.){
	hMuEn2_total_PMTs_hits_short_pot_length->Fill(trueKE, total_PMTs_hits);
      }
      if(pot_length>=300.){
	hMuEn2_total_PMTs_hits_large_pot_length->Fill(trueKE, total_PMTs_hits);
      }
      //----
      if(pot_length<300.){
	hNeuEn2_total_PMTs_hits_short_pot_length->Fill(neutrinoE, total_PMTs_hits);
      }
      if(pot_length>=300.){
	hNeuEn2_total_PMTs_hits_large_pot_length->Fill(neutrinoE, total_PMTs_hits);
      }
      //____________________
      double hits_pot_length=0.;
      if(neutrinoE<5000.){
	if(lambda_max<-10.){lambda_max=0;} 
	if(pot_length>0.){ hits_pot_length = total_hits/pot_length;  }

	double recoE_lookup = recoEnergy[0];
	float total_PMTs_hits2 =total_PMTs_hits/4500.;
	float total_hits2      =total_hits/40000.;
	float total_ring_PEs2  =total_ring_PEs/20000.;
	float pot_length2      =pot_length/2500.;
	float hits_pot_length2 =hits_pot_length/150.; 
	float recoDWallR2      =recoDWallR/550.;
	float recoDWallZ2      =recoDWallZ/1200.;
	float lambda_max_2     =lambda_max/2500.;
	float recoDWall_2      =recoDWall/550.;
	float recoToWall_2     =recoToWall/2500.;
        float vtxTrackBias_2  = vtxTrackBias/550.;

	cout<<"check____: diffKE= "<<diffKE<<" recoE_lookup= "<<recoE_lookup<<" trueKE= "<<trueKE<<endl;

	nu_eneNEW->Branch("i", &i, "i/I");
	nu_eneNEW->Branch("neutrinoE", &neutrinoE, "neutrinoE/D");
	nu_eneNEW->Branch("trueKE", &trueKE, "trueKE/D");
	nu_eneNEW->Branch("recoE_lookup", &recoE_lookup, "recoE_lookup/D");
	//nu_eneNEW->Branch("total_PMTs_hits2", &total_PMTs_hits2, "total_PMTs_hits2/F");
	nu_eneNEW->Branch("total_hits2", &total_hits2, "total_hits2/F");
	nu_eneNEW->Branch("total_ring_PEs2", &total_ring_PEs2, "total_ring_PEs2/F");
	//nu_eneNEW->Branch("pot_length2", &pot_length2, "pot_length2/F");
	//nu_eneNEW->Branch("hits_pot_length2", &hits_pot_length2, "hits_pot_length2/F");
	nu_eneNEW->Branch("recoDWallR2", &recoDWallR2, "recoDWallR2/F");
	nu_eneNEW->Branch("recoDWallZ2", &recoDWallZ2, "recoDWallZ2/F");
	nu_eneNEW->Branch("lambda_max_2", &lambda_max_2, "lambda_max_2/F");
	//nu_eneNEW->Branch("recoDWall_2", &recoDWall_2, "recoDWall_2/F");
        //nu_eneNEW->Branch("recoToWall_2", &recoToWall_2, "recoToWall_2/F");
        //nu_eneNEW->Branch("vtxTrackBias_2", &vtxTrackBias_2, "vtxTrackBias_2/F");

	nu_eneNEW->Fill();
      }
      //____________________
    }//totalPEs>=10
  
  }//end of evts

  //---- for energy studies ----
  //TCanvas *MYEnY = new TCanvas("MYEnY");
  TProfile *amuEn= hMuEn->ProfileX();
  amuEn->SetLineColor(1);
  //amuEn->Draw();
  //TCanvas *MYEnY = new TCanvas("MYEnY2");
  TProfile *amuEn2= hMuEn2->ProfileX();
  amuEn2->SetLineColor(1);
  //amuEn2->Draw();
  TProfile *aneuEn= hNeuEn2->ProfileX();
  aneuEn->SetLineColor(1);
  
  TProfile *amuEn_HITS= hMuEn_HITS->ProfileX();
  amuEn_HITS->SetLineColor(1);
  TProfile *amuEn2_HITS= hMuEn2_HITS->ProfileX();
  amuEn2_HITS->SetLineColor(1);
  TProfile *aneuEn_HITS= hNeuEn_HITS->ProfileX();
  aneuEn_HITS->SetLineColor(1);
  TProfile *aneuEn2_HITS= hNeuEn2_HITS->ProfileX();
  aneuEn2_HITS->SetLineColor(1);

  //------- write --------
  nu_eneNEW->Write();
  
 /*plot_theta_muDir_track->Write();
  plot_vertical_dist->Write();
  plot_vertical_dist_sim->Write();
  amuEn->Write();
  hMuEn->Write();
  amuEn2->Write();
  hMuEn2->Write();
  hMuEn2_trueToWall->Write();
  hMuEn2_short_pot_length->Write();
  hMuEn2_large_pot_length->Write();
  
  hNeuEn2->Write();
  aneuEn->Write();
  hNeuEn2_short_pot_length->Write();
  hNeuEn2_large_pot_length->Write();
  
  hMuEn_HITS->Write();
  hMuEn2_HITS->Write();
  amuEn_HITS->Write();
  amuEn2_HITS->Write();
  hNeuEn_HITS->Write();
  hNeuEn2_HITS->Write();
  aneuEn_HITS->Write();
  aneuEn2_HITS->Write();
  hMuEn2_HITS_short_pot_length->Write();
  hMuEn2_HITS_large_pot_length->Write();
  hNeuEn2_HITS_short_pot_length->Write();
  hNeuEn2_HITS_large_pot_length->Write();
  hNeuEn_HITS_vs_pot_length->Write();
  hNeuEn2_HITS_vs_pot_length->Write();

  hNeuEn2_HITS_short_lambda_max->Write();
  hNeuEn2_HITS_large_lambda_max->Write();

  hMuEn_ring_PEs->Write();
  hMuEn2_ring_PEs->Write();
  //amuEn_ring_PEs->Write();
  //amuEn2_ring_PEs->Write();
  hNeuEn_ring_PEs->Write();
  hNeuEn2_ring_PEs->Write();
  //aneuEn_ring_PEs->Write();
  //aneuEn2_ring_PEs->Write();
  hMuEn2_ring_PEs_short_pot_length->Write();
  hMuEn2_ring_PEs_large_pot_length->Write();
  hNeuEn2_ring_PEs_short_pot_length->Write();
  hNeuEn2_ring_PEs_large_pot_length->Write();

  hMuEn_total_PMTs_hits->Write();
  hMuEn2_total_PMTs_hits->Write();
  //amuEn_ring_PEs->Write();
  //amuEn2_ring_PEs->Write();
  hNeuEn_total_PMTs_hits->Write();
  hNeuEn2_total_PMTs_hits->Write();
  //aneuEn_ring_PEs->Write();
  //aneuEn2_ring_PEs->Write();
  hMuEn2_total_PMTs_hits_short_pot_length->Write();
  hMuEn2_total_PMTs_hits_large_pot_length->Write();
  hNeuEn2_total_PMTs_hits_short_pot_length->Write();
  hNeuEn2_total_PMTs_hits_large_pot_length->Write();

  plot_trueToWall->Write();
  plot_trueDWall->Write();
  plot_trueToWallR->Write();
  plot_trueDWallR->Write();
  plot_trueToWallZ->Write();
  plot_trueDWallZ->Write();
  
  plot_pot_length->Write();
  */

}
