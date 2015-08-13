#ifndef _TMBResidualClass_HH
#define _TMBResidualClass_HH

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <math.h>
#include "TCut.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TChain.h"
#include "TRint.h"
#include "TArrayF.h"
#include "TMath.h"
#include "TString.h"
#include "TNtuple.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TPad.h"
#include "TStyle.h"
#include "TBox.h"

using namespace std;

struct TPCDSTEvent
{
  TPCDSTEvent():
    runID(-1),tpc_eventID(-1),tpc_event_type(-1),tpc_gps_fine(0),tpc_gps_coarse(0),
    tpc_s1_start_time(0),tpc_total_s1(0),tpc_total_f90(0),tpc_t_drift(0),tpc_s1_late(0),
    tpc_npulses(0),tpc_timestamp(0)
  {}

  Int_t    runID;
  Int_t    tpc_eventID;
  Int_t    tpc_event_type; // 0 gamma, 1 neutron, -1 else
  Double_t tpc_gps_fine; //clock cycles
  Double_t tpc_gps_coarse; //[s]
  Double_t tpc_s1_start_time; //[us]
  Double_t tpc_total_s1; //[PE]
  Double_t tpc_total_f90;
  Double_t tpc_t_drift; //[us]
  Double_t tpc_s1_late; //[PE]
  Int_t    tpc_npulses;
  Double_t tpc_timestamp;  //[us]      

};

struct ODDSTEvent
{
  ODDSTEvent():
    od_eventID(-1),od_nclusters(0),od_gps_fine(0),od_gps_coarse(0),od_timestamp(0),
    od_wt_charge(0),od_cluster_charge(0),od_cluster_start(0),od_cluster_height(0),
    od_cluster_multiplicity(0),od_cluster_pass_multcut(0),od_cluster_dtprompt(0)
  {}
  
  Int_t od_eventID;
  Int_t od_nclusters;
  Double_t od_gps_fine;
  Double_t od_gps_coarse;
  Double_t od_timestamp;
  Double_t od_wt_charge;
  vector<double>  *od_cluster_charge;
  vector<double>  *od_cluster_start;
  vector<double>  *od_cluster_height;
  vector<double>  *od_cluster_multiplicity;
  vector<int>     *od_cluster_pass_multcut;
  vector<double>  *od_cluster_dtprompt;

};

struct Plots{
  Plots():
    volume(""),particle(""),time(""),lowtime(0.),uptime(0.),
    ntuple(0),f1(0)
  {}
  
  TString volume;
  TString particle;
  TString time;
  double lowtime;
  double uptime;
  vector<TString> syntax1d;
  vector<TString> syntax2d;
  
  TNtuple *ntuple;
  TF1 *f1;
  vector<TH1F*> hist1d;
  vector<TH2F*> hist2d;  

  void BookHists()
  {
    ntuple = new TNtuple(Form("%s_%s_%s_ntuple",volume.Data(),particle.Data(),time.Data()),"","time:charge:height:multcut");    
    hist1d.push_back(new TH1F(Form("%s_%s_%s_%s_hist",volume.Data(),particle.Data(),time.Data(),syntax1d.at(0).Data())," ;time[ns]",100,lowtime,uptime));
    hist1d.push_back(new TH1F(Form("%s_%s_%s_%s_hist",volume.Data(),particle.Data(),time.Data(),syntax1d.at(1).Data())," ;od_cluster_charge[PE]",500,0,3500));
    hist1d.push_back(new TH1F(Form("%s_%s_%s_%s_hist",volume.Data(),particle.Data(),time.Data(),syntax1d.at(2).Data())," ;od_cluster_charge[PE]",125,0,500));
    
    hist2d.push_back(new TH2F(Form("%s_%s_%s_%s_hist",volume.Data(),particle.Data(),time.Data(),syntax2d.at(0).Data()),
			      " ;od_cluster_charge[PE];time[ns]",500,0,3500,100,lowtime,uptime));  
    hist2d.push_back(new TH2F(Form("%s_%s_%s_%s_hist",volume.Data(),particle.Data(),time.Data(),syntax2d.at(1).Data()),
			      " ;od_cluster_charge[PE];od_cluster_height/od_cluster_multcut",500,0,3500,1000,0,1.e+9));      
    hist2d.push_back(new TH2F(Form("%s_%s_%s_%s_hist",volume.Data(),particle.Data(),time.Data(),syntax2d.at(2).Data()),
			      " ;od_cluster_charge[PE];tpc_total_s1[PE]",500,0,3500,1000,0,40000));      
      
  }
  
};

class TMBResidualClass:public TObject
{
public:
  TMBResidualClass():
    DSTtree(0),tpc_ntuple(0)
  {}
  
  virtual ~TMBResidualClass(){}
  
  TPCDSTEvent t;
  ODDSTEvent e;
  vector<Plots*> p;

  void   SetInputdir(string val) { inputdir=val; }
  string GetInputdir() { return inputdir; }
  void   SetBGInputdir(string val) { bginputdir=val; }
  string GetBGInputdir() { return bginputdir; }
  void   SetOutputdir(string val) { outputdir=val; }
  string GetOutputdir() { return outputdir; }
  void   SetOutFile(string val) { outfile=val; }
  string GetOutFile() { return outfile; }
  void   SetBGFile(string val) { bgfile=val; }
  string GetBGFile() { return bgfile; }
  bool VerifyDataFile(TString);
  int Colors(int);

  void Init();
  void ReadDataFile(int, int);
  void Load_TPCDSTTree(TPCDSTEvent &);
  void Load_ODDSTTree(ODDSTEvent &);
  void BookHistograms();
  void LoopOverEvent(int,int);
  void SliceHistograms();
  void SaveHistograms();
  double ExpFit(double*,double*);
  double FitFunc(double*,double*);
  double TMB_Concentration(double);
  double Attenutation_Coeff(double, vector<double>);  
  
private:
  string inputdir;
  string bginputdir;
  string outputdir;
  string outfile;
  string bgfile;
  TChain *DSTtree;
  TString volume;
  vector<TString> particle;
  vector<TString> time;
  vector<double> uptime;
  vector<double> lowtime;
  vector<TString> syntax1d;
  vector<TString> syntax2d;
  vector<TCanvas*> canv;
  TNtuple *tpc_ntuple;
  vector<TH1F*> tpc_hist1d;
  vector<TH2F*> tpc_hist2d;
  vector<TH1F*> nv_mult1d;
  vector<TH2F*> nv_mult2d;
  vector<TString> multcut;

  ClassDef(TMBResidualClass,0);
};

#endif /* _TMBResidualClass_HH */
