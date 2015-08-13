#include "TMBResidualClass.hh"

using namespace std;

void TMBResidualClass::Init()
{
  volume.Form("nv");

  particle.push_back(Form("gamma"));
  particle.push_back(Form("neutron"));
  
  time.push_back(Form("full"));
  time.push_back(Form("coin"));
  time.push_back(Form("prompt"));
  time.push_back(Form("delay"));
  time.push_back(Form("first"));
  time.push_back(Form("late"));

  lowtime.push_back(-2.e+4);
  lowtime.push_back(1.e+4);
  lowtime.push_back(0);
  lowtime.push_back(2.e+4);
  lowtime.push_back(2.e+4);
  lowtime.push_back(10.e+4);
  
  uptime.push_back(13.e+4);
  uptime.push_back(13.e+4);
  uptime.push_back(2.e+4);
  uptime.push_back(10.e+4);
  uptime.push_back(5.e+4);
  uptime.push_back(13.e+4);

  syntax1d.push_back(Form("gps"));
  syntax1d.push_back(Form("ene"));
  syntax1d.push_back(Form("small_ene"));

  syntax2d.push_back(Form("gps_ene"));
  syntax2d.push_back(Form("charge_vs_height_over_mult"));
  syntax2d.push_back(Form("cluster_charge_vs_total_s1"));

  multcut.push_back(Form("all custers"));
  multcut.push_back(Form("fail mult"));
  multcut.push_back(Form("pass multcut"));

}

void TMBResidualClass::ReadDataFile(int start, int end)
{
  string indir = GetInputdir();
  string middle="DST_Run";
  string last=".root";
  for(int i=start; i<=end; i++)
    {
      TString filename = Form("%s%06d%s",middle.c_str(),i,last.c_str());
      filename.Prepend(indir.c_str());
      if(!VerifyDataFile(filename))
	continue;
      else
	DSTtree->Add(filename);
    }
}

void TMBResidualClass::Load_TPCDSTTree(TPCDSTEvent &h)
{
  DSTtree->SetBranchAddress("runID",             &h.runID);
  DSTtree->SetBranchAddress("tpc_eventID",       &h.tpc_eventID);
  DSTtree->SetBranchAddress("tpc_event_type",    &h.tpc_event_type);
  DSTtree->SetBranchAddress("tpc_gps_fine",      &h.tpc_gps_fine);
  DSTtree->SetBranchAddress("tpc_gps_coarse",    &h.tpc_gps_coarse);
  DSTtree->SetBranchAddress("tpc_s1_start_time", &h.tpc_s1_start_time);
  DSTtree->SetBranchAddress("tpc_total_s1",      &h.tpc_total_s1);
  DSTtree->SetBranchAddress("tpc_total_f90",     &h.tpc_total_f90);
  DSTtree->SetBranchAddress("tpc_t_drift",       &h.tpc_t_drift);
  DSTtree->SetBranchAddress("tpc_s1_late",       &h.tpc_s1_late);
  DSTtree->SetBranchAddress("tpc_npulses",       &h.tpc_npulses);
  DSTtree->SetBranchAddress("tpc_timestamp",     &h.tpc_timestamp);

}

void TMBResidualClass::Load_ODDSTTree(ODDSTEvent &f)
{
  DSTtree->SetBranchAddress("od_eventID",&f.od_eventID);
  DSTtree->SetBranchAddress("od_gps_fine",&f.od_gps_fine);
  DSTtree->SetBranchAddress("od_gps_coarse",&f.od_gps_coarse);
  DSTtree->SetBranchAddress("od_timestamp", &f.od_timestamp);
  DSTtree->SetBranchAddress("od_nclusters", &f.od_nclusters);
  DSTtree->SetBranchAddress("od_wt_charge", &f.od_wt_charge);
  DSTtree->SetBranchAddress("od_cluster_charge", &f.od_cluster_charge);
  DSTtree->SetBranchAddress("od_cluster_start", &f.od_cluster_start); // time in ns
  DSTtree->SetBranchAddress("od_cluster_height", &f.od_cluster_height);
  DSTtree->SetBranchAddress("od_cluster_multiplicity", &f.od_cluster_multiplicity);
  DSTtree->SetBranchAddress("od_cluster_pass_multcut", &f.od_cluster_pass_multcut);
  DSTtree->SetBranchAddress("od_cluster_dtprompt", &f.od_cluster_dtprompt);
}

void TMBResidualClass::BookHistograms()
{
  tpc_ntuple = new TNtuple("tpc_ntuple","","time:s1:f90:tdrift:s1_late");
  tpc_hist2d.push_back(new TH2F("tpc_s1_f90_hist",";tpc_total_s1 [PE];tpc_total_f90",500,0,3500,500,0,1));
  
  for(size_t j=0; j<multcut.size(); j++){
    nv_mult1d.push_back(new TH1F(Form("%s_%s_charge_hist",volume.Data(),multcut.at(j).Data()),";od_cluster_charge[PE]",500,0,500));
    nv_mult2d.push_back(new TH2F(Form("%s_%s_charge_height_multcut_hist",volume.Data(),multcut.at(j).Data()),
				 ";od_cluster_charge[PE];od_cluster_height/od_cluster_multcut",500,0,3500,1000,0,1.e+9));
  }

  //define OD plots
  for(size_t i=0; i<particle.size(); i++){
    for(size_t j=0; j<time.size(); j++){
      p.push_back(new Plots);  
      p.back()->volume = volume;
      p.back()->particle = particle.at(i);
      p.back()->time = time.at(j);
      p.back()->lowtime = lowtime.at(j);
      p.back()->uptime = uptime.at(j);
      p.back()->syntax1d = syntax1d;
      p.back()->syntax2d = syntax2d;
      p.back()->BookHists();

      /*
      p.back()->ntuple = new TNtuple(Form("nv_%s_%s_ntuple",particle.at(i),time.at(j))," ","time:charge:height:multcut");
      for(size_t k=0; k<syntax1d.size(); k++){
	p.back()->syntax1d.push_back(syntax1d.at(k));
	p.back()->hist1d.push_back(new TH1F(Form("nv_%s_%s_%s_hist",particle.at(i),time.at(j),syntax1d.at(k))," ",500,0,3500));	
      }
      for(size_t k=0; k<syntax2d.size(); k++)
	{
	  p.back()->syntax2d.push_back(syntax2d.at(k));
	  p.back()->hist2d.push_back(new TH2F(Form("nv_%s_%s_%s_hist",particle.at(i),time.at(j),syntax2d.at(k))," ",500,0,3500,500,0,3500));	
	}
      */
    }
  }

}

bool SortFunction(double i, double j)
{
  return (TMath::Abs(i) < TMath::Abs(j));
}

void TMBResidualClass::LoopOverEvent(int start, int end)
{
  Init();
  size_t size = time.size();
  DSTtree = new TChain("DSTtree");
  ReadDataFile(start,end);
  int Entries = DSTtree->GetEntries();
  cout<<"Entries: "<<Entries<<endl;
  Load_TPCDSTTree(t);
  Load_ODDSTTree(e);
  BookHistograms();
  
  double afterpulse_cut = 50.;
  
  for(int i=0; i<Entries; i++){
    DSTtree->GetEntry(i);
    tpc_ntuple->Fill(t.tpc_s1_start_time,t.tpc_total_s1,t.tpc_total_f90,t.tpc_t_drift,t.tpc_s1_late);
    tpc_hist2d.at(0)->Fill(t.tpc_total_s1,t.tpc_total_f90);
    vector<bool> afterpulse_tag (particle.size()*time.size(),false);
    vector<int> counts (particle.size()*time.size(),0);
    int h=t.tpc_event_type;

    for(int j=0; j<e.od_nclusters; j++){
      if(e.od_cluster_charge->at(j)>0 && e.od_wt_charge<500){
	double gps = e.od_cluster_dtprompt->at(j)*1.e+3;//ns
	nv_mult1d.at(0)->Fill(e.od_cluster_charge->at(j));
	nv_mult2d.at(0)->Fill(e.od_cluster_charge->at(j),e.od_cluster_height->at(j)/e.od_cluster_multiplicity->at(j)*1.);	
	if(e.od_cluster_pass_multcut->at(j)==0){
	  nv_mult1d.at(1)->Fill(e.od_cluster_charge->at(j));
	  nv_mult2d.at(1)->Fill(e.od_cluster_charge->at(j),e.od_cluster_height->at(j)/e.od_cluster_multiplicity->at(j)*1.);		
	}
	if(e.od_cluster_pass_multcut->at(j)==1){
	  nv_mult1d.at(2)->Fill(e.od_cluster_charge->at(j));
	  nv_mult2d.at(2)->Fill(e.od_cluster_charge->at(j),e.od_cluster_height->at(j)/e.od_cluster_multiplicity->at(j)*1.);	       	
	  p.at(h*size)->ntuple->Fill(e.od_cluster_dtprompt->at(j)*1.e+3,e.od_cluster_charge->at(j),e.od_cluster_height->at(j),e.od_cluster_multiplicity->at(j));
	  p.at(h*size)->hist1d.at(0)->Fill(e.od_cluster_dtprompt->at(j)*1.e+3);
	  p.at(h*size)->hist1d.at(1)->Fill(e.od_cluster_charge->at(j));
	  p.at(h*size)->hist1d.at(2)->Fill(e.od_cluster_charge->at(j));
	  
	  p.at(h*size)->hist2d.at(0)->Fill(e.od_cluster_charge->at(j),e.od_cluster_dtprompt->at(j)*1.e+3);
	  p.at(h*size)->hist2d.at(1)->Fill(e.od_cluster_charge->at(j),e.od_cluster_height->at(j)/e.od_cluster_multiplicity->at(j)*1.0);
	  p.at(h*size)->hist2d.at(2)->Fill(e.od_cluster_charge->at(j),t.tpc_total_s1);
	  
	  for(size_t k=1; k<time.size(); k++){
	    if(gps>p.at(k+h*size)->lowtime-1.e+4 && gps<p.at(k+h*size)->lowtime && e.od_cluster_charge->at(j)>afterpulse_cut)
	      afterpulse_tag.at(k+h*size)=true;      
	  }
	}
      }
    }

    for(int j=0; j<e.od_nclusters; j++){
      if(e.od_cluster_charge->at(j)>0 && e.od_wt_charge<500 && e.od_cluster_pass_multcut->at(j)==1){
      	double gps = e.od_cluster_dtprompt->at(j)*1.e+3;//ns	
	for(size_t k=1; k<time.size(); k++){
	  if(gps>p.at(k+h*size)->lowtime && gps<p.at(k+h*size)->uptime && !afterpulse_tag.at(k+h*size)){
	    ++counts.at(k+h*size);
	    if(counts.at(k+h*size)==1){
	    p.at(k+h*size)->ntuple->Fill(e.od_cluster_dtprompt->at(j)*1.e+3,e.od_cluster_charge->at(j),
					 e.od_cluster_height->at(j),e.od_cluster_multiplicity->at(j));
	    p.at(k+h*size)->hist1d.at(0)->Fill(e.od_cluster_dtprompt->at(j)*1.e+3);
	    p.at(k+h*size)->hist1d.at(1)->Fill(e.od_cluster_charge->at(j));
	    p.at(k+h*size)->hist1d.at(2)->Fill(e.od_cluster_charge->at(j));
	
	    p.at(k+h*size)->hist2d.at(0)->Fill(e.od_cluster_charge->at(j),e.od_cluster_dtprompt->at(j)*1.e+3);
	    p.at(k+h*size)->hist2d.at(1)->Fill(e.od_cluster_charge->at(j),e.od_cluster_height->at(j)/e.od_cluster_multiplicity->at(j)*1.0);
	    p.at(k+h*size)->hist2d.at(2)->Fill(e.od_cluster_charge->at(j),t.tpc_total_s1);	    
	  }
	  }	  	  	  
	}
      }
    }      

  }//end of events

  SliceHistograms();
  SaveHistograms();

}

void TMBResidualClass::SliceHistograms()
{
  int Entries = DSTtree->GetEntries();
  canv.push_back(new TCanvas("mult1","Veto height/multiplicity Vs Charge",600,400));
  nv_mult2d.at(2)->SetMarkerColor(Colors(6));
  nv_mult2d.at(0)->Draw("colz");
  nv_mult2d.at(2)->Draw("same");

  TLegend *leg = new TLegend(0.47,0.7,0.87,0.9);
  canv.push_back(new TCanvas("mult2","Low Energy Spectrum",600,400));
  gPad->SetLogy();
  for(size_t j=0; j<nv_mult1d.size(); j++){
    nv_mult1d.at(j)->SetLineColor(Colors((int)j));
    leg->AddEntry(nv_mult1d.at(j),multcut.at(j).Data(),"l");
  }
  nv_mult1d.at(0)->Draw();
  nv_mult1d.at(1)->Draw("same");
  nv_mult1d.at(2)->Draw("same");
  leg->Draw();

  double pos1 = 0.;
  double pos2 = 50.;
  double pos3 = 150.;
  double pos4 = 310.;
  double height = p.at(0)->hist1d.at(2)->GetMaximum();
  TBox *box1 = new TBox(pos1,0,pos2,height);
  box1->SetLineColor(2);
  box1->SetLineWidth(1);

  TBox *box2 = new TBox(pos3,0,pos4,height);
  box2->SetLineColor(3);
  box2->SetLineWidth(1);

  for(size_t i=0; i<p.size(); i++){   
    for(size_t j=0; j<p.at(i)->hist1d.size(); j++){
      p.at(i)->hist1d.at(j)->Sumw2();
      p.at(i)->hist1d.at(j)->Scale(1./Entries);
      p.at(i)->hist1d.at(j)->SetLineColor(Colors(i));
    }
  }
  vector<double> capture_time;
  vector<double> tmb_fraction;
  for(size_t k=0; k<particle.size(); k++){   
    for(size_t j=0; j<time.size()-1; j++){   
      size_t i = k*time.size() + j;
      p.at(i)->hist1d.push_back(new TH1F(Form("%s_%s_%s_sub_hist",p.at(i)->volume.Data(),p.at(i)->particle.Data(),p.at(i)->time.Data()),
					 " ;Energy [PE]",125,0,500));  
      p.at(i)->hist1d.back()->Add(p.at(i)->hist1d.at(2),p.at((k+1)*time.size()-1)->hist1d.at(2),1,-1);
      canv.push_back(new TCanvas(Form("ene%d",i),"",1000,400));
      canv.back()->Divide(2,1);
      canv.back()->cd(1);
      gPad->SetLogy();
      p.at(i)->hist1d.at(2)->Draw();
      p.at((k+1)*time.size()-1)->hist1d.at(2)->Draw("same");
      box1->Draw();
      box2->Draw();
      canv.back()->cd(2);
      gPad->SetLogy();
      p.at(i)->hist1d.back()->Draw();
      box1->Draw();
      box2->Draw();    
    
      canv.push_back(new TCanvas(Form("fit%d",i),"",1000,400));
      p.at(i)->f1 = new TF1(Form("f%d",(int)i),this,&TMBResidualClass::ExpFit,p.at(i)->lowtime,p.at(i)->uptime,3);
      p.at(i)->f1->SetParNames("Time Offset[ns]","Capture Time[ns]","Constant");
      p.at(i)->f1->SetParameter(0,p.at(i)->lowtime);
      p.at(i)->f1->SetParameter(1,2.e+4);
      p.at(i)->f1->SetParameter(2,0.);
      p.at(i)->f1->SetNpx(200);
      p.at(i)->hist1d.at(0)->Fit(p.at(i)->f1,"RVEM");

      capture_time.push_back(p.at(i)->f1->GetParameter(1)/1.e+3);//us
      tmb_fraction.push_back(TMB_Concentration(capture_time.back()));
    }
  }
  
  for(size_t i=0; i<capture_time.size(); i++)
    cout<<"Capture Time: "<<capture_time.at(i)<<" [us] \t TMB Concentration: "<<tmb_fraction.at(i)<<endl;      
  
}

double TMBResidualClass::ExpFit(double *x,double *params)
{
  double result = TMath::Exp(-(x[0]-params[0])/params[1])+params[2];
  return result;
}

double TMBResidualClass::TMB_Concentration(double capture_time)
{
  vector<double> density{0.8761,0.932};
  vector<double> pc{9.,12.,0,0,0};
  vector<double> tmb{3.,9.,3.,1*0.199,1*0.801};
  vector<double> att_coeff;

  att_coeff.push_back(Attenutation_Coeff(density.at(0),pc));
  att_coeff.push_back(Attenutation_Coeff(density.at(1),tmb));

  TF1* TMB_Concen = new TF1("TMB_Concen",this,&TMBResidualClass::FitFunc,0,1,2);
  TMB_Concen->SetParNames("PC_Attenutation","TMB_Attenutation");
  TMB_Concen->SetParameter(0,att_coeff[0]);
  TMB_Concen->SetParameter(1,att_coeff[1]);
  TMB_Concen->SetNpx(900000);
  TMB_Concen->SetTitle("Thermal Neutron Capture Time vs TMB fraction");
  TMB_Concen->GetXaxis()->SetTitle("TMB fraction");
  TMB_Concen->GetYaxis()->SetTitle("Capture Time [us]");

  double result = TMB_Concen->GetX(capture_time);//207
  TCanvas *c100 = new TCanvas("c100","TMB concentrations",600,400);
  c100->SetLogx();
  TMB_Concen->Draw();
  c100->SaveAs("TMB_concentrations.png");

  const int entry = 11;
  //  double fraction[]={0.,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5};
  for(int i=0; i<entry; i++)
    {
      double fraction = i*0.05;
      cout<<"TMB fraction: "<<fraction<<"  Free Path: "<<TMB_Concen->Eval(fraction)<<endl;
    }
  return result;
}

double TMBResidualClass::Attenutation_Coeff(double density, vector<double> num_isotope)
{
  vector<double> cross_section{0.00339,0.334,1.9e-4,3837,0};
  vector<double> mass_isotope{12,1,16,10.,11.};
  vector<double> number_density;

  double NA = 6.022e-1;
  double sum = 0;
  float mass = 0;
  
  for(size_t i=0; i<mass_isotope.size(); i++){
    sum += num_isotope[i]*cross_section[i];
    mass += num_isotope[i]*mass_isotope[i];
    number_density.push_back(num_isotope[i]*density*6.022e+23);
  }
  double result = NA*density*sum/mass;
  cout<<result<<endl;
  return result;
}

double TMBResidualClass::FitFunc(double *x,double *params)
{  
  double frac = x[0];
  double velocity = 2.22e-1; // cm/s
  double path =  1./(params[0]*(1-frac)+params[1]*frac);
  double time =  path/velocity;
  return time;
}

void TMBResidualClass::SaveHistograms()
{
  string outputname = GetOutputdir() + GetOutFile();
  TFile f2D(outputname.c_str(), "RECREATE");
  for(size_t i=0; i<tpc_hist1d.size(); i++)
    tpc_hist1d.at(i)->Write();
  for(size_t i=0; i<tpc_hist2d.size(); i++)
    tpc_hist2d.at(i)->Write();
  tpc_ntuple->Write();
  for(size_t i=0; i<p.size(); i++){
    p.at(i)->ntuple->Write();
    for(size_t j=0; j<p.at(i)->hist1d.size(); j++)
      p.at(i)->hist1d.at(j)->Write();
    for(size_t j=0; j<p.at(i)->hist2d.size(); j++)
      p.at(i)->hist2d.at(j)->Write();
  }
  for(size_t i=0; i<canv.size(); i++)
    canv.at(i)->Write();

  f2D.Write();
  f2D.Close();
  cout<<"Successfully Saved the Histograms"<<endl;
}

int TMBResidualClass::Colors(int pos)
{
  vector<int> color;
  color.push_back(TColor::GetColor("#FF2007")); //red         0
  color.push_back(TColor::GetColor("#5A1DE8")); //violet      1
  color.push_back(TColor::GetColor("#000000"));
  color.push_back(TColor::GetColor("#F73CFF")); //pink        5
  color.push_back(TColor::GetColor("#1CFFDF")); //low green   7
  color.push_back(TColor::GetColor("#1485CC")); //blue        4
  color.push_back(TColor::GetColor("#FF791F")); //orange      6
  color.push_back(TColor::GetColor("#AF2FCC")); //dark pink   8
  color.push_back(TColor::GetColor("#E8A60C"));
  color.push_back(TColor::GetColor("#B26618"));
  color.push_back(TColor::GetColor("#79FFFF"));
  color.push_back(TColor::GetColor("#11FF8F"));
  color.push_back(TColor::GetColor("#59FF49")); //green       12
  
  int size = (int) color.size();
  return color.at(pos%size);
}

bool TMBResidualClass::VerifyDataFile(TString mainfile)
{
  ifstream NameCheck;
  NameCheck.open(mainfile.Data());
  if(!NameCheck.good())
    { NameCheck.close();
      return false;
    }
  else{
    TFile *f = new TFile(mainfile);
    if(!f->IsOpen() || f->IsZombie())
      {
	NameCheck.close();
	f->Close();
	return false;
      }
    else{
      cout<<"Processing the data file: "<<mainfile<<endl;
      NameCheck.close();
      f->Close();
      return true;
    }
  }
}

