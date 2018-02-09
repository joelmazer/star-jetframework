#include "TMath.h"
#include "TFile.h"
#define PI 3.1415926

int tpc_recenter_getNP()
{
  // pt bin methods
  int ptbin = 4;  // use -99 for standard, and 0, 1, 2, 3, 4, (5,6,7) for pt assoc bin
  bool appendDate = kTRUE;

  // append bin for pt assoc bins to file name
  const char *funcBin = "";
  if(ptbin > -1) { funcBin = Form("_bin%i", ptbin); }

  // append date to file name:  TODO - update every time!
  const char *funcDate = "";
  //if(appendDate) { funcDate = Form("_STEP1_Jan20"); }
  if(appendDate) { funcDate = Form("_STEP1_Feb1"); }

  cout<<"funcbin = "<<funcBin<<endl;
  cout<<"funcDate = "<<funcDate<<endl;

  // method of jet removal
  // "Method1": kRemoveEtaStrip
  // "Method2": kRemoveEtaPhiCone
  // "Method3": kRemoveLeadingJetConstituents
  // "Method5": kRemoveEtaPhiConeLeadSub
  // "Method6": kRemoveLeadingSubJetConstituents
  const char *JetMethod = "";
  JetMethod = "_Method1";

  // input calibration file with histograms
  //	TFile *f = new TFile("./Q_recentered.root");
  //TFile *f = new TFile("./after_recenter.root");
  //TFile *f = new TFile("./before_shift.root");
  TFile *f = new TFile(Form("./recenter_calib_file%s%s.root", funcBin, funcDate)); 

  // output file
  ofstream out(Form("tpc_recenter_data%s%s.h", funcBin, JetMethod));

  // A = m, B = p: for TPC
  double Nx[9][20] = {0};
  double Px[9][20] = {0};
  double Ny[9][20] = {0};
  double Py[9][20] = {0};
  TString *Tn;
  TString *Tp;

  // loop over histograms
  for(int i=0; i<9; i++) {
    for(int j=0; j<20; j++) {
      Tn=new TString("Q2_m"); // m (minus = neg, N)
      Tp=new TString("Q2_p"); // p (positive)

      *Tn += i;
      *Tn +="_";
      *Tn += j;
      *Tp += i;
      *Tp +="_";
      *Tp += j;

      Nx[i][j]=((TProfile *)f->Get(Tn->Data()))->GetBinContent(1);
      Ny[i][j]=((TProfile *)f->Get(Tn->Data()))->GetBinContent(2);
      Px[i][j]=((TProfile *)f->Get(Tp->Data()))->GetBinContent(1);
      Py[i][j]=((TProfile *)f->Get(Tp->Data()))->GetBinContent(2);
      cout<<*Tn<<endl;
      cout<<Nx[i][j]<<endl;

      delete Tn;
      delete Tp;
    }
  }

  // write function for one side:  center_Qnx
  out<<Form("double tpc_center_Qnx%s%s[9][20]=", funcBin, JetMethod)<<endl;
  out<<"{"<<endl;
  for(int i=0; i<9; i++) {
    for(int j=0; j<20; j++) {
      if(i==8 && j==19) { // no comma after last entry
        out<<Nx[i][j]<<" "<<endl; 
      } else {  out<<Nx[i][j]<<","<<endl; }
    }
  }
  out<<"};"<<endl;

  // write function for one side:  center_Qny
  out<<Form("double tpc_center_Qny%s%s[9][20]=", funcBin, JetMethod)<<endl;
  out<<"{"<<endl;
  for(int i=0; i<9; i++) {
    for(int j=0; j<20; j++) {
      if(i==8 && j==19) { // no comma after last entry
        out<<Ny[i][j]<<" "<<endl;     
      } else {  out<<Ny[i][j]<<","<<endl; }
    }
  }
  out<<"};"<<endl;

  // write function for one side:  center_Qpx
  out<<Form("double tpc_center_Qpx%s%s[9][20]=", funcBin, JetMethod)<<endl;
  out<<"{"<<endl;
  for(int i=0; i<9; i++) {
    for(int j=0; j<20; j++) {
      if(i==8 && j==19) { // no comma after last entry
        out<<Px[i][j]<<" "<<endl;     
      } else {  out<<Px[i][j]<<","<<endl; }
    }
  }
  out<<"};"<<endl;

  // write function for one side:  center_Qpy
  out<<Form("double tpc_center_Qpy%s%s[9][20]=", funcBin, JetMethod)<<endl;
  out<<"{"<<endl;
  for(int i=0; i<9; i++) {
    for(int j=0; j<20; j++) {
      if(i==8 && j==19) { // no comma after last entry
        out<<Py[i][j]<<" "<<endl;     
      } else {  out<<Py[i][j]<<","<<endl; }
    }
  }
  out<<"};"<<endl;

  // outfile.close();
  return 0;
}
