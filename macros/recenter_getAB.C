#include "TMath.h"
#include "TFile.h"
#define PI 3.1415926

int recenter_getAB()
{
  // switch
  bool doBBC = kFALSE; // ZDC
//  bool doBBC = kTRUE;  // BBC

  // input calibration file with histograms
  //	TFile *f = new TFile("./Q_recentered.root");
  //TFile *f = new TFile("./after_recenter.root");
  //TFile *f = new TFile("./before_shift.root");
//  TFile *f = new TFile("./recenter_calib_file.root"); // DEFAULT
  // the bin# does not matter as it is only for TPC tracks not hit-based detectors like BBC and ZDC
  TFile *f = new TFile("./recenter_calib_file_bin0_STEP1_Jan20.root");

  // output file
  if(doBBC){
    ofstream out("bbc_recenter_data.h");          // BBC
  } else { ofstream out("zdc_recenter_data.h"); } // ZDC

  double center_ex[1359] = {0};
  double center_ey[1359] = {0};
  double center_wx[1359] = {0};
  double center_wy[1359] = {0};
  TString *Tex;
  TString *Tey;
  TString *Twx;
  TString *Twy; 

  // BBC
  if(doBBC) {
    Tex=new TString("hBBC_center_ex");
    Tey=new TString("hBBC_center_ey");
    Twx=new TString("hBBC_center_wx");
    Twy=new TString("hBBC_center_wy");
  } else { // ZDC
    Tex=new TString("hZDC_center_ex");
    Tey=new TString("hZDC_center_ey");
    Twx=new TString("hZDC_center_wx");
    Twy=new TString("hZDC_center_wy");
  }

  // loop over histograms
  for(int i=0; i<1359; i++) {
    center_ex[i] = ((TProfile *)f->Get(Tex->Data()))->GetBinContent(i+1);
    center_ey[i] = ((TProfile *)f->Get(Tey->Data()))->GetBinContent(i+1);
    center_wx[i] = ((TProfile *)f->Get(Twx->Data()))->GetBinContent(i+1);
    center_wy[i] = ((TProfile *)f->Get(Twy->Data()))->GetBinContent(i+1);
  }

  // write function for one side
  if(doBBC){
    out<<"double bbc_center_ex[1359]="<<endl;          // BBC
  } else { out<<"double zdc_center_ex[1359]="<<endl; } // ZDC
  out<<"{"<<endl;
  for(int i=0; i<1359; i++) {
    if(i==1358) { // no comma after last entry
      out<<center_ex[i]<<" "<<endl; 
    } else {  out<<center_ex[i]<<","<<endl; }
  }
  out<<"};"<<endl;

  // write function for one side
  if(doBBC){
    out<<"double bbc_center_ey[1359]="<<endl;          // BBC
  } else { out<<"double zdc_center_ey[1359]="<<endl; } // ZDC
  out<<"{"<<endl;
  for(int i=0; i<1359; i++) {
    if(i==1358) { // no comma after last entry
      out<<center_ey[i]<<" "<<endl;    
    } else {  out<<center_ey[i]<<","<<endl; }
  }
  out<<"};"<<endl;

  // write function for one side
  if(doBBC){
    out<<"double bbc_center_wx[1359]="<<endl;          // BBC
  } else { out<<"double zdc_center_wx[1359]="<<endl; } // ZDC
  out<<"{"<<endl;
  for(int i=0; i<1359; i++) {
    if(i==1358) { // no comma after last entry
      out<<center_wx[i]<<" "<<endl;    
    } else {  out<<center_wx[i]<<","<<endl; }
  }
  out<<"};"<<endl;

  // write function for one side
  if(doBBC){
    out<<"double bbc_center_wy[1359]="<<endl;          // BBC
  } else { out<<"double zdc_center_wy[1359]="<<endl; } // ZDC
  out<<"{"<<endl;
  for(int i=0; i<1359; i++) {
    if(i==1358) { // no comma after last entry
      out<<center_wy[i]<<" "<<endl;    
    } else {  out<<center_wy[i]<<","<<endl; }
  }
  out<<"};"<<endl;

  // outfile.close();
  return 0;
}
