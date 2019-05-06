#include "TMath.h"
#include "TFile.h"
#define PI 3.1415926

// added "static" to function declaration on March31

int recenter_getAB()
{
  // switch
//  bool doBBC = kFALSE; // ZDC
  bool doBBC = kTRUE;  // BBC
  const char *method = "_Run14";
  
  // input calibration file with histograms
  //	TFile *f = new TFile("./Q_recentered.root");
  //TFile *f = new TFile("./after_recenter.root");
  //TFile *f = new TFile("./before_shift.root");
//  TFile *f = new TFile("./recenter_calib_file.root"); // DEFAULT
  // the bin# does not matter as it is only for TPC tracks not hit-based detectors like BBC and ZDC
//  TFile *f = new TFile("./recenter_calib_file_bin0_STEP1_Jan20.root");
//  TFile *f = new TFile("./recenter_calib_file_bin0_STEP1_Jan20.root");
//  TFile *f = new TFile("/star/u/jmazer19/Y2017/STAR/temp/recenter_calib_file_bin0_STEP1_April28.root");
//  TFile *f = new TFile("/star/u/jmazer19/Y2017/STAR/temp/recenter_calib_file_bin0_STEP1_May26.root");
  TFile *f = new TFile("/star/u/jmazer19/Y2017/STAR/temp/recenter_calib_file_bin0_STEP1_May29.root");
//  TFile *f = new TFile("/star/u/jmazer19/Y2017/STAR/temp/recenter_calib_file_STEP1_June19_R02.root");

  // output file
  if(doBBC){
    ofstream out(Form("bbc_recenter_data%s.h", method));          // BBC
  } else { ofstream out(Form("zdc_recenter_data%s.h", method)); } // ZDC

  //const int nRuns = 1359; // Run16 AuAu
  const int nRuns =  830; // Run14 AuAu
  double center_ex[nRuns] = {0};
  double center_ey[nRuns] = {0};
  double center_wx[nRuns] = {0};
  double center_wy[nRuns] = {0};
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
  for(int i=0; i<nRuns; i++) {
    center_ex[i] = ((TProfile *)f->Get(Tex->Data()))->GetBinContent(i+1);
    center_ey[i] = ((TProfile *)f->Get(Tey->Data()))->GetBinContent(i+1);
    center_wx[i] = ((TProfile *)f->Get(Twx->Data()))->GetBinContent(i+1);
    center_wy[i] = ((TProfile *)f->Get(Twy->Data()))->GetBinContent(i+1);
  }

  // add text to beginning of header to avoid linker problems
  if(doBBC){
    out<<(Form("#ifndef bbc_recenter_data%s_h", method))<<endl;   // BBC
    out<<(Form("#define bbc_recenter_data%s_h", method))<<endl;
  } else { 
    out<<(Form("#ifndef zdc_recenter_data%s_h", method))<<endl;   // ZDC
    out<<(Form("#define zdc_recenter_data%s_h", method))<<endl;
  }
  out<<endl;

  // write function for one side
  if(doBBC){
    out<<Form("static double bbc_center_ex%s[%i]=", method, nRuns)<<endl;          // BBC
  } else { out<<Form("static double zdc_center_ex%s[%i]=", method, nRuns)<<endl; } // ZDC
  out<<"{"<<endl;
  for(int i=0; i<nRuns; i++) {
    if(i==(nRuns-1)) { // no comma after last entry
      out<<center_ex[i]<<" "<<endl; 
    } else {  out<<center_ex[i]<<","<<endl; }
  }
  out<<"};"<<endl;

  // write function for one side
  if(doBBC){
    out<<Form("static double bbc_center_ey%s[%i]=", method, nRuns)<<endl;          // BBC
  } else { out<<Form("static double zdc_center_ey%s[%i]=", method, nRuns)<<endl; } // ZDC
  out<<"{"<<endl;
  for(int i=0; i<nRuns; i++) {
    if(i==(nRuns-1)) { // no comma after last entry
      out<<center_ey[i]<<" "<<endl;    
    } else {  out<<center_ey[i]<<","<<endl; }
  }
  out<<"};"<<endl;

  // write function for one side
  if(doBBC){
    out<<Form("static double bbc_center_wx%s[%i]=", method, nRuns)<<endl;          // BBC
  } else { out<<Form("static double zdc_center_wx%s[%i]=", method, nRuns)<<endl; } // ZDC
  out<<"{"<<endl;
  for(int i=0; i<nRuns; i++) {
    if(i==(nRuns-1)) { // no comma after last entry
      out<<center_wx[i]<<" "<<endl;    
    } else {  out<<center_wx[i]<<","<<endl; }
  }
  out<<"};"<<endl;

  // write function for one side
  if(doBBC){
    out<<Form("static double bbc_center_wy%s[%i]=", method, nRuns)<<endl;          // BBC
  } else { out<<Form("static double zdc_center_wy%s[%i]=", method, nRuns)<<endl; } // ZDC
  out<<"{"<<endl;
  for(int i=0; i<nRuns; i++) {
    if(i==(nRuns-1)) { // no comma after last entry
      out<<center_wy[i]<<" "<<endl;    
    } else {  out<<center_wy[i]<<","<<endl; }
  }
  out<<"};"<<endl;
  out<<"#endif"<<endl;

  cout<<"Finished writing file..."<<endl;
  // outfile.close();
  return 0;
}
