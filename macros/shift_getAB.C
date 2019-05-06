#include "TMath.h"
#include "TFile.h"
#define PI 3.1415926

// March31:  check function names... found incomplete editing from weeks prior

int shift_getAB()
{
  // methods
  int Method = 2;
  //Method = 2;  // TPC shift
  //Method = 3;  // BBC shift
  //Method = 4;  // ZDC shift

  // pt bin methods
  int ptbin = 4;  // use -99 for standard, and 0, 1, 2, 3, 4, (5,6,7) for pt assoc bin
  bool appendDate = kTRUE;

  // append bin for pt assoc bins to file name
  const char *funcBin = "";
  if(ptbin > -1) { funcBin = Form("_bin%i", ptbin); }

  // append date to file name
  const char *funcDate = "";
  //if(appendDate) { funcDate = Form("_STEP2_Jan20"); } // CHANGE NAME case-by-case
  //if(appendDate) { funcDate = Form("_STEP2_Feb5"); } // CHANGE NAME case-by-case
  //if(appendDate) { funcDate = Form("_STEP2_March31"); } // CHANGE NAME case-by-case
  //if(appendDate) { funcDate = Form("_STEP2_April28"); } // CHANGE NAME case-by-case
  //if(appendDate) { funcDate = Form("_STEP2_May15"); } // CHANGE NAME case-by-case
  //if(appendDate) { funcDate = Form("_STEP2_May26"); } // CHANGE NAME case-by-case
  //if(appendDate) { funcDate = Form("_STEP2_May30"); } // CHANGE NAME case-by-case
  //if(appendDate) { funcDate = Form("_STEP2_June7ch"); } // CHANGE NAME case-by-case // charged jets
  //if(appendDate) { funcDate = Form("_STEP2_June22_R02"); } // CHANGE NAME case-by-case // R=0.2 jets
  //if(appendDate) { funcDate = Form("_STEP2_August29"); } // CHANGE NAME case-by-case
  //if(appendDate) { funcDate = Form("_STEP2_Oct15_R03"); } // CHANGE NAME case-by-case // R=0.3 jets
  if(appendDate) { funcDate = Form("_STEP2_Dec5_R03"); } // CHANGE NAME case-by-case // R=0.3 jets

  // method of jet removal
  // "Method1": kRemoveEtaStrip
  // "Method2": kRemoveEtaPhiCone
  // "Method3": kRemoveLeadingJetConstituents
  // "Method5": kRemoveEtaPhiConeLeadSub
  // "Method6": kRemoveLeadingSubJetConstituents
  //JetMethod = "_Method2ch"; // charged jets
  const char *JetMethod = "_Method1_R03";
  //const char *JetMethod = "_Method1";
  const char *RunName = "_Run14_NEW";

  // print append strings
  cout<<"funcbin = "<<funcBin<<endl;
  cout<<"funcDate = "<<funcDate<<endl;
  cout<<"JetMethod = "<<JetMethod<<endl;
  cout<<"RunName = "<<RunName<<endl;

  // input calibration file with histograms
  //	TFile *f = new TFile("./Q_recentered.root");
  //TFile *f = new TFile("./after_recenter.root");
  //TFile *f = new TFile("./before_shift.root");
  //TFile *f = new TFile("./shift_calib_file.root"); // default
  //TFile *f = new TFile(Form("./shift_calib_file%s%s.root", funcBin, funcDate));
  TFile *f = new TFile(Form("/star/u/jmazer19/Y2017/STAR/temp/shift_calib_file%s%s.root", funcBin, funcDate));

  // output file
  if(Method == 2) { ofstream out(Form("tpc_shift_data%s%s%s.h", funcBin, JetMethod, RunName)); }
  if(Method == 3) { ofstream out(Form("bbc_shift_data%s.h", RunName)); } // default
  if(Method == 4) { ofstream out(Form("zdc_shift_data%s.h", RunName)); } // default

  // A = m, B = p: for TPC
  const int N_c =20;
  double An[9][20][20] = {0};
  double Bn[9][20][20] = {0};
  TString *TA;
  TString *TB;

  // loop over histograms
  for(int i=0; i<9; i++) {
    for(int j=0; j<20; j++) {
      if(Method == 2) {
        TA=new TString("hTPC_shift_N");  // N (minus = neg, N)
        TB=new TString("hTPC_shift_P");
        //TA=new TString("hTPC_shift_A"); // new
        //TB=new TString("hTPC_shift_B"); // new
      } else if(Method == 3) {
        TA=new TString("hBBC_shift_A");
        TB=new TString("hBBC_shift_B");
      } else if(Method == 4) {
        TA=new TString("hZDC_shift_A");
        TB=new TString("hZDC_shift_B");
      }

      *TA += i;
      *TA +="_";
      *TA += j;
      *TB += i;
      *TB +="_";
      *TB += j;
      //cout<<*TA<<endl;
      for(int l=0; l<N_c; l++) {
        An[i][j][l]=((TProfile *)f->Get(TA->Data()))->GetBinContent(l+1);
        Bn[i][j][l]=((TProfile *)f->Get(TB->Data()))->GetBinContent(l+1);
        //cout<<*TA<<endl;
        cout<<An[i][j][l]<<endl;
      }	
      delete TA;
      delete TB;
    }
  }

  // =======================================
  // added Feb6, edited April28
  if(Method == 2) { 
    out<<(Form("#ifndef tpc_shift_data%s%s%s_h", funcBin, JetMethod, RunName))<<endl; 
    out<<(Form("#define tpc_shift_data%s%s%s_h", funcBin, JetMethod, RunName))<<endl;
  }
  if(Method == 3) { 
    out<<(Form("#ifndef bbc_shift_data%s_h", RunName))<<endl;
    out<<(Form("#define bbc_shift_data%s_h", RunName))<<endl;
  } 
  if(Method == 4) { 
    out<<(Form("#ifndef zdc_shift_data%s_h", RunName))<<endl;
    out<<(Form("#define zdc_shift_data%s_h", RunName))<<endl;
  } 
  out<<endl;

/*
  out<<"class X {"<<endl;
  out<<endl;
  out<<"public:"<<endl;

  // write function for one side  shift_Qpx
  if(Method == 2) { out<<Form("static constexpr double tpc_shift_N%s%s[9][20][20]=", funcBin, JetMethod)<<endl; }
  if(Method == 3) { out<<"static constexpr double bbc_shift_A[9][20][20]="<<endl; }
  if(Method == 4) { out<<"static constexpr double zdc_shift_A[9][20][20]="<<endl; }
*/
  // ======================================

  // write function for one side  shift_Qpx
  if(Method == 2) { out<<Form("static double tpc_shift_N%s%s%s[9][20][20]=", funcBin, JetMethod, RunName)<<endl; }
  //if(Method == 2) { out<<Form("static double tpc_shift_A%s%s%s[9][20][20]=", funcBin, JetMethod, RunName)<<endl; } // new
  if(Method == 3) { out<<Form("static double bbc_shift_A%s[9][20][20]=", RunName)<<endl; }
  if(Method == 4) { out<<Form("static double zdc_shift_A%s[9][20][20]=", RunName)<<endl; }

  out<<"{"<<endl;
  for(int i=0; i<9; i++) {
    for(int j=0; j<20; j++) {
      for(int l=0; l<N_c; l++) {
        if(i==8 && j==19 && l==(N_c-1)) { // no comma after last entry
          out<<An[i][j][l]<<" "<<endl; 
        } else {  out<<An[i][j][l]<<","<<endl; }
      }
    }
  }
  out<<"};"<<endl;

/*
  // write function for other side
  if(Method == 2) { out<<Form("static double tpc_shift_P%s%s%s[9][20][20]=", funcBin, JetMethod, RunName)<<endl; }
  //if(Method == 2) { out<<Form("static double tpc_shift_B%s%s%s[9][20][20]=", funcBin, JetMethod, RunName)<<endl; } // new
  if(Method == 3) { out<<Form("static double bbc_shift_B%s[9][20][20]=", RunName)<<endl; }
  if(Method == 4) { out<<Form("static double zdc_shift_B%s[9][20][20]=", RunName)<<endl; }
*/

//  if(Method == 2) { out<<Form("static constexpr double tpc_shift_P%s%s%s[9][20][20]=", funcBin, JetMethod, RunName)<<endl; }
//  if(Method == 3) { out<<Form("static constexpr double bbc_shift_B%s[9][20][20]=", RunName)<<endl; }
//  if(Method == 4) { out<<Form("static constexpr double zdc_shift_B%s[9][20][20]=", RunName)<<endl; }

  if(Method == 2) { out<<Form("static double tpc_shift_P%s%s%s[9][20][20]=", funcBin, JetMethod, RunName)<<endl; }
  if(Method == 3) { out<<Form("static double bbc_shift_B%s[9][20][20]=", RunName)<<endl; }
  if(Method == 4) { out<<Form("static double zdc_shift_B%s[9][20][20]=", RunName)<<endl; }

  out<<"{"<<endl;
  for(int i=0; i<9; i++) {
    for(int j=0; j<20; j++) {
      for(int l=0; l<N_c; l++) {
        if(i==8 && j==19 && l==(N_c-1)) { // no comma after last entry
          out<<Bn[i][j][l]<<" "<<endl; 
        } else {  out<<Bn[i][j][l]<<","<<endl; }
      }
    }
  }
  out<<"};"<<endl;
  out<<"#endif"<<endl;

  // outfile.close();
  return 0;
}
