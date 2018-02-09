#include "TMath.h"
#include "TFile.h"
#define PI 3.1415926

int bbc_shift_getAB_orig()
{
	//	TFile *f = new TFile("./Q_recentered.root");
	//TFile *f = new TFile("./after_recenter.root");
	TFile *f = new TFile("./before_shift.root");
	ofstream out("shift_data.h");
	const int N_c =20;
	double An[9][10][20] ={0};
	double Bn[9][10][20] ={0};
	TString *TA;
	TString *TB;

	for(int i=0;i<9;i++)
	{
		for(int j=0;j<10;j++)
		{
			TA=new TString("ShiftA");
			TB=new TString("ShiftB");
			*TA += i;
			*TA +="_";
			*TA += j;
			*TB += i;
			*TB +="_";
			*TB += j;
			cout<<*TA<<endl;
			for(int l=0;l<N_c;l++){
				An[i][j][l]=((TProfile *)f->Get(TA->Data()))->GetBinContent(l+1);
				Bn[i][j][l]=((TProfile *)f->Get(TB->Data()))->GetBinContent(l+1);
				cout<<*TA<<endl;
				cout<<An[i][j][l]<<endl;
			}	
			delete TA;
			delete TB;
		}
	}


	out<<"double shift_A[9][10][20]="<<endl;
//	out<<"double shift_A[9][10][10]="<<endl;
	out<<"{"<<endl;
	for(int i=0;i<9 ;i++){
		for(int j=0;j<10;j++){
			for(int l=0;l<N_c;l++){
				if(i==8 && j==9 && l==(N_c-1))
				{ out<<An[i][j][l]<<" "<<endl;}
				else{  out<<An[i][j][l]<<","<<endl;}
			}
		}
	}
	out<<"};"<<endl;

	out<<"double shift_B[9][10][20]="<<endl;
//	out<<"double shift_B[9][10][10]="<<endl;
	out<<"{"<<endl;
	for(int i=0;i<9 ;i++){
		for(int j=0;j<10;j++){
			for(int l=0;l<N_c;l++){
				if(i==8&& j==9 && l==(N_c-1) )
				{ out<<Bn[i][j][l]<<" "<<endl;}
				else{  out<<Bn[i][j][l]<<","<<endl;}
			}
		}
	}
	out<<"};"<<endl;
	// outfile.close();
	return 0;
}

