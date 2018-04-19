// ************************************** //
// Sept28, 2017
// some current notes:
// Run16 AuAu 200 GeV:  P16ij on RCAS has no EmcalTriggers and no clusters
//     ::P16ij has triggers that need to be accessed a different way
//
//

#include <TSystem>

// basic STAR classes
class StMemStat;
class StMaker;
class StChain;
class StPicoDstMaker;
class StRefMultCorr;
// my added STAR classes
class StJetMakerTask;
//class StRho;
//class StRhoBase;
//class StRhoSparse;
class StJetFrameworkPicoBase;
class StMyAnalysisMaker;

// library and macro loading function
void LoadLibs();
void LoadMacros();

// constants
const double pi = 1.0*TMath::Pi();
Double_t ZVtxMin = -40.0;
Double_t ZVtxMax = 40.0;

bool usePrimaryTracks = kTRUE;
bool doCentSelection = kFALSE; //kTRUE;

StChain *chain;

//void readPicoDstQA(const Char_t *inputFile="test2.list", const Char_t *outputFile="test2.root", Int_t nEv = 10)
//void readPicoDstQA(const Char_t *inputFile="newPicoDsts.list", const Char_t *outputFile="newtest.root", Int_t nEv = 10)
void readPicoDstQA(const Char_t *inputFile="testLIST_Run14.list", const Char_t *outputFile="test.root", Int_t nEv = 10)
{
        Int_t nEvents = 10;
//        Int_t nEvents = 10000;
        if(nEv > 100) nEvents = 100000000;

        // Load necessary libraries and macros
        LoadLibs();
        LoadMacros();

        Int_t centralitySelection = StJetFrameworkPicoBase::kCent2050;
	
        // create chain
        StChain* chain = new StChain();

	// create the picoMaker maker
	//StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");
        StPicoDstMaker *picoMaker = new StPicoDstMaker(2,inputFile,"picoDst"); // updated Aug6th
        picoMaker->setVtxMode((int)(StPicoDstMaker::PicoVtxMode::Default));

        // QA task
        StPicoTrackClusterQA *Task = new StPicoTrackClusterQA("TrackClusterQA", kTRUE, outputFile);
        //StTEST *Task = new StTEST("StTEST", kTRUE, outputFile);
        Task->SetMinTrackPt(0.2);
        Task->SetMaxTrackPt(20.0);
        Task->SetMinClusterPt(0.2);
        Task->SetMaxClusterPt(500.0);
        Task->SetTrackPhiRange(0., 2.*pi);
        Task->SetTrackEtaRange(-1.0, 1.0);
        Task->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task->SetUsePrimaryTracks(usePrimaryTracks);
        //Task->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT1);  // kIsHT1 or kIsHT2
        Task->SetRunFlag(StJetFrameworkPicoBase::Run16_AuAu200);         // FIXME
        Task->SetCentralityDef(StJetFrameworkPicoBase::kgrefmult_P16id); // FIXME
        Task->SetTurnOnCentSelection(doCentSelection);  // run analysis for specific centrality
        Task->SetCentralityBinCut(centralitySelection); // specific centrality range to run
        Task->SetDebugLevel(8);

        // initialize chain
        chain->Init();
        cout<<"chain->Init();"<<endl;
        int total = picoMaker->chain()->GetEntries();
        cout << " Total entries = " << total << endl;
        if(nEvents>total) nEvents = total;
  
        for (Int_t i=0; i<nEvents; i++){
          if(i%100==0) cout << "Working on eventNumber " << i << endl;

          chain->Clear();
          int iret = chain->Make(i);	
          if (iret) { cout << "Bad return code!" << iret << endl; break;}

          total++;		
	}
	
	cout << "****************************************** " << endl;
	cout << "Work done... now its time to close up shop!"<< endl;
	cout << "****************************************** " << endl;
	chain->Finish();
	cout << "****************************************** " << endl;
	cout << "total number of events  " << nEvents << endl;
	cout << "****************************************** " << endl;
	
	delete chain;	

}

void LoadLibs()
{
  // load fastjet libraries 3.x
  gSystem->Load("$FASTJET/lib/libfastjet");
  gSystem->Load("$FASTJET/lib/libsiscone");
  gSystem->Load("$FASTJET/lib/libsiscone_spherical");
  gSystem->Load("$FASTJET/lib/libfastjetplugins");
  gSystem->Load("$FASTJET/lib/libfastjettools");
  gSystem->Load("$FASTJET/lib/libfastjetcontribfragile");

  // add include path to use its functionality
  gSystem->AddIncludePath("-I/global/homes/j/jmazer/STAR/Y2017/mytestinstalldir/FastJet/fastjet-install/include");

  // load the system libraries - these were defaults
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  // these are needed for new / additional classes
  gSystem->Load("libStPicoEvent");
  gSystem->Load("libStPicoDstMaker");

  // my libraries
  gSystem->Load("StRefMultCorr");
  //gSystem->Load("StJetFrameworkPicoBase");
  gSystem->Load("StMyAnalysisMaker");

  gSystem->ListLibraries();
} 

void LoadMacros()
{
}
