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
class StJetFrameworkPicoBase;
class StMyAnalysisMaker;

// library and macro loading function
void LoadLibs();
void LoadMacros();

// constants
const double pi = 1.0*TMath::Pi();
Double_t ZVtxMin = -40.0;
Double_t ZVtxMax = 40.0;

bool usePrimaryTracks = kTRUE; //kTRUE;
bool doCentSelection = kFALSE; //kTRUE;

int RunYear = 4;
// kTRUE for local tests, kFALSE for job submission
bool doTEST = kFALSE;  //FIXME FIXME!!!! be aware before submission

StChain *chain;

//void readPicoDstQA(const Char_t *inputFile="test2.list", const Char_t *outputFile="test2.root", Int_t nEv = 10)
//void readPicoDstQA(const Char_t *inputFile="newPicoDsts.list", const Char_t *outputFile="newtest.root", Int_t nEv = 10)
void readPicoDstQA(const Char_t *inputFile="testLIST_Run14.list", const Char_t *outputFile="towerQA.root", Int_t nEv = 10, const Char_t *fEPoutJobappend="_this_is_a_test")
{
        Int_t nEvents = 10000;
//        Int_t nEvents = 100000;        
        if(nEv > 100) nEvents = 100000000;

        // run enumerators
        enum RunType_t {
          mJobSubmission = 0, // 0th spot - default unless set
          mRun11 = 1, mRun12 = 2, mRun13 = 3,
          mRun14 = 4, mRun15 = 5, mRun16 = 6,
          mRun17 = 7, mRun18 = 8
        };

        // Load necessary libraries and macros
        LoadLibs();
        LoadMacros();

        // =============================================================================== //
        //Int_t RunYear = mRun14;
        if((RunYear == mRun14) && doTEST) inputFile = "testLIST_Run14.list";
        if((RunYear == mRun16) && doTEST) inputFile = "test_run17061011.list";
        cout<<"inputFileName = "<<inputFile<<endl;

        // centrality global flags
        Int_t CentralitySelection = StJetFrameworkPicoBase::kCent2050;
        //Int_t CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P16id; // Run16
        Int_t CentralityDefinition = StJetFrameworkPicoBase::kgrefmult;  // Run14

        // Run/Event Flag
        //Int_t RunFlag = StJetFrameworkPicoBase::Run16_AuAu200;
        Int_t RunFlag = StJetFrameworkPicoBase::Run14_AuAu200;
        Int_t EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2;  // kIsHT1 or kIsHT2 or kIsHT3
        Int_t MBEventType = StJetFrameworkPicoBase::kVPDMB5; // this is default

        // event plane type
        //anaMaker->SetTPCEventPlaneMethod(StMyAnalysisMaker::kRemoveEtaPhiCone); // kRemoveEtaStrip is default
        //anaMaker->SetTPCEventPlaneMethod(StMyAnalysisMaker::kRemoveEtaPhiConeLeadSub);
        //anaMaker->SetTPCEventPlaneMethod(StMyAnalysisMaker::kRemoveEtaStrip); (train2)
        Int_t TPCEPSelectionType = StMyAnalysisMaker::kRemoveEtaPhiCone;
        Int_t EventPlaneTrackWeightMethod = StJetFrameworkPicoBase::kPtLinear2Const5Weight;

        // =============================================================================== //
        // open and close output .root file (so it exist and can be updated by Analysis Tasks)
        TFile *fout = new TFile(outputFile, "RECREATE");
        //fout->cd();
        fout->Close();
	
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
        Task->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        Task->SetRunFlag(RunFlag);                      // RunFlag
        Task->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task->SetTurnOnCentSelection(doCentSelection);  // run analysis for specific centrality
        Task->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        //Task->SetDebugLevel(8);
        //Task->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Task->SetHadronicCorrFrac(1.0);
        Task->SetDoTowerQAforHT(kFALSE);

        // QA task
        StPicoTrackClusterQA *Task2 = new StPicoTrackClusterQA("TrackClusterQAHT", kTRUE, outputFile);
        Task2->SetMinTrackPt(0.2);
        Task2->SetMaxTrackPt(20.0);
        Task2->SetMinClusterPt(0.2);
        Task2->SetMaxClusterPt(500.0);
        Task2->SetTrackPhiRange(0., 2.*pi);
        Task2->SetTrackEtaRange(-1.0, 1.0);
        Task2->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task2->SetUsePrimaryTracks(usePrimaryTracks);
        Task2->SetEmcTriggerEventType(EmcTriggerEventType);    // kIsHT1 or kIsHT2 or kIsHT3
        Task2->SetRunFlag(RunFlag);                      // RunFlag
        Task2->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task2->SetTurnOnCentSelection(doCentSelection);  // run analysis for specific centrality
        Task2->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        //Task2->SetDebugLevel(8);
        //Task2->SetDebugLevel(StPicoTrackClusterQA::kDebugEmcTrigger);
        Task2->SetHadronicCorrFrac(1.0);
        Task2->SetDoTowerQAforHT(kTRUE);


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

        // close output file if open
        if(fout->IsOpen()) fout->Close();
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
