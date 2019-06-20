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

bool doCentSelection = kFALSE; //kTRUE;
bool dopp = kFALSE; // FIXME kTRUE for pp data

int RunYear = 14;
// kTRUE for local tests, kFALSE for job submission
bool doTEST = kFALSE;  //FIXME FIXME!!!! be aware before submission

StChain *chain;

//__________________________________________________________________________________________________________
void readPicoDstQA(const Char_t *inputFile="Run_15164046_files.list", const Char_t *outputFile="towerQA.root", Int_t nEv = 10, const Char_t *fEPoutJobappend="_this_is_a_test")
{
//        Int_t nEvents = 10;
//        Int_t nEvents = 1000;
        Int_t nEvents = 10000;
//        Int_t nEvents = 50000;
//        Int_t nEvents = 100000;        
        if(nEv > 100) nEvents = 100000000;

        // run enumerators
        enum RunType_t {
          mJobSubmission = 0, // 0th spot - default unless set
          mRun07 =  7, mRun08 = 8,  mRun09 = 9, mRun10 = 10,
          mRun11 = 11, mRun12 = 12, mRun13 = 13,
          mRun14 = 14,
          mRun15 = 15, mRun16 = 16, mRun17 = 17, mRun18 = 18
        };

        // Load necessary libraries and macros
        LoadLibs();
        LoadMacros();

        // =============================================================================== //
        // =============================================================================== //
        // =============================================================================== //
        // over-ride functions
        if(dopp) {
          doCentSelection = kFALSE;  // can't ask for a particular centrality if requesting pp collisions
        }

        // input file for tests (based on Run) - updated for new Runs as needed
        if((RunYear == mRun12) && doTEST && dopp) inputFile = "testLIST_Run12pp.list";
        if((RunYear == mRun14) && doTEST) inputFile = "Run_15164046_files.list"; //"Run_15151042_files.list"; //"testLIST_Run14.list";
        if((RunYear == mRun16) && doTEST) inputFile = "test_run17124003_files.list";
        if((RunYear == mRun17) && doTEST && dopp) inputFile = "Run17pp_510GeV.list"; // "filelist_pp2017.list";
        cout<<"inputFileName = "<<inputFile<<endl;

        // centrality global flags - no centrality for pp collisions
        Int_t CentralitySelection = StJetFrameworkPicoBase::kCent2050;
        Int_t CentralityDefinition;
        if(RunYear == mRun12) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30;  // no centrality defintion for Run 12 pp 
        //if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult;                 // Run14 P16id ??FIXME???
        //if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30; // Run14 P17id (NEW - from Nick Oct 23)
        if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30; // Run14 P18ih (NEW - from Nick June10, 2019)
        if(RunYear == mRun16) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P16id;           // Run16
        //if(RunYear == mRun16) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_VpdMBnoVtx;
        cout<<"Centrality definition: "<<CentralityDefinition<<endl;

        // Run/Event Flag
        Int_t RunFlag;
        if(RunYear == mRun14) RunFlag = StJetFrameworkPicoBase::Run14_AuAu200;
        if(RunYear == mRun16) RunFlag = StJetFrameworkPicoBase::Run16_AuAu200;
        if(RunYear == mRun11 && dopp) RunFlag = StJetFrameworkPicoBase::Run11_pp500;
        if(RunYear == mRun12 && dopp) RunFlag = StJetFrameworkPicoBase::Run12_pp200;
        if(RunYear == mRun13 && dopp) RunFlag = StJetFrameworkPicoBase::Run13_pp510;
        if(RunYear == mRun17 && dopp) RunFlag = StJetFrameworkPicoBase::Run17_pp510;
        Bool_t RejectBadRuns = kFALSE; // switch to load and than omit bad runs

        // trigger flags - update default
        Int_t EmcTriggerEventType; // kIsHT1 or kIsHT2 or kIsHT3
        if(RunYear == mRun12) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2;
        if(RunYear == mRun14) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2; // kIsHT2 Run14
        if(RunYear == mRun16) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT1; // kIsHT1 Run16
        if(RunYear == mRun17) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT3;
        Int_t MBEventType = StJetFrameworkPicoBase::kVPDMB5;        // this is default
        if(RunYear == mRun12) MBEventType = StJetFrameworkPicoBase::kRun12main; // default for Run12 pp
        if(RunYear == mRun17) MBEventType = StJetFrameworkPicoBase::kVPDMB; // default for Run17 pp
        //Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerHT;    // kTriggerANY, kTriggerMB, kTriggerHT
        Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerANY;    // kTriggerANY, kTriggerMB, kTriggerHT     FIXME
        Int_t TowerListToUse = 136; // Run14 136-122: jet-hadron, early jet shape // been using 51,  3, 79   - doesn't matter for charged jets
        if(dopp) TowerListToUse = 169;
        // Run12: 1 - Raghav's list, 102 - my initial list, 169 - new list
        // Run14: 136 - main list (updated version for AuAu 200 GeV Run14), 122 - past used list
        // Run14 P18ih: 999 (initial) 

        // track flags
        bool usePrimaryTracks;
        if(RunYear == mRun12) usePrimaryTracks = kTRUE;
        if(RunYear == mRun14) usePrimaryTracks = kTRUE;  // = kTRUE for Run14, kFALSE for Run16
        if(RunYear == mRun16) usePrimaryTracks = kFALSE;
        if(RunYear == mRun17) usePrimaryTracks = kTRUE;

        // event plane type
        //anaMaker->SetTPCEventPlaneMethod(StMyAnalysisMaker::kRemoveEtaPhiCone); // kRemoveEtaStrip is default
        //anaMaker->SetTPCEventPlaneMethod(StMyAnalysisMaker::kRemoveEtaPhiConeLeadSub);
        //anaMaker->SetTPCEventPlaneMethod(StMyAnalysisMaker::kRemoveEtaStrip); (train2)
        Int_t TPCEPSelectionType = StMyAnalysisMaker::kRemoveEtaPhiCone;
        Int_t EventPlaneTrackWeightMethod = StJetFrameworkPicoBase::kPtLinear2Const5Weight;

        // update settings for new centrality definitions
        if(CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30 ||
          CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30
        ) { ZVtxMin = -30.0; ZVtxMax = 30.0; }

        // =============================================================================== //
        // =============================================================================== //
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

        // QA task - MB events
        StPicoTrackClusterQA *Task = new StPicoTrackClusterQA("TrackClusterQA", kTRUE, outputFile);
        Task->SetTrackPtRange(0.2, 30.0);
        Task->SetTrackPhiRange(0., 2.0*pi);
        Task->SetTrackEtaRange(-1.0, 1.0);
        Task->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task->SetClusterPtRange(0.2, 100.0);
        Task->SetTowerERange(0.2, 100.0);
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
        Task->SetdoppAnalysis(dopp);
        Task->SetUseBBCCoincidenceRate(kFALSE);
        Task->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs

        // QA task - HT events
        StPicoTrackClusterQA *Task2 = new StPicoTrackClusterQA("TrackClusterQAHT", kTRUE, outputFile);
        Task2->SetTrackPtRange(0.2, 30.0);
        Task2->SetTrackPhiRange(0., 2.0*pi);
        Task2->SetTrackEtaRange(-1.0, 1.0);
        Task2->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task2->SetClusterPtRange(0.2, 100.0);
        Task2->SetTowerERange(0.2, 100.0);
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
        Task2->SetdoppAnalysis(dopp);
        Task2->SetUseBBCCoincidenceRate(kFALSE);
        Task2->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs

/*
        // QA task - HT1
        StPicoTrackClusterQA *Task1 = new StPicoTrackClusterQA("TrackClusterQAHT1", kTRUE, outputFile);
        Task1->SetTrackPtRange(0.2, 30.0);
        Task1->SetTrackPhiRange(0., 2.0*pi);
        Task1->SetTrackEtaRange(-1.0, 1.0);
        Task1->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task1->SetClusterPtRange(0.2, 100.0);
        Task1->SetTowerERange(0.2, 100.0);
        Task1->SetUsePrimaryTracks(usePrimaryTracks);
        Task1->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT1);    // kIsHT1 or kIsHT2 or kIsHT3
        Task1->SetRunFlag(RunFlag);                      // RunFlag
        Task1->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task1->SetTurnOnCentSelection(doCentSelection);  // run analysis for specific centrality
        Task1->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        Task1->SetHadronicCorrFrac(1.0);
        Task1->SetDoTowerQAforHT(kTRUE);
        Task1->SetdoppAnalysis(dopp);
        Task1->SetUseBBCCoincidenceRate(kFALSE);
        Task1->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs

        // QA task - HT2
        StPicoTrackClusterQA *Task2 = new StPicoTrackClusterQA("TrackClusterQAHT2", kTRUE, outputFile);
        Task2->SetTrackPtRange(0.2, 30.0);
        Task2->SetTrackPhiRange(0., 2.0*pi);
        Task2->SetTrackEtaRange(-1.0, 1.0);
        Task2->SetEventZVtxRange(ZVtxMin, ZVtxMax);      // can be tighter for Run16 (-20,20)
        Task2->SetClusterPtRange(0.2, 100.0);
        Task2->SetTowerERange(0.2, 100.0);
        Task2->SetUsePrimaryTracks(usePrimaryTracks);
        Task2->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT2);    // kIsHT1 or kIsHT2 or kIsHT3
        Task2->SetRunFlag(RunFlag);                      // RunFlag
        Task2->SetCentralityDef(CentralityDefinition);   // FIXME - not needed for Run14 - why?
        Task2->SetTurnOnCentSelection(doCentSelection);  // run analysis for specific centrality
        Task2->SetCentralityBinCut(CentralitySelection); // specific centrality range to run
        Task2->SetHadronicCorrFrac(1.0);
        Task2->SetDoTowerQAforHT(kTRUE);
        Task2->SetdoppAnalysis(dopp);
        Task2->SetUseBBCCoincidenceRate(kFALSE);
        Task2->SetRejectBadRuns(RejectBadRuns);          // switch to load and than omit bad runs
*/

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
  //gSystem->Load("libCGAL"); - not installed 
  gSystem->Load("$FASTJET/lib/libfastjet");
  gSystem->Load("$FASTJET/lib/libsiscone");
  gSystem->Load("$FASTJET/lib/libsiscone_spherical");
  gSystem->Load("$FASTJET/lib/libfastjetplugins");
  gSystem->Load("$FASTJET/lib/libfastjettools");
  gSystem->Load("$FASTJET/lib/libfastjetcontribfragile");

  // add include path to use its functionality
  gSystem->AddIncludePath("-I/star/u/jmazer19/Y2017/STAR/FastJet/fastjet-install/include");

  // load the system libraries - these were defaults
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  // these are needed for new / additional classes
  gSystem->Load("libStPicoEvent");
  gSystem->Load("libStPicoDstMaker");

  // my libraries
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StMyAnalysisMaker");

  gSystem->ListLibraries();
} 

void LoadMacros()
{
}
