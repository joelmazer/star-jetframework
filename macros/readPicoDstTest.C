// ************************************** //
// Author: Joel Mazer
// Institution: Rutgers University
// STAR Collaboration
//
//
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
class StRho;
class StRhoBase;
class StRhoSparse;
class StMyAnalysisMaker;

// library and macro loading function
void LoadLibs();
void LoadMacros();

// constants
const double pi = 1.0*TMath::Pi();

// find kt jets and perform rho subtraction - not needed for constituents 2.0+ GeV or certain analyses
bool doBackgroundJets = kTRUE;

// run analysis for specific centrality bin
bool doCentSelection = kFALSE; // keep false to run over all centralities

// see enumerator below, 14 = Run14
Int_t RunYear = 14; 
// kFALSE when submitting jobs, kTRUE for tests
bool doTEST = kTRUE;     // FIXME double check before submitting!
bool dopp = kFALSE;      // kTRUE for pp data

// z-vertex cuts, go down to {-30, 30} for new centrality definitons - done below in setup section
Double_t ZVtxMin = -40.0;
Double_t ZVtxMax = 40.0;

void readPicoDstTest(const Char_t *inputFile="", const Char_t *outputFile="test.root", Int_t nEv = 10, const Char_t *fOutJobappend="")
{
        Int_t nEvents = 1000;
        if(nEv > 100) nEvents = 100000000;

        // run enumerators
        enum RunType_t {
          mJobSubmission = 0, // 0th spot
          mRun07 =  7, mRun08 = 8,  mRun09 = 9, mRun10 = 10,
          mRun11 = 11, mRun12 = 12, mRun13 = 13,
          mRun14 = 14,
          mRun15 = 15, mRun16 = 16, mRun17 = 17, mRun18 = 18
        };

        // set up Jet enumerators:
        // jet type
        enum EJetType_t {
          kFullJet,  // tracks + clusters
          kChargedJet,
          kNeutralJet
        };

        // jet algorithm
        enum EJetAlgo_t {
          kt_algorithm                    = 0, // background jets
          antikt_algorithm                = 1, // signal jets
          cambridge_algorithm             = 2,
          genkt_algorithm                 = 3,
          cambridge_for_passive_algorithm = 11,
          genkt_for_passive_algorithm     = 13,
          plugin_algorithm                = 99,
          undefined_jet_algorithm         = 999
        };

        // jet recombination scheme
        enum ERecoScheme_t {
          E_scheme        = 0,
          pt_scheme       = 1, pt2_scheme      = 2,
          Et_scheme       = 3, Et2_scheme      = 4,
          BIpt_scheme     = 5, BIpt2_scheme    = 6,
          WTA_pt_scheme   = 7, WTA_modp_scheme = 8,
          external_scheme = 99
        };

        // jet shape jet type
        enum EJetShapeJetType_t {
          kInclusiveJets,
          kLeadingJets,
          kSubLeadingJets
        };

        // Load necessary libraries and macros
        LoadLibs();
        LoadMacros();

        // =============================================================================== //
        // =============================================================================== //
        // over-ride functions
        if(dopp) {
          doCentSelection = kFALSE;  // can't ask for a particular centrality if requesting pp collisions
          doBackgroundJets = kFALSE; // don't do a rho background calculation when using pp collisions   
        }

        // input file for tests (based on Run) - updated for new Runs as needed - FIXME update for your own lists
        if((RunYear == mRun12) && doTEST && dopp) inputFile = "testLIST_Run12pp.list";
        if((RunYear == mRun14) && doTEST) inputFile = "Run14_P18ih_HPSS_15164046.list"; //"Run_15164046_files.list"; 
        if((RunYear == mRun16) && doTEST) inputFile = "test_run17124003_files.list";
        if((RunYear == mRun17) && doTEST && dopp) inputFile = "filelist_pp2017.list";
        cout<<"inputFileName = "<<inputFile<<endl;

        // centrality global flags: see StCentMaker
        // Centrality Selection can be set up in certain classes, to only run for a given centrality range
        Int_t CentralitySelection = StJetFrameworkPicoBase::kCent2050; // set as an example
        Int_t CentralityDefinition;
        if(RunYear == mRun12) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30;  // no centrality defintion for Run 12 pp
        if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30_AllLumi; // (NEW - from Nick: Aug16, 2019 set for all lumi)
        //if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30; // Run14 P18ih (NEW - from Nick June10, 2019)
        if(RunYear == mRun16) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P16id; // Run16
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
        Int_t fBadRunListVers = StJetFrameworkPicoBase::fBadRuns_w_missing_HT; // fBadRuns_w_missing_HT, fBadRuns_wo_missing_HT,

        // trigger flags - set as such
        // EmcTriggerEventType: selects which HT trigger to use
        Int_t EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2; // kIsHT1 or kIsHT2 or kIsHT3 (Set for Run14)
        if(RunYear == mRun12) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2; // StJetFrameworkPicoBase::kIsHT2;
        if(RunYear == mRun14) EmcTriggerEventTYpe = StJetFrameworkPicoBase::kIsHT2; // 
        if(RunYear == mRun16) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT1; // kIsHT1 Run16
        if(RunYear == mRun17) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT3; 
        Int_t MBEventType = StJetFrameworkPicoBase::kVPDMB5;        // this is default, want kVPDMB30 for new centrality definitions
        if(RunYear == mRun12) MBEventType = StJetFrameworkPicoBase::kRun12main; // default for Run12 pp
        if(RunYear == mRun17) MBEventType = StJetFrameworkPicoBase::kVPDMB; // default for Run17 pp
        Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerANY;   // kTriggerANY, kTriggerMB, kTriggerHT - only used by JetMaker and EPMaker (set to HT when doing EP corrections)

        // track flags
        bool usePrimaryTracks;
        if(RunYear == mRun12) usePrimaryTracks = kTRUE;
        if(RunYear == mRun14) usePrimaryTracks = kTRUE;
        if(RunYear == mRun16) usePrimaryTracks = kFALSE; // don't have primary tracks (at least mid 2018 with that current production)
        if(RunYear == mRun17) usePrimaryTracks = kTRUE;  

        // jet type flags
        Int_t fJetType;
        if(RunYear == mRun12) fJetType = kFullJet;
        if(RunYear == mRun14) fJetType = kFullJet; // kChargedJet
        if(RunYear == mRun16) fJetType = kChargedJet;
        if(RunYear == mRun17) fJetType = kFullJet;
        double fJetRadius = 0.3;  // 0.4, 0.3, 0.2
        double fJetConstituentCut = 0.2; // correlation analysis: 2.0, jet shape analysis: 2.0 (been using 2.0 for corrections)
        Int_t fJetAnalysisJetType = kLeadingJets;  // Jet analysis jet types - options: kInclusiveJets, kLeadingJets, kSubLeadingJets (need to set up in analysis when you want to use)
        cout<<"fJetType: "<<fJetType<<endl;

        // FIXME - be aware of which list is used! 
        // tower flags - lists to load for bad towers, see StJetFrameworkPicoBase and below
        Int_t TowerListToUse = 136; // doesn't matter for charged jets - Run14 136-122: jet-hadron, early jet shape - Pt dep lists set below
        if(dopp) TowerListToUse = 169;
        // see StJetFrameworkPicoBase:   9992000 - 2 GeV, 9991000 - 1 GeV, 9990200 - 0.2 GeV  (applicable currently for Run12 pp and Run14 AuAu)
        if(fJetConstituentCut == 2.0) TowerListToUse = 9992000;
        if(fJetConstituentCut == 1.0) TowerListToUse = 9991000;
        if(fJetConstituentCut == 0.2) TowerListToUse = 9990200;
        // Run12: 1 - Raghav's list, 102 - my initial list, 169 - new list
        // Run14: 136 - main list (updated version for AuAu 200 GeV Run14), 122 - past used list
        // Run14 P18ih: 999 (initial) 
        cout<<"TowerListUsed: "<<TowerListToUse<<endl;

        // event plane type and selection
        Int_t TPCEPSelectionType = StJetFrameworkPicoBase::kRemoveEtaStrip;
        Int_t EventPlaneTrackWeightMethod = StJetFrameworkPicoBase::kPtLinear2Const5Weight;

        // update settings for new centrality definitions
        if(CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30 ||
           CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30 ||
           CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30_AllLumi
        ) { ZVtxMin = -30.0; ZVtxMax = 30.0; }

        // =============================================================================== //
        // =============================================================================== //
        // =============================================================================== //

        // open and close output .root file (so it exist and can be updated by Analysis Tasks)
        TFile *fout = new TFile(outputFile, "RECREATE");
        fout->Close();

        // create chain
        StChain* chain = new StChain();

	// create the picoMaker maker
	//StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst"); // writing Pico's
        StPicoDstMaker *picoMaker = new StPicoDstMaker(2,inputFile,"picoDst");   // reading Pico's
        picoMaker->setVtxMode((int)(StPicoDstMaker::PicoVtxMode::Default));

        // create base class maker pointer
        StJetFrameworkPicoBase *baseMaker = new StJetFrameworkPicoBase("baseClassMaker");
        baseMaker->SetRunFlag(RunFlag);                  // run flag (year)
        baseMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
        baseMaker->SetBadRunListVers(fBadRunListVers);          // switch to select specific bad run version file
        baseMaker->SetBadTowerListVers(TowerListToUse);
        cout<<baseMaker->GetName()<<endl;  // print name of class instance

        // create centrality class maker pointer
        StCentMaker *CentMaker = new StCentMaker("CentMaker", picoMaker, outputFile, kFALSE);
        CentMaker->SetUsePrimaryTracks(usePrimaryTracks);       // use primary tracks
        CentMaker->SetEventZVtxRange(ZVtxMin, ZVtxMax);         // can be tighter for Run16 (-20,20)
        CentMaker->SetRunFlag(RunFlag);                         // Run Flag
        CentMaker->SetdoppAnalysis(dopp);                       // pp-analysis switch
        CentMaker->SetCentralityDef(CentralityDefinition);      // centrality definition
        CentMaker->SetUseBBCCoincidenceRate(kFALSE);            // BBC or ZDC (default) rate used?
        CentMaker->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        CentMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
        cout<<CentMaker->GetName()<<endl;  // print name of class instance
 
        // create JetFinder first (JetMaker)
        // task Name, track constituent cut, doHistos, output file name
        StJetMakerTask *jetTask = new StJetMakerTask("JetMaker", fJetConstituentCut, kTRUE, outputFile);
        jetTask->SetJetType(fJetType);          // jetType
        jetTask->SetJetAlgo(antikt_algorithm);  // jetAlgo
        jetTask->SetRecombScheme(BIpt2_scheme); // recombination scheme
        jetTask->SetRadius(fJetRadius);         // jet radius
        jetTask->SetJetsName("Jets");
        jetTask->SetMinJetPt(10.0);             // 15.0 signal jets
        jetTask->SetMaxJetTrackPt(30.0);        // max track constituent
        jetTask->SetMinJetTowerE(fJetConstituentCut);  // 2.0 correlations
        jetTask->SetHadronicCorrFrac(1.0);      // fractional hadronic correction
        jetTask->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // set for default options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
        jetTask->SetGhostArea(0.005);           // ghost area
        jetTask->SetMinJetArea(0.0);            // minimum jet area
        jetTask->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // fiducial eta acceptance
        jetTask->SetJetPhiRange(0,2.0*pi);      // phi acceptance
        jetTask->SetUsePrimaryTracks(usePrimaryTracks);
        jetTask->SetRunFlag(RunFlag);           // run flag      
        jetTask->SetdoppAnalysis(dopp);         // pp switch
        jetTask->SetEventZVtxRange(ZVtxMin, ZVtxMax);     // can be tighter for Run16 (-20,20)
        jetTask->SetTurnOnCentSelection(doCentSelection);
        jetTask->SetCentralityBinCut(CentralitySelection);
        jetTask->SetEmcTriggerEventType(EmcTriggerEventType);  // kIsHT1 or kIsHT2 or kIsHT3
        jetTask->SetTriggerToUse(TriggerToUse);
        jetTask->SetdoConstituentSubtr(kFALSE);      // implement constituent subtractor, if TRUE, don't want to subtract underlying event Rho
        jetTask->SetRejectBadRuns(RejectBadRuns);    // switch to load and than omit bad runs

        // create JetFinder for background now (JetMakerBG) - background jets to be used in Rho Maker
        StJetMakerTask *jetTaskBG;
        if(doBackgroundJets) {
          //StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", 0.2, kTRUE, outputFile); // all inclusive
          // task Name, track constituent cut, doHistos, output file name
          jetTaskBG = new StJetMakerTask("JetMakerBG", fJetConstituentCut, kTRUE, outputFile); 
          jetTaskBG->SetJetType(fJetType);          // jetType
          jetTaskBG->SetJetAlgo(kt_algorithm);      // jetAlgo
          jetTaskBG->SetRecombScheme(BIpt2_scheme); // recombination scheme
          jetTaskBG->SetRadius(fJetRadius);         // jet radius
          jetTaskBG->SetJetsName("JetsBG");
          jetTaskBG->SetMinJetPt(0.0);
          jetTaskBG->SetMaxJetTrackPt(30.0);
          jetTaskBG->SetMinJetTowerE(fJetConstituentCut);     // inclusive: 0.2
          jetTaskBG->SetHadronicCorrFrac(1.0);      // hadronic correlation fraction 0-1
          jetTaskBG->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // set for default options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
          jetTaskBG->SetGhostArea(0.005);           // ghost area
          jetTaskBG->SetMinJetArea(0.0);            // minimum jet area
          jetTaskBG->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // -0.5,0.5
          jetTaskBG->SetJetPhiRange(0, 2.0*pi);   // 0,pi
          jetTaskBG->SetUsePrimaryTracks(usePrimaryTracks);
          jetTaskBG->SetRunFlag(RunFlag);
          jetTaskBG->SetdoppAnalysis(dopp);
          jetTaskBG->SetEventZVtxRange(ZVtxMin, ZVtxMax);        // can be tighter for Run16 (-20,20)
          jetTaskBG->SetTurnOnCentSelection(doCentSelection);
          jetTaskBG->SetCentralityBinCut(CentralitySelection);
          jetTaskBG->SetEmcTriggerEventType(EmcTriggerEventType);// kIsHT1 or kIsHT2 or kIsHT3
          jetTaskBG->SetTriggerToUse(TriggerToUse);
          jetTaskBG->SetdoConstituentSubtr(kFALSE);
          jetTaskBG->SetRejectBadRuns(RejectBadRuns);            // switch to load and than omit bad runs
        }

        bool dohisto = kFALSE;  // histogram switch for Rho Maker

        // Rho task, and scale it up to include neutral constituents
        StRho *rhoTask;
        // RhoMaker name, doHisto switch, output filename, Background JetMaker name
        if(doBackgroundJets) { rhoTask = new StRho("StRho_JetsBG", dohisto, outputFile, "JetMakerBG");
        } else { rhoTask = new StRho("StRho_JetsBG", dohisto, outputFile, "JetMaker"); }
        rhoTask->SetExcludeLeadJets(2);
        rhoTask->SetOutRhoName("OutRho");
        rhoTask->SetRunFlag(RunFlag);
        rhoTask->SetdoppAnalysis(dopp);
        rhoTask->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
        rhoTask->SetTurnOnCentSelection(doCentSelection);
        rhoTask->SetCentralityBinCut(CentralitySelection);
        rhoTask->SetRejectBadRuns(RejectBadRuns);     // switch to load and than omit bad runs

        // create the analysis maker!
        bool doCorrJetPt = kFALSE;
        if(doBackgroundJets) doCorrJetPt = kTRUE;

        // MakerName, picoMakerPtr, outputFileName, JetMakerName, RhoMakerName
        StAnMaker *anaMaker = new StAnMaker("AnMaker", picoMaker, outputFile, "JetMaker", "StRho_JetsBG");
        anaMaker->SetUsePrimaryTracks(usePrimaryTracks); // use primary tracks
        anaMaker->SetCorrectJetPt(doCorrJetPt);          // subtract Rho BG
        anaMaker->SetJetRad(fJetRadius);                 // jet radius
        anaMaker->SetMinJetPt(10.0);                     // minimum jet pt cut
        anaMaker->SetJetMaxTrackPt(0.0);                 // jet track bias
        anaMaker->SetJetMaxTowerE(0.0);                  // jet tower bias
        anaMaker->SetMinTrackPt(0.2);                    // track quality cut (not related to constituents!)
        anaMaker->SetEventZVtxRange(ZVtxMin, ZVtxMax);   // can be tighter for Run16 (-20,20)
        anaMaker->SetTrackPhiRange(0.0, 2.0*TMath::Pi());// track phi acceptance
        anaMaker->SetTrackEtaRange(-1.0, 1.0);           // track eta acceptance
        anaMaker->SetJetConstituentCut(fJetConstituentCut);   // 0.2 for inclusive 
        anaMaker->SetRunFlag(RunFlag);                        // run flag: i.e. - Run14, Run16...
        anaMaker->SetdoppAnalysis(dopp);                      // running analysis over pp?
        anaMaker->SetTurnOnCentSelection(doCentSelection);    // run analysis for specific centrality: BOOLEAN
        anaMaker->SetCentralityBinCut(CentralitySelection);   // specific centrality range to run: if above is FALSE, this doesn't matter
        anaMaker->SetEmcTriggerEventType(EmcTriggerEventType);// kIsHT1 or kIsHT2 or kIsHT3
        anaMaker->SetRejectBadRuns(RejectBadRuns);            // switch to load and than omit bad runs
        cout<<anaMaker->GetName()<<endl;  // print name of class instance

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
        if(fout->IsOpen())   fout->Close();
}

void LoadLibs()
{
  // load fastjet libraries 3.x
  //gSystem->Load("libCGAL"); - not installed 
  // FIXME FIXME! make sure if you get a crash from running this macro the first time, that the $FASTJET variable points to the right place
  // else you will not find the libraries
  gSystem->Load("$FASTJET/lib/libfastjet");
  gSystem->Load("$FASTJET/lib/libsiscone");
  gSystem->Load("$FASTJET/lib/libsiscone_spherical");
  gSystem->Load("$FASTJET/lib/libfastjetplugins");
  gSystem->Load("$FASTJET/lib/libfastjettools");
  gSystem->Load("$FASTJET/lib/libfastjetcontribfragile");

  // add include path to use its functionality - this needs to be updated for different users FIXME
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
