// ************************************** //
// Sept28, 2017
// some current notes:
// Run16 AuAu 200 GeV:  P16ij on RCAS has no EmcalTriggers and no clusters
//     ::P16ij has triggers that need to be accessed a different way
//
// April, 2019
// Run14 AuAu 200 GeV: P18ih is new production
//

#include <TSystem>

// basic STAR classes
class StMemStat;
class StMaker;
class StChain;
class StPicoDstMaker;
class StRefMultCorr;

// jet-framework classes
class StJetMakerTask;
class StRho;
class StRhoBase;
class StRhoSparse;
class StMyAnalysisMaker;

// library and macro loading function
void LoadLibs();
void LoadMacros();
void PrintSetup();

// constants
const double pi = 1.0*TMath::Pi();

// find kt jets and perform rho subtraction: not used for main analyses, as pt constit >= 2.0 GeV
bool doBackgroundJets = kFALSE;

// run analysis for specific centrality bin
bool doCentSelection = kFALSE; //kTRUE; // FIXME

// select run year
Int_t RunYear = 14;      // (17) 14 for 2014 AuAu //mJobSubmission; (0)
// kFALSE when submitting jobs, kTRUE for tests
bool doTEST = kTRUE;     // FIXME double check before submitting!
bool dopp = kFALSE;       // FIXME kTRUE for pp data
bool doSetupQA = kFALSE; //kTRUE; // FIXME
bool doCentralityQA = kFALSE; 
bool doAnalysisQAoutputFile = kFALSE;
// bool dohisto = kFALSE;  // histogram switch for Rho Maker (SOME PROBLEM WITH THIS - FIXME)

// additional switches - Event plane calculations and corrections
bool doEventPlaneCorrections = kTRUE;   // needs to be set to kTRUE to run corrections
bool phi_shift_switch = kFALSE;         // not used! (TPC)
bool tpc_recenter_read_switch = kFALSE; // STEP1: kTRUE = write output for recentering (shift)
bool tpc_shift_read_switch = kFALSE;    // STEP2: kTRUE = write output for shift (flattening)
bool tpc_apply_corr_switch = kFALSE;    // STEP3: kTRUE = read shift info for final TPC calculation
bool zdc_recenter_read_switch = kFALSE; // STEP1: kTRUE = write output for recentering of zdc (shift)
bool zdc_shift_read_switch = kFALSE;    // STEP2: kTRUE = write output for zdc shift (flattening)
bool zdc_apply_corr_switch = kFALSE;    // STEP3: kTRUE = read shift info for final ZDC calculation
bool bbc_recenter_read_switch = kFALSE; // STEP1: kTRUE = write output for recentering of bbc (shift)
bool bbc_shift_read_switch = kFALSE;    // STEP2: kTRUE = write output for bbc shift (flattening)
bool bbc_apply_corr_switch = kFALSE;    // STEP3: kTRUE = read shift info for final BBC calculation
bool doSTEP1 = kFALSE;
bool doSTEP2 = kFALSE;
bool doSTEP3 = kTRUE;
bool doEPresolutions = kTRUE;

// z-vertex cuts (tighter cuts below, based on centrality definitions and cuts used to create them)
Double_t ZVtxMin = -40.0;
Double_t ZVtxMax = 40.0;

// testLIST_Run14.list
//void readPicoDst(const Char_t *inputFile="test2.list", const Char_t *outputFile="test2.root", Int_t nEv = 10)
//void readPicoDst(const Char_t *inputFile="newPicoDsts.list", const Char_t *outputFile="newtest.root", Int_t nEv = 10)
void readPicoDstMultPtBins(const Char_t *inputFile="Run_15164046_files.list", const Char_t *outputFile="test.root", Int_t nEv = 10, const Char_t *fEPoutJobappend="")
{
//  Int_t nEvents = 100000; // pp
//  Int_t nEvents = 50000;
//  Int_t nEvents = 10000;
//  Int_t nEvents = 4000;
  Int_t nEvents = 1000;
//  Int_t nEvents = 10;
  if(nEv > 100) nEvents = 100000000;

  // event plane correction switches: BBC / ZDC / TPC
  if(doSTEP1) { bbc_recenter_read_switch = kTRUE; }
  if(doSTEP2) { bbc_shift_read_switch = kTRUE; }
  if(doSTEP3) { bbc_apply_corr_switch = kTRUE; bbc_shift_read_switch = kTRUE; bbc_recenter_read_switch = kTRUE; } 

  if(doSTEP1) { zdc_recenter_read_switch = kTRUE; }
  if(doSTEP2) { zdc_shift_read_switch = kTRUE; }
  if(doSTEP3) { zdc_apply_corr_switch = kTRUE; zdc_shift_read_switch = kTRUE; zdc_recenter_read_switch = kTRUE; }

  if(doSTEP1) { tpc_recenter_read_switch = kTRUE; }
  if(doSTEP2) { tpc_shift_read_switch = kTRUE; }
  if(doSTEP3) { tpc_apply_corr_switch = kTRUE; tpc_shift_read_switch = kTRUE; tpc_recenter_read_switch = kTRUE; }

  // run enumerators
  enum RunType_t {
    mJobSubmission = 0, // 0th spot - default unless set
    mRun07 =  7, mRun08 = 8,  mRun09 = 9, mRun10 = 10,
    mRun11 = 11, mRun12 = 12, mRun13 = 13,
    mRun14 = 14,
    mRun15 = 15, mRun16 = 16, mRun17 = 17, mRun18 = 18
  };

  // set up Jet enumerators:
  // jet type
  enum EJetType_t {
    kFullJet,    // tracks + clusters
    kChargedJet, // tracks only
    kNeutralJet  // towers only
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
  // =============================================================================== //
  // over-ride functions
  if(dopp) {
    doCentSelection = kFALSE;  // can't ask for a particular centrality if requesting pp collisions
    doBackgroundJets = kFALSE; // don't do a rho background calculation when using pp collisions   
  }

  // input file for tests (based on Run) - updated for new Runs as needed
  if((RunYear == mRun12) && doTEST && dopp) inputFile = "testLIST_Run12pp.list"; 
  if((RunYear == mRun14) && doTEST) inputFile = "Run_15164046_files.list"; //Run_15151042_files.list"; //"testLIST_Run14.list";
  if((RunYear == mRun16) && doTEST) inputFile = "test_run17124003_files.list";
  if((RunYear == mRun17) && doTEST && dopp) inputFile = "filelist_pp2017.list";
  cout<<"inputFileName = "<<inputFile<<endl;

  // centrality global flags - no centrality for pp collisions
  Int_t CentralitySelection = StJetFrameworkPicoBase::kCent2050;
  Int_t CentralityDefinition;
  if(RunYear == mRun12) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30;  // no centrality defintion for Run 12 pp
  //if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult;                 // Run14 P16id ??FIXME???
  //if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30; // Run14 P17id (NEW - from Nick Oct 23)
  if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30; // Run14 P18ih (NEW - from Nick June 10, 2019)
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
  Int_t EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2; // kIsHT1 or kIsHT2 or kIsHT3 (set for Run14)
  if(RunYear == mRun12) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2; // StJetFrameworkPicoBase::kIsHT2;
  if(RunYear == mRun16) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT1; // kIsHT1 Run16
  if(RunYear == mRun17) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT3; 
  Int_t MBEventType = StJetFrameworkPicoBase::kVPDMB5;        // this is default (Run14)
  if(RunYear == mRun12) MBEventType = StJetFrameworkPicoBase::kRun12main; // default for Run12 pp
  if(RunYear == mRun17) MBEventType = StJetFrameworkPicoBase::kVPDMB; // default for Run17 pp
  Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerHT;    // kTriggerANY, kTriggerMB, kTriggerHT
//  Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerANY;    // kTriggerANY, kTriggerMB, kTriggerHT     FIXME
  Int_t TowerListToUse = 136; // Run14 136-122: jet-hadron, early jet shape // been using 51,  3, 79   - doesn't matter for charged jets
  TowerListToUse = 1; //10000;
  //Int_t TowerListToUse = 999; // LIST to use for initial Run14 config 
  if(dopp) TowerListToUse = 169;
  // Run12: 1 - Raghav's list, 102 - my initial list, 169 - new list
  // Run14: 136 - main list, 122 - past used list
  // Run14 P18ih: 999 (initial) 

  // track flags
  bool usePrimaryTracks;
  if(RunYear == mRun12) usePrimaryTracks = kTRUE;
  if(RunYear == mRun14) usePrimaryTracks = kTRUE;
  if(RunYear == mRun16) usePrimaryTracks = kFALSE; // doesn't have primary tracks 
  if(RunYear == mRun17) usePrimaryTracks = kTRUE;  

  // jet type flags
  Int_t fJetType;
  if(RunYear == mRun12) fJetType = kFullJet;
  if(RunYear == mRun14) fJetType = kFullJet;
  if(RunYear == mRun16) fJetType = kChargedJet; // not towers, only has kChargedJet
  if(RunYear == mRun17) fJetType = kFullJet;
  double fJetRadius = 0.3;  // 0.4, 0.3, 0.2
  double fJetConstituentCut = 1.0; // FIXME  correlation analysis: 2.0, jet shape analysis: 1.0 (been using 2.0 for corrections)
  Int_t fJetAnalysisJetType = kLeadingJets; // Jet analysis jet types - options: kInclusiveJets, kLeadingJets, kSubLeadingJets
  bool fDoJetShapeAnalysis = kTRUE;              // switch for doing jet shape analysis
  bool fDoJetHadronCorrelationAnalysis = kFALSE; // switch for doing jet-hadron correlation analysis
  cout<<"JetType: "<<fJetType<<"     JetRad: "<<fJetRadius<<"     JetConstit Cut: "<<fJetConstituentCut<<endl;

  // event plane type / configuration flags - make sure this is set properly!!!
  // StMyAnalysisMaker::kRemoveEtaPhiCone, StMyAnalysisMaker::kRemoveEtaPhiConeLeadSub
  Int_t TPCEPSelectionType = StJetFrameworkPicoBase::kRemoveEtaStrip;
  Int_t EventPlaneTrackWeightMethod = StJetFrameworkPicoBase::kPtLinear2Const5Weight;
  bool doUseMainEPAngle = kFALSE;  // FIXME     kFALSE: pt-bin approach, kTRUE: use 0.2-2.0 GeV charged tracks used for event plane reconstruction

  // update settings for new centrality definitions
  if(CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30 ||
     CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30
  ) { ZVtxMin = -28.0; ZVtxMax = 28.0; }

  // =============================================================================== //
  // =============================================================================== //
  // =============================================================================== //

  // open and close output .root file (so it exist and can be updated by Analysis Tasks)
  TFile *fout = new TFile(outputFile, "RECREATE");
  fout->Close();

  // don't need event plane files for pp collisions
  if(!dopp && !doSetupQA) {
    TFile *fEPout;
    if(doSTEP1) fEPout = new TFile(Form("recenter_calib_file%s.root", fEPoutJobappend), "RECREATE");
    if(doSTEP2) fEPout = new TFile(Form("shift_calib_file%s.root", fEPoutJobappend), "RECREATE"); 
    if(doSTEP3) fEPout = new TFile(Form("final_eventplane%s.root", fEPoutJobappend), "RECREATE");
    fEPout->Close();        
  }

  // don't create empty files if your not running QA below for analysis task
  if(doAnalysisQAoutputFile){
    TFile *fQAout = new TFile(Form("QAtracksANDjets%s.root", fEPoutJobappend), "RECREATE");
    fQAout->Close();
  }

  // create chain
  StChain* chain = new StChain();

  // create the picoMaker maker
  //StPicoDstMaker *picoMaker = new StPicoDstMaker(0,inputFile,"picoDst");
  StPicoDstMaker *picoMaker = new StPicoDstMaker(2,inputFile,"picoDst"); // updated Aug6th
  picoMaker->setVtxMode((int)(StPicoDstMaker::PicoVtxMode::Default));

  // create base class maker pointer
  StJetFrameworkPicoBase *baseMaker = new StJetFrameworkPicoBase("baseClassMaker");
  baseMaker->SetRunFlag(RunFlag);           // run flag (year)
  baseMaker->SetBadTowerListVers(TowerListToUse);

  // this doesn't work
  //StJetPicoDefinitions::SetDebugLevel(9);

  // create JetFinder first (JetMaker)
  StJetMakerTask *jetTask = new StJetMakerTask("JetMaker", fJetConstituentCut, kTRUE, outputFile);
  jetTask->SetJetType(fJetType);          // jetType
  jetTask->SetJetAlgo(antikt_algorithm);  // jetAlgo
  jetTask->SetRecombScheme(BIpt2_scheme); // recomb
  jetTask->SetRadius(fJetRadius);         // jet radius
  jetTask->SetJetsName("Jets");           // jets name
  jetTask->SetMinJetPt(10.0);             // 15.0(?) signal jets
  jetTask->SetMaxJetTrackPt(30.0);        // max track constituent
  jetTask->SetMinJetTowerE(fJetConstituentCut);  // 2.0 correlations
  jetTask->SetHadronicCorrFrac(1.0);      // fractional hadronic correction
  jetTask->SetGhostArea(0.005);
  jetTask->SetMinJetArea(0.0);
  jetTask->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // fiducial eta acceptance
  jetTask->SetJetPhiRange(0,2.0*pi);        // phi acceptance
  jetTask->SetUsePrimaryTracks(usePrimaryTracks);
  jetTask->SetRunFlag(RunFlag);           // run flag (year)
  jetTask->SetdoppAnalysis(dopp);
  jetTask->SetCentralityDef(CentralityDefinition); // run based centrality definition
  jetTask->SetEventZVtxRange(ZVtxMin, ZVtxMax);    // can be tighter for Run16 (-20,20)
  jetTask->SetTurnOnCentSelection(doCentSelection);
  jetTask->SetCentralityBinCut(CentralitySelection);
  jetTask->SetUseBBCCoincidenceRate(kFALSE);
  jetTask->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
  jetTask->SetTriggerToUse(TriggerToUse);
  jetTask->SetBadTowerListVers(TowerListToUse);        
  jetTask->SetdoConstituentSubtr(kFALSE); // FIXME
  jetTask->SetRejectBadRuns(RejectBadRuns);        // switch to load and than omit bad runs
  //jetTask->SetDebugLevel(2); // 8 spits out cluster/tower stuff

/*
  if (bFillGhosts) jetTask->SetFillGhost();
  if (lockTask) jetTask->SetLocked();
*/

  // background jets to be used in Rho Maker
  StJetMakerTask *jetTaskBG;
  if(doBackgroundJets) {
    // create JetFinder for background now (JetMakerBG)
    //StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", 0.2, kTRUE, outputFile); // all inclusive
    ////StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", fJetConstituentCut, kTRUE, outputFile);
    jetTaskBG = new StJetMakerTask("JetMakerBG", fJetConstituentCut, kTRUE, outputFile);
    jetTaskBG->SetJetType(fJetType);          // jetType
    jetTaskBG->SetJetAlgo(kt_algorithm);      // jetAlgo
    jetTaskBG->SetRecombScheme(BIpt2_scheme); // recombination scheme
    jetTaskBG->SetRadius(fJetRadius);         // 0.4
    jetTaskBG->SetJetsName("JetsBG");
    jetTaskBG->SetMinJetPt(0.0);
    jetTaskBG->SetMaxJetTrackPt(30.0);
    jetTaskBG->SetMinJetTowerE(fJetConstituentCut); // inclusive: 0.2
    jetTaskBG->SetHadronicCorrFrac(1.0);
    jetTaskBG->SetGhostArea(0.005);
    jetTaskBG->SetMinJetArea(0.0);
    jetTaskBG->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // -0.5,0.5
    jetTaskBG->SetJetPhiRange(0, 2.0*pi);        // 0,pi
    jetTaskBG->SetUsePrimaryTracks(usePrimaryTracks);
    jetTaskBG->SetRunFlag(RunFlag);
    jetTaskBG->SetdoppAnalysis(dopp);
    jetTaskBG->SetCentralityDef(CentralityDefinition);
    jetTaskBG->SetEventZVtxRange(ZVtxMin, ZVtxMax);         // can be tighter for Run16 (-20,20)
    jetTaskBG->SetTurnOnCentSelection(doCentSelection);
    jetTaskBG->SetCentralityBinCut(CentralitySelection);
    jetTaskBG->SetUseBBCCoincidenceRate(kFALSE);
    jetTaskBG->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
    jetTaskBG->SetTriggerToUse(TriggerToUse);
    jetTaskBG->SetBadTowerListVers(TowerListToUse);
    jetTaskBG->SetdoConstituentSubtr(kFALSE);
    jetTaskBG->SetRejectBadRuns(RejectBadRuns);        // switch to load and than omit bad runs
  }

  // this is the centrality dependent scaling for RHO as used in ALICE - don't use right now
  // DON'T need this for STAR due to the same acceptance range of towers and tracks
  // s(Centrality) = 0.00015 ×Centrality^2 ? 0.016 ×Centrality + 1.91
  // scaled function for rho and setting parameters
  TF1* sfunc = new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
  sfunc->SetParameter(2,1.80642);
  sfunc->SetParameter(1,-0.0112331);
  sfunc->SetParameter(0,0.0001139561);
  
  // histogram switch
  bool dohisto = kFALSE;  // histogram switch for Rho Maker

  // run the below tasks when running an Analysis and not strictly QA
  if(!doSetupQA) {
    // Rho task, and scale it up to include neutral constituents
    StRho *rhoTask;
    if(doBackgroundJets) { rhoTask = new StRho("StRho_JetsBG", dohisto, outputFile, "JetMakerBG"); // kt jets, for background
    } else { rhoTask = new StRho("StRho_JetsBG", dohisto, outputFile, "JetMaker"); } // signal jets, bc not doing a rho subtraction
    rhoTask->SetExcludeLeadJets(2);
    rhoTask->SetOutRhoName("OutRho");
    rhoTask->SetRunFlag(RunFlag);
    rhoTask->SetdoppAnalysis(dopp);
    rhoTask->SetCentralityDef(CentralityDefinition);
    rhoTask->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
    rhoTask->SetTurnOnCentSelection(doCentSelection);
    rhoTask->SetCentralityBinCut(CentralitySelection);
    rhoTask->SetUseBBCCoincidenceRate(kFALSE);
    //rhoTask->SetRejectBadRuns(RejectBadRuns);        // switch to load and than omit bad runs
    //rhoTask->SetScaleFunction(sfunc); // don't NEED

    //if(doTPCptassocBin) {
    //  if(ptbin == 0) { if((pt > 0.20) && (pt <= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
    //  if(ptbin == 1) { if((pt > 0.50) && (pt <= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
    //  if(ptbin == 2) { if((pt > 1.00) && (pt <= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
    //  if(ptbin == 3) { if((pt > 1.50) && (pt <= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
    //  if(ptbin == 4) { if((pt > 2.00) && (pt <= 20.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations
    //}

    if(!dopp) {
      //int ptAssocBins[8] = {0,1,2,3,4,5,6,7};
      //int ptAssocBins[5] = {0, 1, 2, 3, 4};         
      StEventPlaneMaker *EPMaker[5];

      // May 6th, changed named from EventPlaneMaker%i to EventPlaneMaker_bin%i
      for(int i = 0; i < 5; i++) {
        //if(i < 4) continue; // test
        ///EPMaker[i] = new StEventPlaneMaker(Form("EventPlaneMaker%i", i), picoMaker, "JetMaker", "StRho_JetsBG");
        EPMaker[i] = new StEventPlaneMaker(Form("EventPlaneMaker_bin%i", i), picoMaker, "JetMaker", "StRho_JetsBG"); //EPMaker[i]->SetJetMaxTrackPt(4.0);
        EPMaker[i]->SetMinJetPt(10.0);                     // perhaps lower this TODO
        EPMaker[i]->SetUsePrimaryTracks(usePrimaryTracks); // use primary tracks
        EPMaker[i]->SetCorrectJetPt(kFALSE);               // subtract Rho BG
        EPMaker[i]->SetMinTrackPt(0.2);                    // may be able to just set this to 0? - TODO
        EPMaker[i]->SetEventZVtxRange(ZVtxMin, ZVtxMax);   // can be tighter for Run16 (-20,20)
        EPMaker[i]->SetTrackPhiRange(0.0, 2.0*TMath::Pi());
        EPMaker[i]->SetTrackEtaRange(-1.0, 1.0);
        EPMaker[i]->SetEventPlaneMaxTrackPtCut(2.0);   // 5.0 is default (can use < 2 to avoid autocorrelation issues)
        EPMaker[i]->SetExcludeLeadingJetsFromFit(1.0); // 1.0 is default
        EPMaker[i]->SetEventPlaneTrackWeight(EventPlaneTrackWeightMethod); // type of track weighting selection
        EPMaker[i]->SetTPCEventPlaneMethod(TPCEPSelectionType);            // TPC type method
        EPMaker[i]->SetdoEPTPCptAssocMethod(kTRUE);    // calculate TPC event plane / RES on pt assoc bin basis
        EPMaker[i]->SetEPTPCptAssocBin(i);             // pt assoc bin to use
        EPMaker[i]->SetRunFlag(RunFlag);               // run flag (year)
        EPMaker[i]->SetdoppAnalysis(dopp);
        EPMaker[i]->SetJetType(fJetType);              // jet type (full, charged, neutral)
        EPMaker[i]->SetJetRad(fJetRadius);             // jet radius
        EPMaker[i]->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        EPMaker[i]->SetTriggerToUse(TriggerToUse);               // switch to require MB or HT to run analysis
        EPMaker[i]->SetCentralityDef(CentralityDefinition);      // centrality definition
        EPMaker[i]->SetTurnOnCentSelection(doCentSelection);     // run analysis for specific centrality
        EPMaker[i]->SetCentralityBinCut(CentralitySelection);    // specific centrality range to run 
        EPMaker[i]->SetUseBBCCoincidenceRate(kFALSE);            // only use ZDC coincidence rate for refMult2
        EPMaker[i]->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
        if(doSTEP1) EPMaker[i]->SetOutFileNameEP(Form("recenter_calib_file%s.root", fEPoutJobappend)); // set output file for recentering calibration of event plane
        if(doSTEP2) EPMaker[i]->SetOutFileNameEP(Form("shift_calib_file%s.root", fEPoutJobappend));    // set output file for shift calibration of event plane
        if(doSTEP3) EPMaker[i]->SetOutFileNameEP(Form("final_eventplane%s.root", fEPoutJobappend));    // set output file for final event plane calculation
        // switches for event plane correction
        if(doEventPlaneCorrections) {
          EPMaker[i]->SetTPCRecenterRead(tpc_recenter_read_switch);   // TPC - STEP1
          EPMaker[i]->SetTPCShiftRead(tpc_shift_read_switch);         // TPC - STEP2
          EPMaker[i]->SetTPCApplyCorrections(tpc_apply_corr_switch);  // TPC - STEP3
          EPMaker[i]->SetZDCrecentering(zdc_recenter_read_switch);    // ZDC - STEP1
          EPMaker[i]->SetZDCShiftRead(zdc_shift_read_switch);         // ZDC - STEP2 
          EPMaker[i]->SetZDCApplyCorrections(zdc_apply_corr_switch);  // ZDC - STEP3
          EPMaker[i]->SetBBCrecentering(bbc_recenter_read_switch);    // BBC - STEP1
          EPMaker[i]->SetBBCShiftRead(bbc_shift_read_switch);         // BBC - STEP2
          EPMaker[i]->SetBBCApplyCorrections(bbc_apply_corr_switch);  // BBC - STEP3

          EPMaker[i]->SetdoEventPlaneRes(doEPresolutions);            // Event Plane resolutions
        }
      } // loop over pt bins
    }   // only run event plane maker for non-pp datasets

    // create the analysis maker!
    bool doComments = kFALSE;
    bool doCorrJetPt = kFALSE;
    if(doBackgroundJets) doCorrJetPt = kTRUE;
    //StMyAnalysisMaker *anaMaker = new StMyAnalysisMaker("AnalysisMaker", picoMaker, outputFile, doComments, 5.0, 0.15, "JetMaker", "StRho_JetsBG");
    StMyAnalysisMaker3 *anaMaker[9]; // was [5] for correlation analysis
    int ptAssocBins[9] = {0, 1, 2, 3, 4, 4,4,4,4};

    for(int i = 0; i < 9; i++) { // jet shape analysis
    //for(int i = 0; i < 5; i++) { // correlation analysis
      //if(i < 10) continue; // for jet shape analysis
      if(i<8) continue; // TEST for just inclusive bin
      //if(i < 6) continue;
      //if(i < 4) continue; // cut for tests
      // update the below when running analysis - minjet pt and bias requirement
      anaMaker[i] = new StMyAnalysisMaker3(Form("AnalysisMaker_bin%i", i), picoMaker, outputFile, doComments, 10.0, 4.0, "JetMaker", "StRho_JetsBG");
      //anaMaker[i] = new StMyAnalysisMaker3(Form("AnalysisMaker_bin%i", i), picoMaker, outputFile, doComments, 15.0, 4.0, "JetMaker", "StRho_JetsBG");
      anaMaker[i]->SetEventMixing(kTRUE);   // event mixing switch
      anaMaker[i]->SetCentBinSize(5);       // FIXME centrality bin size for mixed events - NOT currently used
      anaMaker[i]->SetMixingTracks(25000);  // 50,000 in ALICE (2.76 TeV)
      anaMaker[i]->SetNMixedTr(1500);       // was 2500
      anaMaker[i]->SetNMixedEvt(1);         // was 5, was 2 before June11, 2018
      anaMaker[i]->SetDoUseMultBins(kTRUE); // use multiplicity bins (hand defined) instead of cent bins - Jet Shape Analysis
      anaMaker[i]->SetdoUseEPBins(kFALSE);  // use event plane bins // TODO 
      anaMaker[i]->SetDoFilterPtMixEvents(kFALSE);            // DONT USE, filter mixed event pool by pt cut switch
      anaMaker[i]->SetUsePrimaryTracks(usePrimaryTracks);     // use primary tracks
      anaMaker[i]->SetMinTrackPt(0.2);      // track quality cut (not related to constituents!)
      anaMaker[i]->SetCorrectJetPt(doCorrJetPt); // subtract Rho BG
      //anaMaker[i]->SetJetMaxTrackPt(4.0); // jet track bias (set in constructor)
      anaMaker[i]->SetJetMaxTowerE(4.0);    // jet tower bias
      anaMaker[i]->SetJetRad(fJetRadius);   // jet radius
      anaMaker[i]->SetJetConstituentCut(fJetConstituentCut);  // 2.0 is default 
      anaMaker[i]->SetJetLJSubLJPtThresholds(20.0, 10.0);     // LJ pt > 20.0, SubLJ pt > 10.0 GeV
      anaMaker[i]->SetdoRequireAjSelection(kFALSE);           // requirement of Aj selection on jets for jet shape analysis
      anaMaker[i]->SetEventZVtxRange(ZVtxMin, ZVtxMax);       // can be tighter for Run16 (-20,20)
      anaMaker[i]->SetTrackPhiRange(0.0, 2.0*TMath::Pi());
      anaMaker[i]->SetTrackEtaRange(-1.0, 1.0);
      anaMaker[i]->SetEventPlaneMaxTrackPtCut(2.0);   // 5.0 is default (can use < 2 to avoid autocorrelation issues)
      anaMaker[i]->SetExcludeLeadingJetsFromFit(1.0); // 1.0 is default
      anaMaker[i]->SetEventPlaneTrackWeight(EventPlaneTrackWeightMethod);
      anaMaker[i]->SetTPCEventPlaneMethod(TPCEPSelectionType);
      anaMaker[i]->SetdoEPTPCptAssocMethod(kTRUE);    // calculate TPC event plane / RES on pt assoc bin basis - TODO keep updated
      anaMaker[i]->SetEPTPCptAssocBin(ptAssocBins[i]);//i);     // pt assoc bin to use ************
      anaMaker[i]->SetEventPlaneMakerName("EventPlaneMaker_bin");// event plane maker name
      anaMaker[i]->SetdoUseMainEPAngle(doUseMainEPAngle);       // sets EP angle to use 0.2-2.0 GeV tracks for reconstruction instead of pt-bin approach
      anaMaker[i]->SetRunFlag(RunFlag);                         // run flag: i.e. - Run14, Run16...
      anaMaker[i]->SetdoppAnalysis(dopp);                       // running analysis over pp?
      anaMaker[i]->SetCentralityDef(CentralityDefinition);      // centrality definition
      anaMaker[i]->SetTurnOnCentSelection(doCentSelection);     // run analysis for specific centrality
      anaMaker[i]->SetCentralityBinCut(CentralitySelection);    // specific centrality range to run 
      anaMaker[i]->SetUseBBCCoincidenceRate(kFALSE);            // only use ZDC coincidence rate for refMult2
      anaMaker[i]->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
      anaMaker[i]->SetPrintEventCounter(kFALSE); 
      anaMaker[i]->SetdoEventPlaneRes(kFALSE);                  // write EP resolution profiles
      anaMaker[i]->SetdoJetHadronCorrelationAnalysis(fDoJetHadronCorrelationAnalysis); // do jet-hadron correlation analysis
      anaMaker[i]->SetdoJetShapeAnalysis(fDoJetShapeAnalysis);  // do jet shape analysis 
      anaMaker[i]->SetJetAnalysisJetType(fJetAnalysisJetType);  // jet type for use with jet analysis
      anaMaker[i]->SetJetShapePtAssocBin(i);                    // pt associated bin for jet shape analysis
      anaMaker[i]->SetJetShapeTrackPtRange(0.2, 30.0);          // jet shape analysis - track pt range
      anaMaker[i]->SetCentBinSizeJS(5);                         // cent bin size for jet shape analysis mixed events, change to 5
      anaMaker[i]->SetBadTowerListVers(TowerListToUse);         // bad tower list version
      anaMaker[i]->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
      cout<<anaMaker[i]->GetName()<<endl;                       // print name of class instance

      anaMaker[i]->SetOutFileNameQA(Form("QAtracksANDjets%s.root", fEPoutJobappend));
      anaMaker[i]->SetWriteTrackQAHistograms(kFALSE);
      anaMaker[i]->SetWriteJetQAHistograms(kFALSE);

      ////anaMaker[i]->SetDebugLevel(StMyAnalysisMaker3::kDebugEventPlaneCalc);  // 0 turned off
      ////anaMaker[i]->SetDebugLevel(StMyAnalysisMaker3::kDebugGeneralEvt);
      ////anaMaker[i]->SetDebugLevel(StMyAnalysisMaker3::kDebugJetvsEPtype);
      //anaMaker[i]->SetDebugLevel(StMyAnalysisMaker3::kDebugEmcTrigger);
    }

    if(doCentralityQA) {
      // update the below when running analysis - minjet pt and bias requirement
      StCentralityQA *cenMaker = new StCentralityQA("CentralityQA", picoMaker, outputFile, kFALSE);
      cenMaker->SetUsePrimaryTracks(usePrimaryTracks); // use primary tracks
      cenMaker->SetMinTrackPt(0.2);      // track quality cut (not related to constituents!)
      cenMaker->SetEventZVtxRange(ZVtxMin, ZVtxMax);         // can be tighter for Run16 (-20,20)
      cenMaker->SetTrackPhiRange(0.0, 2.0*TMath::Pi());        // track phi acceptance
      cenMaker->SetTrackEtaRange(-1.0, 1.0);                 // track eta acceptance
      ////cenMaker->SetDebugLevel(StMyAnalysisMaker3::kDebugGeneralEvt);
      ////cenMaker->SetDebugLevel(StMyAnalysisMaker3::kDebugEmcTrigger);
      cenMaker->SetRunFlag(RunFlag);                         // run flag (year)
      cenMaker->SetdoppAnalysis(dopp);
      cenMaker->SetCentralityDef(CentralityDefinition);      // centrality definition
      cenMaker->SetTurnOnCentSelection(doCentSelection);     // run analysis for specific centrality
      //cenMaker->SetCentralityBinCut(CentralitySelection);    // specific centrality range to run 
      cenMaker->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
      cenMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
      cout<<cenMaker->GetName()<<endl;  // print name of class instance
    }

  } // if !doSetupQA

  // initialize chain
  chain->Init();
  cout<<"chain->Init();"<<endl;
  int total = picoMaker->chain()->GetEntries();
  cout << " Total entries = " << total << endl;
  if(nEvents > total) nEvents = total;
  
  for (Int_t i = 0; i < nEvents; i++){
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
  if(fout->IsOpen()) {  fout->Close();  }
  if(!dopp && !doSetupQA)    {  if(fEPout->IsOpen()) fEPout->Close();  }
  if(doAnalysisQAoutputFile) {  if(fQAout->IsOpen()) fQAout->Close();  }

  //StMemStat::PrintMem("load StChain");
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
  //gROOT->LoadMacro("StRoot/StMyAnalysisMaker/ConfigJetFinders.C");
}

void PrintSetup()
{

}
