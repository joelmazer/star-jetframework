// ************************************** //
// Sept28, 2017
// some current notes:
// Run16 AuAu 200 GeV:  P16ij on RCAS has no EmcalTriggers and no clusters
//     ::P16ij has triggers that need to be accessed a different way
//
// April, 2019
// Run14 AuAu 200 GeV: P18ih is new production
//
// March, 2020
// Run14 AuAu 200 GeV: P18ih, SL20a
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
class StMyAnalysisMaker;
class StJetFrameworkPicoBase;

// library and macro loading function
void LoadLibs();
void LoadMacros();
void PrintSetup();

// constants
const double pi = 1.0*TMath::Pi();

// find kt jets and perform rho subtraction: not used for main analyses, as pt constit >= 2.0 GeV
bool doBackgroundJets = kFALSE;

// run analysis for specific centrality bin - useful when analysis produces large output
bool doCentSelection = kTRUE;

// for running on the Run14 MB dataset for Neil
// 	- doMBset = kTRUE, and local: doTEST = kTRUE
// select run year
Int_t RunYear = 14;      // (17) 14 for 2014 AuAu
// kFALSE when submitting jobs, kTRUE for tests
bool doMBset = kFALSE;  // HFT dataset (MB triggers for Run14 AuAu)
bool doTEST = kFALSE;     // FIXME double check before submitting!
bool dopp = kFALSE;       // FIXME kTRUE for pp data
bool doSetupQA = kFALSE; //kFALSE; //kTRUE; // FIXME
bool doEPseparate = kFALSE; // added this Sept 14- due to new setup and thus different 'output files' resulting from a submission
bool doAnalysisQAoutputFile = kFALSE;
// bool dohisto = kFALSE;  // histogram switch for Rho Maker (SOME PROBLEM WITH THIS - FIXME)

double fJetRadius = 0.4;  // 0.4, 0.3, 0.2

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
// keep at 40 when generating event plane corrections - over-written below
Double_t ZVtxMin = -40.0;
Double_t ZVtxMax = 40.0;

// "Run14_P18ih_HPSS_15164046.list" - old file
//void readPicoDstMultPtBins(const Char_t *inputFile="Run_15164046_files.list", const Char_t *outputFile="test.root", Int_t nEv = 10, const Char_t *fEPoutJobappend="")
void readPicoDstMultPtBins(const Char_t *inputFile="Run14_P18ih_SL20a_15110029.list", const Char_t *outputFile="test.root", Int_t nEv = 10, const Char_t *fEPoutJobappend="")
{
//  Int_t nEvents = 100000; // pp
//  Int_t nEvents = 50000;
//  Int_t nEvents = 10000;
  Int_t nEvents = 2500;
//  Int_t nEvents = 1000;
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
    mRun07 =  7, mRun08 = 8,  mRun09 = 9,  mRun10 = 10, mRun11 = 11, mRun12 = 12, mRun13 = 13, 
    mRun14 = 14, mRun15 = 15, mRun16 = 16, mRun17 = 17, mRun18 = 18, mRun19 = 19, mRun20 = 20
  };

  // set up Jet enumerators:
  // jet type
  enum EJetType_t {
    kFullJet,    // tracks + towers
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
  // E_scheme - for jet mass measurements
  // BIpt2_scheme - standard scheme for jet reconstruction
  enum ERecoScheme_t {
    E_scheme        = 0, // used when calculating jet mass
    pt_scheme       = 1, pt2_scheme      = 2,
    Et_scheme       = 3, Et2_scheme      = 4,
    BIpt_scheme     = 5, BIpt2_scheme    = 6,
    WTA_pt_scheme   = 7, WTA_modp_scheme = 8,
    external_scheme = 99
  };

  // jet type - need to configure in your own analysis class
  enum EJetShapeJetType_t {
    kInclusiveJets, kLeadingJets, kSubLeadingJets
  };

  // enumerator for tracking efficiency systematic uncertainty type
  enum ESystematicUncType_t {
    kDoNothing, kTrkEffMin, kTrkEffMax
  };

  // enumerator for hadronic correction systematic uncertainty type
  enum ESystematicUncHadCorrType_t {
    kHadCorrOFF, kHadCorrON, kHadCorrVAR
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

  // input file for tests (based on Run) - update for new Runs as needed - FIXME
  if((RunYear == mRun12) && doTEST && dopp)  inputFile = "Run12pp_P12id_SL18f_Run13066104_2020-testFiles.list"; // Sept 2020 new re-production from tape
// inputFile = "Run12pp_P12id_SL18f_xrootd.list"; // "testLIST_Run12pp.list"; 
  if((RunYear == mRun14) && doTEST) inputFile = "Run14_P18ih_SL20a_15110029.list";
  if((RunYear == mRun14) && doMBset && doTEST) inputFile = "Run14_P16id_SL18f_MB_test.list"; // test files for HFT dataset - (MB triggers)
  cout<<"inputFileName = "<<inputFile<<endl;

  // centrality global flags - no centrality for pp collisions
  Int_t CentralitySelection = StJetFrameworkPicoBase::kCent2050;  // ::kCent010;
  Int_t CentralityDefinition;
  if(RunYear == mRun12) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30; // no centrality defintion for Run 12 pp
  if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30_AllLumi; // (NEW - from Nick: Aug16, 2019 set for all lumi)
  //if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30; // Run14 P18ih (Nick: June 10, 2019)
  cout<<"Centrality definition: "<<CentralityDefinition<<endl;

  // Run/Event Flag
  Int_t RunFlag;
  if(RunYear == mRun14 && !doMBset) RunFlag = StJetFrameworkPicoBase::Run14_AuAu200;
  if(RunYear == mRun14 &&  doMBset) RunFlag = StJetFrameworkPicoBase::Run14_AuAu200_MB;
  if(RunYear == mRun12 && dopp) RunFlag = StJetFrameworkPicoBase::Run12_pp200;
  Bool_t RejectBadRuns = kFALSE; // switch to load and than omit bad runs
  Int_t fBadRunListVers = StJetFrameworkPicoBase::fBadRuns_w_missing_HT; // options: fBadRuns_w_missing_HT, fBadRuns_wo_missing_HT,

  // trigger flags - update default
  Int_t EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2; // kIsHT1 or kIsHT2 or kIsHT3 (set for Run14)
  if(RunYear == mRun12) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2; 
  if(RunYear == mRun14) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2; // kIsHT2 Run14
  Int_t MBEventType = StJetFrameworkPicoBase::kVPDMB5;        // this is default (Run14)  - THIS variable is *NOT* used in this MACRO
  if(RunYear == mRun12) MBEventType = StJetFrameworkPicoBase::kRun12main; // default for Run12 pp
  Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerANY;   // kTriggerANY, kTriggerMB, kTriggerHT FIXME - runs over all triggers, USER needs to code specific cuts in their class
  if(doSTEP1 || doSTEP2) TriggerToUse = StJetFrameworkPicoBase::kTriggerHT;  // kTriggerANY, kTriggerMB, kTriggerHT,    (EP corrections to be done with kTriggerHT)
  cout<<"TriggerToUse: "<<TriggerToUse<<endl;  // for STEP1 and STEP2, since jetmaker and epmaker are all that is ran, use only HT triggered events

  // track flags
  bool usePrimaryTracks; // primary track switch
  if(RunYear == mRun12) usePrimaryTracks = kTRUE;
  if(RunYear == mRun14) usePrimaryTracks = kTRUE;
  bool doTrkEff = kTRUE; // track efficiency switch
  bool doCorrectTracksforEffBeforeJetReco = kFALSE; // THIS should only be turned on to CORRECT charged tracks for efficiency before giving to FastJet for jet reconstruction
  Int_t effType = StJetFrameworkPicoBase::kNormalPtEtaBased; // options: kNormalPtEtaBased (DEFAULT), kPtBased, kEtaBased, kHeaderArray

  // HadCorrON, uncType = doNothing  - DEFAULT!
  // systematic flags: tracking efficiency
  Int_t kSystematicUncType = kDoNothing;            // kDoNothing, kTrkEffMin, kTrkEffMax - available options 
  Int_t kSystematicUncHadCorrType = kHadCorrON;    // kHadCorrOFF, kHadCorrON, kHadCorrVAR - available options

  // mixed event flags
  Bool_t doGenBadMixEventBGcone = kFALSE; // kTRUE changes normalization and removes high frac BGcones

  // jet flags
  Int_t fJetType = kFullJet; //kFullJet, kChargedJet, kNeutralJet
  bool doCorrJetPt = (doBackgroundJets) ? kTRUE : kFALSE; // used in EP and AN makers to correct jet for underlying event (when using constituents < 2.0 GeV)
  bool doSkip1ParticleJets = kTRUE;            // only relevent to jet shape analyses
  bool doRequireJetTowFireTrig = kTRUE;        // require jet contain tower firing trigger
  bool doBiasJet = kTRUE;                      // jet shape / jet-had jets: require leading bias (requiring tow firing trig kinda overrides this)
  double fJetTrackBias = 4.0;  // leading jet constituent track bias - doBiasJet must be kTRUE 
  double fJetTowerBias = 4.0;  // leading jet constituent tower bias - doBiasJet must be kTRUE
  //double fJetRadius = 0.4;  // 0.4, 0.3, 0.2  // set up top, as this is often overlooked when submitting multiple submission in a row
  double fJetConstituentCut = 2.0;              // correlation analysis: 2.0, jet shape analysis: 1.0 (been using 2.0 for corrections)
  Int_t fJetAnalysisJetType = kInclusiveJets;     // Jet analysis jet types - options: kInclusiveJets, kLeadingJets, kSubLeadingJets, USER needs to define and use w/ setter in their class
  cout<<"JetType: "<<fJetType<<"     JetRad: "<<fJetRadius<<"     JetConstit Cut: "<<fJetConstituentCut<<endl;

  // FIXME - be aware of which list is used! 
  // tower flags - lists to load for bad towers, see StJetFrameworkPicoBase and below
  Int_t TowerListToUse = 136; // doesn't matter for charged jets - Run14 136-122: jet-hadron, early jet shape - Pt dep lists set below
  if(dopp) TowerListToUse = 169;
  // see StJetFrameworkPicoBase:   9992000 - 2 GeV, 9991000 - 1 GeV, 9990200 - 0.2 GeV  (applicable currently for Run12 pp and Run14 AuAu)
  // as of February 13, 2020, any jet constituent analysis below 0.7 is bad due to track-tower matching, need updated Pico Production - update came in April 2020
  if(fJetConstituentCut == 2.0) TowerListToUse = 9992000;
  if(fJetConstituentCut == 1.0) TowerListToUse = 9991000;
  if(fJetConstituentCut == 0.2) TowerListToUse = 9990200;
  // Run12: 1 - Raghav's list, 102 - my initial list, 169 - new list
  // Run14: 136 - main list (updated version for AuAu 200 GeV Run14), 122 - past used list
  // Run14 P18ih: 999 (initial) 
  cout<<"TowerListUsed: "<<TowerListToUse<<endl;

  // analysis type flags (booleans): TODO
  bool fDoJetShapeAnalysis = kFALSE;              // switch for doing jet shape analysis - one or the other true
  bool fDoJetHadronCorrelationAnalysis = kTRUE; // switch for doing jet-hadron correlation analysis - one or the other true

  // event plane type / configuration flags - make sure this is set properly!!!
  Int_t TPCEPSelectionType = StJetFrameworkPicoBase::kRemoveEtaStrip;
  Int_t EventPlaneTrackWeightMethod = StJetFrameworkPicoBase::kPtLinear2Const5Weight; // cuts off at 2.0 though from other cut
  bool doUseMainEPAngle = kFALSE;  // kFALSE: pt-bin approach, kTRUE: use 0.2-2.0 GeV charged tracks used for event plane reconstruction

  // update settings for new centrality definitions - certain productions had settings for z-vertex < 30 when calculating centrality definitions, etc..
  if((CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30  ||
      CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30  ||
      CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P18ih_VpdMB30_AllLumi) && (!doSTEP1 && !doSTEP2) // don't restrict if calculating event plane corrections 1&2
  ) { ZVtxMin = -30.0; ZVtxMax = 30.0; }
  cout<<"zVtxMin = "<<ZVtxMin<<"   zVtxMax = "<<ZVtxMax<<endl;

  // =============================================================================== //
  // =============================================================================== //

  // open and close output .root file (so it exist and can be updated by Analysis Tasks)
  TFile *fout = new TFile(outputFile, "RECREATE");
  fout->Close();

  // don't need event plane files for pp collisions
  if(!dopp && !doSetupQA && doEPseparate) {
    // for now not writing to these files - perhaps write to one and split
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

  // create the picoMaker maker:  (PicoIoMode, inputFile, name="picoDst")
  // - Write PicoDst's: PicoIoMode::IoWrite -> StPicoDstMaker::IoWrite
  // - Read  PicoDst's: PicoIoMode::IoRead  -> StPicoDstMaker::IoRead
  StPicoDstMaker *picoMaker = new StPicoDstMaker(StPicoDstMaker::IoRead, inputFile, "picoDst");
  picoMaker->setVtxMode((int)(StPicoDstMaker::PicoVtxMode::Default));  // PicoVtxMode:PicoVtxDefault

  // create base class maker pointer
  StJetFrameworkPicoBase *baseMaker = new StJetFrameworkPicoBase("baseClassMaker");
  baseMaker->SetRunFlag(RunFlag);                         // run flag (year)
  baseMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
  baseMaker->SetBadRunListVers(fBadRunListVers);          // switch to select specific bad run version file
  baseMaker->SetBadTowerListVers(TowerListToUse);
  cout<<baseMaker->GetName()<<endl;  // print name of class instance

  // create centrality class maker pointer
  StCentMaker *CentMaker = new StCentMaker("CentMaker", picoMaker, outputFile, kFALSE); //kTRUE);
  CentMaker->SetUsePrimaryTracks(usePrimaryTracks);       // use primary tracks
  CentMaker->SetEventZVtxRange(ZVtxMin, ZVtxMax);         // can be tighter for Run16 (-20,20)
  CentMaker->SetRunFlag(RunFlag);                         // Run Flag
  CentMaker->SetdoppAnalysis(dopp);                       // pp-analysis switch
  CentMaker->SetCentralityDef(CentralityDefinition);      // centrality definition
  CentMaker->SetUseBBCCoincidenceRate(kFALSE);            // BBC or ZDC (default) rate used? - kFALSE is default, uses ZDC coincidence
  CentMaker->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3 - *NO*T used in StCentMaker
  CentMaker->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
  cout<<CentMaker->GetName()<<endl;  // print name of class instance

  // create JetFinder first (JetMaker)
  StJetMakerTask *jetTask = new StJetMakerTask("JetMaker", fJetConstituentCut, kTRUE, outputFile);
  jetTask->SetJetType(fJetType);          // jetType
  jetTask->SetJetAlgo(antikt_algorithm);  // jetAlgo
  jetTask->SetRecombScheme(BIpt2_scheme); // recomb
  //jetTask->SetRecombScheme(E_scheme); // recomb - this scheme actually doesn't pre-process the 4-vectors during the recombination scheme to set mass to 0 - USED for jet mass
  //jetTask->SetRecombScheme(WTA_pt_scheme); // Winner-Takes-All scheme
  jetTask->SetRadius(fJetRadius);         // jet radius
  jetTask->SetJetsName("Jets");           // jets name
  jetTask->SetMinJetPt(10.0);             // signal jets
  jetTask->SetMaxJetTrackPt(30.0);        // max track constituent
  jetTask->SetMinJetTowerE(fJetConstituentCut);  // 2.0 correlations
  jetTask->SetHadronicCorrFrac(1.0);      // fractional hadronic correction
  if(kSystematicUncHadCorrType = kHadCorrOFF) jetTask->SetHadronicCorrFrac(0.0); // OFF
  jetTask->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
  if(kSystematicUncHadCorrType = kHadCorrVAR) jetTask->SetJetHadCorrType(StJetFrameworkPicoBase::kHighestEMatchedTrack);
  jetTask->SetGhostArea(0.005);
  jetTask->SetMinJetArea(0.0);
  jetTask->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // fiducial eta acceptance
  jetTask->SetJetPhiRange(0, 2.0*pi);     // phi acceptance
  jetTask->SetUsePrimaryTracks(usePrimaryTracks);
  jetTask->SetRunFlag(RunFlag);           // run flag (year)
  jetTask->SetdoppAnalysis(dopp);
  jetTask->SetEventZVtxRange(ZVtxMin, ZVtxMax);    // can be tighter for Run16 (-20,20)
  jetTask->SetTurnOnCentSelection(doCentSelection);
  jetTask->SetCentralityBinCut(CentralitySelection);
  jetTask->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
  jetTask->SetTriggerToUse(TriggerToUse);
  jetTask->SetdoConstituentSubtr(kFALSE); // constituent subtractor switch: kTRUE produces a second fJets array containing jets whom underwent constituent subtraction
  jetTask->SetRejectBadRuns(RejectBadRuns);        // switch to load and than omit bad runs
  jetTask->SetDoEffCorr(doTrkEff);       // Loads efficiency file, tells call to efficiency function to use or not use correction
  jetTask->SetDoCorrectTracksforEffBeforeJetReco(doCorrectTracksforEffBeforeJetReco); // set above, only use to correct charged tracks before jet reconstruction for efficiency
  //jetTask->SetDebugLevel(2); // 8 spits out cluster/tower stuff
  //////////// not using:  if (bFillGhosts) jetTask->SetFillGhost();

  // create JetFinder for background now (JetMakerBG) - background jets to be used in Rho Maker
  StJetMakerTask *jetTaskBG;
  if(doBackgroundJets) {
    //StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", 0.2, kTRUE, outputFile); // all inclusive
    jetTaskBG = new StJetMakerTask("JetMakerBG", fJetConstituentCut, kTRUE, outputFile);
    jetTaskBG->SetJetType(fJetType);          // jetType
    jetTaskBG->SetJetAlgo(kt_algorithm);      // jetAlgo
    jetTaskBG->SetRecombScheme(BIpt2_scheme); // recombination scheme
    //jetTaskBG->SetRecombScheme(E_scheme); // recomb - this scheme actually doesn't pre-process the 4-vectors during the recombination scheme to set mass to 0 - USED for jet mass
    //jetTaskBG->SetRecombScheme(WTA_pt_scheme); // Winner-Takes-All scheme
    jetTaskBG->SetRadius(fJetRadius);         // 0.4
    jetTaskBG->SetJetsName("JetsBG");
    jetTaskBG->SetMinJetPt(0.0);
    jetTaskBG->SetMaxJetTrackPt(30.0);
    jetTaskBG->SetMinJetTowerE(fJetConstituentCut); // inclusive: 0.2
    jetTaskBG->SetHadronicCorrFrac(1.0);
    if(kSystematicUncHadCorrType = kHadCorrOFF) jetTaskBG->SetHadronicCorrFrac(0.0); // OFF
    jetTaskBG->SetJetHadCorrType(StJetFrameworkPicoBase::kAllMatchedTracks); // options:  kLastMatchedTrack, kHighestEMatchedTrack, kAllMatchedTracks
    if(kSystematicUncHadCorrType = kHadCorrVAR) jetTaskBG->SetJetHadCorrType(StJetFrameworkPicoBase::kHighestEMatchedTrack);
    jetTaskBG->SetGhostArea(0.005);
    jetTaskBG->SetMinJetArea(0.0);
    jetTaskBG->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // -0.5,0.5
    jetTaskBG->SetJetPhiRange(0, 2.0*pi);        // 0,pi
    jetTaskBG->SetUsePrimaryTracks(usePrimaryTracks);
    jetTaskBG->SetRunFlag(RunFlag);
    jetTaskBG->SetdoppAnalysis(dopp);
    jetTaskBG->SetEventZVtxRange(ZVtxMin, ZVtxMax);         // can be tighter for Run16 (-20,20)
    jetTaskBG->SetTurnOnCentSelection(doCentSelection);
    jetTaskBG->SetCentralityBinCut(CentralitySelection);
    jetTaskBG->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
    jetTaskBG->SetTriggerToUse(TriggerToUse);
    jetTaskBG->SetdoConstituentSubtr(kFALSE);
    jetTaskBG->SetRejectBadRuns(RejectBadRuns);        // switch to load and than omit bad runs
    // TODO make sure the next two lines make sense when using this to calculate background kt jets
    jetTaskBG->SetDoEffCorr(doTrkEff);           // Loads efficiency file, tells call to efficiency function to use or not use correction
    jetTaskBG->SetDoCorrectTracksforEffBeforeJetReco(doCorrectTracksforEffBeforeJetReco); // set above, only use to correct charged tracks before jet reconstruction for efficiency
  }

  // histogram switch
  bool dohisto = kFALSE;  // histogram switch for Rho Maker

  // Rho task - for underlying event background subtraction
  // run the below tasks when running an Analysis and not strictly QA
  if(!doSetupQA) {
    StRho *rhoTask;
    if(doBackgroundJets) { rhoTask = new StRho("StRho_JetsBG", dohisto, outputFile, "JetMakerBG"); // kt jets, for background
    } else { rhoTask = new StRho("StRho_JetsBG", dohisto, outputFile, "JetMaker"); } // signal jets, bc not doing a rho subtraction
    rhoTask->SetExcludeLeadJets(2);
    rhoTask->SetOutRhoName("OutRho");
    rhoTask->SetRunFlag(RunFlag);
    rhoTask->SetdoppAnalysis(dopp);                  // pp switch
    rhoTask->SetEventZVtxRange(ZVtxMin, ZVtxMax);    // can be tighter for Run16 (-20,20)
    rhoTask->SetTurnOnCentSelection(doCentSelection);
    rhoTask->SetCentralityBinCut(CentralitySelection);
    rhoTask->SetRejectBadRuns(RejectBadRuns);        // switch to load and than omit bad runs

    //  if(ptbin == 0) { if((pt > 0.20) && (pt <= 0.5)) continue; }  // 0.20 - 0.5 GeV assoc bin used for correlations
    //  if(ptbin == 1) { if((pt > 0.50) && (pt <= 1.0)) continue; }  // 0.50 - 1.0 GeV assoc bin used for correlations
    //  if(ptbin == 2) { if((pt > 1.00) && (pt <= 1.5)) continue; }  // 1.00 - 1.5 GeV assoc bin used for correlations
    //  if(ptbin == 3) { if((pt > 1.50) && (pt <= 2.0)) continue; }  // 1.50 - 2.0 GeV assoc bin used for correlations
    //  if(ptbin == 4) { if((pt > 2.00) && (pt <= 20.)) continue; }  // 2.00 - MAX GeV assoc bin used for correlations

    if(!dopp) {
      StEventPlaneMaker *EPMaker[5];
      for(int i = 0; i < 5; i++) {
        //if(i < 4) continue; // test
        EPMaker[i] = new StEventPlaneMaker(Form("EventPlaneMaker_bin%i", i), picoMaker, "JetMaker", "StRho_JetsBG");
        EPMaker[i]->SetMinJetPt(10.0);                     // perhaps lower this TODO, but *CAN* only go as low as whats in StJetMakerTask
        EPMaker[i]->SetUsePrimaryTracks(usePrimaryTracks); // use primary tracks
        EPMaker[i]->SetCorrectJetPt(doCorrJetPt);          // subtract Rho BG from jet pt
        EPMaker[i]->SetMinTrackPt(0.2);                    // may be able to just set this to 0? - TODO
        EPMaker[i]->SetEventZVtxRange(ZVtxMin, ZVtxMax);   // can be tighter for Run16 (-20,20)
        EPMaker[i]->SetTrackPhiRange(0.0, 2.0*TMath::Pi());
        EPMaker[i]->SetTrackEtaRange(-1.0, 1.0);
        EPMaker[i]->SetDoEffCorr(doTrkEff);            // tracking efficiency switch: NOT USED currently in this class 
        EPMaker[i]->SetEventPlaneMaxTrackPtCut(2.0);   // 5.0 is default (can use < 2 to avoid autocorrelation issues)
        EPMaker[i]->SetExcludeLeadingJetsFromFit(1.0); // 1.0 is default
        EPMaker[i]->SetEventPlaneTrackWeight(EventPlaneTrackWeightMethod); // type of track weighting selection
        EPMaker[i]->SetTPCEventPlaneMethod(TPCEPSelectionType);            // TPC type method
        EPMaker[i]->SetdoEPTPCptAssocMethod(kTRUE);    // calculate TPC event plane / RES on pt assoc bin basis
        EPMaker[i]->SetEPTPCptAssocBin(i);             // pt assoc bin to use
        EPMaker[i]->SetRunFlag(RunFlag);               // run flag (year)
        EPMaker[i]->SetdoppAnalysis(dopp);             // pp analysis switch kTRUE exits as there is no event plane in pp
        EPMaker[i]->SetJetType(fJetType);              // jet type (full, charged, neutral)
        EPMaker[i]->SetJetRad(fJetRadius);             // jet radius
        EPMaker[i]->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3
        EPMaker[i]->SetTriggerToUse(TriggerToUse);               // switch to require MB or HT to run analysis
        EPMaker[i]->SetTurnOnCentSelection(doCentSelection);     // run analysis for specific centrality
        EPMaker[i]->SetCentralityBinCut(CentralitySelection);    // specific centrality range to run 
        EPMaker[i]->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
        if(doSTEP1) EPMaker[i]->SetOutFileNameEP(Form("%s", outputFile)); // set output file for recentering calibration of event plane
        if(doSTEP2) EPMaker[i]->SetOutFileNameEP(Form("%s", outputFile)); // set output file for shift calibration of event plane
        if(doSTEP3) EPMaker[i]->SetOutFileNameEP(Form("%s", outputFile)); // set output file for final event plane calculation
        EPMaker[i]->SetdoReadCalibFile(kFALSE);        // switch to read from input root file for event plane corrections
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
    StMyAnalysisMaker3 *anaMaker[9]; // was [5] for correlation analysis
    int ptAssocBins[9] = {0, 1, 2, 3, 4, 4,4,4,4};

    for(int i = 0; i < 9; i++) { // jet shape analysis
      if(doSTEP1 || doSTEP2) continue;
      if(i<8 && fDoJetShapeAnalysis) continue;             // - used for jet shape 
      if(i!=4 && fDoJetHadronCorrelationAnalysis) continue; // - used for jet-hadron correlation analysis
      // MinJetPt = 10 GeV for Jet Shape, MinJetPt = 15 GeV for jet-correlations
      anaMaker[i] = new StMyAnalysisMaker3(Form("AnalysisMaker_bin%i", i), picoMaker, outputFile, doComments, 10.0, "JetMaker", "StRho_JetsBG");
      anaMaker[i]->SetEventMixing(kTRUE);   // event mixing switch
      anaMaker[i]->SetCentBinSize(5);       // centrality bin size for mixed events
      anaMaker[i]->SetMixingTracks(25000);  // 50,000 in ALICE (2.76 TeV)
      anaMaker[i]->SetNMixedTr(1500);       // was 2500
      anaMaker[i]->SetNMixedEvt(1);         // was 5, was 2 before June11, 2018
      anaMaker[i]->SetdoGenerateBadMixEventBGcone(doGenBadMixEventBGcone); // kFALSE doesn't run cut
      anaMaker[i]->SetBGConeFractionCut(0.7); // BG cone fraction to cut on and exclude mixed event in question, 0.3 was default
      anaMaker[i]->SetDoUseMultBins(kFALSE);  // use multiplicity bins (hand defined) instead of cent bins - Jet Shape Analysis
      anaMaker[i]->SetdoUseEPBins(kTRUE);     // use event plane bins
      anaMaker[i]->SetnEPBins(4);             // number of event plane bins to use for event mixing (0, pi) range, default = 4
      if(dopp) { 
        anaMaker[i]->SetDoUseMultBins(kTRUE); // for pp
        anaMaker[i]->SetdoUseEPBins(kFALSE);  // for pp
      }
      anaMaker[i]->SetDoFilterPtMixEvents(kFALSE);              // DONT USE, filter mixed event pool by pt cut switch

      anaMaker[i]->SetCorrectJetPt(doCorrJetPt);                // subtract Rho BG, when constituents < 2.0 GeV
      anaMaker[i]->SetJetMaxTrackPt(fJetTrackBias);             // jet track bias
      anaMaker[i]->SetJetMaxTowerEt(fJetTowerBias);             // jet tower bias
      anaMaker[i]->SetJetRad(fJetRadius);                       // jet radius
      anaMaker[i]->SetJetConstituentCut(fJetConstituentCut);    // 2.0 is default 
      anaMaker[i]->SetJetLJSubLJPtThresholds(20.0, 10.0);       // LJ pt > 20.0, SubLJ pt > 10.0 GeV
      anaMaker[i]->SetdoRequireAjSelection(kFALSE);             // requirement of Aj selection on jets for jet shape analysis
      anaMaker[i]->SetEventZVtxRange(ZVtxMin, ZVtxMax);         // can be tighter for Run16 (-20,20)
      anaMaker[i]->SetUsePrimaryTracks(usePrimaryTracks);       // use primary tracks
      anaMaker[i]->SetMinTrackPt(0.2);                          // track quality cut (not related to constituents!)
      anaMaker[i]->SetTrackPhiRange(0.0, 2.0*TMath::Pi());      // track phi range
      anaMaker[i]->SetTrackEtaRange(-1.0, 1.0);                 // track eta range
      anaMaker[i]->SetDoEffCorr(doTrkEff);                      // track reco efficiency switch
      anaMaker[i]->SetTrackEfficiencyType(effType);             // tracking efficiency type:  pt-eta, pt-based, eta-based

      anaMaker[i]->SetEventPlaneMaxTrackPtCut(2.0);             // EP max track pt cut: 5.0 is default (can use < 2 to avoid autocorrelation issues)
      anaMaker[i]->SetExcludeLeadingJetsFromFit(1.0);           // 1.0 is default
      anaMaker[i]->SetEventPlaneTrackWeight(EventPlaneTrackWeightMethod);
      anaMaker[i]->SetTPCEventPlaneMethod(TPCEPSelectionType);
      anaMaker[i]->SetdoEPTPCptAssocMethod(kTRUE);              // calculate TPC event plane / RES on pt assoc bin basis - TODO keep updated
      anaMaker[i]->SetEPTPCptAssocBin(ptAssocBins[i]);//i);     // pt assoc bin to use ************
      anaMaker[i]->SetEventPlaneMakerName("EventPlaneMaker_bin");// event plane maker name - doesn't matter in current setup: the name is looped over the 5 pt bins in this AN-class
      anaMaker[i]->SetdoUseMainEPAngle(doUseMainEPAngle);       // sets EP angle to use 0.2-2.0 GeV tracks for reconstruction instead of pt-bin approach
      anaMaker[i]->SetdoEventPlaneRes(kTRUE);                   // write EP resolution profiles -- TODO updated

      anaMaker[i]->SetRunFlag(RunFlag);                         // run flag: i.e. - Run14, Run16...
      anaMaker[i]->SetdoppAnalysis(dopp);                       // running analysis over pp?
      anaMaker[i]->SetTurnOnCentSelection(doCentSelection);     // run analysis for specific centrality
      anaMaker[i]->SetCentralityBinCut(CentralitySelection);    // specific centrality range to run 
      anaMaker[i]->SetEmcTriggerEventType(EmcTriggerEventType); // kIsHT1 or kIsHT2 or kIsHT3 
      anaMaker[i]->SetdoJetHadronCorrelationAnalysis(fDoJetHadronCorrelationAnalysis); // do jet-hadron correlation analysis
      anaMaker[i]->SetdoJetShapeAnalysis(fDoJetShapeAnalysis);  // do jet shape analysis 
      anaMaker[i]->SetJetShapeTrackPtRange(0.2, 30.0);          // jet shape analysis - track pt range
      anaMaker[i]->SetJetAnalysisJetType(fJetAnalysisJetType);  // jet type for use with jet analysis
      anaMaker[i]->SetRejectBadRuns(RejectBadRuns);             // switch to load and than omit bad runs
      anaMaker[i]->SetPrintEventCounter(kFALSE);                // print event #, useful for tests or finding which event caused problem
      anaMaker[i]->SetSystematicUncType(kSystematicUncType);    // systematic uncertainty type to run/toggle

      anaMaker[i]->SetdoSkip1ParticleJets(doSkip1ParticleJets);         // do skip 1 particle jets 
      anaMaker[i]->SetdoBiasJetLeadConstituent(doBiasJet);              // bias jet with leading constituent requirement
      anaMaker[i]->SetdoRequireJetTowFireTrig(doRequireJetTowFireTrig); // require jets contain constituent tower which fired the events HT trigger

      anaMaker[i]->SetOutFileNameQA(Form("QAtracksANDjets%s.root", fEPoutJobappend));
      anaMaker[i]->SetWriteTrackQAHistograms(kFALSE);
      anaMaker[i]->SetWriteJetQAHistograms(kFALSE);

      //anaMaker[i]->SetdoRunAnalysis(kFALSE); // TODO FIXME - this is only tmp - March 2021
      //anaMaker[i]->SetEventMixing(kFALSE);   // event mixing switch

      // when running jet-hadron correlation analysis, overwrite some specific setters
      if(fDoJetHadronCorrelationAnalysis) {
        anaMaker[i]->SetDoUseMultBins(kFALSE);               // use multiplicity bins (hand defined) instead of cent bins - Jet Shape Analysis
        anaMaker[i]->SetdoUseEPBins(kFALSE);               // use event plane bins
        //anaMaker[i]->SetDoUseMultBins(kTRUE);             // use multiplicity bins (hand defined) instead of cent bins - Jet Shape Analysis
        //anaMaker[i]->SetdoUseEPBins(kTRUE);                  // use event plane bins FIXME - not used for prelim
        anaMaker[i]->SetCentBinSize(10);                      // centrality bin size for mixed events - 
        anaMaker[i]->SetJetAnalysisJetType(kInclusiveJets);  // jet type for use with jet analysis
      }

/*
      // when running jet shape analysis, overwrite some specific setters
      if(fDoJetShapeAnalysis) {
        anaMaker[i]->SetDoUseMultBins(kFALSE);               // use multiplicity bins (hand defined) instead of cent bins - Jet Shape Analysis
        anaMaker[i]->SetdoUseEPBins(kTRUE);                  // use event plane bins
        anaMaker[i]->SetCentBinSize(5);                      // centrality bin size for mixed events - NOT currently used
        anaMaker[i]->SetJetAnalysisJetType(fJetAnalysisJetType);  // jet type for use with jet analysis
      }
*/

      //anaMaker[i]->SetDebugLevel(StMyAnalysisMaker3::kDebugTowersFiringTriggers);
      cout<<anaMaker[i]->GetName()<<endl;                       // print name of class instance
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
    if (iret) { cout << "Bad return code!" << iret << endl; break; }

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
  if(!dopp && !doSetupQA && doEPseparate)    {  if(fEPout->IsOpen()) fEPout->Close();  }
  if(doAnalysisQAoutputFile)                 {  if(fQAout->IsOpen()) fQAout->Close();  }

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
  gSystem->AddIncludePath("-I$FASTJET/include");

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

void PrintSetup()
{
}
