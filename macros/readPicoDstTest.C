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

// see enumerator below, 4 = Run14
Int_t RunYear = 4; 
// kFALSE when submitting jobs, kTRUE for tests
bool doTEST = kTRUE;     // FIXME double check before submitting!
bool dopp = kFALSE;      // kTRUE for pp data

// z-vertex cuts, go down to {-30, 30} for new centrality definitons - done below in setup section
Double_t ZVtxMin = -40.0;
Double_t ZVtxMax = 40.0;

StChain *chain;

void readPicoDstTest(const Char_t *inputFile="", const Char_t *outputFile="test.root", Int_t nEv = 10, const Char_t *fOutJobappend="")
{
        Int_t nEvents = 1000;
        if(nEv > 100) nEvents = 100000000;

        // run enumerators
        enum RunType_t {
          mJobSubmission = 0, // 0th spot - default unless set
          mRun11 = 1, mRun12 = 2, mRun13 = 3, mRun14 = 4, 
          mRun15 = 5, mRun16 = 6, mRun17 = 7, mRun18 = 8
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
        // over-ride functions
        if(dopp) doCentSelection = kFALSE; // can't ask for a particular centrality if requesting pp collisions

        // input file for tests (based on Run)
        if((RunYear == mRun14) && doTEST) inputFile = "Run_15151042_files.list";;
        if((RunYear == mRun16) && doTEST) inputFile = "test_run17124003_files.list";
        if((RunYear == mRun17) && doTEST && dopp) inputFile = "filelist_pp2017.list";
        cout<<"inputFileName = "<<inputFile<<endl;

        // centrality global flags
        // Centrality Selection can be set up in certain classes, to only run for a given centrality range
        Int_t CentralitySelection = StJetFrameworkPicoBase::kCent2050;
        Int_t CentralityDefinition;
        //if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult;       // Run14
        if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30; // Run14 P17id (NEW - from Nick Oct 23)
        if(RunYear == mRun16) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P16id; // Run16
        cout<<"Centrality definition: "<<CentralityDefinition<<endl;

        // Run/Event Flag
        Int_t RunFlag;
        if(RunYear == mRun14) RunFlag = StJetFrameworkPicoBase::Run14_AuAu200;
        if(RunYear == mRun16) RunFlag = StJetFrameworkPicoBase::Run16_AuAu200;
        if(RunYear == mRun11 && dopp) RunFlag = StJetFrameworkPicoBase::Run11_pp500;
        if(RunYear == mRun13 && dopp) RunFlag = StJetFrameworkPicoBase::Run13_pp510;
        if(RunYear == mRun17 && dopp) RunFlag = StJetFrameworkPicoBase::Run17_pp510;

        // trigger flags
        // EmcTriggerEventType: selects which HT trigger to use
        Int_t EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2; // kIsHT1 or kIsHT2 or kIsHT3 (Set for Run14)
        if(RunYear == mRun16) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT1; // kIsHT1 Run16
        if(RunYear == mRun17) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT3; 
        Int_t MBEventType = StJetFrameworkPicoBase::kVPDMB5;        // this is default, want kVPDMB30 for new centrality definitions
        if(RunYear == mRun17) MBEventType = StJetFrameworkPicoBase::kVPDMB; // default for Run17 pp
        Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerHT;    // kTriggerANY, kTriggerMB, kTriggerHT
        Int_t TowerListToUse = 136; //122: - doesn't matter for charged jets

        // track flags
        bool usePrimaryTracks;
        if(RunYear == mRun14) usePrimaryTracks = kTRUE;  // = kTRUE for Run14, kFALSE for Run16
        if(RunYear == mRun16) usePrimaryTracks = kFALSE; 
        if(RunYear == mRun17) usePrimaryTracks = kTRUE;  

        // jet type flags
        Int_t fJetType;
        if(RunYear == mRun14) fJetType = kFullJet; // kChargedJet
        //fJetType = kChargedJet; // running Run14 analaysis for charged jets
        if(RunYear == mRun16) fJetType = kChargedJet;
        if(RunYear == mRun17) fJetType = kFullJet;
        double fJetRadius = 0.3;  // 0.4, 0.3, 0.2
        double fJetConstituentCut = 0.2; // correlation analysis: 2.0, jet shape analysis: 1.0 (been using 2.0 for corrections)
        Int_t fJetShapeAnalysisJetType = kLeadingJets;  // Jet shape jet types - options: kInclusiveJets, kLeadingJets, kSubLeadingJets

        // event plane type and selection
        Int_t TPCEPSelectionType = StJetFrameworkPicoBase::kRemoveEtaStrip;
        Int_t EventPlaneTrackWeightMethod = StJetFrameworkPicoBase::kPtLinear2Const5Weight;

        // update settings for new centrality definitions
        if(CentralityDefinition == StJetFrameworkPicoBase::kgrefmult_P17id_VpdMB30) { ZVtxMin = -30.0; ZVtxMax = 30.0; }

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

        // create JetFinder first (JetMaker)
        // task Name, track constituent cut, doHistos, output file name
        StJetMakerTask *jetTask = new StJetMakerTask("JetMaker", fJetConstituentCut, kTRUE, outputFile);
        jetTask->SetJetType(fJetType);          //jetType
        jetTask->SetJetAlgo(antikt_algorithm);  //jetAlgo
        jetTask->SetRecombScheme(BIpt2_scheme); //recomb
        jetTask->SetRadius(fJetRadius);
        jetTask->SetJetsName("Jets");
        jetTask->SetMinJetPt(10.0);             // 15.0 signal jets
        jetTask->SetMaxJetTrackPt(30.0);        // max track constituent
        jetTask->SetMinJetTowerE(fJetConstituentCut);  // 2.0 correlations
        jetTask->SetHadronicCorrFrac(1.0);      // fractional hadronic correction
        jetTask->SetGhostArea(0.005);
        jetTask->SetMinJetArea(0.0);
        jetTask->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // fiducial eta acceptance
        jetTask->SetJetPhiRange(0,2*pi);        // phi acceptance
        jetTask->SetUsePrimaryTracks(usePrimaryTracks);
        jetTask->SetRunFlag(RunFlag);
        jetTask->SetCentralityDef(CentralityDefinition);  // run based centrality definition
        jetTask->SetUseBBCCoincidenceRate(kFALSE);
        jetTask->SetEventZVtxRange(ZVtxMin, ZVtxMax);     // can be tighter for Run16 (-20,20)
        jetTask->SetTurnOnCentSelection(doCentSelection);
        jetTask->SetCentralityBinCut(CentralitySelection);
        jetTask->SetUseBBCCoincidenceRate(kFALSE);
        jetTask->SetEmcTriggerEventType(EmcTriggerEventType);  // kIsHT1 or kIsHT2 or kIsHT3
        jetTask->SetTriggerToUse(TriggerToUse);
        jetTask->SetBadTowerListVers(TowerListToUse);
        jetTask->SetdoConstituentSubtr(kFALSE); // implement constituent subtractor, if TRUE, don't want to subtract underlying event Rho

        // create JetFinder for background now (JetMakerBG) - background jets to be used in Rho Maker
        StJetMakerTask *jetTaskBG;
        if(doBackgroundJets) {
          //StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", 0.2, kTRUE, outputFile); // all inclusive
          // task Name, track constituent cut, doHistos, output file name
          jetTaskBG = new StJetMakerTask("JetMakerBG", fJetConstituentCut, kTRUE, outputFile); 
          jetTaskBG->SetJetType(fJetType);          // jetType
          jetTaskBG->SetJetAlgo(kt_algorithm);      // jetAlgo
          jetTaskBG->SetRecombScheme(BIpt2_scheme); // recombination scheme
          jetTaskBG->SetRadius(fJetRadius);         // 0.4
          jetTaskBG->SetJetsName("JetsBG");
          jetTaskBG->SetMinJetPt(0.0);
          jetTaskBG->SetMaxJetTrackPt(30.0);
          jetTaskBG->SetMinJetTowerE(fJetConstituentCut);     // inclusive: 0.2
          jetTaskBG->SetHadronicCorrFrac(1.0);
          jetTaskBG->SetGhostArea(0.005);
          jetTaskBG->SetMinJetArea(0.0);
          jetTaskBG->SetJetEtaRange(-1.0 + fJetRadius, 1.0 - fJetRadius); // -0.5,0.5
          jetTaskBG->SetJetPhiRange(0, 2*pi);   // 0,pi
          jetTaskBG->SetUsePrimaryTracks(usePrimaryTracks);
          jetTaskBG->SetRunFlag(RunFlag);
          jetTaskBG->SetdoppAnalysis(dopp);
          jetTaskBG->SetCentralityDef(CentralityDefinition);
          jetTaskBG->SetEventZVtxRange(ZVtxMin, ZVtxMax);        // can be tighter for Run16 (-20,20)
          jetTaskBG->SetTurnOnCentSelection(doCentSelection);
          jetTaskBG->SetCentralityBinCut(CentralitySelection);
          jetTaskBG->SetUseBBCCoincidenceRate(kFALSE);
          jetTaskBG->SetEmcTriggerEventType(EmcTriggerEventType);  // kIsHT1 or kIsHT2 or kIsHT3
          jetTaskBG->SetTriggerToUse(TriggerToUse);
          jetTaskBG->SetBadTowerListVers(TowerListToUse);
          jetTaskBG->SetdoConstituentSubtr(kFALSE);
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
        rhoTask->SetCentralityDef(CentralityDefinition);
        rhoTask->SetUseBBCCoincidenceRate(kFALSE);
        rhoTask->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
        rhoTask->SetTurnOnCentSelection(doCentSelection);
        rhoTask->SetCentralityBinCut(CentralitySelection);
        rhoTask->SetUseBBCCoincidenceRate(kFALSE);

/*
        // Rho Sparse - for Ultra-peripheral centrality events
        StRhoSparse *rhoTaskSparse = new StRhoSparse("StRhoSparse_JetsBG", kTRUE, outputFile, Form("%s", BGjetsName));
        rhoTaskSparse->SetRunFlag(RunFlag);
        rhoTaskSparse->SetdoppAnalysis(dopp);
        rhoTaskSparse->SetCentralityDef(CentralityDefinition);
        rhoTaskSparse->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
        rhoTaskSparse->SetExcludeLeadJets(2);
        rhoTaskSparse->SetJetMakerName("JetMaker");
        rhoTaskSparse->SetJetBGMakerName("JetMakerBG");
        rhoTaskSparse->SetDebugLevel(0); // 0 does nothing
        rhoTaskSparse->SetTurnOnCentSelection(doCentSelection);
        rhoTaskSparse->SetCentralityBinCut(CentralitySelection);
        rhoTaskSparse->SetUseBBCCoincidenceRate(kFALSE);
        //rhoTaskSparse->SetScaleFunction(sfunc); // don't NEED
*/

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
        anaMaker->SetTrackPhiRange(0.0, 2*TMath::Pi());  // track phi acceptance
        anaMaker->SetTrackEtaRange(-1.0, 1.0);           // track eta acceptance
        anaMaker->SetJetConstituentCut(fJetConstituentCut);   // 0.2 for inclusive 
        anaMaker->SetRunFlag(RunFlag);                        // run flag: i.e. - Run14, Run16...
        anaMaker->SetdoppAnalysis(dopp);                      // running analysis over pp?
        anaMaker->SetCentralityDef(CentralityDefinition);     // Centrality Definition
        anaMaker->SetUseBBCCoincidenceRate(kFALSE);           // use ZDC default - KEEP false
        anaMaker->SetTurnOnCentSelection(doCentSelection);    // run analysis for specific centrality: BOOLEAN
        anaMaker->SetCentralityBinCut(CentralitySelection);   // specific centrality range to run: if above is FALSE, this doesn't matter
        anaMaker->SetEmcTriggerEventType(EmcTriggerEventType);// kIsHT1 or kIsHT2 or kIsHT3
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
//  gSystem->Load("StPicoDstMaker");
  //gSystem->Load("StJetMakerTask");
  //gSystem->Load("StRho");
  gSystem->Load("libStPicoEvent");
  gSystem->Load("libStPicoDstMaker");

  // my libraries
  //gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  //gSystem->Load("StJetFrameworkPicoBase");
  gSystem->Load("StMyAnalysisMaker");

  gSystem->ListLibraries();
} 

void LoadMacros()
{
}
