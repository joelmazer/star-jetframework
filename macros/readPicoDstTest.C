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
class StRho;
class StRhoBase;
class StRhoSparse;
class StMyAnalysisMaker;

// library and macro loading function
void LoadLibs();
void LoadMacros();

// constants
const double pi = 1.0*TMath::Pi();

// run analysis for specific centrality bin
bool doCentSelection = kFALSE; //kTRUE; // keep false to run over all centralities

// see enumerator below, 4 = Run14
Int_t RunYear = 4; 
// kFALSE when submitting jobs, kTRUE for tests
bool doTEST = kTRUE;   // FIXME double check before submitting!

Double_t ZVtxMin = -40.0;
Double_t ZVtxMax = 40.0;

StChain *chain;

void readPicoDstTest(const Char_t *inputFile="", const Char_t *outputFile="test.root", Int_t nEv = 10, const Char_t *fEPoutJobappend="")
{
        Int_t nEvents = 1000;
        if(nEv > 100) nEvents = 100000000;

        // run enumerators
        enum RunType_t {
          mJobSubmission = 0, // 0th spot - default unless set
          mRun11 = 1,
          mRun12 = 2,
          mRun13 = 3,
          mRun14 = 4,
          mRun15 = 5,
          mRun16 = 6,
          mRun17 = 7,
          mRun18 = 8
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
          pt_scheme       = 1,
          pt2_scheme      = 2,
          Et_scheme       = 3,
          Et2_scheme      = 4,
          BIpt_scheme     = 5,
          BIpt2_scheme    = 6,
          WTA_pt_scheme   = 7,
          WTA_modp_scheme = 8,
          external_scheme = 99
        };

        // Load necessary libraries and macros
        LoadLibs();
        LoadMacros();

        // =============================================================================== //
        // input file for tests (based on Run)
        if((RunYear == mRun14) && doTEST) inputFile = "Run_15151042_files.list";;
        if((RunYear == mRun16) && doTEST) inputFile = "test_run17124003_files.list";
        cout<<"inputFileName = "<<inputFile<<endl;

        // centrality global flags
        Int_t CentralitySelection = StJetFrameworkPicoBase::kCent2050;
        Int_t CentralityDefinition;
        if(RunYear == mRun16) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_P16id; // Run16
        //if(RunYear == mRun16) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult_VpdMBnoVtx;
        if(RunYear == mRun14) CentralityDefinition = StJetFrameworkPicoBase::kgrefmult;       // Run14
        cout<<"Centrality definition: "<<CentralityDefinition<<endl;

        // Run/Event Flag
        Int_t RunFlag;
        if(RunYear == mRun16) RunFlag = StJetFrameworkPicoBase::Run16_AuAu200;
        if(RunYear == mRun14) RunFlag = StJetFrameworkPicoBase::Run14_AuAu200;

        // trigger flags
        Int_t EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT2; // kIsHT1 or kIsHT2 or kIsHT3
        if(RunYear == mRun16) EmcTriggerEventType = StJetFrameworkPicoBase::kIsHT1; // kIsHT1 Run16
        Int_t MBEventType = StJetFrameworkPicoBase::kVPDMB5;        // this is default
        Int_t TriggerToUse = StJetFrameworkPicoBase::kTriggerHT;    // kTriggerANY, kTriggerMB, kTriggerHT
        Int_t TowerListToUse = 122; //3; // been using 51,  3, 79   - doesn't matter for charged jets

        // track flags
        bool usePrimaryTracks;
        if(RunYear == mRun14) usePrimaryTracks = kTRUE;  // = kTRUE for Run14, kFALSE for Run16
        if(RunYear == mRun16) usePrimaryTracks = kFALSE; 

        // jet type flags
        Int_t fJetType;
        if(RunYear == mRun14) fJetType = kFullJet; // kChargedJet
        //fJetType = kChargedJet; // running Run14 analaysis for charged jets
        if(RunYear == mRun16) fJetType = kChargedJet;

        // event plane type
        //StJetFrameworkPicoBase::kRemoveEtaPhiCone); // kRemoveEtaStrip is default
        //StJetFrameworkPicoBase::kRemoveEtaPhiConeLeadSub);
        //StJetFrameworkPicoBase::kRemoveEtaStrip);
        //Int_t TPCEPSelectionType = StJetFrameworkPicoBase::kRemoveEtaStrip;
        Int_t TPCEPSelectionType = StJetFrameworkPicoBase::kRemoveEtaPhiCone;
        Int_t EventPlaneTrackWeightMethod = StJetFrameworkPicoBase::kPtLinear2Const5Weight;
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
        StJetMakerTask *jetTask = new StJetMakerTask("JetMaker", 0.2, kTRUE, outputFile);
        jetTask->SetJetType(fJetType);          //jetType
        jetTask->SetJetAlgo(antikt_algorithm);  //jetAlgo
        jetTask->SetRecombScheme(BIpt2_scheme); //recomb
        jetTask->SetRadius(0.4);
        jetTask->SetJetsName("Jets");
        jetTask->SetMinJetPt(10.0);
        jetTask->SetMaxJetTrackPt(20.0);
        jetTask->SetMinJetTowerE(2.0);
        jetTask->SetHadronicCorrFrac(1.0);
        jetTask->SetGhostArea(0.005);
        jetTask->SetMinJetArea(0.0);
        jetTask->SetJetEtaRange(-0.6,0.6);
        jetTask->SetJetPhiRange(0,2*pi); 
        jetTask->SetUsePrimaryTracks(usePrimaryTracks);
        //jetTask->SetDebugLevel(2); // 8 spits out cluster/tower stuff
        jetTask->SetRunFlag(RunFlag);
        jetTask->SetCentralityDef(CentralityDefinition);
        jetTask->SetEventZVtxRange(ZVtxMin, ZVtxMax);          // can be tighter for Run16 (-20,20)
        jetTask->SetTurnOnCentSelection(doCentSelection);
        jetTask->SetCentralityBinCut(CentralitySelection);
        jetTask->SetEmcTriggerEventType(EmcTriggerEventType);  // kIsHT1 or kIsHT2 or kIsHT3
        jetTask->SetTriggerToUse(TriggerToUse);
        jetTask->SetBadTowerListVers(TowerListToUse);


        // background jets to be used in Rho Maker
        // create JetFinder for background now (JetMakerBG)
        StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", 0.2, kTRUE, outputFile);   // all inclusive
        jetTaskBG->SetJetType(fJetType); //jetType
        jetTaskBG->SetJetAlgo(kt_algorithm); //jetAlgo
        jetTaskBG->SetRecombScheme(BIpt2_scheme); //reco
        jetTaskBG->SetRadius(0.4);
        jetTaskBG->SetJetsName("JetsBG");
        jetTaskBG->SetMinJetPt(0.0);
        jetTaskBG->SetMaxJetTrackPt(20.0);
        jetTaskBG->SetMinJetTowerE(0.2);
        jetTaskBG->SetHadronicCorrFrac(1.0);
        jetTaskBG->SetGhostArea(0.005);
        jetTaskBG->SetMinJetArea(0.0);
        jetTaskBG->SetJetEtaRange(-0.6,0.6);
        jetTaskBG->SetJetPhiRange(0,2*pi);
        jetTaskBG->SetUsePrimaryTracks(usePrimaryTracks);
        jetTaskBG->SetRunFlag(RunFlag);
        jetTaskBG->SetCentralityDef(CentralityDefinition);
        jetTaskBG->SetEventZVtxRange(ZVtxMin, ZVtxMax);        // can be tighter for Run16 (-20,20)
        jetTaskBG->SetTurnOnCentSelection(doCentSelection);
        jetTaskBG->SetCentralityBinCut(CentralitySelection);
        jetTaskBG->SetEmcTriggerEventType(EmcTriggerEventType);  // kIsHT1 or kIsHT2 or kIsHT3
        jetTaskBG->SetTriggerToUse(TriggerToUse);
        jetTaskBG->SetBadTowerListVers(TowerListToUse);

        // this is the centrality dependent scaling for RHO as used in ALICE - don't use right now
        // DON'T need this for STAR due to the same acceptance range of towers and tracks
        // s(Centrality) = 0.00015 ×Centrality^2 ? 0.016 ×Centrality + 1.91
        // scaled function for rho and setting parameters
        TF1* sfunc = new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
        sfunc->SetParameter(2,1.80642);
        sfunc->SetParameter(1,-0.0112331);
        sfunc->SetParameter(0,0.0001139561);
  
        bool dohisto = kFALSE;  // histogram switch for Rho Maker

        // Rho task, and scale it up to include neutral constituents
        StRho *rhoTask = new StRho("StRho_JetsBG", dohisto, outputFile, "JetMakerBG"); // background jets
        //StRho *rhoTask = new StRho("StRho_JetsBG", dohisto, outputFile, "JetMaker");   // signal jets, bc not using bg jets
        rhoTask->SetExcludeLeadJets(2);
        rhoTask->SetOutRhoName("OutRho");
        rhoTask->SetRunFlag(RunFlag);
        rhoTask->SetCentralityDef(CentralityDefinition);
        rhoTask->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
        rhoTask->SetTurnOnCentSelection(doCentSelection);
        rhoTask->SetCentralityBinCut(CentralitySelection);
        //rhoTask->SetScaleFunction(sfunc); // don't NEED

/*
        // Rho Sparse - for Ultra-peripheral centrality events
        StRhoSparse *rhoTaskSparse = new StRhoSparse("StRhoSparse_JetsBG", kTRUE, outputFile, Form("%s", BGjetsName));
        rhoTaskSparse->SetRunFlag(RunFlag);
        rhoTaskSparse->SetCentralityDef(CentralityDefinition);
        rhoTaskSparse->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
        rhoTaskSparse->SetExcludeLeadJets(2);
        rhoTaskSparse->SetJetMakerName("JetMaker");
        rhoTaskSparse->SetJetBGMakerName("JetMakerBG");
        rhoTaskSparse->SetDebugLevel(0); // 0 does nothing
        rhoTaskSparse->SetTurnOnCentSelection(doCentSelection);
        rhoTaskSparse->SetCentralityBinCut(CentralitySelection);
        //rhoTaskSparse->SetScaleFunction(sfunc); // don't NEED
*/

        // create the analysis maker!
        // MakerName, picoMakerPtr, outputFileName, JetMakerName, RhoMakerName
        StAnMaker *anaMaker = new StAnMaker("AnMaker", picoMaker, outputFile, "JetMaker", "StRho_JetsBG");
        anaMaker->SetUsePrimaryTracks(usePrimaryTracks); // use primary tracks
        anaMaker->SetCorrectJetPt(kFALSE);               // subtract Rho BG
        anaMaker->SetMinJetPt(10.0);
        anaMaker->SetJetMaxTrackPt(0.0);                 // jet track bias (probably don't want to use for trees?, so use 0.0)
        anaMaker->SetJetMaxTowerE(0.0);                  // jet tower bias
        anaMaker->SetMinTrackPt(0.2);                    // track quality cut (not related to constituents!)
        anaMaker->SetEventZVtxRange(ZVtxMin, ZVtxMax);   // can be tighter for Run16 (-20,20)
        anaMaker->SetTrackPhiRange(0.0, 2*TMath::Pi());
        anaMaker->SetTrackEtaRange(-1.0, 1.0);
        anaMaker->SetJetConstituentCut(2.0);             // 2.0 is default 
        anaMaker->SetRunFlag(RunFlag);
        anaMaker->SetCentralityDef(CentralityDefinition);     // Centrality Definition
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
  gSystem->AddIncludePath("-I/global/homes/j/jmazer/STAR/Y2017/mytestinstalldir/FastJet/fastjet-install/include");

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
