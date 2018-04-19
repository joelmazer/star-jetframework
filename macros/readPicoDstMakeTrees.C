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
//class StJetFrameworkPicoBase;
class StMyAnalysisMaker;

// library and macro loading function
void LoadLibs();
void LoadMacros();

// ======= constants and switches ============
const double pi = 1.0*TMath::Pi();

// kTRUE for Run14, kFALSE for Run16
bool usePrimaryTracks = kFALSE;

// keep kFALSE for all centralities, kTRUE runs analysis for specific cent range defined below
bool doCentSelection = kFALSE;

// bool dohisto = kFALSE;  // histogram switch for Rho Maker (SOME PROBLEM WITH THIS - FIXME)

Double_t ZVtxMin = -20.0;
Double_t ZVtxMax = 20.0;

StChain *chain;

void readPicoDstMakeTrees(const Char_t *inputFile="test_run17061011.list", const Char_t *outputFile="test_CreateTree.root")
{
        Int_t nEvents = 100;

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

        Int_t centralitySelection = StJetFrameworkPicoBase::kCent2050;
	
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

        // create JetFinder first (JetMaker)
        // makerName, constituent cut, doHistos?, outputFileName
        StJetMakerTask *jetTask = new StJetMakerTask("JetMaker", 2.0, kTRUE, outputFile);
        jetTask->SetJetType(kFullJet);          //kChargedJet); //jetType
        jetTask->SetJetAlgo(antikt_algorithm);  //jetAlgo
        jetTask->SetRecombScheme(BIpt2_scheme); //recomb
        jetTask->SetRadius(0.4);
        jetTask->SetJetsName("Jets");
        jetTask->SetMinJetPt(10.0); // 15.0 signal jets
        jetTask->SetMaxJetTrackPt(20.0);
        jetTask->SetGhostArea(0.005);
        jetTask->SetMinJetArea(0.0);
        jetTask->SetJetEtaRange(-0.6,0.6);
        jetTask->SetJetPhiRange(0,2*pi); 
        jetTask->SetUsePrimaryTracks(usePrimaryTracks);
        jetTask->SetRunFlag(StJetFrameworkPicoBase::Run16_AuAu200);
        jetTask->SetCentralityDef(StJetFrameworkPicoBase::kgrefmult_P16id);
        jetTask->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
        jetTask->SetTurnOnCentSelection(doCentSelection);
        jetTask->SetCentralityBinCut(centralitySelection);
        //if (bFillGhosts) jetTask->SetFillGhost();

        // create JetFinder for background now (JetMakerBG) - to be used in Rho Maker
        // makerName, constituent cut, doHistos?, outputFileName
        //StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", 0.2, kTRUE, outputFile); // all inclusive
        StJetMakerTask *jetTaskBG = new StJetMakerTask("JetMakerBG", 2.0, kTRUE, outputFile); // TEST
        jetTaskBG->SetJetType(kChargedJet);       // jetType
        jetTaskBG->SetJetAlgo(kt_algorithm);      // jetAlgo
        jetTaskBG->SetRecombScheme(BIpt2_scheme); // reco
        jetTaskBG->SetRadius(0.4);                // 0.4
        jetTaskBG->SetJetsName("JetsBG");
        jetTaskBG->SetMinJetPt(0.0);
        jetTaskBG->SetMaxJetTrackPt(20.0);
        jetTaskBG->SetGhostArea(0.005);
        jetTaskBG->SetMinJetArea(0.0);
        jetTaskBG->SetJetEtaRange(-0.6,0.6); 
        jetTaskBG->SetJetPhiRange(0,2*pi);       
        jetTaskBG->SetUsePrimaryTracks(usePrimaryTracks);
        jetTaskBG->SetRunFlag(StJetFrameworkPicoBase::Run16_AuAu200);
        jetTaskBG->SetCentralityDef(StJetFrameworkPicoBase::kgrefmult_P16id);
        jetTaskBG->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
        jetTaskBG->SetTurnOnCentSelection(doCentSelection);
        jetTaskBG->SetCentralityBinCut(centralitySelection);

        // this is the centrality dependent scaling for RHO as used in ALICE - don't use right now
        // s(Centrality) = 0.00015 ×Centrality^2 ? 0.016 ×Centrality + 1.91
        // scaled function for rho and setting parameters
        TF1* sfunc = new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
        sfunc->SetParameter(2,1.80642);
        sfunc->SetParameter(1,-0.0112331);
        sfunc->SetParameter(0,0.0001139561);
  
        // write Rho histograms
        bool dohisto = kFALSE;  // histogram switch for Rho Maker

        // Rho task, and scale it up to include neutral constituents
        // MakerName, doHistograms?, outputFileName, JetMakerBGName
        StRho *rhoTask = new StRho("StRho_JetsBG", dohisto, outputFile, "JetMakerBG"); // background jets
        rhoTask->SetExcludeLeadJets(2);
        rhoTask->SetOutRhoName("OutRho");
        rhoTask->SetRunFlag(StJetFrameworkPicoBase::Run16_AuAu200);
        rhoTask->SetCentralityDef(StJetFrameworkPicoBase::kgrefmult_P16id);
        rhoTask->SetEventZVtxRange(ZVtxMin, ZVtxMax); // can be tighter for Run16 (-20,20)
        rhoTask->SetTurnOnCentSelection(doCentSelection);
        rhoTask->SetCentralityBinCut(centralitySelection);
        //rhoTask->SetScaleFunction(sfunc); // don't NEED for STAR

	// create the analysis maker!
	// MakerName, picoMakerPtr, outputFileName, JetMakerName, RhoMakerName
	StAnMaker *anaMaker = new StAnMaker("Trees", picoMaker, outputFile, "JetMaker", "StRho_JetsBG");
        anaMaker->SetUsePrimaryTracks(usePrimaryTracks); // use primary tracks
        anaMaker->SetCorrectJetPt(kFALSE);               // subtract Rho BG
        anaMaker->SetMinJetPt(10.0);
        anaMaker->SetJetMaxTrackPt(0.0);                 // track bias (probably don't want to use for trees?, so use 0.0)
        anaMaker->SetMinTrackPt(0.2);
        anaMaker->SetEventZVtxRange(ZVtxMin, ZVtxMax);   // can be tighter for Run16 (-20,20)
        anaMaker->SetTrackPhiRange(0.0, 2*TMath::Pi());
        anaMaker->SetTrackEtaRange(-1.0, 1.0);
        anaMaker->SetJetConstituentCut(2.0);             // 2.0 is default 
        anaMaker->SetRunFlag(StJetFrameworkPicoBase::Run16_AuAu200);
        anaMaker->SetCentralityDef(StJetFrameworkPicoBase::kgrefmult_P16id);
        anaMaker->SetEmcTriggerEventType(StJetFrameworkPicoBase::kIsHT1);  // kIsHT1 or kIsHT2
        anaMaker->SetTurnOnCentSelection(doCentSelection);              // run analysis for specific centrality
        anaMaker->SetCentralityBinCut(centralitySelection);             // specific centrality range to run 


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
  gSystem->AddIncludePath("-I/global/homes/j/jmazer/STAR/Y2017/mytestinstalldir/FastJet/fastjet-install/include");

  // load the system libraries - these were defaults
  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  // these are needed for new / additional classes
  //gSystem->Load("StPicoDstMaker");
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
