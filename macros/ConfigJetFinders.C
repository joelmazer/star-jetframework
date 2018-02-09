StJetMakerTask* ConfigJetFinders(
  const char *nTracks                        = "usedefault",
  const char *nClusters                      = "usedefault",
//  const StMyAnalysisMaker::EJetAlgo_t jetAlgo  = StMyAnalysisMaker::antikt_algorithm,
  const EJetAlgo_t jetAlgo  = antikt_algorithm,
  const Double_t radius                      = 0.4,
  const EJetType_t jetType  = kFullJet,
  const Double_t minTrPt                     = 0.15,
  const Double_t minClPt                     = 0.30,
  const Double_t ghostArea                   = 0.005,
  const ERecoScheme_t reco  = pt_scheme,
  const char *tag                            = "Jet",
  const Double_t minJetPt                    = 0.,
  const Bool_t lockTask                      = kTRUE,
  const Bool_t bFillGhosts                   = kFALSE
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------
  TString trackName(nTracks);
  TString clusName(nClusters);

  TString name = StMyAnalysisMaker::GenerateJetName(jetType, jetAlgo, reco, radius, partCont, clusCont, tag);
  Printf("Jet task name: %s", name.Data());
 
  //StJetMakerTask* jetTask = new StJetMakerTask(name);
  StJetMakerTask* jetTask = new StJetMakerTask();
  jetTask->SetJetType(jetType);
  jetTask->SetJetAlgo(jetAlgo);
  jetTask->SetRecombScheme(reco);
  jetTask->SetRadius(radius);
//  if (partCont) jetTask->AdoptParticleContainer(partCont);
//  if (clusCont) jetTask->AdoptClusterContainer(clusCont);
  jetTask->SetJetsName(tag);
  jetTask->SetMinJetPt(minJetPt);
  jetTask->SetGhostArea(ghostArea);

//  if (bFillGhosts) jetTask->SetFillGhost();
//  if (lockTask) jetTask->SetLocked();

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------
/*
  mgr->AddTask(jetTask);
  
  // Create containers for input/output
  AliAnalysisDataContainer* cinput = mgr->GetCommonInputContainer();
  mgr->ConnectInput(jetTask, 0, cinput);

  TObjArray* cnt = mgr->GetContainers();
*/
  return jetTask;
}
