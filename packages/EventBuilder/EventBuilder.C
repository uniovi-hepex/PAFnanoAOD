//////////////////////////////////////////////////////////////////////////////////////
//
//  EventBuilder: compute all important variables related with the event
//  Triggers, PU reweighting, global SFs and systematics...
//
//  To do: Add trigger SFs when become available
//         Add LHE weights when become available
//         Update PU weights
//
//////////////////////////////////////////////////////////////////////////////////////


#include "EventBuilder.h"

ClassImp(EventBuilder);
EventBuilder::EventBuilder() : PAFChainItemSelector(),
			       METfilters(false),
			       passTrigger(false),
			       isSS(false),
			       gIsFastSim(false),
			       TriggerSF(0),
			       TriggerSF_Up(0),
			       TriggerSF_Down(0),
			       TriggerSF_err(0),
			       PUSF(0),
			       PUSF_Up(0),
			       PUSF_Down(0),
			       NormWeight(0),
			       Weight(0),
			       genWeight(0),
			       nTrueInt(0),
			       gChannel(0),
			       //fPUWeight(0),
			       //fPUWeightUp(0),
			       //fPUWeightDown(0),
			       gIsSingleMuon(false),
			       gIsSingleElec(false),
			       gIsDoubleMuon(false),
			       gIsDoubleElec(false),
			       gIsMuonEG(false),
			       gIsData(false),
			       run(-1),
			       gSelection(-1),
			       gSampleName(""),
			       gPathToHeppyTrees(""),
			       gXSec(0),
			       gCount(0),
			       gIsMCatNLO(false),
			       gNEntries(0),
			       gSumOfWeights(0),
			       nEntries(0),
			       Count(0),
			       xsec(0),
			       nProcessedEvents(0),
			       gOptions(""),
			       gIsData2017(0)
{}



EventBuilder::~EventBuilder() {
  //delete fPUWeight;
  //delete fPUWeightUp;
  //delete fPUWeightDown;
}

void EventBuilder::Initialise(){
  year    = GetParam<Int_t>("Year");
  gIsData = GetParam<Bool_t>("IsData");
  gSelection = GetParam<Int_t>("iSelection");
  gSampleName  = GetParam<TString>("sampleName");
  gIsMCatNLO   = GetParam<Bool_t>("IsMCatNLO");
  gIsFastSim   = GetParam<Bool_t>("IsFastSim");
  gXSec        = GetParam<Float_t>("xsec");
  gCount       = GetParam<Int_t>("Count");
  gNEntries    = GetParam<Int_t>("nEntries");
  gSumOfWeights= GetParam<Double_t>("SumOfWeights");
  gOptions     = GetParam<TString>("_options");
  gIs2017 = false;
  if(gOptions.Contains("2017")) gIs2017 = true;
  gChannel = -1;
  nProcessedEvents = 0; 
  //if(gSelection == iTopSelec) gIsFastSim = true;

  selLeptons = std::vector<Lepton>();
  vetoLeptons = std::vector<Lepton>();

  gIsDoubleElec = false; gIsDoubleMuon = false; gIsSingleElec = false;
  gIsSingleMuon = false; gIsMuonEG = false;
  if(gSampleName.Contains("DoubleEG")) gIsDoubleElec = true;
  else if(gSampleName.Contains("DoubleMuon")) gIsDoubleMuon = true;
  else if(gSampleName.Contains("SingleElec")) gIsSingleElec = true;
  else if(gSampleName.Contains("SingleMuon")) gIsSingleMuon = true;
  else if(gSampleName.Contains("MuonEG"))     gIsMuonEG     = true;

/*
  fPUWeight     = new PUWeight(19468.3, Moriond17MC_PoissonOOTPU, "2016_Moriond17");
  if (!gIsData) {
    fPUWeightUp   = new PUWeight(18494.9,  Moriond17MC_PoissonOOTPU, "2016_Moriond17"); //  18494.9
    fPUWeightDown = new PUWeight(20441.7,  Moriond17MC_PoissonOOTPU, "2016_Moriond17"); //  20441.7
  }
*/

  Weight = GetParam<Float_t>("weight");

  passTrigger = 1;
  isSS = 0;
  nTrueInt = 0;

  PUSF = 1;
  PUSF_Up = 1;
  PUSF_Down = 1;
  gIsData2017 = gOptions.Contains("Data2017")? true : false;

}

void EventBuilder::InsideLoop(){
  nProcessedEvents++;
  // >>>>>>>>>>>>>> Get selected leptons:
  selLeptons = GetParam<std::vector<Lepton>>("selLeptons");
  vetoLeptons = GetParam<std::vector<Lepton>>("vetoLeptons");

  // Set channel
  if(selLeptons.size() >= 2){ // Dilepton Channels
    if     (selLeptons.at(0).isElec && selLeptons.at(1).isMuon) gChannel = iElMu;
    else if(selLeptons.at(0).isMuon && selLeptons.at(1).isElec) gChannel = iElMu;
    else if(selLeptons.at(0).isMuon && selLeptons.at(1).isMuon) gChannel = iMuon;
    else if(selLeptons.at(0).isElec && selLeptons.at(1).isElec) gChannel = iElec;
    isSS = (selLeptons[0].charge*selLeptons[1].charge) > 0;
  }
  else{
    isSS = false;
    gChannel = -1;
  }
  passTrigger  = false;

  if      (gChannel == iElMu && TrigElMu()) passTrigger = true;
  else if (gChannel == iMuon && TrigMuMu()) passTrigger = true;
  else if (gChannel == iElec && TrigElEl()) passTrigger = true;

  METfilters = PassesMETfilters();

  // >>>>>>>>> Calculate norm weight
  if(gIsMCatNLO) genWeight = Get<Float_t>("genWeight");
  else           genWeight = 1;
  NormWeight = Weight*genWeight;

  // >>>>>>>>> Calculate PU weight and variations
  if(!gIsData){
    //nTrueInt = Get<Float_t>("nTrueInt");
    //PUSF      = fPUWeight    ->GetWeight(nTrueInt);
    //PUSF_Up   = fPUWeightUp  ->GetWeight(nTrueInt);
    //PUSF_Down = fPUWeightDown->GetWeight(nTrueInt);
      PUSF = 1; PUSF_Up = 1; PUSF_Down = 1;
  } 
  else{
    PUSF      = 1;
    PUSF_Up   = 1;
    PUSF_Down = 1;
  }

  // Set Params to pass all the info...
  SetParam("PUSF",      PUSF);
  SetParam("PUSF_Up",   PUSF_Up);
  SetParam("PUSF_Down", PUSF_Down);

  SetParam("gChannel",        gChannel);
  SetParam("NormWeight",      NormWeight);
  SetParam("passTrigger",     passTrigger);
  SetParam("isSS",            isSS);
  SetParam("METfilters",      METfilters);
}



void EventBuilder::Summary(){
  cout << endl << endl << " ========================================================= " << endl;
  cout << " ====== Sample: \"" << gSampleName <<"\"" << flush;
  cout << " for selection: \"" << LabSelection[gSelection] << "\" ======\n";
  cout << " -----------> Is data?.......... "; if(gIsData)    cout << "YES\n"; else cout << "NO\n";
  cout << " -----------> Is aMCatNLO?...... "; if(gIsMCatNLO) cout << "YES\n"; else cout << "NO\n";
  cout << " -----------> Is FastSim?-...... "; if(gIsFastSim) cout << "YES\n"; else cout << "NO\n";
  cout << " >>> Number of gen events      : " << gCount           << endl;
  cout << " >>> Number of events in sample: " << gNEntries        << endl;
  cout << " >>> Number of processed events: " << nProcessedEvents << endl;
  cout << " >>> Sum of gen weights        : " << gSumOfWeights    << endl;
  cout << " >>> Cross section for norm    : " << gXSec            << endl;
  cout << " >>> Total weight  for norm    : " << Weight           << endl;
  cout << " >>> Processed events          : " << nProcessedEvents << endl;
  cout << " ========================================================= " << endl;
}








// Compute single and double lepton triggers for each year (and analysis)
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Bool_t EventBuilder::PassesDoubleElecTrigger(){
  if(gIsFastSim) return true; // no triger in FastSim samples
  Bool_t pass = false;
  gIsData = GetParam<Bool_t>("IsData");
  if(gIsData) run = Get<UInt_t>("run");

  // SAME FOR 2017... if(gIsData2017)  pass = Get<Int_t>("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  if     (year == 2016){
    pass = (Get<Int_t>("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"));
  }
  else if(year == 2017){
    pass = Get<Int_t>("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") ||
      Get<Int_t>("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  }
  else if(year == 2018){
    pass = false;
  }
  else{
    cout << "[EventBuilder] Wrong selection for checking trigger requirements!!" << endl;
    return false;
  }
  return pass;
}

Bool_t EventBuilder::PassesDoubleMuonTrigger(){
  if(gIsFastSim) return true; // no triger in FastSim samples
  Bool_t pass = false;
  if (gIsData) run     = Get<UInt_t>("run");

  if     (year == 2016){
   if ( (gIsData && run <= 280385) || (!gIsData)){
      pass = (Get<Int_t>("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL")  ||
	      Get<Int_t>("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL"));
    }
    else{
      pass = ( Get<Int_t>("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ") ||
	       Get<Int_t>("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"));
    }
  }
  else if(year == 2017){
    pass = Get<Int_t>("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
  }
  else if(year == 2018){
    pass = false;
  }
  else{
    cout << "[EventBuilder] Wrong selection for checking trigger requirements!!" << endl;
  }
  return pass;
}

Bool_t EventBuilder::PassesElMuTrigger(){
  if(gIsFastSim) return true; // no triger in FastSim samples
  Bool_t pass = false;
  if (gIsData) run     = Get<UInt_t>("run");

  if     (year == 2016){
    if ( (gIsData && run <= 280385) || (!gIsData)){
      pass = ( Get<Int_t>("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL")  ||
          Get<Int_t>("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") );
    }
    else{
      pass = ( Get<Int_t>("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")||
          Get<Int_t>("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") );
    }
  }
  else if(year == 2017){
    pass = Get<Int_t>("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") ||
      Get<Int_t>("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") ||
      Get<Int_t>("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  }
  else if(year == 2018){
    pass = false;
  }
  else{
    cout << "[EventBuilder] Wrong selection for checking trigger requirements!!" << endl;
    return false;
  }
  return pass;
}

Bool_t EventBuilder::PassesSingleElecTrigger(){
  if(gIsFastSim) return true; // no triger in FastSim samples
  Bool_t pass = false;
  if     (year == 2016){
    pass =  Get<Int_t>("HLT_Ele27_WPTight_Gsf");
  }
  else if(year == 2017){
    pass = Get<Int_t>("HLT_Ele32_WPTight_Gsf") || Get<Int_t>("HLT_Ele35_WPLoose_Gsf");
  }
  else if(year == 2018){
    pass = false;
  }
  else{
    cout << "[EventBuilder] Wrong selection for checking trigger requirements!!" << endl;
    return false;
  }
  return pass;
}

Bool_t EventBuilder::PassesSingleMuonTrigger(){
  if(gIsFastSim) return true; // no triger in FastSim samples
  Bool_t pass = false;
  if     (year == 2016){
    pass = Get<Int_t>("HLT_IsoTkMu24") || Get<Int_t>("HLT_IsoMu24") ;
  }
  else if(year == 2017){
    pass = Get<Int_t>("HLT_IsoMu24") || Get<Int_t>("HLT_IsoMu27");
  }
  else if(year == 2018){
    pass = false;
  }
  return pass;
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// In principle, you don't want to change the functions below...

Bool_t EventBuilder::TrigElEl(){
  Bool_t pass = false;
  if(gIsData){
    if     (gIsDoubleElec) pass =  PassesDoubleElecTrigger();
    else if(gIsSingleElec) pass = !PassesDoubleElecTrigger() && PassesSingleElecTrigger();
  }
  else pass = PassesDoubleElecTrigger() || PassesSingleElecTrigger();
  return pass;
}

Bool_t EventBuilder::TrigMuMu(){
  Bool_t pass = false;
  if(gIsData){
    if     (gIsDoubleMuon) pass =  PassesDoubleMuonTrigger();
    else if(gIsSingleMuon) pass = !PassesDoubleMuonTrigger() && PassesSingleMuonTrigger();
  }
  else pass = PassesDoubleMuonTrigger() || PassesSingleMuonTrigger();
  return pass;
}

Bool_t EventBuilder::TrigElMu(){
  Bool_t pass = false;
  if(gIsData){
    if(gIsMuonEG    ) pass =  PassesElMuTrigger();
    else if(gIsSingleMuon) pass = !PassesElMuTrigger() &&  PassesSingleMuonTrigger();
    else if(gIsSingleElec) pass = !PassesElMuTrigger() && !PassesSingleMuonTrigger() && PassesSingleElecTrigger();
  }
  else pass = PassesElMuTrigger() || PassesSingleMuonTrigger() || PassesSingleElecTrigger();
  return pass;
}

// ########################### MET FILTERS

Bool_t EventBuilder::PassesMETfilters(){
  // https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2#Moriond_2017
  if(gIsFastSim) return true;
  if( (Get<Int_t>("Flag_HBHENoiseFilter") &&        // MET filters for data and MC
        Get<Int_t>("Flag_HBHENoiseIsoFilter") &&
        Get<Int_t>("Flag_EcalDeadCellTriggerPrimitiveFilter") &&
        Get<Int_t>("Flag_goodVertices") ) //&&
        //Get<Int_t>("Flag_badMuonFilter") &&
       // Get<Int_t>("Flag_badChargedHadronFilter"))
      && (
        gIsFastSim || // no more MET filters for Fast Sim
        (!gIsData && Get<Int_t>("Flag_globalTightHalo2016Filter")) || // for MC
        // --> This is the right thing!: //
        ( gIsData && Get<Int_t>("Flag_globalTightHalo2016Filter") && Get<Int_t>("Flag_eeBadScFilter")) ) // for Data
        //( Get<Int_t>("Flag_eeBadScFilter")) ) // for Data
        //( gIsData && Get<Int_t>("Flag_globalTightHalo2016Filter")) ) // for Data
    ) return true;
  else return false;
}
