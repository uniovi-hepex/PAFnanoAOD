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
                               gIsMET(false),
                               gIsData(false),
                               run(-1),
                               gSelection(-1),
                               gSampleName(""),
                               gPathToHeppyTrees(""),
                               gXSec(0),
                               gCount(0),
                               gIsMCatNLO(false),
                               gNEntries(0),
                               xsec(0),
                               nProcessedEvents(0),
                               gOptions("")
{}



EventBuilder::~EventBuilder() {
  //delete fPUWeight;
  //delete fPUWeightUp;
  //delete fPUWeightDown;
}

void EventBuilder::Initialise(){
  year         = GetParam<TString>("year").Atoi();
  gIsData      = GetParam<Bool_t>("IsData");
  selection    = GetParam<TString>("selection");
  gSampleName  = GetParam<TString>("sampleName");
  gIsMCatNLO   = GetParam<Bool_t>("IsMCatNLO");
  gXSec        = GetParam<Float_t>("xsec");
  gOptions     = GetParam<TString>("_options");
  gSelection   = GetSelection(selection);
  gChannel = -1;
  nProcessedEvents = 0;
  
  gIs2017 = false; gIs2016 = false; gIs2018 = false;
  if     (year == 2017) gIs2017 = true;
  else if(year == 2018) gIs2018 = true;
  else if(year == 2016) gIs2016 = true;

  gSelection     = GetSelection(selection);

  Float_t binsEta[] = {0, 1, 2.1, 2.4};     Int_t nbinsEta = 3;
  Float_t binsPt[]  = {20, 30, 40, 50, 60}; Int_t nbinsPt  = 4;
  
  if (gSelection ==  itt) makeeffhistos = true;
  else                    makeeffhistos = false;

  if (makeeffhistos) {
    ElecTrigEffNum = CreateH2F("ElecTrigEffNum","", nbinsEta, binsEta, nbinsPt, binsPt);
    ElecTrigEffDen = CreateH2F("ElecTrigEffDen","", nbinsEta, binsEta, nbinsPt, binsPt);
    MuonTrigEffNum = CreateH2F("MuonTrigEffNum","", nbinsEta, binsEta, nbinsPt, binsPt);
    MuonTrigEffDen = CreateH2F("MuonTrigEffDen","", nbinsEta, binsEta, nbinsPt, binsPt);
    ElElTrigEffNum = CreateH2F("ElElTrigEffNum", "", nbinsEta, binsEta, nbinsPt, binsPt);
    MuMuTrigEffNum = CreateH2F("MuMuTrigEffNum", "", nbinsEta, binsEta, nbinsPt, binsPt);
    ElMuTrigEffNum = CreateH2F("ElMuTrigEffNum", "", nbinsEta, binsEta, nbinsPt, binsPt);
    ElElTrigEffDen = CreateH2F("ElElTrigEffDen", "", nbinsEta, binsEta, nbinsPt, binsPt);
    MuMuTrigEffDen = CreateH2F("MuMuTrigEffDen", "", nbinsEta, binsEta, nbinsPt, binsPt);
    ElMuTrigEffDen = CreateH2F("ElMuTrigEffDen", "", nbinsEta, binsEta, nbinsPt, binsPt);
  }
  

  selLeptons = std::vector<Lepton>();
  vetoLeptons = std::vector<Lepton>();

  gIsDoubleElec = false; gIsDoubleMuon = false; gIsSingleElec = false;
  gIsSingleMuon = false; gIsMuonEG = false; gIsMET = false; gIsEGamma = false;
  if(gSampleName.Contains("DoubleEG")) gIsDoubleElec = true;
  else if(gSampleName.Contains("DoubleMuon")) gIsDoubleMuon = true;
  else if(gSampleName.Contains("SingleElec")) gIsSingleElec = true;
  else if(gSampleName.Contains("SingleMuon")) gIsSingleMuon = true;
  else if(gSampleName.Contains("MuonEG"))     gIsMuonEG     = true;
  else if(gSampleName.Contains("HighEGJet"))  gIsSingleElec = true;
  else if(gSampleName.Contains("MET"))        gIsMET = true;
  else if(gSampleName.Contains("EGamma"))     gIsEGamma = true;
  if(gOptions.Contains("DoubleEG")) gIsDoubleElec = true;
  else if(gOptions.Contains("SingleElec")) gIsSingleElec = true;

/*
  fPUWeight     = new PUWeight(19468.3, Moriond17MC_PoissonOOTPU, "2016_Moriond17");
  if (!gIsData) {
    fPUWeightUp   = new PUWeight(18494.9,  Moriond17MC_PoissonOOTPU, "2016_Moriond17"); //  18494.9
    fPUWeightDown = new PUWeight(20441.7,  Moriond17MC_PoissonOOTPU, "2016_Moriond17"); //  20441.7
  }
*/

  Weight = GetParam<Double_t>("weight");

  passTrigger = 1;
  isSS = 0;
  nTrueInt = 0;

  PUSF = 1;
  PUSF_Up = 1;
  PUSF_Down = 1;
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

  if(gChannel == iElMu){
    Lepton muon = selLeptons.at(0).isMuon ? selLeptons.at(0) : selLeptons.at(1);
    Lepton elec = selLeptons.at(0).isMuon ? selLeptons.at(1) : selLeptons.at(0);
    float mupt = muon.Pt() < 200 ? muon.Pt() : 199;
    float elpt = elec.Pt() < 200 ? elec.Pt() : 199;
    if (makeeffhistos) {
      if(muon.Pt() > 20 && PassesSingleMuonTrigger())                              ElecTrigEffDen->Fill(elec.Eta(), elpt);
      if(muon.Pt() > 20 && PassesSingleMuonTrigger() && PassesSingleElecTrigger()) ElecTrigEffNum->Fill(elec.Eta(), elpt);
      if(elec.Pt() > 20 && PassesSingleElecTrigger())                              MuonTrigEffDen->Fill(muon.Eta(), mupt);
      if(elec.Pt() > 20 && PassesSingleElecTrigger() && PassesSingleMuonTrigger()) MuonTrigEffNum->Fill(muon.Eta(), mupt);
    }
  }

  if(selLeptons.size() >= 2){ // Dilepton Channels
    float pt = selLeptons.at(0).Pt() < 200 ? selLeptons.at(0).Pt() : 199;
    float pt2= selLeptons.at(1).Pt();
    float TWeight = 1;
    //if(!gIsData){
    //  Float_t lepSF = selLeptons.at(0).GetSF( 0)*selLeptons.at(1).GetSF( 0);
    //  Float_t PUSF  = Get<Float_t>("puWeight");
    //  TWeight       = NormWeight*lepSF*PUSF;
    //}
    Float_t MET = gIs2017? Get<Float_t>("METFixEE2017_pt") : Get<Float_t>("MET_pt");

    if(pt > 25 && pt2 > 20 && MET > 120){
      float eta = TMath::Abs(selLeptons.at(0).Eta());
      float eta2 = selLeptons.at(1).Eta();
      
      if (makeeffhistos) {
        if (gChannel == iMuon && PassesMETtrigger()) {
          MuMuTrigEffDen->Fill(eta,pt);
          if (TrigMuMu()) MuMuTrigEffNum->Fill(eta,pt);
        }
        else if (gChannel==iElec && PassesMETtrigger()) {
          ElElTrigEffDen->Fill(eta,pt);
          if (TrigElEl()) ElElTrigEffNum->Fill(eta,pt);
        }
        else if (gChannel==iElMu && PassesMETtrigger()) {
          ElMuTrigEffDen->Fill(eta,pt);
          if (TrigElMu()) ElMuTrigEffNum->Fill(eta,pt);
        }
      }
    }
  }

  if      (gChannel == iElMu && TrigElMu()) passTrigger = true;
  else if (gChannel == iMuon && TrigMuMu()) passTrigger = true;
  else if (gChannel == iElec && TrigElEl()) passTrigger = true;

  METfilters = PassesMETfilters();

  // >>>>>>>>> Calculate norm weight
  if(gIsMCatNLO) genWeight = (Double_t)Get<Float_t>("genWeight");
  else           genWeight = 1;
  NormWeight = Weight*genWeight;
  if(gIsData) NormWeight = 1;
 
  SetParam("gChannel",        gChannel);
  SetParam("NormWeight",      NormWeight);
  SetParam("passTrigger",     passTrigger);
  SetParam("isSS",            isSS);
  SetParam("METfilters",      METfilters);
}



void EventBuilder::Summary(){
}





// Compute single and double lepton triggers for each year (and analysis)
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Bool_t EventBuilder::PassesDoubleElecTrigger(){
  int era = -1; if(gIsData) era = GetRunEra(Get<Int_t>("run"), year);
  Bool_t pass = false;
  if(gIsData) run = Get<UInt_t>("run");

  if     (gIs2016){
    pass = (Get<Bool_t>("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ"));
  }
  else if(gIs2017){
    pass = Get<Bool_t>("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") ||
      Get<Bool_t>("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  }
  else if(gIs2018){
    pass = Get<Bool_t>("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") ||
      Get<Bool_t>("HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL");
  }
  return pass;
}

Bool_t EventBuilder::PassesDoubleMuonTrigger(){
  int era = -1; if(gIsData) era = GetRunEra(Get<Int_t>("run"), year);
  Bool_t pass = false;
  if (gIsData) run     = Get<UInt_t>("run");

  if     (gIs2016){
   if ( (gIsData && run <= 280385) || (!gIsData)){
      pass = (Get<Bool_t>("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL")  ||
	      Get<Bool_t>("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL"));
    }
    else{
      pass = ( Get<Bool_t>("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ") ||
	       Get<Bool_t>("HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ"));
    }
  }
  else if(gIs2017){
    if(era == runA || era == runB){
      pass = Get<Bool_t>("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ");
    }
    else if(era == runC || era == runD || era == runE || era == runF || era == runG){
      pass = Get<Bool_t>("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") ||  Get<Bool_t>("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
    }
    else{
      pass =  Get<Bool_t>("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ") || Get<Bool_t>("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") ||  Get<Bool_t>("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
    }
  }
  else if(gIs2018){
    pass = Get<Bool_t>("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8") ||
      Get<Bool_t>("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8");
  }
  return pass;
}

Bool_t EventBuilder::PassesElMuTrigger(){
  int era = -1; if(gIsData) era = GetRunEra(Get<Int_t>("run"), year);
  Bool_t pass = false;
  if (gIsData) run     = Get<UInt_t>("run");

  if     (gIs2016){
    if ( (gIsData && run <= 280385) || (!gIsData)){
      pass = ( Get<Bool_t>("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL")  ||
          Get<Bool_t>("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") );
    }
    else{
      pass = ( Get<Bool_t>("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ")||
          Get<Bool_t>("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ") );
    }
  }
  else if(gIs2017){
    if(era == runB){
      pass =
        Get<Bool_t>("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") ||
        Get<Bool_t>("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") ||
        Get<Bool_t>("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
    }
    else{
      pass =
        Get<Bool_t>("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") ||
        Get<Bool_t>("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL") ||
        Get<Bool_t>("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL") ||
        Get<Bool_t>("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") ||
        Get<Bool_t>("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") ||
        Get<Bool_t>("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
    }
  }
  else if(gIs2018){
    pass =
      Get<Bool_t>("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") ||
      Get<Bool_t>("HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL") ||
      Get<Bool_t>("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL") ||
      Get<Bool_t>("HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ") ||
      Get<Bool_t>("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL") ||
      Get<Bool_t>("HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ");
  }
  else{
    cout << "[EventBuilder] Wrong selection for checking trigger requirements!!" << endl;
    return false;
  }
  return pass;
}

Bool_t EventBuilder::PassesSingleElecTrigger(){
  int era = -1; if(gIsData) era = GetRunEra(Get<Int_t>("run"), year);
  Bool_t pass = false;
  if     (gIs2016){
    if(era == runH) pass =  Get<Bool_t>("HLT_Ele27_WPTight_Gsf");
    else pass =  Get<Bool_t>("HLT_Ele27_WPTight_Gsf") || Get<Bool_t>("HLT_Ele23_WPLoose_Gsf");
    //if(era == runH) pass =  1;
    //else pass = Get<Bool_t>("HLT_Ele23_WPLoose_Gsf");
  }
  else if(gIs2017){
    if(gIsData){
      pass = Get<Bool_t>("HLT_Ele32_WPTight_Gsf_L1DoubleEG") || Get<Bool_t>("HLT_Ele35_WPTight_Gsf");// || Get<Bool_t>("HLT_Ele38_WPTight_Gsf") || Get<Bool_t>("HLT_Ele40_WPTight_Gsf");
    }
    else{
      pass = Get<Bool_t>("HLT_Ele32_WPTight_Gsf_L1DoubleEG") || Get<Bool_t>("HLT_Ele35_WPTight_Gsf");// || Get<Bool_t>("HLT_Ele38_WPTight_Gsf") || Get<Bool_t>("HLT_Ele40_WPTight_Gsf");
    }
  }

  else if(gIs2018){
      pass = Get<Bool_t>("HLT_Ele32_WPTight_Gsf_L1DoubleEG") || Get<Bool_t>("HLT_Ele35_WPTight_Gsf") || Get<Bool_t>("HLT_Ele38_WPTight_Gsf");
  }
  return pass;
}

Bool_t EventBuilder::PassesSingleMuonTrigger(){
  int era = -1; if(gIsData) era = GetRunEra(Get<Int_t>("run"), year);
  Bool_t pass = false;
  if     (gIs2016){
    pass = Get<Bool_t>("HLT_IsoTkMu24") || Get<Bool_t>("HLT_IsoMu24") || Get<Bool_t>("HLT_IsoTkMu20") || Get<Bool_t>("HLT_IsoMu20") ;
  }
  else if(gIs2017){
    pass = Get<Bool_t>("HLT_IsoMu24") || Get<Bool_t>("HLT_IsoMu27") || Get<Bool_t>("HLT_IsoMu24_eta2p1");
  }
  else if(gIs2018){
    pass = Get<Bool_t>("HLT_IsoMu24");
  }
  return pass;
}

Bool_t EventBuilder::PassesMETtrigger(){ // THIS IS WROOOOOOOOOOOOOOONG (at least for 2017's nanoAODv4)
  Bool_t pass = false;
  int era = -1; if(gIsData) era = GetRunEra(Get<Int_t>("run"), year);
  if(gIsData){
    if(era == runB) pass = Get<Bool_t>("HLT_PFMET120_PFMHT120_IDTight");
    else pass = Get<Bool_t>("HLT_PFMET120_PFMHT120_IDTight") || Get<Bool_t>("HLT_PFMET120_PFMHT120_IDTight_PFHT60");
  }
  else pass = Get<Bool_t>("HLT_PFMET120_PFMHT120_IDTight") || Get<Bool_t>("HLT_PFMET120_PFMHT120_IDTight_PFHT60");
  return pass;
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
// In principle, you don't want to change the functions below...

Bool_t EventBuilder::TrigElEl(){
  Bool_t pass = false;
  if(gIsData){
    if     (gIsDoubleElec) pass =  PassesDoubleElecTrigger();
    else if(gIsSingleElec) pass = !PassesDoubleElecTrigger() && PassesSingleElecTrigger();
    else if(gIsMET)        pass = PassesDoubleElecTrigger() || PassesSingleElecTrigger();
    else if(gIsEGamma)     pass = PassesDoubleElecTrigger() || PassesSingleElecTrigger();
    else                   pass = false; //PassesDoubleElecTrigger() || PassesSingleElecTrigger();
  }
  else pass = PassesDoubleElecTrigger() || PassesSingleElecTrigger();
  return pass;
}

Bool_t EventBuilder::TrigMuMu(){
  Bool_t pass = false;
  if(gIsData){
    if     (gIsDoubleMuon) pass =  PassesDoubleMuonTrigger();
    else if(gIsSingleMuon) pass = !PassesDoubleMuonTrigger() && PassesSingleMuonTrigger();
    else if(gIsMET)        pass = PassesDoubleMuonTrigger() || PassesSingleMuonTrigger();
    else                   pass = false; //PassesDoubleMuonTrigger() || PassesSingleMuonTrigger();
  }
  else pass = PassesDoubleMuonTrigger() || PassesSingleMuonTrigger();
  return pass;
}

Bool_t EventBuilder::TrigElMu(){
  Bool_t pass = false;
  if(gIsData){
    if     (gIsMuonEG    ) pass =  PassesElMuTrigger();
    else if(gIsSingleMuon) pass = !PassesElMuTrigger() &&  PassesSingleMuonTrigger();
    else if(gIsSingleElec) pass = !PassesElMuTrigger() && !PassesSingleMuonTrigger() && PassesSingleElecTrigger();
    else if(gIsMET)        pass = PassesElMuTrigger() || PassesSingleMuonTrigger() || PassesSingleElecTrigger();
    else                   pass = false; //PassesElMuTrigger() || PassesSingleMuonTrigger() | PassesSingleElecTrigger();
  }
  else pass = PassesElMuTrigger() || PassesSingleMuonTrigger() | PassesSingleElecTrigger();
  return pass;
}



// ########################### MET FILTERS
Bool_t EventBuilder::PassesMETfilters() { // Updated on 2019-06-14 for both data and MC
  if (gIsData) {
    if      (gIs2016) {
      if (Get<Bool_t>("Flag_goodVertices")                       &&
          Get<Bool_t>("Flag_globalSuperTightHalo2016Filter")     &&
          Get<Bool_t>("Flag_HBHENoiseFilter")                    &&
          Get<Bool_t>("Flag_HBHENoiseIsoFilter")                 &&
          Get<Bool_t>("Flag_EcalDeadCellTriggerPrimitiveFilter") &&
          Get<Bool_t>("Flag_BadPFMuonFilter")                    &&
//           Get<Bool_t>("Flag_BadChargedCandidateFilter")          && // Currently NOT RECOMMENDED
          Get<Bool_t>("Flag_eeBadScFilter")
         ) return true;
      else return false;
    }
    else if (gIs2017) {
      if (Get<Bool_t>("Flag_goodVertices")                       &&
          Get<Bool_t>("Flag_globalSuperTightHalo2016Filter")     &&
          Get<Bool_t>("Flag_HBHENoiseFilter")                    &&
          Get<Bool_t>("Flag_HBHENoiseIsoFilter")                 &&
          Get<Bool_t>("Flag_EcalDeadCellTriggerPrimitiveFilter") &&
          Get<Bool_t>("Flag_BadPFMuonFilter")                    &&
//           Get<Bool_t>("Flag_BadChargedCandidateFilter")          && // Currently NOT RECOMMENDED
          Get<Bool_t>("Flag_eeBadScFilter")                      &&
          Get<Bool_t>("Flag_ecalBadCalibFilterV2")                 //&& // WE ARE APPLYING THIS ONE INSTEAD OF THE NEXT ONE THAT DOES NOT EXIST EXACTLY
          //Get<Bool_t>("ecalBadCalibReducedMINIAODFilter")             // WE DO NOT HAVE THIS ONE
        ) return true;
      else return false;
    }
    else if (gIs2018) {
      if (Get<Bool_t>("Flag_goodVertices")                       &&
          Get<Bool_t>("Flag_globalSuperTightHalo2016Filter")     &&
          Get<Bool_t>("Flag_HBHENoiseFilter")                    &&
          Get<Bool_t>("Flag_HBHENoiseIsoFilter")                 &&
          Get<Bool_t>("Flag_EcalDeadCellTriggerPrimitiveFilter") &&
          Get<Bool_t>("Flag_BadPFMuonFilter")                    &&
//           Get<Bool_t>("Flag_BadChargedCandidateFilter")          && // Currently NOT RECOMMENDED
          Get<Bool_t>("Flag_eeBadScFilter")                      &&
          Get<Bool_t>("Flag_ecalBadCalibFilterV2")                 //&& // WE ARE APPLYING THIS ONE INSTEAD OF THE NEXT ONE THAT DOES NOT EXIST EXACTLY
          //Get<Bool_t>("ecalBadCalibReducedMINIAODFilter")             // WE DO NOT HAVE THIS ONE
        ) return true;
      else return false;
    }
  }
  else { // ONLY FULLSIM
    if      (gIs2016) {
      if (Get<Bool_t>("Flag_goodVertices")                       &&
          Get<Bool_t>("Flag_globalSuperTightHalo2016Filter")     &&
          Get<Bool_t>("Flag_HBHENoiseFilter")                    &&
          Get<Bool_t>("Flag_HBHENoiseIsoFilter")                 &&
          Get<Bool_t>("Flag_EcalDeadCellTriggerPrimitiveFilter") &&
          Get<Bool_t>("Flag_BadPFMuonFilter")
//           Get<Bool_t>("Flag_BadChargedCandidateFilter") // Currently NOT RECOMMENDED
         ) return true;
      else return false;
    }
    else if (gIs2017) {
      if (Get<Bool_t>("Flag_goodVertices")                       &&
          Get<Bool_t>("Flag_globalSuperTightHalo2016Filter")     &&
          Get<Bool_t>("Flag_HBHENoiseFilter")                    &&
          Get<Bool_t>("Flag_HBHENoiseIsoFilter")                 &&
          Get<Bool_t>("Flag_EcalDeadCellTriggerPrimitiveFilter") &&
          Get<Bool_t>("Flag_BadPFMuonFilter")                    &&
//           Get<Bool_t>("Flag_BadChargedCandidateFilter")          && // Currently NOT RECOMMENDED
          Get<Bool_t>("Flag_ecalBadCalibFilterV2")                 //&& // WE ARE APPLYING THIS ONE INSTEAD OF THE NEXT ONE THAT DOES NOT EXIST EXACTLY
          //Get<Bool_t>("ecalBadCalibReducedMINIAODFilter")             // WE DO NOT HAVE THIS ONE
        ) return true;
      else return false;
    }
    else if (gIs2018) {
      if (Get<Bool_t>("Flag_goodVertices")                       &&
          Get<Bool_t>("Flag_globalSuperTightHalo2016Filter")     &&
          Get<Bool_t>("Flag_HBHENoiseFilter")                    &&
          Get<Bool_t>("Flag_HBHENoiseIsoFilter")                 &&
          Get<Bool_t>("Flag_EcalDeadCellTriggerPrimitiveFilter") &&
          Get<Bool_t>("Flag_BadPFMuonFilter")                    &&
//           Get<Bool_t>("Flag_BadChargedCandidateFilter")         && // Currently NOT RECOMMENDED
          Get<Bool_t>("Flag_ecalBadCalibFilterV2")                 //&& // WE ARE APPLYING THIS ONE INSTEAD OF THE NEXT ONE THAT DOES NOT EXIST EXACTLY
          //Get<Bool_t>("ecalBadCalibReducedMINIAODFilter")             // WE DO NOT HAVE THIS ONE
        ) return true;
      else return false;
    }
  }
  return false;
}
