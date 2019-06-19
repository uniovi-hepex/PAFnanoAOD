#include "TWTTbarAnalysis.h"
ClassImp(TWTTbarAnalysis);


//#####################################################################
// Core PAF methods
//---------------------------------------------------------------------
TWTTbarAnalysis::TWTTbarAnalysis() : PAFChainItemSelector() {
  fhDummy        = 0;
  passMETfilters = 0;
  passTrigger    = 0;
  isSS           = 0;
  
  for (UShort_t i = 0; i < 254; i++) TLHEWeight[i] = 0;
}



void TWTTbarAnalysis::Initialise() {
  gIsData     = GetParam<Bool_t>("IsData");
  gSampleName = GetParam<TString>("sampleName");
  gOptions    = GetParam<TString>("_options");
  if (gOptions.Contains("Semi")) {
    cout << "> Running the semileptonic ttbar sample" << endl;
  }
  gPUWeight   = gOptions.Contains("PUweight")? true : false;
  gIsTTbar    = false;
  gIsLHE      = false;

  if (gSampleName.Contains("TTbar") || gSampleName.Contains("TTJets")) gIsTTbar = true;
  if (gSampleName == "TTbar_Powheg")   gIsLHE = true;
  
  if (gSampleName.Contains("_")) {
    TObjArray *tx = gSampleName.Tokenize("_");
    if (((TObjString*) tx->Last())->GetString().IsDigit()) gIsLHE = true;
  }
  
  fhDummy = CreateH1F("fhDummy", "fhDummy", 1, 0, 1);
  
  fMiniTree = CreateTree("fMiniTree", "MiniTree");
  SetTWTTbarVariables();
  
  DressLeptons          = std::vector<Lepton>();
  selLeptons            = std::vector<Lepton>();
  selJets               = std::vector<Jet>();
  selJetsJecUp          = std::vector<Jet>();
  selJetsJecDown        = std::vector<Jet>();
  selJetsJER            = std::vector<Jet>();
  genJets               = std::vector<Jet>();
  vetoJets              = std::vector<Jet>();
  
  DressJets             = std::vector<Jet>();
  DressLooseCentralJets = std::vector<Jet>();
  DressLooseFwdJets     = std::vector<Jet>();
}



void TWTTbarAnalysis::InsideLoop() {
  ResetTWTTbarVariables();

  TEvent = Get<ULong64_t>("event");

  // Vectors with the objects
  DressLeptons      = GetParam<vector<Lepton>>("genLeptons");
  selLeptons        = GetParam<vector<Lepton>>("selLeptons");
  selJets           = GetParam<vector<Jet>>("selJets");
  selJetsJecUp      = GetParam<vector<Jet>>("selJetsJecUp");
  selJetsJecDown    = GetParam<vector<Jet>>("selJetsJecDown");
  selJetsJER        = GetParam<vector<Jet>>("selJetsJER");
  TNBJets           = (UShort_t)GetParam<Int_t>("nSelBJets");
  vetoJets          = GetParam<vector<Jet>>("vetoJets");
  genJets           = GetParam<vector<Jet>>("genJets");
  
  // Weights and SFs
  NormWeight        = GetParam<Float_t>("NormWeight");
  TrigSF            = GetParam<Float_t>("TriggerSF");
  TrigSFerr         = GetParam<Float_t>("TriggerSFerr");

  if (!gIsData && gPUWeight) {
    PUSF         = Get<Float_t>("puWeight");
    PUSF_Up      = Get<Float_t>("puWeightUp");
    PUSF_Down    = Get<Float_t>("puWeightDown");
  }
  else {PUSF = 1; PUSF_Up = 1; PUSF_Down = 1;}

  BtagSF            = GetParam<Float_t>("BtagSF");
  BtagSFBtagUp      = GetParam<Float_t>("BtagSFBtagUp");
  BtagSFBtagDown    = GetParam<Float_t>("BtagSFBtagDown");
  BtagSFMistagUp    = GetParam<Float_t>("BtagSFMistagUp");
  BtagSFMistagDown  = GetParam<Float_t>("BtagSFMistagDown");

  // Event variables
  passMETfilters    = GetParam<Bool_t>("METfilters");
  passTrigger       = GetParam<Bool_t>("passTrigger");
  isSS              = GetParam<Bool_t>("isSS");
  year              = (UShort_t)GetParam<Int_t>("year");
  
  // Leptons and Jets
  GetLeptonVariables();
  GetGenLepVariables();
  
  if (gOptions.Contains("Semi")) {
    if (gIsTTbar && DressNLeps > 1 ) return;
  } else {
    if (gIsTTbar && DressNLeps < 2 ) return; // Dilepton selection for ttbar!
  }
  
  GetJetVariables();
  GetGenJetVariables();
  GetMETandGenMET();
  
  TWeight_normal = NormWeight;
  fhDummy->Fill(0.5);
  
  // Particle level selection
  if ((DressNLeps >= 2) && (DressNJets == 2) && (DressNBJets == 2) && ((DressLeptons.at(0).p + DressLeptons.at(1).p).M() > 20) &&
      (DressLeptons.at(0).p.Pt() > 25) && !TDressIsSS && (DressNLooseCentral == 2)) {
    
    CalculateDressTWTTbarVariables();
    DoesItReallyPassDress();
  }
  
  // Detector level selection
  if ((TNSelLeps >= 2) && passTrigger && passMETfilters && ((selLeptons.at(0).p + selLeptons.at(1).p).M() > 20) &&
      (selLeptons.at(0).p.Pt() > 25) && (TChannel == iElMu || TChannel == iElec || TChannel == iMuon) && (!isSS)) {
    
    CalculateSFAndWeights();
    CalculateTWTTbarVariables();
    DoesItReallyPassReco();
  }
  
  // Filling choice
  if (TPassPart || TPassDress || TPassReco || TPassRecoJESUp || TPassRecoJESDown || TPassRecoJERUp) { // If needed, filling.
    SetMinimaAndMaxima();
    fMiniTree->Fill();
  }
}



void TWTTbarAnalysis::Summary(){}




//#####################################################################
// Functions
//---------------------------------------------------------------------
void TWTTbarAnalysis::SetTWTTbarVariables() {
  // Detector level variables
  fMiniTree->Branch("TEvent",                &TEvent,                "TEvent/l");
  fMiniTree->Branch("TChannel",              &TChannel,              "TChannel/S");
  fMiniTree->Branch("TIsSS",                 &TIsSS,                 "TIsSS/O");
  fMiniTree->Branch("TWeight",               &TWeight,               "TWeight/F");
  fMiniTree->Branch("TWeight_ElecEffUp",     &TWeight_ElecEffUp,     "TWeight_ElecEffUp/F");
  fMiniTree->Branch("TWeight_ElecEffDown",   &TWeight_ElecEffDown,   "TWeight_ElecEffDown/F");
  fMiniTree->Branch("TWeight_MuonEffUp",     &TWeight_MuonEffUp,     "TWeight_MuonEffUp/F");
  fMiniTree->Branch("TWeight_MuonEffDown",   &TWeight_MuonEffDown,   "TWeight_MuonEffDown/F");
  fMiniTree->Branch("TWeight_TrigUp",        &TWeight_TrigUp,        "TWeight_TrigUp/F");
  fMiniTree->Branch("TWeight_TrigDown",      &TWeight_TrigDown,      "TWeight_TrigDown/F");
  fMiniTree->Branch("TWeight_PUUp",          &TWeight_PUUp,          "TWeight_PUUp/F");
  fMiniTree->Branch("TWeight_PUDown",        &TWeight_PUDown,        "TWeight_PUDown/F");
  fMiniTree->Branch("TWeight_BtagUp",        &TWeight_BtagUp,        "TWeight_BtagUp/F");
  fMiniTree->Branch("TWeight_BtagDown",      &TWeight_BtagDown,      "TWeight_BtagDown/F");
  fMiniTree->Branch("TWeight_MistagUp",      &TWeight_MistagUp,      "TWeight_MistagUp/F");
  fMiniTree->Branch("TWeight_MistagDown",    &TWeight_MistagDown,    "TWeight_MistagDown/F");
  fMiniTree->Branch("TWeight_normal",        &TWeight_normal,        "TWeight_normal/F");
  fMiniTree->Branch("TLHEWeight",            TLHEWeight,             "TLHEWeight[254]/F");
  
  
  fMiniTree->Branch("TPassReco",             &TPassReco,             "TPassReco/O");
  fMiniTree->Branch("TPassRecoJESUp",        &TPassRecoJESUp,        "TPassRecoJESUp/O");
  fMiniTree->Branch("TPassRecoJESDown",      &TPassRecoJESDown,      "TPassRecoJESDown/O");
  fMiniTree->Branch("TPassRecoJERUp",        &TPassRecoJERUp,        "TPassRecoJERUp/O");
  fMiniTree->Branch("TNJets"       ,         &TNJets,                "TNJets/s");
  fMiniTree->Branch("TNJetsJESUp"  ,         &TNJetsJESUp,           "TNJetsJESUp/s");
  fMiniTree->Branch("TNJetsJESDown",         &TNJetsJESDown,         "TNJetsJESDown/s");
  fMiniTree->Branch("TNJetsJERUp",           &TNJetsJERUp,           "TNJetsJERUp/s");
  fMiniTree->Branch("TNBJets"       ,        &TNBJets,               "TNBJets/s");
  fMiniTree->Branch("TNBJetsJESUp"  ,        &TNBJetsJESUp,          "TNBJetsJESUp/s");
  fMiniTree->Branch("TNBJetsJESDown",        &TNBJetsJESDown,        "TNBJetsJESDown/s");
  fMiniTree->Branch("TNBJetsJERUp",          &TNBJetsJERUp,          "TNBJetsJERUp/s");
  fMiniTree->Branch("TNLooseCentral"        ,&NLooseCentral      ,   "TNLooseCentral/s");
  fMiniTree->Branch("TNLooseCentralJESUp"   ,&NLooseCentralJESUp  ,  "TNLooseCentralJESUp/s");
  fMiniTree->Branch("TNLooseCentralJESDown" ,&NLooseCentralJESDown,  "TNLooseCentralJESDown/s");
  fMiniTree->Branch("TNLooseCentralJERUp"   ,&NLooseCentralJERUp  ,  "TNLooseCentralJERUp/s");
  fMiniTree->Branch("TNBLooseCentral"       ,&NBLooseCentral      ,  "TNBLooseCentral/s");
  fMiniTree->Branch("TNBLooseCentralJESUp"  ,&NBLooseCentralJESUp  , "TNBLooseCentralJESUp/s");
  fMiniTree->Branch("TNBLooseCentralJESDown",&NBLooseCentralJESDown, "TNBLooseCentralJESDown/s");
  fMiniTree->Branch("TNBLooseCentralJERUp"  ,&NBLooseCentralJERUp  , "TNBLooseCentralJERUp/s");
  fMiniTree->Branch("TNLooseFwd"            ,&NLooseFwd            , "TNLooseFwd/s");
  fMiniTree->Branch("TNLooseFwdJESUp"       ,&NLooseFwdJESUp        ,"TNLooseFwdJESUp/s");
  fMiniTree->Branch("TNLooseFwdJESDown"     ,&NLooseFwdJESDown      ,"TNLooseFwdJESDown/s");
  fMiniTree->Branch("TNLooseFwdJERUp"       ,&NLooseFwdJERUp        ,"TNLooseFwdJERUp/s");
  
  fMiniTree->Branch("TLep1_Pt",              &TLep1_Pt,              "TLep1_Pt/F");
  fMiniTree->Branch("TLep1_E",               &TLep1_E,               "TLep1_E/F");
  fMiniTree->Branch("TLep1_Eta",             &TLep1_Eta,             "TLep1_Eta/F");
  fMiniTree->Branch("TLep2_Pt",              &TLep2_Pt,              "TLep2_Pt/F");
  fMiniTree->Branch("TLep2_E",               &TLep2_E,               "TLep2_E/F");
  fMiniTree->Branch("TLep2_Eta",             &TLep2_Eta,             "TLep2_Eta/F");
  fMiniTree->Branch("TJet1_Pt",              &TJet1_Pt,              "TJet1_Pt/F");
  fMiniTree->Branch("TJet1_E",               &TJet1_E,               "TJet1_E/F");
  fMiniTree->Branch("TJet1_Eta",             &TJet1_Eta,             "TJet1_Eta/F");
  fMiniTree->Branch("TJet2_Pt",              &TJet2_Pt,              "TJet2_Pt/F");
  fMiniTree->Branch("TJet2_E",               &TJet2_E,               "TJet2_E/F");
  fMiniTree->Branch("TJet2_Eta",             &TJet2_Eta,             "TJet2_Eta/F");
  
  fMiniTree->Branch("TLep1Lep2_Pt",          &TLep1Lep2_Pt,          "TLep1Lep2_Pt/F");
  fMiniTree->Branch("TLep1Lep2_M",           &TLep1Lep2_M,           "TLep1Lep2_M/F");
  fMiniTree->Branch("TLep1Lep2_DPhi",        &TLep1Lep2_DPhi,        "TLep1Lep2_DPhi/F");
  fMiniTree->Branch("TLep1Jet1_M" ,          &TLep1Jet1_M,           "TLep1Jet1_M/F");
  fMiniTree->Branch("TLep1Jet1_DPhi",        &TLep1Jet1_DPhi,        "TLep1Jet1_DPhi/F");
  fMiniTree->Branch("TLep1Jet2_M",           &TLep1Jet2_M,           "TLep1Jet2_M/F");
  fMiniTree->Branch("TLep1Jet2_DPhi",        &TLep1Jet2_DPhi,        "TLep1Jet2_DPhi/F");
  fMiniTree->Branch("TLep2Jet1_M",           &TLep2Jet1_M,           "TLep2Jet1_M/F");
  fMiniTree->Branch("TLep2Jet1_DPhi",        &TLep2Jet1_DPhi,        "TLep2Jet1_DPhi/F");
  fMiniTree->Branch("TLep2Jet2_M",           &TLep2Jet2_M,           "TLep2Jet2_M/F");
  fMiniTree->Branch("TLep2Jet2_DPhi",        &TLep2Jet2_DPhi,        "TLep2Jet2_DPhi/F");
  
  fMiniTree->Branch("TSys_Pt",               &Sys_Pt,                "TSys_Pt/F");
  fMiniTree->Branch("TSys_E",                &Sys_E,                 "TSys_E/F");
  fMiniTree->Branch("TSys_Eta",              &Sys_Eta,               "TSys_Eta/F");
  fMiniTree->Branch("TSys_M",                &Sys_M,                 "TSys_M/F");
  fMiniTree->Branch("TSys_Pz",               &Sys_Pz,                "TSys_Pz/F");
  fMiniTree->Branch("TMiniMax",              &TMiniMax,              "TMiniMax/F");
  fMiniTree->Branch("THT",                   &THT,                   "THT/F");
  fMiniTree->Branch("TMET",                  &TMET,                  "TMET/F");
  
  
  // JESUp
  fMiniTree->Branch("TLep1_PtJESUp",         &TLep1_PtJESUp,         "TLep1_PtJESUp/F");
  fMiniTree->Branch("TLep1_EJESUp",          &TLep1_EJESUp,          "TLep1_EJESUp/F");
  fMiniTree->Branch("TLep1_EtaJESUp",        &TLep1_EtaJESUp,        "TLep1_EtaJESUp/F");
  fMiniTree->Branch("TLep2_PtJESUp",         &TLep2_PtJESUp,         "TLep2_PtJESUp/F");
  fMiniTree->Branch("TLep2_EJESUp",          &TLep2_EJESUp,          "TLep2_EJESUp/F");
  fMiniTree->Branch("TLep2_EtaJESUp",        &TLep2_EtaJESUp,        "TLep2_EtaJESUp/F");
  fMiniTree->Branch("TJet1_PtJESUp",         &TJet1_PtJESUp,         "TJet1_PtJESUp/F");
  fMiniTree->Branch("TJet1_EJESUp",          &TJet1_EJESUp,          "TJet1_EJESUp/F");
  fMiniTree->Branch("TJet1_EtaJESUp",        &TJet1_EtaJESUp,        "TJet1_EtaJESUp/F");
  fMiniTree->Branch("TJet2_PtJESUp",         &TJet2_PtJESUp,         "TJet2_PtJESUp/F");
  fMiniTree->Branch("TJet2_EJESUp",          &TJet2_EJESUp,          "TJet2_EJESUp/F");
  fMiniTree->Branch("TJet2_EtaJESUp",        &TJet2_EtaJESUp,        "TJet2_EtaJESUp/F");
  
  fMiniTree->Branch("TLep1Lep2_PtJESUp",     &TLep1Lep2_PtJESUp,     "TLep1Lep2_PtJESUp/F");
  fMiniTree->Branch("TLep1Lep2_MJESUp",      &TLep1Lep2_MJESUp,      "TLep1Lep2_MJESUp/F");
  fMiniTree->Branch("TLep1Lep2_DPhiJESUp",   &TLep1Lep2_DPhiJESUp,   "TLep1Lep2_DPhiJESUp/F");
  fMiniTree->Branch("TLep1Jet1_MJESUp" ,     &TLep1Jet1_MJESUp,      "TLep1Jet1_MJESUp/F");
  fMiniTree->Branch("TLep1Jet1_DPhiJESUp",   &TLep1Jet1_DPhiJESUp,   "TLep1Jet1_DPhiJESUp/F");
  fMiniTree->Branch("TLep1Jet2_MJESUp",      &TLep1Jet2_MJESUp,      "TLep1Jet2_MJESUp/F");
  fMiniTree->Branch("TLep1Jet2_DPhiJESUp",   &TLep1Jet2_DPhiJESUp,   "TLep1Jet2_DPhiJESUp/F");
  fMiniTree->Branch("TLep2Jet1_MJESUp",      &TLep2Jet1_MJESUp,      "TLep2Jet1_MJESUp/F");
  fMiniTree->Branch("TLep2Jet1_DPhiJESUp",   &TLep2Jet1_DPhiJESUp,   "TLep2Jet1_DPhiJESUp/F");
  fMiniTree->Branch("TLep2Jet2_MJESUp",      &TLep2Jet2_MJESUp,      "TLep2Jet2_MJESUp/F");
  fMiniTree->Branch("TLep2Jet2_DPhiJESUp",   &TLep2Jet2_DPhiJESUp,   "TLep2Jet2_DPhiJESUp/F");
  
  fMiniTree->Branch("TSys_PtJESUp",          &Sys_PtJESUp,           "TSys_PtJESUp/F");
  fMiniTree->Branch("TSys_EJESUp",           &Sys_EJESUp,            "TSys_EJESUp/F");
  fMiniTree->Branch("TSys_EtaJESUp",         &Sys_EtaJESUp,          "TSys_EtaJESUp/F");
  fMiniTree->Branch("TSys_MJESUp",           &Sys_MJESUp,            "TSys_MJESUp/F");
  fMiniTree->Branch("TSys_PzJESUp",          &Sys_PzJESUp,           "TSys_PzJESUp/F");
  fMiniTree->Branch("TMiniMaxJESUp",         &TMiniMaxJESUp,         "TMiniMaxJESUp/F");
  fMiniTree->Branch("THTJESUp",              &THTJESUp,              "THTJESUp/F");
  fMiniTree->Branch("TMETJESUp",             &TMETJESUp,             "TMETJESUp/F");
  
  
  // JESDown
  fMiniTree->Branch("TLep1_PtJESDown",       &TLep1_PtJESDown,       "TLep1_PtJESDown/F");
  fMiniTree->Branch("TLep1_EJESDown",        &TLep1_EJESDown,        "TLep1_EJESDown/F");
  fMiniTree->Branch("TLep1_EtaJESDown",      &TLep1_EtaJESDown,      "TLep1_EtaJESDown/F");
  fMiniTree->Branch("TLep2_PtJESDown",       &TLep2_PtJESDown,       "TLep2_PtJESDown/F");
  fMiniTree->Branch("TLep2_EJESDown",        &TLep2_EJESDown,        "TLep2_EJESDown/F");
  fMiniTree->Branch("TLep2_EtaJESDown",      &TLep2_EtaJESDown,      "TLep2_EtaJESDown/F");
  fMiniTree->Branch("TJet1_PtJESDown",       &TJet1_PtJESDown,       "TJet1_PtJESDown/F");
  fMiniTree->Branch("TJet1_EJESDown",        &TJet1_EJESDown,        "TJet1_EJESDown/F");
  fMiniTree->Branch("TJet1_EtaJESDown",      &TJet1_EtaJESDown,      "TJet1_EtaJESDown/F");
  fMiniTree->Branch("TJet2_PtJESDown",       &TJet2_PtJESDown,       "TJet2_PtJESDown/F");
  fMiniTree->Branch("TJet2_EJESDown",        &TJet2_EJESDown,        "TJet2_EJESDown/F");
  fMiniTree->Branch("TJet2_EtaJESDown",      &TJet2_EtaJESDown,      "TJet2_EtaJESDown/F");
  
  fMiniTree->Branch("TLep1Lep2_PtJESDown",   &TLep1Lep2_PtJESDown,   "TLep1Lep2_PtJESDown/F");
  fMiniTree->Branch("TLep1Lep2_MJESDown",    &TLep1Lep2_MJESDown,    "TLep1Lep2_MJESDown/F");
  fMiniTree->Branch("TLep1Lep2_DPhiJESDown", &TLep1Lep2_DPhiJESDown, "TLep1Lep2_DPhiJESDown/F");
  fMiniTree->Branch("TLep1Jet1_MJESDown" ,   &TLep1Jet1_MJESDown,    "TLep1Jet1_MJESDown/F");
  fMiniTree->Branch("TLep1Jet1_DPhiJESDown", &TLep1Jet1_DPhiJESDown, "TLep1Jet1_DPhiJESDown/F");
  fMiniTree->Branch("TLep1Jet2_MJESDown",    &TLep1Jet2_MJESDown,    "TLep1Jet2_MJESDown/F");
  fMiniTree->Branch("TLep1Jet2_DPhiJESDown", &TLep1Jet2_DPhiJESDown, "TLep1Jet2_DPhiJESDown/F");
  fMiniTree->Branch("TLep2Jet1_MJESDown",    &TLep2Jet1_MJESDown,    "TLep2Jet1_MJESDown/F");
  fMiniTree->Branch("TLep2Jet1_DPhiJESDown", &TLep2Jet1_DPhiJESDown, "TLep2Jet1_DPhiJESDown/F");
  fMiniTree->Branch("TLep2Jet2_MJESDown",    &TLep2Jet2_MJESDown,    "TLep2Jet2_MJESDown/F");
  fMiniTree->Branch("TLep2Jet2_DPhiJESDown", &TLep2Jet2_DPhiJESDown, "TLep2Jet2_DPhiJESDown/F");
  
  fMiniTree->Branch("TSys_PtJESDown",        &Sys_PtJESDown,         "TSys_PtJESDown/F");
  fMiniTree->Branch("TSys_EJESDown",         &Sys_EJESDown,          "TSys_EJESDown/F");
  fMiniTree->Branch("TSys_EtaJESDown",       &Sys_EtaJESDown,        "TSys_EtaJESDown/F");
  fMiniTree->Branch("TSys_MJESDown",         &Sys_MJESDown,          "TSys_MJESDown/F");
  fMiniTree->Branch("TSys_PzJESDown",        &Sys_PzJESDown,         "TSys_PzJESDown/F");
  fMiniTree->Branch("TMiniMaxJESDown",       &TMiniMaxJESDown,       "TMiniMaxJESDown/F");
  fMiniTree->Branch("THTJESDown",            &THTJESDown,            "THTJESDown/F");
  fMiniTree->Branch("TMETJESDown",           &TMETJESDown,           "TMETJESDown/F");
  
  
  // JERUp
  fMiniTree->Branch("TLep1_PtJERUp",         &TLep1_PtJERUp,         "TLep1_PtJERUp/F");
  fMiniTree->Branch("TLep1_EJERUp",          &TLep1_EJERUp,          "TLep1_EJERUp/F");
  fMiniTree->Branch("TLep1_EtaJERUp",        &TLep1_EtaJERUp,        "TLep1_EtaJERUp/F");
  fMiniTree->Branch("TLep2_PtJERUp",         &TLep2_PtJERUp,         "TLep2_PtJERUp/F");
  fMiniTree->Branch("TLep2_EJERUp",          &TLep2_EJERUp,          "TLep2_EJERUp/F");
  fMiniTree->Branch("TLep2_EtaJERUp",        &TLep2_EtaJERUp,        "TLep2_EtaJERUp/F");
  fMiniTree->Branch("TJet1_PtJERUp",         &TJet1_PtJERUp,         "TJet1_PtJERUp/F");
  fMiniTree->Branch("TJet1_EJERUp",          &TJet1_EJERUp,          "TJet1_EJERUp/F");
  fMiniTree->Branch("TJet1_EtaJERUp",        &TJet1_EtaJERUp,        "TJet1_EtaJERUp/F");
  fMiniTree->Branch("TJet2_PtJERUp",         &TJet2_PtJERUp,         "TJet2_PtJERUp/F");
  fMiniTree->Branch("TJet2_EJERUp",          &TJet2_EJERUp,          "TJet2_EJERUp/F");
  fMiniTree->Branch("TJet2_EtaJERUp",        &TJet2_EtaJERUp,        "TJet2_EtaJERUp/F");
  
  fMiniTree->Branch("TLep1Lep2_PtJERUp",     &TLep1Lep2_PtJERUp,     "TLep1Lep2_PtJERUp/F");
  fMiniTree->Branch("TLep1Lep2_MJERUp",      &TLep1Lep2_MJERUp,      "TLep1Lep2_MJERUp/F");
  fMiniTree->Branch("TLep1Lep2_DPhiJERUp",   &TLep1Lep2_DPhiJERUp,   "TLep1Lep2_DPhiJERUp/F");
  fMiniTree->Branch("TLep1Jet1_MJERUp" ,     &TLep1Jet1_MJERUp,      "TLep1Jet1_MJERUp/F");
  fMiniTree->Branch("TLep1Jet1_DPhiJERUp",   &TLep1Jet1_DPhiJERUp,   "TLep1Jet1_DPhiJERUp/F");
  fMiniTree->Branch("TLep1Jet2_MJERUp",      &TLep1Jet2_MJERUp,      "TLep1Jet2_MJERUp/F");
  fMiniTree->Branch("TLep1Jet2_DPhiJERUp",   &TLep1Jet2_DPhiJERUp,   "TLep1Jet2_DPhiJERUp/F");
  fMiniTree->Branch("TLep2Jet1_MJERUp",      &TLep2Jet1_MJERUp,      "TLep2Jet1_MJERUp/F");
  fMiniTree->Branch("TLep2Jet1_DPhiJERUp",   &TLep2Jet1_DPhiJERUp,   "TLep2Jet1_DPhiJERUp/F");
  fMiniTree->Branch("TLep2Jet2_MJERUp",      &TLep2Jet2_MJERUp,      "TLep2Jet2_MJERUp/F");
  fMiniTree->Branch("TLep2Jet2_DPhiJERUp",   &TLep2Jet2_DPhiJERUp,   "TLep2Jet2_DPhiJERUp/F");
  
  fMiniTree->Branch("TSys_PtJERUp",          &Sys_PtJERUp,           "TSys_PtJERUp/F");
  fMiniTree->Branch("TSys_EJERUp",           &Sys_EJERUp,            "TSys_EJERUp/F");
  fMiniTree->Branch("TSys_EtaJERUp",         &Sys_EtaJERUp,          "TSys_EtaJERUp/F");
  fMiniTree->Branch("TSys_MJERUp",           &Sys_MJERUp,            "TSys_MJERUp/F");
  fMiniTree->Branch("TSys_PzJERUp",          &Sys_PzJERUp,           "TSys_PzJERUp/F");
  fMiniTree->Branch("TMiniMaxJERUp",         &TMiniMaxJERUp,         "TMiniMaxJERUp/F");
  fMiniTree->Branch("THTJERUp",              &THTJERUp,              "THTJERUp/F");
  fMiniTree->Branch("TMETJERUp",             &TMETJERUp,             "TMETJERUp/F");
  
  
  // VAINAS DE GENERACION
  // PARTICLE LEVEL: DRESS
  // PARTON LEVEL: PART
  // SHARED: GEN
  fMiniTree->Branch("TGenChannel",           &GenChannel,            "GenChannel/S");
  
  // Particle level variables
  fMiniTree->Branch("TPassDress",            &TPassDress,            "TPassDress/O");
  fMiniTree->Branch("TDressIsSS",            &TDressIsSS,            "TDressIsSS/O");
  fMiniTree->Branch("TDressNJets",           &DressNJets,            "TDressNJets/s");
  fMiniTree->Branch("TDressNBJets",          &DressNBJets,           "TDressNBJets/s");
  fMiniTree->Branch("TDressNLooseCentral",   &DressNLooseCentral,    "TDressNLooseCentral/s");
  fMiniTree->Branch("TDressNBLooseCentral",  &DressNBLooseCentral,   "TDressNBLooseCentral/s");
  fMiniTree->Branch("TDressNLooseFwd",       &DressNLooseFwd,        "TDressNLooseFwd/s");
  
  fMiniTree->Branch("TDressLep1_Pt",         &TDressLep1_Pt,         "TDressLep1_Pt/F");
  fMiniTree->Branch("TDressLep1_E",          &TDressLep1_E,          "TDressLep1_E/F");
  fMiniTree->Branch("TDressLep1_Eta",        &TDressLep1_Eta,        "TDressLep1_Eta/F");
  fMiniTree->Branch("TDressLep2_Pt",         &TDressLep2_Pt,         "TDressLep2_Pt/F");
  fMiniTree->Branch("TDressLep2_E",          &TDressLep2_E,          "TDressLep2_E/F");
  fMiniTree->Branch("TDressLep2_Eta",        &TDressLep2_Eta,        "TDressLep2_Eta/F");
  fMiniTree->Branch("TDressJet1_Pt",         &TDressJet1_Pt,         "TDressJet1_Pt/F");
  fMiniTree->Branch("TDressJet1_E",          &TDressJet1_E,          "TDressJet1_E/F");
  fMiniTree->Branch("TDressJet1_Eta",        &TDressJet1_Eta,        "TDressJet1_Eta/F");
  fMiniTree->Branch("TDressJet2_Pt",         &TDressJet2_Pt,         "TDressJet2_Pt/F");
  fMiniTree->Branch("TDressJet2_E",          &TDressJet2_E,          "TDressJet2_E/F");
  fMiniTree->Branch("TDressJet2_Eta",        &TDressJet2_Eta,        "TDressJet2_Eta/F");
  
  fMiniTree->Branch("TDressLep1Lep2_Pt",     &TDressLep1Lep2_Pt,     "TDressLep1Lep2_Pt/F");
  fMiniTree->Branch("TDressLep1Lep2_M",      &TDressLep1Lep2_M,      "TDressLep1Lep2_M/F");
  fMiniTree->Branch("TDressLep1Lep2_DPhi",   &TDressLep1Lep2_DPhi,   "TDressLep1Lep2_DPhi/F");
  fMiniTree->Branch("TDressLep1Jet1_M" ,     &TDressLep1Jet1_M,      "TDressLep1Jet1_M/F");
  fMiniTree->Branch("TDressLep1Jet1_DPhi",   &TDressLep1Jet1_DPhi,   "TDressLep1Jet1_DPhi/F");
  fMiniTree->Branch("TDressLep1Jet2_M",      &TDressLep1Jet2_M,      "TDressLep1Jet2_M/F");
  fMiniTree->Branch("TDressLep1Jet2_DPhi",   &TDressLep1Jet2_DPhi,   "TDressLep1Jet2_DPhi/F");
  fMiniTree->Branch("TDressLep2Jet1_M",      &TDressLep2Jet1_M,      "TDressLep2Jet1_M/F");
  fMiniTree->Branch("TDressLep2Jet1_DPhi",   &TDressLep2Jet1_DPhi,   "TDressLep2Jet1_DPhi/F");
  fMiniTree->Branch("TDressLep2Jet2_M",      &TDressLep2Jet2_M,      "TDressLep2Jet2_M/F");
  fMiniTree->Branch("TDressLep2Jet2_DPhi",   &TDressLep2Jet2_DPhi,   "TDressLep2Jet2_DPhi/F");
  
  fMiniTree->Branch("TDressSys_Pt",          &DressSys_Pt,           "TDressSys_Pt/F");
  fMiniTree->Branch("TDressSys_E",           &DressSys_E,            "TDressSys_E/F");
  fMiniTree->Branch("TDressSys_Eta",         &DressSys_Eta,          "TDressSys_Eta/F");
  fMiniTree->Branch("TDressSys_M",           &DressSys_M,            "TDressSys_M/F");
  fMiniTree->Branch("TDressSys_Pz",          &DressSys_Pz,           "TDressSys_Pz/F");
  fMiniTree->Branch("TDressMiniMax",         &TDressMiniMax,         "TDressMiniMax/F");
  fMiniTree->Branch("TDressHT",              &TDressHT,              "TDressHT/F");
  fMiniTree->Branch("TDressMET",             &TDressMET,             "TDressMET/F");
  
  // Parton level variables
  fMiniTree->Branch("TPassPart",             &TPassPart,             "TPassPart/O");
  fMiniTree->Branch("TPartMET",              &TPartMET,              "TPartMET/F");
}


void TWTTbarAnalysis::GetLeptonVariables() {
  TNSelLeps = selLeptons.size();

  if (TNSelLeps < 2) TChannel = -1;
  else if (selLeptons.at(0).isMuon && selLeptons.at(1).isElec) TChannel = iElMu;
  else if (selLeptons.at(0).isElec && selLeptons.at(1).isMuon) TChannel = iElMu;
  else if (selLeptons.at(0).isMuon && selLeptons.at(1).isMuon) TChannel = iMuon;
  else if (selLeptons.at(0).isElec && selLeptons.at(1).isElec) TChannel = iElec;
  
  TIsSS = isSS;
  
  if (TNSelLeps >= 1) {
    TLep1_Pt          = selLeptons.at(0).Pt();
    TLep1_E           = selLeptons.at(0).E();
    TLep1_Phi         = selLeptons.at(0).Phi();
    TLep1_Eta         = selLeptons.at(0).Eta();
    TLep1_PtJESUp     = selLeptons.at(0).Pt();
    TLep1_EJESUp      = selLeptons.at(0).E();
    TLep1_PhiJESUp    = selLeptons.at(0).Phi();
    TLep1_EtaJESUp    = selLeptons.at(0).Eta();
    TLep1_PtJESDown   = selLeptons.at(0).Pt();
    TLep1_EJESDown    = selLeptons.at(0).E();
    TLep1_PhiJESDown  = selLeptons.at(0).Phi();
    TLep1_EtaJESDown  = selLeptons.at(0).Eta();
    TLep1_PtJERUp     = selLeptons.at(0).Pt();
    TLep1_EJERUp      = selLeptons.at(0).E();
    TLep1_PhiJERUp    = selLeptons.at(0).Phi();
    TLep1_EtaJERUp    = selLeptons.at(0).Eta();
    if (TNSelLeps >= 2) {
      TLep2_Pt          = selLeptons.at(1).Pt();
      TLep2_E           = selLeptons.at(1).E();
      TLep2_Phi         = selLeptons.at(1).Phi();
      TLep2_Eta         = selLeptons.at(1).Eta();
      TLep2_PtJESUp     = selLeptons.at(1).Pt();
      TLep2_EJESUp      = selLeptons.at(1).E();
      TLep2_PhiJESUp    = selLeptons.at(1).Phi();
      TLep2_EtaJESUp    = selLeptons.at(1).Eta();
      TLep2_PtJESDown   = selLeptons.at(1).Pt();
      TLep2_EJESDown    = selLeptons.at(1).E();
      TLep2_PhiJESDown  = selLeptons.at(1).Phi();
      TLep2_EtaJESDown  = selLeptons.at(1).Eta();
      TLep2_PtJERUp     = selLeptons.at(1).Pt();
      TLep2_EJERUp      = selLeptons.at(1).E();
      TLep2_PhiJERUp    = selLeptons.at(1).Phi();
      TLep2_EtaJERUp    = selLeptons.at(1).Eta();
    }
  }
}


void TWTTbarAnalysis::GetGenLepVariables() {
  if (gIsData) return;
  DressNLeps = DressLeptons.size();
  
  if (DressNLeps >= 1) {
    TDressLep1_Pt   = DressLeptons.at(0).Pt();
    TDressLep1_E    = DressLeptons.at(0).E();
    TDressLep1_Phi  = DressLeptons.at(0).Phi();
    TDressLep1_Eta  = DressLeptons.at(0).Eta();
    if (DressNLeps >= 2) {
      TDressLep2_Pt   = DressLeptons.at(1).Pt();
      TDressLep2_E    = DressLeptons.at(1).E();
      TDressLep2_Phi  = DressLeptons.at(1).Phi();
      TDressLep2_Eta  = DressLeptons.at(1).Eta();
      
      if (DressLeptons.at(0).isElec && DressLeptons.at(1).isMuon) GenChannel = iElMu;
      if (DressLeptons.at(0).isMuon && DressLeptons.at(1).isElec) GenChannel = iElMu;
      if (DressLeptons.at(0).isMuon && DressLeptons.at(1).isMuon) GenChannel = iMuon;
      if (DressLeptons.at(0).isElec && DressLeptons.at(1).isElec) GenChannel = iElec;
      TDressIsSS = (DressLeptons.at(0).charge * DressLeptons.at(1).charge) > 0;
    }
  }
}


void TWTTbarAnalysis::GetJetVariables() {
  TNJets        = selJets.size();
  TNJetsJESUp   = selJetsJecUp.size();
  TNJetsJESDown = selJetsJecDown.size();
  TNJetsJERUp   = selJetsJER.size();
  
  for (UShort_t i = 0; i < TNJets; i++) THT += selJets.at(i).Pt();
  
  if (TNJets > 0) {
    TJet1_Pt  = selJets.at(0).Pt();
    TJet1_E   = selJets.at(0).E();
    TJet1_Phi = selJets.at(0).Phi();
    TJet1_Eta = selJets.at(0).Eta();
    if (TNJets > 1) {
      TJet2_Pt  = selJets.at(1).Pt();
      TJet2_E   = selJets.at(1).E();
      TJet2_Eta = selJets.at(1).Eta();
    }
  }
  
  for (UShort_t j = 0; j < vetoJets.size(); ++j) {
    if (TMath::Abs(vetoJets.at(j).p.Eta()) < 2.4) {
      LooseCentralJets.push_back(vetoJets.at(j));
      if (vetoJets.at(j).isBtag) NBLooseCentral++;
    }
    else LooseFwdJets.push_back(vetoJets.at(j));
  }
  
  NLooseCentral = LooseCentralJets.size();
  NLooseFwd     = LooseFwdJets.size();
  
  
  // For systematics
  if (gIsData) return;
  
  for (auto& jet : selJetsJecUp) {
    THTJESUp += jet.pTJESUp;
    if (jet.isBtag) TNBJetsJESUp++;
  }
  for (auto& jet : selJetsJecDown) {
    THTJESDown += jet.pTJESDown;
    if (jet.isBtag) TNBJetsJESDown++;
  }
  for (auto& jet : selJetsJER) {
    if (jet.isBtag) TNBJetsJERUp++;
  }
  
  if (TNJetsJESUp >= 1) {
    TJet1_PtJESUp  = selJetsJecUp.at(0).Pt();
    TJet1_EJESUp   = selJetsJecUp.at(0).E();
    TJet1_PhiJESUp = selJetsJecUp.at(0).Phi();
    TJet1_EtaJESUp = selJetsJecUp.at(0).Eta();
  }
  if (TNJetsJESDown >= 1) {
    TJet1_PtJESDown  = selJetsJecDown.at(0).Pt();
    TJet1_EJESDown   = selJetsJecDown.at(0).E();
    TJet1_PhiJESDown = selJetsJecDown.at(0).Phi();
    TJet1_EtaJESDown = selJetsJecDown.at(0).Eta();
  }
  if (TNJetsJERUp >= 1) {
    TJet1_PtJERUp  = selJetsJER.at(0).Pt();
    TJet1_EJERUp   = selJetsJER.at(0).E();
    TJet1_PhiJERUp = selJetsJER.at(0).Phi();
    TJet1_EtaJERUp = selJetsJER.at(0).Eta();
  }
  
  
  for (UShort_t j = 0; j < vetoJets.size(); ++j) {
    if (vetoJets.at(j).pTJESUp > 20.) {
      if (TMath::Abs(vetoJets.at(j).p.Eta()) < 2.4) {
        LooseCentralJetsJESUp.push_back(vetoJets.at(j));
        if (vetoJets.at(j).isBtag) NBLooseCentralJESUp++;
      }
      else LooseFwdJetsJESUp.push_back(vetoJets.at(j));
    }

    if (vetoJets.at(j).pTJESDown > 20.) {
      if (TMath::Abs(vetoJets.at(j).p.Eta()) < 2.4) {
        LooseCentralJetsJESDown.push_back(vetoJets.at(j));
        if (vetoJets.at(j).isBtag) NBLooseCentralJESDown++;
      }
      else LooseFwdJetsJESDown.push_back(vetoJets.at(j));
    }

    if (vetoJets.at(j).pTJERUp > 20.) {
      if (TMath::Abs(vetoJets.at(j).p.Eta()) < 2.4) {
        LooseCentralJetsJERUp.push_back(vetoJets.at(j));
        if (vetoJets.at(j).isBtag) NBLooseCentralJERUp++;
      }
      else LooseFwdJetsJERUp.push_back(vetoJets.at(j));
    }
  }
  NLooseCentralJESUp   = LooseCentralJetsJESUp.size();
  NLooseFwdJESUp       = LooseFwdJetsJESUp.size();
  NLooseCentralJESDown = LooseCentralJetsJESDown.size();
  NLooseFwdJESDown     = LooseFwdJetsJESDown.size();
  NLooseCentralJERUp   = LooseCentralJetsJERUp.size();
  NLooseFwdJERUp       = LooseFwdJetsJERUp.size();
}


void TWTTbarAnalysis::GetGenJetVariables() {  // TERMINAR DE REHACER DESDE JET SELECTOR
  if (gIsData) return;
  DressNJets = genJets.size();
  
  for (UShort_t i = 0; i < (UShort_t)DressNJets; i++) {
    if (TMath::Abs(genJets.at(i).p.Eta()) < 2.4 && Cleaning(genJets.at(i), DressLeptons, 0.4)) {
      DressLooseCentralJets.push_back(genJets.at(i));
      if (genJets.at(i).isBtag) DressNBLooseCentral++;
      if (genJets.at(i).p.Pt() > 30) {
        DressJets.push_back(genJets.at(i));
        TDressHT += genJets.at(i).Pt();
        if (TMath::Abs(genJets.at(i).flavmc) == 5) DressNBJets++;
      }
    }
    else if (TMath::Abs(genJets.at(i).p.Eta()) < 4.7) {
      DressLooseFwdJets.push_back(genJets.at(i));
    }
  }
  
  DressNJets          = DressJets.size();
  DressNLooseCentral  = DressLooseCentralJets.size();
  DressNLooseFwd      = DressLooseFwdJets.size();
  
  if (DressNJets >= 1) {
    TDressJet1_Pt   = DressJets.at(0).Pt();
    TDressJet1_E    = DressJets.at(0).E();
    TDressJet1_Phi  = DressJets.at(0).Phi();
    TDressJet1_Eta  = DressJets.at(0).Eta();
    if (DressNJets >= 2) {
      TDressJet1_Pt   = DressJets.at(1).Pt();
      TDressJet1_E    = DressJets.at(1).E();
      TDressJet1_Eta  = DressJets.at(1).Eta();
    }
  }
}


void TWTTbarAnalysis::GetMETandGenMET() {
  if      (year == 2017) {
    TMET        = Get<Float_t>("METFixEE2017_pt");
    TMET_Phi    = Get<Float_t>("METFixEE2017_phi");
  }
  else if (year == 2018) { // CAMBIAR PA QUE SEA LISTO Y DETECTE EL NOM CUANDO LO HAYA
    TMET        = Get<Float_t>("MET_pt_nom");
    TMET_Phi    = Get<Float_t>("MET_phi_nom");
  }
  
  if (gIsData) return;
//   TMETJESUp   = Get<Float_t>("met_jecUp_pt"  );
//   TMETJESDown = Get<Float_t>("met_jecDown_pt");
//   TMET_PhiJESUp   = Get<Float_t>("met_jecUp_phi"  );
//   TMET_PhiJESDown = Get<Float_t>("met_jecDown_phi");
//   Float_t  diff_MET_JER_phi = GetParam<Float_t>("diff_MET_JER_phi");
//   Float_t  diff_MET_JER_pt  = GetParam<Float_t>("diff_MET_JER_pt");
// 
//   TLorentzVector diff_MET_JER; diff_MET_JER.SetPtEtaPhiM(diff_MET_JER_pt, 0.,diff_MET_JER_phi, 0.);
//   TLorentzVector vMET; vMET.SetPtEtaPhiM(TMET, 0., TMET_Phi, 0);
// 
//   TMET_PhiJERUp     = (vMET + diff_MET_JER).Phi();
//   TMETJERUp         = (vMET + diff_MET_JER).Pt();

  TPartMET      = Get<Float_t>("GenMET_pt");
  TPartMET_Phi  = Get<Float_t>("GenMET_phi");
  TDressMET     = Get<Float_t>("MET_fiducialGenPt");
  TDressMET_Phi = Get<Float_t>("MET_fiducialGenPhi");
  if (gIsLHE)  for(UShort_t i = 0; i < Get<Int_t>("nLHEweight"); i++) TLHEWeight[i] = Get<Float_t>("LHEweight_wgt", i);
}


void TWTTbarAnalysis::ResetTWTTbarVariables() {
  ElecSF    = 1;  MuonSF   = 1; ElecSFUp  = 1;  ElecSFDo = 1;  MuonSFUp = 1;  MuonSFDo = 1; lepSF = 1;
  
  TPassReco              = false;
  TPassRecoJESUp         = false;
  TPassRecoJESDown       = false;
  TPassRecoJERUp         = false;
  TPassDress             = false;
  TPassPart              = false;
  TDressIsSS             = false;
  GenChannel             = -1;
  
  TNBJetsJESUp = 0; TNBJetsJESDown = 0; TNBJetsJERUp = 0;
  NBLooseCentral = 0; NBLooseCentralJESUp = 0; NBLooseCentralJESDown = 0; NBLooseCentralJERUp = 0;
  LooseCentralJets.clear(); LooseFwdJets.clear();
  LooseCentralJetsJESUp.clear(); LooseFwdJetsJESUp.clear();
  LooseCentralJetsJESDown.clear(); LooseFwdJetsJESDown.clear();
  LooseCentralJetsJERUp.clear(); LooseFwdJetsJERUp.clear();
  
  DressNBJets = 0; DressNBLooseCentral = 0;
  DressJets.clear();
  DressLooseCentralJets.clear();
  DressLooseFwdJets.clear();
  
  TLep1_Pt               = -99;
  TLep1_E                = -99;
  TLep1_Phi              = -99;
  TLep1_Eta              = -99;
  TLep2_Pt               = -99;
  TLep2_E                = -99;
  TLep2_Phi              = -99;
  TLep2_Eta              = -99;
  TJet1_Pt               = -99;
  TJet1_E                = -99;
  TJet1_Eta              = -99;
  TJet2_Pt               = -99;
  TJet2_E                = -99;
  TJet2_Eta              = -99;
  TLep1Lep2_Pt           = -99;
  TLep1Lep2_M            = -99;
  TLep1Lep2_DPhi         = -99;
  TLep1Jet1_M            = -99;
  TLep1Jet1_DPhi         = -99;
  TLep1Jet2_M            = -99;
  TLep1Jet2_DPhi         = -99;
  TLep2Jet1_M            = -99;
  TLep2Jet1_DPhi         = -99;
  TLep2Jet2_M            = -99;
  TLep2Jet2_DPhi         = -99;
  Sys_Pt                 = -99;
  Sys_E                  = -99;
  Sys_Eta                = -99;
  Sys_M                  = -99;
  Sys_Pz                 = -99;
  TMiniMax               = -99;
  THT                    = -99;
  TMET                   = -99;
  TMET_Phi               = -99;
  
  TLep1_PtJESUp          = -99;
  TLep1_EJESUp           = -99;
  TLep1_PhiJESUp         = -99;
  TLep1_EtaJESUp         = -99;
  TLep2_PtJESUp          = -99;
  TLep2_EJESUp           = -99;
  TLep2_PhiJESUp         = -99;
  TLep2_EtaJESUp         = -99;
  TJet1_PtJESUp          = -99;
  TJet1_EJESUp           = -99;
  TJet1_EtaJESUp         = -99;
  TJet2_PtJESUp          = -99;
  TJet2_EJESUp           = -99;
  TJet2_EtaJESUp         = -99;
  TLep1Lep2_PtJESUp      = -99;
  TLep1Lep2_MJESUp       = -99;
  TLep1Lep2_DPhiJESUp    = -99;
  TLep1Jet1_MJESUp       = -99;
  TLep1Jet1_DPhiJESUp    = -99;
  TLep1Jet2_MJESUp       = -99;
  TLep1Jet2_DPhiJESUp    = -99;
  TLep2Jet1_MJESUp       = -99;
  TLep2Jet1_DPhiJESUp    = -99;
  TLep2Jet2_MJESUp       = -99;
  TLep2Jet2_DPhiJESUp    = -99;
  Sys_PtJESUp            = -99;
  Sys_EJESUp             = -99;
  Sys_EtaJESUp           = -99;
  Sys_MJESUp             = -99;
  Sys_PzJESUp            = -99;
  TMiniMaxJESUp          = -99;
  THTJESUp               = -99;
  TMETJESUp              = -99;
  TMET_PhiJESUp          = -99;
  
  TLep1_PtJESDown        = -99;
  TLep1_EJESDown         = -99;
  TLep1_PhiJESDown       = -99;
  TLep1_EtaJESDown       = -99;
  TLep2_PtJESDown        = -99;
  TLep2_EJESDown         = -99;
  TLep2_PhiJESDown       = -99;
  TLep2_EtaJESDown       = -99;
  TJet1_PtJESDown        = -99;
  TJet1_EJESDown         = -99;
  TJet1_EtaJESDown       = -99;
  TJet2_PtJESDown        = -99;
  TJet2_EJESDown         = -99;
  TJet2_EtaJESDown       = -99;
  TLep1Lep2_PtJESDown    = -99;
  TLep1Lep2_MJESDown     = -99;
  TLep1Lep2_DPhiJESDown  = -99;
  TLep1Jet1_MJESDown     = -99;
  TLep1Jet1_DPhiJESDown  = -99;
  TLep1Jet2_MJESDown     = -99;
  TLep1Jet2_DPhiJESDown  = -99;
  TLep2Jet1_MJESDown     = -99;
  TLep2Jet1_DPhiJESDown  = -99;
  TLep2Jet2_MJESDown     = -99;
  TLep2Jet2_DPhiJESDown  = -99;
  Sys_PtJESDown          = -99;
  Sys_EJESDown           = -99;
  Sys_EtaJESDown         = -99;
  Sys_MJESDown           = -99;
  Sys_PzJESDown          = -99;
  TMiniMaxJESDown        = -99;
  THTJESDown             = -99;
  TMETJESDown            = -99;
  TMET_PhiJESDown        = -99;
  
  TLep1_PtJERUp          = -99;
  TLep1_EJERUp           = -99;
  TLep1_PhiJERUp         = -99;
  TLep1_EtaJERUp         = -99;
  TLep2_PtJERUp          = -99;
  TLep2_EJERUp           = -99;
  TLep2_PhiJERUp         = -99;
  TLep2_EtaJERUp         = -99;
  TJet1_PtJERUp          = -99;
  TJet1_EJERUp           = -99;
  TJet1_EtaJERUp         = -99;
  TJet2_PtJERUp          = -99;
  TJet2_EJERUp           = -99;
  TJet2_EtaJERUp         = -99;
  TLep1Lep2_PtJERUp      = -99;
  TLep1Lep2_MJERUp       = -99;
  TLep1Lep2_DPhiJERUp    = -99;
  TLep1Jet1_MJERUp       = -99;
  TLep1Jet1_DPhiJERUp    = -99;
  TLep1Jet2_MJERUp       = -99;
  TLep1Jet2_DPhiJERUp    = -99;
  TLep2Jet1_MJERUp       = -99;
  TLep2Jet1_DPhiJERUp    = -99;
  TLep2Jet2_MJERUp       = -99;
  TLep2Jet2_DPhiJERUp    = -99;
  Sys_PtJERUp            = -99;
  Sys_EJERUp             = -99;
  Sys_EtaJERUp           = -99;
  Sys_MJERUp             = -99;
  Sys_PzJERUp            = -99;
  TMiniMaxJERUp          = -99;
  THTJERUp               = -99;
  TMETJERUp              = -99;
  TMET_PhiJERUp          = -99;
  
  TDressLep1_Pt          = -99;
  TDressLep1_E           = -99;
  TDressLep1_Phi         = -99;
  TDressLep1_Eta         = -99;
  TDressLep2_Pt          = -99;
  TDressLep2_E           = -99;
  TDressLep2_Phi         = -99;
  TDressLep2_Eta         = -99;
  TDressJet1_Pt          = -99;
  TDressJet1_E           = -99;
  TDressJet1_Eta         = -99;
  TDressJet2_Pt          = -99;
  TDressJet2_E           = -99;
  TDressJet2_Eta         = -99;
  TDressLep1Lep2_Pt      = -99;
  TDressLep1Lep2_M       = -99;
  TDressLep1Lep2_DPhi    = -99;
  TDressLep1Jet1_M       = -99;
  TDressLep1Jet1_DPhi    = -99;
  TDressLep1Jet2_M       = -99;
  TDressLep1Jet2_DPhi    = -99;
  TDressLep2Jet1_M       = -99;
  TDressLep2Jet1_DPhi    = -99;
  TDressLep2Jet2_M       = -99;
  TDressLep2Jet2_DPhi    = -99;
  DressSys_Pt            = -99;
  DressSys_E             = -99;
  DressSys_Eta           = -99;
  DressSys_M             = -99;
  DressSys_Pz            = -99;
  TDressMiniMax          = -99;
  TDressHT               = -99;
  TDressMET              = -99;
  TDressMET_Phi          = -99;
  
  
  TPartMET               = -99;
  TPartMET_Phi           = -99;
}


void TWTTbarAnalysis::CalculateTWTTbarVariables() {  // VERY IMPORTAAAAAAAANT        Now ONLY for 2j2b
  TLorentzVector tmpsys;
  if (TNJets == 2 && TNBJets == 2) {
    TLep1Lep2_Pt   = (selLeptons.at(0).p + selLeptons.at(1).p).Pt();
    TLep1Lep2_M    = (selLeptons.at(0).p + selLeptons.at(1).p).M();
    TLep1Lep2_DPhi = selLeptons.at(0).p.DeltaPhi(selLeptons.at(1).p);
    TLep1Jet1_M    = (selLeptons.at(0).p + selJets.at(0).p).M();
    TLep1Jet1_DPhi = selLeptons.at(0).p.DeltaPhi(selJets.at(0).p);
    TLep1Jet2_M    = (selLeptons.at(0).p + selJets.at(1).p).M();
    TLep1Jet2_DPhi = selLeptons.at(0).p.DeltaPhi(selJets.at(1).p);
    TLep2Jet1_M    = (selLeptons.at(1).p + selJets.at(0).p).M();
    TLep2Jet1_DPhi = selLeptons.at(1).p.DeltaPhi(selJets.at(0).p);
    TLep2Jet2_M    = (selLeptons.at(1).p + selJets.at(1).p).M();
    TLep2Jet2_DPhi = selLeptons.at(1).p.DeltaPhi(selJets.at(1).p);
    tmpsys         = getSysVector();
    Sys_Pt         = tmpsys.Pt();
    Sys_E          = tmpsys.E();
    Sys_Eta        = tmpsys.Eta();
    Sys_M          = tmpsys.M();
    Sys_Pz         = tmpsys.Pz();
    TMiniMax       = getMiniMax(TLep1Jet1_M, TLep1Jet2_M, TLep2Jet1_M, TLep2Jet2_M);
  }
  
  // For systematic calculations
  if (!gIsData) {
    if (TNJetsJESUp == 2 && TNBJetsJESUp == 2) {
      TLep1Lep2_PtJESUp   = (selLeptons.at(0).p + selLeptons.at(1).p).Pt();
      TLep1Lep2_MJESUp    = (selLeptons.at(0).p + selLeptons.at(1).p).M();
      TLep1Lep2_DPhiJESUp = selLeptons.at(0).p.DeltaPhi(selLeptons.at(1).p);
      TLep1Jet1_MJESUp    = (selLeptons.at(0).p + selJetsJecUp.at(0).p).M();
      TLep1Jet1_DPhiJESUp = selLeptons.at(0).p.DeltaPhi(selJetsJecUp.at(0).p);
      TLep1Jet2_MJESUp    = (selLeptons.at(0).p + selJetsJecUp.at(1).p).M();
      TLep1Jet2_DPhiJESUp = selLeptons.at(0).p.DeltaPhi(selJetsJecUp.at(1).p);
      TLep2Jet1_MJESUp    = (selLeptons.at(1).p + selJetsJecUp.at(0).p).M();
      TLep2Jet1_DPhiJESUp = selLeptons.at(1).p.DeltaPhi(selJetsJecUp.at(0).p);
      TLep2Jet2_MJESUp    = (selLeptons.at(1).p + selJetsJecUp.at(1).p).M();
      TLep2Jet2_DPhiJESUp = selLeptons.at(1).p.DeltaPhi(selJetsJecUp.at(1).p);
      tmpsys              = getSysVector("JESUp");
      Sys_PtJESUp         = tmpsys.Pt();
      Sys_EJESUp          = tmpsys.E();
      Sys_EtaJESUp        = tmpsys.Eta();
      Sys_MJESUp          = tmpsys.M();
      Sys_PzJESUp         = tmpsys.Pz();
      TMiniMaxJESUp       = getMiniMax(TLep1Jet1_M, TLep1Jet2_M, TLep2Jet1_M, TLep2Jet2_M);
    }
    
    if (TNJetsJESDown == 2 && TNBJetsJESDown == 2) {
      TLep1Lep2_PtJESDown   = (selLeptons.at(0).p + selLeptons.at(1).p).Pt();
      TLep1Lep2_MJESDown    = (selLeptons.at(0).p + selLeptons.at(1).p).M();
      TLep1Lep2_DPhiJESDown = selLeptons.at(0).p.DeltaPhi(selLeptons.at(1).p);
      TLep1Jet1_MJESDown    = (selLeptons.at(0).p + selJetsJecDown.at(0).p).M();
      TLep1Jet1_DPhiJESDown = selLeptons.at(0).p.DeltaPhi(selJetsJecDown.at(0).p);
      TLep1Jet2_MJESDown    = (selLeptons.at(0).p + selJetsJecDown.at(1).p).M();
      TLep1Jet2_DPhiJESDown = selLeptons.at(0).p.DeltaPhi(selJetsJecDown.at(1).p);
      TLep2Jet1_MJESDown    = (selLeptons.at(1).p + selJetsJecDown.at(0).p).M();
      TLep2Jet1_DPhiJESDown = selLeptons.at(1).p.DeltaPhi(selJetsJecDown.at(0).p);
      TLep2Jet2_MJESDown    = (selLeptons.at(1).p + selJetsJecDown.at(1).p).M();
      TLep2Jet2_DPhiJESDown = selLeptons.at(1).p.DeltaPhi(selJetsJecDown.at(1).p);
      tmpsys                = getSysVector("JESDown");
      Sys_PtJESDown         = tmpsys.Pt();
      Sys_EJESDown          = tmpsys.E();
      Sys_EtaJESDown        = tmpsys.Eta();
      Sys_MJESDown          = tmpsys.M();
      Sys_PzJESDown         = tmpsys.Pz();
      TMiniMaxJESDown       = getMiniMax(TLep1Jet1_M, TLep1Jet2_M, TLep2Jet1_M, TLep2Jet2_M);
    }
    
    if (TNJetsJERUp == 2 && TNBJetsJERUp == 2) {
      TLep1Lep2_PtJERUp   = (selLeptons.at(0).p + selLeptons.at(1).p).Pt();
      TLep1Lep2_MJERUp    = (selLeptons.at(0).p + selLeptons.at(1).p).M();
      TLep1Lep2_DPhiJERUp = selLeptons.at(0).p.DeltaPhi(selLeptons.at(1).p);
      TLep1Jet1_MJERUp    = (selLeptons.at(0).p + selJetsJER.at(0).p).M();
      TLep1Jet1_DPhiJERUp = selLeptons.at(0).p.DeltaPhi(selJetsJER.at(0).p);
      TLep1Jet2_MJERUp    = (selLeptons.at(0).p + selJetsJER.at(1).p).M();
      TLep1Jet2_DPhiJERUp = selLeptons.at(0).p.DeltaPhi(selJetsJER.at(1).p);
      TLep2Jet1_MJERUp    = (selLeptons.at(1).p + selJetsJER.at(0).p).M();
      TLep2Jet1_DPhiJERUp = selLeptons.at(1).p.DeltaPhi(selJetsJER.at(0).p);
      TLep2Jet2_MJERUp    = (selLeptons.at(1).p + selJetsJER.at(1).p).M();
      TLep2Jet2_DPhiJERUp = selLeptons.at(1).p.DeltaPhi(selJetsJER.at(1).p);
      tmpsys              = getSysVector("JERUp");
      Sys_PtJERUp         = tmpsys.Pt();
      Sys_EJERUp          = tmpsys.E();
      Sys_EtaJERUp        = tmpsys.Eta();
      Sys_MJERUp          = tmpsys.M();
      Sys_PzJERUp         = tmpsys.Pz();
      TMiniMaxJERUp       = getMiniMax(TLep1Jet1_M, TLep1Jet2_M, TLep2Jet1_M, TLep2Jet2_M);
    }
  }
}


void TWTTbarAnalysis::CalculateDressTWTTbarVariables() {  // VERY IMPORTAAAAAAAANT        Now ONLY for 2j2b
  TLorentzVector tmpsys;
  if (DressNJets == 2 && DressNBJets == 2) {
    TDressLep1Lep2_Pt   = (DressLeptons.at(0).p + DressLeptons.at(1).p).Pt();
    TDressLep1Lep2_M    = (DressLeptons.at(0).p + DressLeptons.at(1).p).M();
    TDressLep1Lep2_DPhi = DressLeptons.at(0).p.DeltaPhi(DressLeptons.at(1).p);
    TDressLep1Jet1_M    = (DressLeptons.at(0).p + DressJets.at(0).p).M();
    TDressLep1Jet1_DPhi = DressLeptons.at(0).p.DeltaPhi(DressJets.at(0).p);
    TDressLep1Jet2_M    = (DressLeptons.at(0).p + DressJets.at(1).p).M();
    TDressLep1Jet2_DPhi = DressLeptons.at(0).p.DeltaPhi(DressJets.at(1).p);
    TDressLep2Jet1_M    = (DressLeptons.at(1).p + DressJets.at(0).p).M();
    TDressLep2Jet1_DPhi = DressLeptons.at(1).p.DeltaPhi(DressJets.at(0).p);
    TDressLep2Jet2_M    = (DressLeptons.at(1).p + DressJets.at(1).p).M();
    TDressLep2Jet2_DPhi = DressLeptons.at(1).p.DeltaPhi(DressJets.at(1).p);
    tmpsys              = getSysVector("Dress");
    DressSys_Pt         = tmpsys.Pt();
    DressSys_E          = tmpsys.E();
    DressSys_Eta        = tmpsys.Eta();
    DressSys_M          = tmpsys.M();
    DressSys_Pz         = tmpsys.Pz();
    TDressMiniMax       = getMiniMax(TDressLep1Jet1_M, TDressLep1Jet2_M, TDressLep2Jet1_M, TDressLep2Jet2_M);
  }
}


Float_t TWTTbarAnalysis::getMiniMax(Float_t ml1j1, Float_t ml1j2, Float_t ml2j1, Float_t ml2j2) {
  return min(max(ml1j1, ml2j2), max(ml1j2, ml2j1));
}


TLorentzVector TWTTbarAnalysis::getSysVector(const TString& sys) {
  vector<TLorentzVector> col;
  TLorentzVector met;
  
  if (sys == "Norm" || sys == "") {
    met.SetPtEtaPhiE(TMET, 0, TMET_Phi, TMET);
    col.push_back(met);
    col.push_back(selLeptons.at(0).p);
    col.push_back(selLeptons.at(1).p);
    col.push_back(selJets.at(0).p);
    col.push_back(selJets.at(1).p);
  }
  else if (sys == "JESUp") {
    met.SetPtEtaPhiE(TMETJESUp, 0, TMET_PhiJESUp, TMETJESUp);
    col.push_back(met);
    col.push_back(selLeptons.at(0).p);
    col.push_back(selLeptons.at(1).p);
    col.push_back(selJetsJecUp.at(0).p);
    col.push_back(selJetsJecUp.at(1).p);
  }
  else if (sys == "JESDown") {
    met.SetPtEtaPhiE(TMETJESDown, 0, TMET_PhiJESDown, TMETJESDown);
    col.push_back(met);
    col.push_back(selLeptons.at(0).p);
    col.push_back(selLeptons.at(1).p);
    col.push_back(selJetsJecDown.at(0).p);
    col.push_back(selJetsJecDown.at(1).p);
  }
  else if (sys == "JER") {
    met.SetPtEtaPhiE(TMETJERUp, 0, TMET_PhiJERUp, TMETJERUp);
    col.push_back(met);
    col.push_back(selLeptons.at(0).p);
    col.push_back(selLeptons.at(1).p);
    col.push_back(selJetsJER.at(0).p);
    col.push_back(selJetsJER.at(1).p);
  }
  else if (sys == "Dress") {
    met.SetPtEtaPhiE(TDressMET, 0, TDressMET_Phi, TDressMET);
    col.push_back(met);
    col.push_back(DressLeptons.at(0).p);
    col.push_back(DressLeptons.at(1).p);
    col.push_back(DressJets.at(0).p);
    col.push_back(DressJets.at(1).p);
  }
  return GetColVector(col);
}


void TWTTbarAnalysis::DoesItReallyPassDress() { // FINAL (not all) requirements for particle level selection
  if (GenChannel == iMuon || GenChannel == iElec) {
    if ((TDressMET > 20) && (abs(TDressLep1Lep2_M - Zm) > 15)) TPassDress = true;
  }
  else TPassDress = true;
}


void TWTTbarAnalysis::DoesItReallyPassReco() { // FINAL (not all) requirements for detector level selection
  if ((TNJets == 2) && (TNBJets == 2) && (NLooseCentral == 2)) {
    if (TChannel == iMuon || TChannel == iElec) {
      if ((TMET > 20) && (abs(TLep1Lep2_M - Zm) > 15)) TPassReco = true;
    }
    else TPassReco = true;
  }
  if ((TNJetsJESUp == 2) && (TNBJetsJESUp == 2) && (NLooseCentralJESUp == 2)) {
    if (TChannel == iMuon || TChannel == iElec) {
      if ((TMETJESUp > 20) && (abs(TLep1Lep2_MJESUp - Zm) > 15)) TPassRecoJESUp = true;
    }
    else TPassRecoJESUp = true;
  }
  if ((TNJetsJESDown == 2) && (TNBJetsJESDown == 2) && (NLooseCentralJESDown == 2)) {
    if (TChannel == iMuon || TChannel == iElec) {
      if ((TMETJESDown > 20) && (abs(TLep1Lep2_MJESDown - Zm) > 15)) TPassRecoJESDown = true;
    }
    else TPassRecoJESDown = true;
  }
  if ((TNJetsJERUp == 2) && (TNBJetsJERUp == 2) && (NLooseCentralJERUp == 2)) {
    if (TChannel == iMuon || TChannel == iElec) {
      if ((TMETJERUp > 20) && (abs(TLep1Lep2_MJERUp - Zm) > 15)) TPassRecoJERUp = true;
    }
    else TPassRecoJERUp = true;
  }
}


void TWTTbarAnalysis::CalculateSFAndWeights() {
  if (gIsData) {
    TWeight = 1;
    TWeight_ElecEffUp = 1; TWeight_ElecEffDown = 1;
    TWeight_MuonEffUp = 1; TWeight_MuonEffDown = 1;
    TWeight_PUUp      = 1; TWeight_PUDown      = 1;
    TWeight_BtagUp    = 1; TWeight_BtagDown    = 1;
    TWeight_MistagUp  = 1; TWeight_MistagDown  = 1;
  }
  else {
    // SF calculations
    lepSF = selLeptons.at(0).GetSF( 0) * selLeptons.at(1).GetSF( 0);
    if (TChannel == iElec) {
      ElecSF   = selLeptons.at(0).GetSF( 0)*selLeptons.at(1).GetSF( 0);
      ElecSFUp = selLeptons.at(0).GetSF( 1)*selLeptons.at(1).GetSF( 1);
      ElecSFDo = selLeptons.at(0).GetSF(-1)*selLeptons.at(1).GetSF(-1);
      MuonSFUp = 1; MuonSFDo = 1; MuonSF = 1;
    }
    else if (TChannel == iMuon) {
      MuonSFUp = selLeptons.at(0).GetSF( 1)*selLeptons.at(1).GetSF( 1);
      MuonSFDo = selLeptons.at(0).GetSF(-1)*selLeptons.at(1).GetSF(-1);
      ElecSFUp = 1; ElecSFDo = 1; ElecSF = 1;
    }
    else {
      if(selLeptons.at(0).isMuon) {
        MuonSF   *= selLeptons.at(0).GetSF( 0);
        MuonSFUp *= selLeptons.at(0).GetSF( 1);
        MuonSFDo *= selLeptons.at(0).GetSF(-1);
      }
      else {
        ElecSF   *= selLeptons.at(0).GetSF( 0);
        ElecSFUp *= selLeptons.at(0).GetSF( 1);
        ElecSFDo *= selLeptons.at(0).GetSF(-1);
      }
      if (selLeptons.at(1).isMuon) {
        MuonSF   *= selLeptons.at(1).GetSF( 0);
        MuonSFUp *= selLeptons.at(1).GetSF( 1);
        MuonSFDo *= selLeptons.at(1).GetSF(-1);
      }
      else {
        ElecSF   *= selLeptons.at(1).GetSF( 0);
        ElecSFUp *= selLeptons.at(1).GetSF( 1);
        ElecSFDo *= selLeptons.at(1).GetSF(-1);
      }
    }
    // Weight calculations
    TWeight               = NormWeight * ElecSF   * MuonSF * TrigSF * PUSF * BtagSF;
    TWeight_ElecEffUp     = NormWeight * ElecSFUp * MuonSF * TrigSF * PUSF * BtagSF;
    TWeight_ElecEffDown   = NormWeight * ElecSFDo * MuonSF * TrigSF * PUSF * BtagSF;
    TWeight_MuonEffUp     = NormWeight * ElecSF   * (MuonSF+TMath::Sqrt(TMath::Power(MuonSFUp-MuonSF,2)+TMath::Power(MuonSF*0.0122,2)))*TrigSF*PUSF*BtagSF;
    TWeight_MuonEffDown   = NormWeight * ElecSF   * (MuonSF-TMath::Sqrt(TMath::Power(MuonSFDo-MuonSF,2)+TMath::Power(MuonSF*0.0122,2)))*TrigSF*PUSF*BtagSF;
    TWeight_TrigUp        = NormWeight * lepSF    * (TrigSF+TrigSFerr) * PUSF   * BtagSF;
    TWeight_TrigDown      = NormWeight * lepSF    * (TrigSF-TrigSFerr) * PUSF   * BtagSF;
    TWeight_PUDown        = NormWeight * lepSF    * TrigSF * PUSF_Up   * BtagSF;
    TWeight_PUUp          = NormWeight * lepSF    * TrigSF * PUSF_Down * BtagSF;
    TWeight_BtagUp        = NormWeight * ElecSF   * MuonSF * TrigSF * PUSF * BtagSFBtagUp    ;
    TWeight_BtagDown      = NormWeight * ElecSF   * MuonSF * TrigSF * PUSF * BtagSFBtagDown  ;
    TWeight_MistagUp      = NormWeight * ElecSF   * MuonSF * TrigSF * PUSF * BtagSFMistagUp  ;
    TWeight_MistagDown    = NormWeight * ElecSF   * MuonSF * TrigSF * PUSF * BtagSFMistagDown;
  }
}


void TWTTbarAnalysis::SetMinimaAndMaxima() {
  //   Setting maximum value of the unfolding candidate variables.
  if (TLep1_Pt               >= 300)         TLep1_Pt             = 299.999;
  if (TLep1_E                >= 350)         TLep1_E              = 349.999;
  if (TLep1_Eta              >= 2.4)         TLep1_Eta            = 2.39999;
  if (TLep2_Pt               >= 150)         TLep2_Pt             = 149.999;
  if (TLep2_E                >= 250)         TLep2_E              = 249.999;
  if (TLep2_Eta              >= 2.4)         TLep2_Eta            = 2.39999;
  if (TJet1_Pt               >= 300)         TJet1_Pt             = 299.999;
  if (TJet1_E                >= 400)         TJet1_E              = 399.999;
  if (TJet1_Eta              >= 2.4)         TJet1_Eta            = 2.39999;
  if (TJet2_Pt               >= 300)         TJet2_Pt             = 299.999;
  if (TJet2_E                >= 400)         TJet2_E              = 399.999;
  if (TJet2_Eta              >= 2.4)         TJet2_Eta            = 2.39999;
  if (TLep1Lep2_Pt           >= 200)         TLep1Lep2_Pt         = 199.999;
  if (TLep1Lep2_M            >= 300)         TLep1Lep2_M          = 299.999;
  if (TLep1Lep2_DPhi         >= TMath::Pi()) TLep1Lep2_DPhi       = 3.14;
  if (TLep1Jet1_M            >= 400)         TLep1Jet1_M          = 399.999;
  if (TLep1Jet1_DPhi         >= TMath::Pi()) TLep1Jet1_DPhi       = 3.14;
  if (TLep1Jet2_M            >= 400)         TLep1Jet2_M          = 399.999;
  if (TLep1Jet2_DPhi         >= TMath::Pi()) TLep1Jet2_DPhi       = 3.14;
  if (TLep2Jet1_M            >= 300)         TLep2Jet1_M          = 299.999;
  if (TLep2Jet1_DPhi         >= TMath::Pi()) TLep2Jet1_DPhi       = 3.14;
  if (TLep2Jet2_M            >= 300)         TLep2Jet2_M          = 299.999;
  if (TLep2Jet2_DPhi         >= TMath::Pi()) TLep2Jet2_DPhi       = 3.14;
  if (Sys_Pt                 >= 700)         Sys_Pt               = 699.999;
  if (Sys_E                  >= 700)         Sys_E                = 699.999;
  if (Sys_Eta                >= 2.4)         Sys_Eta              = 2.39999;
  if (Sys_M                  >= 700)         Sys_M                = 699.999;
  if (Sys_Pz                 >= 450)         Sys_Pz               = 449.999;
  if (TMiniMax               >= 420)         TMiniMax             = 419.999;
  if (THT                    >= 600)         THT                  = 599.999;
  if (TMET                   >= 200)         TMET                 = 199.999;
  
  if (TLep1_Pt               < 0)         TLep1_Pt                = 0;
  if (TLep1_E                < 0)         TLep1_E                 = 0;
  if (TLep1_Eta              <= -2.4)     TLep1_Eta               = -2.39999;
  if (TLep2_Pt               < 0)         TLep2_Pt                = 0;
  if (TLep2_E                < 0)         TLep2_E                 = 0;
  if (TLep2_Eta              <= -2.4)     TLep2_Eta               = -2.39999;
  if (TJet1_Pt               < 0)         TJet1_Pt                = 0;
  if (TJet1_E                < 0)         TJet1_E                 = 0;
  if (TJet1_Eta              <= -2.4)     TJet1_Eta               = -2.39999;
  if (TJet2_Pt               < 0)         TJet2_Pt                = 0;
  if (TJet2_E                < 0)         TJet2_E                 = 0;
  if (TJet2_Eta              <= -2.4)     TJet2_Eta               = -2.39999;
  if (TLep1Lep2_Pt           < 0)         TLep1Lep2_Pt            = 0;
  if (TLep1Lep2_M            < 0)         TLep1Lep2_M             = 0;
  if (TLep1Lep2_DPhi         <= -TMath::Pi()) TLep1Lep2_DPhi      = -3.14;
  if (TLep1Jet1_M            < 0)         TLep1Jet1_M             = 0;
  if (TLep1Jet1_DPhi         <= -TMath::Pi())  TLep1Jet1_DPhi     = -3.14;
  if (TLep1Jet2_M            < 0)         TLep1Jet2_M             = 0;
  if (TLep1Jet2_DPhi         <= -TMath::Pi())  TLep1Jet2_DPhi     = -3.14;
  if (TLep2Jet1_M            < 0)         TLep2Jet1_M             = 0;
  if (TLep2Jet1_DPhi         <= -TMath::Pi()) TLep2Jet1_DPhi      = -3.14;
  if (TLep2Jet2_M            < 0)         TLep2Jet2_M             = 0;
  if (TLep2Jet2_DPhi         <= -TMath::Pi()) TLep2Jet2_DPhi      = -3.14;
  if (Sys_Pt                 < 0)         Sys_Pt                  = 0;
  if (Sys_E                  < 0)         Sys_E                   = 0;
  if (Sys_Eta                <= -2.4)     Sys_Eta                 = -2.39999;
  if (Sys_M                  < 0)         Sys_M                   = 0;
  if (Sys_Pz                 <= -450)     Sys_Pz                  = -449.999;
  if (TMiniMax               < 0)         TMiniMax                = 0;
  if (THT                    < 0)         THT                     = 0;
  if (TMET                   < 0)         TMET                    = 0;
  
  
  if (TLep1_PtJESUp          >= 300)      TLep1_PtJESUp           = 299.999;
  if (TLep1_EJESUp           >= 350)      TLep1_EJESUp            = 349.999;
  if (TLep1_EtaJESUp         >= 2.4)      TLep1_EtaJESUp          = 2.39999;
  if (TLep2_PtJESUp          >= 150)      TLep2_PtJESUp           = 149.999;
  if (TLep2_EJESUp           >= 250)      TLep2_EJESUp            = 249.999;
  if (TLep2_EtaJESUp         >= 2.4)      TLep2_EtaJESUp          = 2.39999;
  if (TJet1_PtJESUp          >= 300)      TJet1_PtJESUp           = 299.999;
  if (TJet1_EJESUp           >= 400)      TJet1_EJESUp            = 399.999;
  if (TJet1_EtaJESUp         >= 2.4)      TJet1_EtaJESUp          = 2.39999;
  if (TJet2_PtJESUp          >= 300)      TJet2_PtJESUp           = 299.999;
  if (TJet2_EJESUp           >= 400)      TJet2_EJESUp            = 399.999;
  if (TJet2_EtaJESUp         >= 2.4)      TJet2_EtaJESUp          = 2.39999;
  if (TLep1Lep2_PtJESUp      >= 200)      TLep1Lep2_PtJESUp       = 199.999;
  if (TLep1Lep2_MJESUp       >= 300)      TLep1Lep2_MJESUp        = 299.999;
  if (TLep1Lep2_DPhiJESUp    >= TMath::Pi()) TLep1Lep2_DPhiJESUp  = 3.14;
  if (TLep1Jet1_MJESUp       >= 400)      TLep1Jet1_MJESUp        = 399.999;
  if (TLep1Jet1_DPhiJESUp    >= TMath::Pi()) TLep1Jet1_DPhiJESUp  = 3.14;
  if (TLep1Jet2_MJESUp       >= 400)      TLep1Jet2_MJESUp        = 399.999;
  if (TLep1Jet2_DPhiJESUp    >= TMath::Pi()) TLep1Jet2_DPhiJESUp  = 3.14;
  if (TLep2Jet1_MJESUp       >= 300)      TLep2Jet1_MJESUp        = 299.999;
  if (TLep2Jet1_DPhiJESUp    >= TMath::Pi()) TLep2Jet1_DPhiJESUp  = 3.14;
  if (TLep2Jet2_MJESUp       >= 300)      TLep2Jet2_MJESUp        = 299.999;
  if (TLep2Jet2_DPhiJESUp    >= TMath::Pi()) TLep2Jet2_DPhiJESUp  = 3.14;
  if (Sys_PtJESUp            >= 700)      Sys_PtJESUp             = 699.999;
  if (Sys_EJESUp             >= 700)      Sys_EJESUp              = 699.999;
  if (Sys_EtaJESUp           >= 2.4)      Sys_EtaJESUp            = 2.39999;
  if (Sys_MJESUp             >= 700)      Sys_MJESUp              = 699.999;
  if (Sys_PzJESUp            >= 450)      Sys_PzJESUp             = 449.999;
  if (TMiniMaxJESUp          >= 420)      TMiniMaxJESUp           = 419.999;
  if (THTJESUp               >= 600)      THTJESUp                = 599.999;
  if (TMETJESUp              >= 200)      TMETJESUp               = 199.999;
  
  if (TLep1_PtJESUp          < 0)          TLep1_PtJESUp          = 0;
  if (TLep1_EJESUp           < 0)          TLep1_EJESUp           = 0;
  if (TLep1_EtaJESUp         <= -2.4)      TLep1_EtaJESUp         = -2.39999;
  if (TLep2_PtJESUp          < 0)          TLep2_PtJESUp          = 0;
  if (TLep2_EJESUp           < 0)          TLep2_EJESUp           = 0;
  if (TLep2_EtaJESUp         <= -2.4)      TLep2_EtaJESUp         = -2.39999;
  if (TJet1_PtJESUp          < 0)          TJet1_PtJESUp          = 0;
  if (TJet1_EJESUp           < 0)          TJet1_EJESUp           = 0;
  if (TJet1_EtaJESUp         <= -2.4)      TJet1_EtaJESUp         = -2.39999;
  if (TJet2_PtJESUp          < 0)          TJet2_PtJESUp          = 0;
  if (TJet2_EJESUp           < 0)          TJet2_EJESUp           = 0;
  if (TJet2_EtaJESUp         <= -2.4)      TJet2_EtaJESUp         = -2.39999;
  if (TLep1Lep2_PtJESUp      < 0)          TLep1Lep2_PtJESUp      = 0;
  if (TLep1Lep2_MJESUp       < 0)          TLep1Lep2_MJESUp       = 0;
  if (TLep1Lep2_DPhiJESUp    <= -TMath::Pi()) TLep1Lep2_DPhiJESUp = -3.14;
  if (TLep1Jet1_MJESUp       < 0)          TLep1Jet1_MJESUp       = 0;
  if (TLep1Jet1_DPhiJESUp    <= -TMath::Pi())  TLep1Jet1_DPhiJESUp= -3.14;
  if (TLep1Jet2_MJESUp       < 0)          TLep1Jet2_MJESUp       = 0;
  if (TLep1Jet2_DPhiJESUp    <= -TMath::Pi())  TLep1Jet2_DPhiJESUp= -3.14;
  if (TLep2Jet1_MJESUp       < 0)          TLep2Jet1_MJESUp       = 0;
  if (TLep2Jet1_DPhiJESUp    <= -TMath::Pi()) TLep2Jet1_DPhiJESUp = -3.14;
  if (TLep2Jet2_MJESUp       < 0)          TLep2Jet2_MJESUp       = 0;
  if (TLep2Jet2_DPhiJESUp    <= -TMath::Pi()) TLep2Jet2_DPhiJESUp = -3.14;
  if (Sys_PtJESUp            < 0)          Sys_PtJESUp            = 0;
  if (Sys_EJESUp             < 0)          Sys_EJESUp             = 0;
  if (Sys_EtaJESUp           <= -2.4)      Sys_EtaJESUp           = -2.39999;
  if (Sys_MJESUp             < 0)          Sys_MJESUp             = 0;
  if (Sys_PzJESUp            <= -450)      Sys_PzJESUp            = -449.999;
  if (TMiniMaxJESUp          < 0)          TMiniMaxJESUp          = 0;
  if (THTJESUp               < 0)          THTJESUp               = 0;
  if (TMETJESUp              < 0)          TMETJESUp              = 0;
  
  
  if (TLep1_PtJESDown        >= 300)       TLep1_PtJESDown        = 299.999;
  if (TLep1_EJESDown         >= 350)       TLep1_EJESDown         = 349.999;
  if (TLep1_EtaJESDown       >= 2.4)       TLep1_EtaJESDown       = 2.39999;
  if (TLep2_PtJESDown        >= 150)       TLep2_PtJESDown        = 149.999;
  if (TLep2_EJESDown         >= 250)       TLep2_EJESDown         = 249.999;
  if (TLep2_EtaJESDown       >= 2.4)       TLep2_EtaJESDown       = 2.39999;
  if (TJet1_PtJESDown        >= 300)       TJet1_PtJESDown        = 299.999;
  if (TJet1_EJESDown         >= 400)       TJet1_EJESDown         = 399.999;
  if (TJet1_EtaJESDown       >= 2.4)       TJet1_EtaJESDown       = 2.39999;
  if (TJet2_PtJESDown        >= 300)       TJet2_PtJESDown        = 299.999;
  if (TJet2_EJESDown         >= 400)       TJet2_EJESDown         = 399.999;
  if (TJet2_EtaJESDown       >= 2.4)       TJet2_EtaJESDown       = 2.39999;
  if (TLep1Lep2_PtJESDown    >= 200)       TLep1Lep2_PtJESDown    = 199.999;
  if (TLep1Lep2_MJESDown     >= 300)       TLep1Lep2_MJESDown     = 299.999;
  if (TLep1Lep2_DPhiJESDown  >= TMath::Pi()) TLep1Lep2_DPhiJESDown= 3.14;
  if (TLep1Jet1_MJESDown     >= 400)       TLep1Jet1_MJESDown     = 399.999;
  if (TLep1Jet1_DPhiJESDown  >= TMath::Pi()) TLep1Jet1_DPhiJESDown= 3.14;
  if (TLep1Jet2_MJESDown     >= 400)       TLep1Jet2_MJESDown     = 399.999;
  if (TLep1Jet2_DPhiJESDown  >= TMath::Pi()) TLep1Jet2_DPhiJESDown= 3.14;
  if (TLep2Jet1_MJESDown     >= 300)       TLep2Jet1_MJESDown     = 299.999;
  if (TLep2Jet1_DPhiJESDown  >= TMath::Pi()) TLep2Jet1_DPhiJESDown= 3.14;
  if (TLep2Jet2_MJESDown     >= 300)       TLep2Jet2_MJESDown     = 299.999;
  if (TLep2Jet2_DPhiJESDown  >= TMath::Pi()) TLep2Jet2_DPhiJESDown= 3.14;
  if (Sys_PtJESDown          >= 700)       Sys_PtJESDown          = 699.999;
  if (Sys_EJESDown           >= 700)       Sys_EJESDown           = 699.999;
  if (Sys_EtaJESDown         >= 2.4)       Sys_EtaJESDown         = 2.39999;
  if (Sys_MJESDown           >= 700)       Sys_MJESDown           = 699.999;
  if (Sys_PzJESDown          >= 450)       Sys_PzJESDown          = 449.999;
  if (TMiniMaxJESDown        >= 420)       TMiniMaxJESDown        = 419.999;
  if (THTJESDown             >= 600)       THTJESDown             = 599.999;
  if (TMETJESDown            >= 200)       TMETJESDown            = 199.999;
  
  if (TLep1_PtJESDown        < 0)          TLep1_PtJESDown        = 0;
  if (TLep1_EJESDown         < 0)          TLep1_EJESDown         = 0;
  if (TLep1_EtaJESDown       <= -2.4)      TLep1_EtaJESDown       = -2.39999;
  if (TLep2_PtJESDown        < 0)          TLep2_PtJESDown        = 0;
  if (TLep2_EJESDown         < 0)          TLep2_EJESDown         = 0;
  if (TLep2_EtaJESDown       <= -2.4)      TLep2_EtaJESDown       = -2.39999;
  if (TJet1_PtJESDown        < 0)          TJet1_PtJESDown        = 0;
  if (TJet1_EJESDown         < 0)          TJet1_EJESDown         = 0;
  if (TJet1_EtaJESDown       <= -2.4)      TJet1_EtaJESDown       = -2.39999;
  if (TJet2_PtJESDown        < 0)          TJet2_PtJESDown        = 0;
  if (TJet2_EJESDown         < 0)          TJet2_EJESDown         = 0;
  if (TJet2_EtaJESDown       <= -2.4)      TJet2_EtaJESDown       = -2.39999;
  if (TLep1Lep2_PtJESDown    < 0)          TLep1Lep2_PtJESDown    = 0;
  if (TLep1Lep2_MJESDown     < 0)          TLep1Lep2_MJESDown     = 0;
  if (TLep1Lep2_DPhiJESDown  <= -TMath::Pi()) TLep1Lep2_DPhiJESDown = -3.14;
  if (TLep1Jet1_MJESDown     < 0)          TLep1Jet1_MJESDown       = 0;
  if (TLep1Jet1_DPhiJESDown  <= -TMath::Pi())  TLep1Jet1_DPhiJESDown= -3.14;
  if (TLep1Jet2_MJESDown     < 0)          TLep1Jet2_MJESDown       = 0;
  if (TLep1Jet2_DPhiJESDown  <= -TMath::Pi())  TLep1Jet2_DPhiJESDown= -3.14;
  if (TLep2Jet1_MJESDown     < 0)          TLep2Jet1_MJESDown       = 0;
  if (TLep2Jet1_DPhiJESDown  <= -TMath::Pi()) TLep2Jet1_DPhiJESDown = -3.14;
  if (TLep2Jet2_MJESDown     < 0)          TLep2Jet2_MJESDown       = 0;
  if (TLep2Jet2_DPhiJESDown  <= -TMath::Pi()) TLep2Jet2_DPhiJESDown = -3.14;
  if (Sys_PtJESDown          < 0)          Sys_PtJESDown         = 0;
  if (Sys_EJESDown           < 0)          Sys_EJESDown          = 0;
  if (Sys_EtaJESDown         <= -2.4)      Sys_EtaJESDown        = -2.39999;
  if (Sys_MJESDown           < 0)          Sys_MJESDown          = 0;
  if (Sys_PzJESDown          <= -450)      Sys_PzJESDown         = -449.999;
  if (TMiniMaxJESDown        < 0)          TMiniMaxJESDown       = 0;
  if (THTJESDown             < 0)          THTJESDown            = 0;
  if (TMETJESDown            < 0)          TMETJESDown           = 0;
  
  if (TLep1_PtJERUp          >= 300)       TLep1_PtJERUp         = 299.999;
  if (TLep1_EJERUp           >= 350)       TLep1_EJERUp          = 349.999;
  if (TLep1_EtaJERUp         >= 2.4)       TLep1_EtaJERUp        = 2.39999;
  if (TLep2_PtJERUp          >= 150)       TLep2_PtJERUp         = 149.999;
  if (TLep2_EJERUp           >= 250)       TLep2_EJERUp          = 249.999;
  if (TLep2_EtaJERUp         >= 2.4)       TLep2_EtaJERUp        = 2.39999;
  if (TJet1_PtJERUp          >= 300)       TJet1_PtJERUp         = 299.999;
  if (TJet1_EJERUp           >= 400)       TJet1_EJERUp          = 399.999;
  if (TJet1_EtaJERUp         >= 2.4)       TJet1_EtaJERUp        = 2.39999;
  if (TJet2_PtJERUp          >= 300)       TJet2_PtJERUp         = 299.999;
  if (TJet2_EJERUp           >= 400)       TJet2_EJERUp          = 399.999;
  if (TJet2_EtaJERUp         >= 2.4)       TJet2_EtaJERUp        = 2.39999;
  if (TLep1Lep2_PtJERUp      >= 200)       TLep1Lep2_PtJERUp     = 199.999;
  if (TLep1Lep2_MJERUp       >= 300)       TLep1Lep2_MJERUp      = 299.999;
  if (TLep1Lep2_DPhiJERUp    >= TMath::Pi()) TLep1Lep2_DPhiJERUp = 3.14;
  if (TLep1Jet1_MJERUp       >= 400)       TLep1Jet1_MJERUp      = 399.999;
  if (TLep1Jet1_DPhiJERUp    >= TMath::Pi()) TLep1Jet1_DPhiJERUp = 3.14;
  if (TLep1Jet2_MJERUp       >= 400)       TLep1Jet2_MJERUp      = 399.999;
  if (TLep1Jet2_DPhiJERUp    >= TMath::Pi()) TLep1Jet2_DPhiJERUp = 3.14;
  if (TLep2Jet1_MJERUp       >= 300)       TLep2Jet1_MJERUp      = 299.999;
  if (TLep2Jet1_DPhiJERUp    >= TMath::Pi()) TLep2Jet1_DPhiJERUp = 3.14;
  if (TLep2Jet2_MJERUp       >= 300)       TLep2Jet2_MJERUp      = 299.999;
  if (TLep2Jet2_DPhiJERUp    >= TMath::Pi()) TLep2Jet2_DPhiJERUp = 3.14;
  if (Sys_PtJERUp            >= 700)       Sys_PtJERUp           = 699.999;
  if (Sys_EJERUp             >= 700)       Sys_EJERUp            = 699.999;
  if (Sys_EtaJERUp           >= 2.4)       Sys_EtaJERUp          = 2.39999;
  if (Sys_MJERUp             >= 700)       Sys_MJERUp            = 699.999;
  if (Sys_PzJERUp            >= 450)       Sys_PzJERUp           = 449.999;
  if (TMiniMaxJERUp          >= 420)       TMiniMaxJERUp         = 419.999;
  if (THTJERUp               >= 600)       THTJERUp              = 599.999;
  if (TMETJERUp              >= 200)       TMETJERUp             = 199.999;
  
  if (TLep1_PtJERUp          < 0)          TLep1_PtJERUp         = 0;
  if (TLep1_EJERUp           < 0)          TLep1_EJERUp          = 0;
  if (TLep1_EtaJERUp         <= -2.4)      TLep1_EtaJERUp        = -2.39999;
  if (TLep2_PtJERUp          < 0)          TLep2_PtJERUp         = 0;
  if (TLep2_EJERUp           < 0)          TLep2_EJERUp          = 0;
  if (TLep2_EtaJERUp         <= -2.4)      TLep2_EtaJERUp        = -2.39999;
  if (TJet1_PtJERUp          < 0)          TJet1_PtJERUp         = 0;
  if (TJet1_EJERUp           < 0)          TJet1_EJERUp          = 0;
  if (TJet1_EtaJERUp         <= -2.4)      TJet1_EtaJERUp        = -2.39999;
  if (TJet2_PtJERUp          < 0)          TJet2_PtJERUp         = 0;
  if (TJet2_EJERUp           < 0)          TJet2_EJERUp          = 0;
  if (TJet2_EtaJERUp         <= -2.4)      TJet2_EtaJERUp        = -2.39999;
  if (TLep1Lep2_PtJERUp      < 0)          TLep1Lep2_PtJERUp     = 0;
  if (TLep1Lep2_MJERUp       < 0)          TLep1Lep2_MJERUp      = 0;
  if (TLep1Lep2_DPhiJERUp    <= -TMath::Pi()) TLep1Lep2_DPhiJERUp = -3.14;
  if (TLep1Jet1_MJERUp       < 0)          TLep1Jet1_MJERUp       = 0;
  if (TLep1Jet1_DPhiJERUp    <= -TMath::Pi())  TLep1Jet1_DPhiJERUp= -3.14;
  if (TLep1Jet2_MJERUp       < 0)          TLep1Jet2_MJERUp       = 0;
  if (TLep1Jet2_DPhiJERUp    <= -TMath::Pi())  TLep1Jet2_DPhiJERUp= -3.14;
  if (TLep2Jet1_MJERUp       < 0)          TLep2Jet1_MJERUp       = 0;
  if (TLep2Jet1_DPhiJERUp    <= -TMath::Pi()) TLep2Jet1_DPhiJERUp = -3.14;
  if (TLep2Jet2_MJERUp       < 0)          TLep2Jet2_MJERUp       = 0;
  if (TLep2Jet2_DPhiJERUp    <= -TMath::Pi()) TLep2Jet2_DPhiJERUp = -3.14;
  if (Sys_PtJERUp            < 0)          Sys_PtJERUp           = 0;
  if (Sys_EJERUp             < 0)          Sys_EJERUp            = 0;
  if (Sys_EtaJERUp           <= -2.4)      Sys_EtaJERUp          = -2.39999;
  if (Sys_MJERUp             < 0)          Sys_MJERUp            = 0;
  if (Sys_PzJERUp            <= -450)      Sys_PzJERUp           = -449.999;
  if (TMiniMaxJERUp          < 0)          TMiniMaxJERUp         = 0;
  if (THTJERUp               < 0)          THTJERUp              = 0;
  if (TMETJERUp              < 0)          TMETJERUp             = 0;
  
  // Particle level variables
  if (TDressLep1_Pt               >= 300)         TDressLep1_Pt             = 299.999;
  if (TDressLep1_E                >= 350)         TDressLep1_E              = 349.999;
  if (TDressLep1_Eta              >= 2.4)         TDressLep1_Eta            = 2.39999;
  if (TDressLep2_Pt               >= 150)         TDressLep2_Pt             = 149.999;
  if (TDressLep2_E                >= 250)         TDressLep2_E              = 249.999;
  if (TDressLep2_Eta              >= 2.4)         TDressLep2_Eta            = 2.39999;
  if (TDressJet1_Pt               >= 300)         TDressJet1_Pt             = 299.999;
  if (TDressJet1_E                >= 400)         TDressJet1_E              = 399.999;
  if (TDressJet1_Eta              >= 2.4)         TDressJet1_Eta            = 2.39999;
  if (TDressJet2_Pt               >= 300)         TDressJet2_Pt             = 299.999;
  if (TDressJet2_E                >= 400)         TDressJet2_E              = 399.999;
  if (TDressJet2_Eta              >= 2.4)         TDressJet2_Eta            = 2.39999;
  if (TDressLep1Lep2_Pt           >= 200)         TDressLep1Lep2_Pt         = 199.999;
  if (TDressLep1Lep2_M            >= 300)         TDressLep1Lep2_M          = 299.999;
  if (TDressLep1Lep2_DPhi         >= TMath::Pi()) TDressLep1Lep2_DPhi       = 3.14;
  if (TDressLep1Jet1_M            >= 400)         TDressLep1Jet1_M          = 399.999;
  if (TDressLep1Jet1_DPhi         >= TMath::Pi()) TDressLep1Jet1_DPhi       = 3.14;
  if (TDressLep1Jet2_M            >= 400)         TDressLep1Jet2_M          = 399.999;
  if (TDressLep1Jet2_DPhi         >= TMath::Pi()) TDressLep1Jet2_DPhi       = 3.14;
  if (TDressLep2Jet1_M            >= 300)         TDressLep2Jet1_M          = 299.999;
  if (TDressLep2Jet1_DPhi         >= TMath::Pi()) TDressLep2Jet1_DPhi       = 3.14;
  if (TDressLep2Jet2_M            >= 300)         TDressLep2Jet2_M          = 299.999;
  if (TDressLep2Jet2_DPhi         >= TMath::Pi()) TDressLep2Jet2_DPhi       = 3.14;
  if (DressSys_Pt                 >= 700)         DressSys_Pt               = 699.999;
  if (DressSys_E                  >= 700)         DressSys_E                = 699.999;
  if (DressSys_Eta                >= 2.4)         DressSys_Eta              = 2.39999;
  if (DressSys_M                  >= 700)         DressSys_M                = 699.999;
  if (DressSys_Pz                 >= 450)         DressSys_Pz               = 449.999;
  if (TDressMiniMax               >= 420)         TDressMiniMax             = 419.999;
  if (TDressHT                    >= 600)         TDressHT                  = 599.999;
  if (TDressMET                   >= 200)         TDressMET                 = 199.999;
  
  if (TDressLep1_Pt               < 0)         TDressLep1_Pt                = 0;
  if (TDressLep1_E                < 0)         TDressLep1_E                 = 0;
  if (TDressLep1_Eta              <= -2.4)     TDressLep1_Eta               = -2.39999;
  if (TDressLep2_Pt               < 0)         TDressLep2_Pt                = 0;
  if (TDressLep2_E                < 0)         TDressLep2_E                 = 0;
  if (TDressLep2_Eta              <= -2.4)     TDressLep2_Eta               = -2.39999;
  if (TDressJet1_Pt               < 0)         TDressJet1_Pt                = 0;
  if (TDressJet1_E                < 0)         TDressJet1_E                 = 0;
  if (TDressJet1_Eta              <= -2.4)     TDressJet1_Eta               = -2.39999;
  if (TDressJet2_Pt               < 0)         TDressJet2_Pt                = 0;
  if (TDressJet2_E                < 0)         TDressJet2_E                 = 0;
  if (TDressJet2_Eta              <= -2.4)     TDressJet2_Eta               = -2.39999;
  if (TDressLep1Lep2_Pt           < 0)         TDressLep1Lep2_Pt            = 0;
  if (TDressLep1Lep2_M            < 0)         TDressLep1Lep2_M             = 0;
  if (TDressLep1Lep2_DPhi         <= -TMath::Pi()) TDressLep1Lep2_DPhi      = -3.14;
  if (TDressLep1Jet1_M            < 0)         TDressLep1Jet1_M             = 0;
  if (TDressLep1Jet1_DPhi         <= -TMath::Pi())  TDressLep1Jet1_DPhi     = -3.14;
  if (TDressLep1Jet2_M            < 0)         TDressLep1Jet2_M             = 0;
  if (TDressLep1Jet2_DPhi         <= -TMath::Pi())  TDressLep1Jet2_DPhi     = -3.14;
  if (TDressLep2Jet1_M            < 0)         TDressLep2Jet1_M             = 0;
  if (TDressLep2Jet1_DPhi         <= -TMath::Pi()) TDressLep2Jet1_DPhi      = -3.14;
  if (TDressLep2Jet2_M            < 0)         TDressLep2Jet2_M             = 0;
  if (TDressLep2Jet2_DPhi         <= -TMath::Pi()) TDressLep2Jet2_DPhi      = -3.14;
  if (DressSys_Pt                 < 0)         DressSys_Pt                  = 0;
  if (DressSys_E                  < 0)         DressSys_E                   = 0;
  if (DressSys_Eta                <= -2.4)     DressSys_Eta                 = -2.39999;
  if (DressSys_M                  < 0)         DressSys_M                   = 0;
  if (DressSys_Pz                 <= -450)     DressSys_Pz                  = -449.999;
  if (TDressMiniMax               < 0)         TDressMiniMax                = 0;
  if (TDressHT                    < 0)         TDressHT                     = 0;
  if (TDressMET                   < 0)         TDressMET                    = 0;
  //   Setting dummy value for gen events that don't pass the reco selection for
  // unfolding procedures.
  if (TPassDress && !TPassReco) {
    TLep1_Pt               = 99999;
    TLep1_E                = 99999;
    TLep1_Eta              = 99999;
    TLep2_Pt               = 99999;
    TLep2_E                = 99999;
    TLep2_Eta              = 99999;
    TJet1_Pt               = 99999;
    TJet1_E                = 99999;
    TJet1_Eta              = 99999;
    TJet2_Pt               = 99999;
    TJet2_E                = 99999;
    TJet2_Eta              = 99999;
    TLep1Lep2_Pt           = 99999;
    TLep1Lep2_M            = 99999;
    TLep1Lep2_DPhi         = 99999;
    TLep1Jet1_M            = 99999;
    TLep1Jet1_DPhi         = 99999;
    TLep1Jet2_M            = 99999;
    TLep1Jet2_DPhi         = 99999;
    TLep2Jet1_M            = 99999;
    TLep2Jet1_DPhi         = 99999;
    TLep2Jet2_M            = 99999;
    TLep2Jet2_DPhi         = 99999;
    Sys_Pt                 = 99999;
    Sys_E                  = 99999;
    Sys_Eta                = 99999;
    Sys_M                  = 99999;
    Sys_Pz                 = 99999;
    TMiniMax               = 99999;
    THT                    = 99999;
    TMET                   = 99999;
  }
  if (TPassDress && !TPassRecoJESUp) {
    TLep1_PtJESUp          = 99999;
    TLep1_EJESUp           = 99999;
    TLep1_EtaJESUp         = 99999;
    TLep2_PtJESUp          = 99999;
    TLep2_EJESUp           = 99999;
    TLep2_EtaJESUp         = 99999;
    TJet1_PtJESUp          = 99999;
    TJet1_EJESUp           = 99999;
    TJet1_EtaJESUp         = 99999;
    TJet2_PtJESUp          = 99999;
    TJet2_EJESUp           = 99999;
    TJet2_EtaJESUp         = 99999;
    TLep1Lep2_PtJESUp      = 99999;
    TLep1Lep2_MJESUp       = 99999;
    TLep1Lep2_DPhiJESUp    = 99999;
    TLep1Jet1_MJESUp       = 99999;
    TLep1Jet1_DPhiJESUp    = 99999;
    TLep1Jet2_MJESUp       = 99999;
    TLep1Jet2_DPhiJESUp    = 99999;
    TLep2Jet1_MJESUp       = 99999;
    TLep2Jet1_DPhiJESUp    = 99999;
    TLep2Jet2_MJESUp       = 99999;
    TLep2Jet2_DPhiJESUp    = 99999;
    Sys_PtJESUp            = 99999;
    Sys_EJESUp             = 99999;
    Sys_EtaJESUp           = 99999;
    Sys_MJESUp             = 99999;
    Sys_PzJESUp            = 99999;
    TMiniMaxJESUp          = 99999;
    THTJESUp               = 99999;
    TMETJESUp              = 99999;
  }
  if (TPassDress && !TPassRecoJESDown) {
    TLep1_PtJESDown        = 99999;
    TLep1_EJESDown         = 99999;
    TLep1_EtaJESDown       = 99999;
    TLep2_PtJESDown        = 99999;
    TLep2_EJESDown         = 99999;
    TLep2_EtaJESDown       = 99999;
    TJet1_PtJESDown        = 99999;
    TJet1_EJESDown         = 99999;
    TJet1_EtaJESDown       = 99999;
    TJet2_PtJESDown        = 99999;
    TJet2_EJESDown         = 99999;
    TJet2_EtaJESDown       = 99999;
    TLep1Lep2_PtJESDown    = 99999;
    TLep1Lep2_MJESDown     = 99999;
    TLep1Lep2_DPhiJESDown  = 99999;
    TLep1Jet1_MJESDown     = 99999;
    TLep1Jet1_DPhiJESDown  = 99999;
    TLep1Jet2_MJESDown     = 99999;
    TLep1Jet2_DPhiJESDown  = 99999;
    TLep2Jet1_MJESDown     = 99999;
    TLep2Jet1_DPhiJESDown  = 99999;
    TLep2Jet2_MJESDown     = 99999;
    TLep2Jet2_DPhiJESDown  = 99999;
    Sys_PtJESDown          = 99999;
    Sys_EJESDown           = 99999;
    Sys_EtaJESDown         = 99999;
    Sys_MJESDown           = 99999;
    Sys_PzJESDown          = 99999;
    TMiniMaxJESDown        = 99999;
    THTJESDown             = 99999;
    TMETJESDown            = 99999;
  }
  if (TPassDress && !TPassRecoJERUp) {
    TLep1_PtJERUp          = 99999;
    TLep1_EJERUp           = 99999;
    TLep1_EtaJERUp         = 99999;
    TLep2_PtJERUp          = 99999;
    TLep2_EJERUp           = 99999;
    TLep2_EtaJERUp         = 99999;
    TJet1_PtJERUp          = 99999;
    TJet1_EJERUp           = 99999;
    TJet1_EtaJERUp         = 99999;
    TJet2_PtJERUp          = 99999;
    TJet2_EJERUp           = 99999;
    TJet2_EtaJERUp         = 99999;
    TLep1Lep2_PtJERUp      = 99999;
    TLep1Lep2_MJERUp       = 99999;
    TLep1Lep2_DPhiJERUp    = 99999;
    TLep1Jet1_MJERUp       = 99999;
    TLep1Jet1_DPhiJERUp    = 99999;
    TLep1Jet2_MJERUp       = 99999;
    TLep1Jet2_DPhiJERUp    = 99999;
    TLep2Jet1_MJERUp       = 99999;
    TLep2Jet1_DPhiJERUp    = 99999;
    TLep2Jet2_MJERUp       = 99999;
    TLep2Jet2_DPhiJERUp    = 99999;
    Sys_PtJERUp            = 99999;
    Sys_EJERUp             = 99999;
    Sys_EtaJERUp           = 99999;
    Sys_MJERUp             = 99999;
    Sys_PzJERUp            = 99999;
    TMiniMaxJERUp          = 99999;
    THTJERUp               = 99999;
    TMETJERUp              = 99999;
  }
}
