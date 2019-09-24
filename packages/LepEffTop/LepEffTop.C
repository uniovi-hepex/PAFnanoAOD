///////////////////////////////////////////////////////////////////////////////////////////
//
//  Lepton Selector: create vectors with selected leptons, veto (or fakeable) leptons...
//
//  All SFs and variables are within the Lepton definition
//
//  To do: definition of veto leptons
//
//
///////////////////////////////////////////////////////////////////////////////////////////


#include "LepEffTop.h"


ClassImp(LepEffTop);
LepEffTop::LepEffTop() : PAFChainItemSelector() {}
void LepEffTop::Summary(){}

void LepEffTop::Initialise(){
  year    = GetParam<TString>("year").Atoi();
  gIsData        = GetParam<Bool_t>("IsData");
  selection      = GetParam<TString>("selection");
  gOptions       = GetParam<TString>("_options");
  localPath      = GetParam<TString>("WorkingDir");
  gSelection     = GetSelection(selection);
  gDoLepGood     = gOptions.Contains("LepGood")? true : false;
  gPUWeigth    = gOptions.Contains("PUweight")? true : false;
  gPrefire     = gOptions.Contains("prefire")? true : false;


  gIs2017 = false; gIs2016 = false; gIs2018 = false;
  if     (year == 2017) gIs2017 = true;
  else if(year == 2018) gIs2018 = true;
  else if(year == 2016) gIs2016 = true;

  leptons   = std::vector<Lepton>();
  InitHistos();
}


TString LepEffTop::GetLabel(Int_t i, Int_t j, Int_t k){
  // [pass][prompt][OS]
  TString l1 = i == kPass  ? "Pass"   : "Fail";
  TString l2 = j == kPrompt? "Prompt" : "Fake";
  TString l3 = k == kOS?     "OS"     : "SS";
  TString lab = Form("%s_%s_%s", l1.Data(), l2.Data(), l3.Data());
  return lab;
}

void LepEffTop::InitHistos(){
  TString lab;
  for(Int_t i = 0; i < 2; i++){
    for(Int_t j = 0; j < 2; j++){
      for(Int_t k = 0; k < 2; k++){
        lab = GetLabel(i, j, k);
        hMuonPt [i][j][k] = CreateH1F("H_Muon_Pt_" +lab, "Muon p_{T} (GeV)"   , 150, 0, 150);
        hMuonEta[i][j][k] = CreateH1F("H_Muon_Eta_"+lab, "Muon #eta"          , 120, -2.5, 2-5);
        hMuonPhi[i][j][k] = CreateH1F("H_Muon_Phi_"+lab, "Muon #phi (rad/#pi)", 120, -1, 1);
        hMuonIso[i][j][k] = CreateH1F("H_Muon_Iso_"+lab, "Muon isolation"     , 500, 0, 1);
        hElecPt [i][j][k] = CreateH1F("H_Elec_Pt_" +lab, "Elec p_{T} (GeV)"   , 150, 0, 150);
        hElecEta[i][j][k] = CreateH1F("H_Elec_Eta_"+lab, "Elec #eta"          , 120, -2.5, 2-5);
        hElecPhi[i][j][k] = CreateH1F("H_Elec_Phi_"+lab, "Elec #phi (rad/#pi)", 120, -1, 1);
        hElecIso[i][j][k] = CreateH1F("H_Elec_Iso_"+lab, "Elec isolation"     , 500, 0, 1);
      }
    }
  }
}

//################################################################
//## InsideLoop
//################################################################
void LepEffTop::InsideLoop(){
  // Clear vectors...
  leptons.clear();
  genLeptons.clear();

  // Get parameters
  leptons = GetParam<vector<Lepton>>("looseLeptons");
  njets   = GetParam<Int_t>("NJets");
  nbtags  = GetParam<Int_t>("NBtags");

  // Weights and SFs
  NormWeight     = GetParam<Double_t>("NormWeight");
  TriggerSF      = GetParam<Float_t>("TriggerSF");
  PUSF = 1; PrefWeight = 1;
  if(!gIsData && gPUWeigth){ PUSF       = Get<Float_t>("puWeight");}
  if(!gIsData && gPrefire ){ PrefWeight = Get<Float_t>("PrefireWeight");}

  // Event variables
  PassMETFilters = GetParam<Bool_t>("METfilters");
  PassTrigger    = GetParam<Bool_t>("passTrigger");
  isSS           = GetParam<Bool_t>("isSS");

  Int_t nlep = leptons.size();

  // 2 Leptons, emu
  if(nlep < 2) return;
  if (! ( (leptons.at(0).IsElec() && leptons.at(1).IsMuon()) || (leptons.at(0).IsMuon() && leptons.at(1).IsElec()) ) ) return; 
  Elec = leptons.at(0).IsElec() ? leptons.at(0) : leptons.at(1);
  Muon = leptons.at(1).IsElec() ? leptons.at(0) : leptons.at(1);
  Bool_t passDilep = (Elec.Pt() > 20 && Muon.Pt() > 20 && (Elec.p + Muon.p).M() > 20 && leptons.at(0).Pt() > 25);

  if(PassTrigger && PassMETFilters && passDilep && njets >= 2 && nbtags >= 1)
    FillHistos();
}

void LepEffTop::FillAll(Int_t i, Int_t j, Int_t k){
  Float_t SF = Elec.GetSF()*Muon.GetSF()*TriggerSF;
  Float_t weight = NormWeight*SF*PUSF*PrefWeight;
  hMuonPt [i][j][k]->Fill(Muon.Pt(),  weight);
  hMuonEta[i][j][k]->Fill(Muon.Eta(), weight);
  hMuonPhi[i][j][k]->Fill(Muon.Phi(), weight);
  hMuonIso[i][j][k]->Fill(Muon.GetIso(), weight);
}

void LepEffTop::FillHistos(){
  //cout << Form("Muon iso: %1.4f", Muon.GetIso());
  if(Muon.GetIso() < 0.15){ // Passing
    //cout << "... pass ISO!!! (passing)"<< endl;
    if(!isSS){ // OS
      if(Muon.IsPrompt()){ // Prompt
        FillAll(kPass, kPrompt, kOS);
      }
      else{ // Nonprompt
        FillAll(kPass, kNonprompt, kOS);
      }
    }
    else{ // SS
      if(Muon.IsPrompt()){ // Prompt
        FillAll(kPass, kPrompt, kSS);
      }
      else{ // Nonprompt
        FillAll(kPass, kNonprompt, kSS);
      }
    }
  }
  else{ // Failing
    cout << "... NO pass ISO!!! (failing)"<< endl;
    if(!isSS){ // OS
      if(Muon.IsPrompt()){ // Prompt
        FillAll(kFail, kPrompt, kOS);
      }
      else{ // Nonprompt
        FillAll(kFail, kNonprompt, kOS);
      }
    }
    else{ // SS
      if(Muon.IsPrompt()){ // Prompt
        FillAll(kFail, kPrompt, kSS);
      }
      else{ // Nonprompt
        FillAll(kFail, kNonprompt, kSS);
      }
    }
  }
}

//################################################################
//################################################################
//################################################################
//## Other
//################################################################
/*
void LepEffTop::GetLeptonVariables(Int_t i, int LepType){ // Once per muon, get all the info
  TString LeptonName = abs(LepType) == 11? "Electron" : "Muon";
  rho = Get<Float_t>("fixedGridRhoFastjetCentralCalo");
  if(LepType == 0) LeptonName = "LepGood";
  tP.SetPtEtaPhiM(Get<Float_t>(LeptonName + "_pt", i), Get<Float_t>(LeptonName + "_eta", i), Get<Float_t>(LeptonName + "_phi", i), Get<Float_t>(LeptonName + "_mass", i));

  // PDG id...
  if(LepType != 0){
    type          = TMath::Abs(LepType) == 11 ? 1 : 0;
    pdgid         = LepType;
  }
  else{
    pdgid = TMath::Abs(Get<Int_t>("LepGood_pdgId",i));
    type  = pdgid == 11 ? 1 : 0;
  }

  pt = tP.Pt(); eta = tP.Eta(); energy = tP.Energy();
  charge        = Get<Int_t>(LeptonName + "_charge", i);
  dxy           = TMath::Abs(Get<Float_t>(LeptonName + "_dxy", i));
  dz            = TMath::Abs(Get<Float_t>(LeptonName + "_dz", i));
  miniIso       = Get<Float_t>(LeptonName + "_miniPFRelIso_all",i);
  RelIso03      = Get<Float_t>(LeptonName + "_pfRelIso03_all",i);
  sip           = Get<Float_t>(LeptonName + "_sip3d",i);
  TightCharge   = Get<Int_t>(LeptonName + "_tightCharge",i);      //

  jetindex      = Get<Int_t>(LeptonName + "_jetIdx", i);  //index of the associated jet (-1 if none)
  if(!gIsData){
    genPartIndex  = Get<Int_t>(LeptonName + "_genPartIdx", i);
    genPartFlav   = Get<Int_t>(LeptonName + "_genPartFlav", i);
  }

  RelIso04 = -1; mediumMuonId = -1; SegComp = -1; dEtaSC = -1; HoE = -1; eImpI = -1; 
  lostHits = -1; convVeto = -1; sigmaIEtaIEta = -1; MVAID = -1; R9 = -1; etaSC = -1;
  SF = 1;
  if(pdgid == 13){  // ONLY MUONS
    RelIso04       = Get<Float_t>(LeptonName + "_pfRelIso04_all",i);
    mediumMuonId   = Get<Bool_t>(LeptonName + "_mediumId",i);
    SegComp        = Get<Float_t>(LeptonName + "_segmentComp", i);
    tightVar       = Get<Bool_t>(LeptonName + "_tightId", i);
  }
  else{              // ONLY ELECTRON
    tightVar       = Get<Int_t>(LeptonName + "_cutBased", i);
    dEtaSC         = Get<Float_t>(LeptonName + "_deltaEtaSC", i);
    etaSC          = Get<Float_t>(LeptonName + "_eta") + dEtaSC;
    HoE            = Get<Float_t>(LeptonName + "_hoe", i);
    eImpI          = Get<Float_t>(LeptonName + "_eInvMinusPInv", i);
    lostHits       = Get<UChar_t>(LeptonName + "_lostHits", i);
    convVeto       = Get<Bool_t>(LeptonName + "_convVeto", i);
    sigmaIEtaIEta  = Get<Float_t>(LeptonName + "_sieie", i); // Es esta???
    MVAID          = 0; //= Get<Float_t>("Electron_mvaFall17Iso",i);  //Electron_mvaSpring16GP_WP80, Electron_mvaSpring16GP_WP90, Electron_mvaSpring16HZZ_WPL
    R9             = Get<Float_t>(LeptonName + "_r9",i);
  }
}

void LepEffTop::GetGenLeptonVariables(Int_t i){
  tP.SetPtEtaPhiM(Get<Float_t>("GenDressedLepton_pt", i), Get<Float_t>("GenDressedLepton_eta", i), Get<Float_t>("GenDressedLepton_phi", i), Get<Float_t>("GenDressedLepton_mass", i));
  pdgid = Get<Int_t>("GenDressedLepton_pdgId",i);
  charge = pdgid < 0 ? 1 : -1;
  type = TMath::Abs(pdgid) == 11 ? 1 : 0;
}

void LepEffTop::IsElecBarrel(){
  return TMath::Abs(etaSC) < 1.479;
}

void LepEffTop::GetEffArea(Float_t eta){
  eta = TMath::Abs(eta);
  if     (eta < 1.000) return 0.1440; 
  else if(eta < 1.479) return 0.1562;
  else if(eta < 2.000) return 0.1032;
  else if(eta < 2.200) return 0.0859;
  else if(eta < 2.300) return 0.1116;
  else if(eta < 2.400) return 0.1321;
  else if(eta < 2.500) return 0.1654;
  cout << "[LepEffTop::GetEffArea] ERROR : eta out of range!" << endl;
  return 0;
}

void LepEffTop::IsElecTight_ISOonly(){
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  //  Offline selection criteria for V2  -- tight (94X)
  Float_t cut;
  cut = IsElecBarrel() ? 0.0287+0.506/pt : 0.0445+0.963/pt;
  return RelIso03 < cut;
}

void LepEffTop::IsElecTight_IDonly(){
  Float_t cut;
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/CutBasedElectronIdentificationRun2
  //  Offline selection criteria for V2  -- tight (94X)
  //
  // full5x5_sigmaIetaIeta
  cut = IsElecBarrel() ? 0.0104 : 0.0353;
  if(sigmaIEtaIEta > cut) return false;
  
  // abs(dEtaSeed)
  cut = IsElecBarrel() ? 0.00255 : 0.00501;

  // abs(dPhiIn)
  cut = IsElecBarrel() ? 0.022 : 0.0236;

  // H/E
  cut = IsElecBarrel() ? 0.026 + 1.15/ESC + 0.0324*rho/ESC : 0.0188 + 2.06/ESC + 0.183*rho/ESC;
  if(HoE > cut) return false;

  // abs(1/E-1/p)
  cut = IsElecBarrel() ? 0.159 : 0.0197;
  if(eImpI > cut) return false;

  // expected missing inner hits
  if(lostHits > 1) return false;

  // pass conversion veto
  if(!convVeto) return false;
}
*/

