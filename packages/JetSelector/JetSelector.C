//////////////////////////////////////////////////////////////////////////
//
//  Jet Selector
//  Produce vectors with: genJets, selJets, vetoJets...
//  also some variables like number of b-jets
//
//  All SFs and variables are within the Jet definition
//
//  ToDo: add JER syst (and JER jet variations to selJets)...
//
/////////////////////////////////////////////////////////////////////////


#include "JetSelector.h"
#include <string>

ClassImp(JetSelector);
JetSelector::JetSelector() : PAFChainItemSelector() {
  fBTagSFnom = 0;
  fBTagSFbUp = 0;
  fBTagSFbDo = 0;
  fBTagSFlUp = 0;
  fBTagSFlDo = 0;
  MeasType = "comb";
  minDR = 0;
  jet_MaxEta = 0;
  jet_MinPt = 0;
  vetoJet_minPt = 0;
  BtagSFFS = 1;
}

JetSelector::~JetSelector() {
  delete fBTagSFnom; 
  delete fBTagSFbUp; 
  delete fBTagSFbDo;
  delete fBTagSFlUp; 
  delete fBTagSFlDo;
}

void JetSelector::Summary(){}

void JetSelector::Initialise(){
  gIsData    = GetParam<Bool_t>("IsData");
  selection = GetParam<Int_t>("iSelection");
  gIsFastSim   = GetParam<Bool_t>("IsFastSim");
  gSampleName  = GetParam<TString>("sampleName");
  gDoSys       = GetParam<Bool_t>("doSyst");
  gOptions     = GetParam<TString>("_options");
  gSelection     = GetSelection(selection);

  gIs2017 = false; gIs2016 = false; gIs2018 = false;
  year = 0;
  if(gOptions.Contains("2017")) gIs2017 = true;
  if(gOptions.Contains("2018")) gIs2018 = true;
  if(gOptions.Contains("2016")) gIs2016 = true;
  if(gIs2018) year = 2018;
  if(gIs2017) year = 2017;
  if(gIs2016) year = 2016;

  gIsFSRUp = false; gIsFSRDown = false;
  if     (gSampleName.Contains("TTbar_Powheg") && gSampleName.Contains("fsrUp"))   gIsFSRUp = true;
  else if(gSampleName.Contains("TTbar_Powheg") && gSampleName.Contains("fsrDown")) gIsFSRDown = true;

  //---- Select your wp for b-tagging and pt, eta for the jets
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  if(gSelection == itt){
    taggerName="DeepFlav";
    stringWP = "Medium";
    jet_MaxEta = 2.4;
    jet_MinPt  = 30;
    vetoJet_minPt = 30;
    vetoJet_maxEta = 5.0;
    minDR = 0.4;
  }
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  MeasType = "mujets";
  TString pwd  = GetParam<TString>("WorkingDir");
  TString BTagSFPath = Form("%s/packages/BTagSFUtil", pwd.Data());
  
  TString FastSimDataset = "";
  if(gIsFastSim) FastSimDataset = "1";
  fBTagSFnom = new BTagSFUtil(MeasType, BTagSFPath, taggerName.Data(), stringWP,  0, year);
  fBTagSFbUp = new BTagSFUtil(MeasType, BTagSFPath, taggerName.Data(), stringWP,  1, year);
  fBTagSFbDo = new BTagSFUtil(MeasType, BTagSFPath, taggerName.Data(), stringWP, -1, year);
  fBTagSFlUp = new BTagSFUtil(MeasType, BTagSFPath, taggerName.Data(), stringWP,  3, year);
  fBTagSFlDo = new BTagSFUtil(MeasType, BTagSFPath, taggerName.Data(), stringWP, -3, year);

  Leptons  = std::vector<Lepton>();
  selJets  = std::vector<Jet>();
  mcJets   = std::vector<Jet>();
  genJets  = std::vector<Jet>();
  vetoJets = std::vector<Jet>();
  selJetsJecUp = std::vector<Jet>();
  selJetsJER   = std::vector<Jet>();
  selJetsJecDown = std::vector<Jet>();
}

void JetSelector::GetJetVariables(Int_t i){
  //tpJ.SetPxPyPzE(Get<Float_t>("Jet"+jec+"_px",i), Get<Float_t>("Jet"+jec+"_py",i), Get<Float_t>("Jet"+jec+"_pz", i), Get<Float_t>("Jet"+jec+"_energy",i));
  Float_t FSRSF = 1;
  //if(gIsFSRUp)   FSRSF = GetFSR_JECSF_Up(  Get<Float_t>("Jet"+jec+"_pt",i));
  //if(gIsFSRDown) FSRSF = GetFSR_JECSF_Down(Get<Float_t>("Jet"+jec+"_pt",i));
  tpJ.SetPtEtaPhiM(Get<Float_t>("Jet_pt",i), Get<Float_t>("Jet_eta",i), Get<Float_t>("Jet_phi", i), Get<Float_t>("Jet_mass",i));
  eta = tpJ.Eta();;
  pt = tpJ.Pt();
  //rawPt       = Get<Float_t>("Jet_rawPt",i);
  //pt_corrUp   = Get<Float_t>("Jet_corr_JECUp",i); 
  //pt_corrDown = Get<Float_t>("Jet_corr_JECDown",i);
  jetId       = Get<Int_t>("Jet_jetId",i); 
  csv         = Get<Float_t>("Jet_btagCSVV2", i);
  deepcsv     = Get<Float_t>("Jet_btagDeepB", i);
  deepcsvC    = Get<Float_t>("Jet_btagDeepC", i);
  deepflav    = Get<Float_t>("Jet_btagDeepFlavB", i);
  flavmc = -999999;
  if(!gIsData){
    flavmc = Get<Int_t>("Jet_hadronFlavour", i);
    //tmcJ.SetPxPyPzE(Get<Float_t>("Jet_mcPx",i), Get<Float_t>("Jet_mcPy",i), Get<Float_t>("Jet_mcPz",i), Get<Float_t>("Jet_mcEnergy",i));
  }
}

void JetSelector::GetGenJetVariables(Int_t i){
  tpJ.SetPtEtaPhiM(Get<Float_t>("GenJet_pt",i), Get<Float_t>("GenJet_eta",i), Get<Float_t>("GenJet_phi", i), Get<Float_t>("GenJet_mass",i));
  eta = Get<Float_t>("GenJet_eta",i);
  pt =  Get<Float_t>("GenJet_pt",i);
  flavmc = Get<Float_t>("GenJet_hadronFlavour", i);
}

void JetSelector::GetmcJetVariables(Int_t i){
  tpJ.SetPxPyPzE(Get<Float_t>("Jet_mcPx",i), Get<Float_t>("Jet_mcPy",i), Get<Float_t>("Jet_mcPz", i), Get<Float_t>("Jet_mcEnergy",i));
  eta = tpJ.Eta();
  pt =  Get<Float_t>("Jet_mcPt",i);
  flavmc = TMath::Abs(Get<Int_t>("Jet_mcFlavour",i));
  csv         = Get<Float_t>("Jet_btagCSV", i);
}

void JetSelector::InsideLoop(){
  // Clear vectors...
  selJets.clear();
  selJetsJecUp.clear();
  selJetsJecDown.clear();
  selJetsJER.clear();
  mcJets.clear();
  genJets.clear();
  vetoJets.clear();
  Jets15.clear();
  nBtagJets = 0;
  nBtagJetsJECUp = 0;
  nBtagJetsJECDown = 0;
  Leptons.clear();
  VetoLeptons.clear();

  Leptons = GetParam<vector<Lepton>>("selLeptons"); 
  VetoLeptons = GetParam<vector<Lepton>>("vetoLeptons");

  //evt = (UInt_t)Get<ULong64_t>("evt");
  //rho = Get<Float_t>("rho");

  // Loop over the jets
  nJet = Get<Int_t>("nJet");

  Float_t mcTag     = 1.;
  Float_t dataTag   = 1.;
  Float_t mcNoTag   = 1.;
  Float_t dataNoTag = 1.;
  Float_t errHup  = 0.;
  Float_t errHdn  = 0.;
  Float_t errLup  = 0.;
  Float_t errLdn  = 0.;

  
  // Loop over all the jets
  for(Int_t i = 0; i < nJet; i++){

    // Get all jet variables
    GetJetVariables(i);
    tJ = Jet(tpJ, csv, jetId, flavmc);
    tJ.SetDeepCSVB(deepcsv);
    tJ.SetDeepCSVC(deepcsvC);
    tJ.SetDeepFlav(deepflav);
    tJ.isBtag = IsBtag(tJ);

    // Check and clean
    if(tJ.id > 1 && Cleaning(tJ, Leptons, minDR)){ // Jet id == 2, 6

      // MC info
      if(!gIsData){
        //GetmcJetVariables(i);
        //tJ.SetMCjet(tpJ);
        GetJetVariables(i);
      }

      // Get systematics
      if(gDoSys) SetSystematics(&tJ);

      // Fill the selected jets 
      tJ.isBtag = IsBtag(tJ);
      if (TMath::Abs(tJ.p.Eta()) < jet_MaxEta){
        if(tJ.p.Pt() > 15 || tJ.pTJESUp > 15 || tJ.pTJESDown > 15 || tJ.pTJERUp > 15 ) Jets15.push_back(tJ);
        if(tJ.p.Pt() > jet_MinPt){
          selJets.push_back(tJ);
          if(tJ.isBtag) nBtagJets++; 
        } 
      }

      // Fill the veto collections
      if (tJ.p.Pt() > vetoJet_minPt && TMath::Abs(tJ.p.Eta()) < vetoJet_maxEta){
        if (gSelection == itt){
          if(TMath::Abs(tJ.p.Eta()) > 2.4)  vetoJets.push_back(tJ);
        }
        else  vetoJets.push_back(tJ);
      }
    }
  }

  // Loop over Gen and MC jets...
  if (!gIsData) {
    ngenJet = Get<Int_t>("nGenJet");
    for(Int_t i = 0; i < ngenJet; i++){
      GetGenJetVariables(i);
      tJ = Jet(tpJ, 0, 1, flavmc);
      if (gSelection == itWtt) {
        if (tJ.p.Pt() > 20) {
          genJets.push_back(tJ);
        }
      }
      else genJets.push_back(tJ);
    }

    /*nJet = Get<Int_t>("nJet");
      for(Int_t i = 0; i < nJet; i++){
      GetmcJetVariables(i);
      tJ = Jet(tpJ, csv, 1, flavmc);
      tJ.isBtag = (flavmc == 5);
      mcJets.push_back(tJ);
      }*/
  }

  selJets        = SortJetsByPt(selJets);
  selJetsJecUp   = SortJetsByPt(selJetsJecUp);
  selJetsJecDown = SortJetsByPt(selJetsJecDown);
  selJetsJER     = SortJetsByPt(selJetsJER);
  Jets15         = SortJetsByPt(Jets15);
  mcJets         = SortJetsByPt(mcJets);
  vetoJets       = SortJetsByPt(vetoJets);

  nSelJets  = selJets.size();
  nJets15   = Jets15.size();
  nVetoJets = vetoJets.size();
  nGenJets  = genJets.size();

  // Set params...
  SetParam("selJets",  selJets);
  SetParam("selJetsJecUp", selJetsJecUp);
  SetParam("selJetsJecDown", selJetsJecDown);
  SetParam("selJetsJER" , selJetsJER);
  SetParam("Jets15",   Jets15);
  SetParam("vetoJets", vetoJets);
  SetParam("genJets",  genJets);
  SetParam("mcJets",  mcJets);
  SetParam("nSelJets",  nSelJets);
  SetParam("nJets15",  nJets15);
  SetParam("nVetoJets",  nVetoJets);
  SetParam("nGenJets",  nGenJets);
  SetParam("nSelBJets",  nBtagJets);

  // Propagate JES to MET
  Float_t met_pt  = Get<Float_t>("MET_pt");
  Float_t met_phi = Get<Float_t>("MET_phi");
  MET_JERUp   = met_pt;
  MET_JESUp   = met_pt;
  MET_JESDown = met_pt;
  if(nSelJets > 0 && gDoSys){
    MET_JERUp   = JERtoMET(selJets, met_pt, met_phi);
    MET_JESUp   = JEStoMET(selJets, met_pt, met_phi,  1);
    MET_JESDown = JEStoMET(selJets, met_pt, met_phi, -1);
  }
  SetParam("MET_JERUp",   MET_JERUp);
  SetParam("MET_JESUp",   MET_JESUp);
  SetParam("MET_JESDown", MET_JESDown);
}

Bool_t JetSelector::IsBtag(Jet j){
  if(j.Pt() < 20) return false;
  Bool_t isbtag;
  // using "weights" as scale factors in the tW analysis :)
  if(gSelection == itWtt) isbtag = fBTagSFnom->IsTagged(j.csv, -999999, j.p.Pt(), j.p.Eta(), (UInt_t)j.p.Pt());

  else                         isbtag = fBTagSFnom->IsTagged(j.csv,j.flavmc, j.p.Pt(), j.p.Eta(), (UInt_t)j.p.Pt());
  if(gIsFastSim && BtagSFFS == 1. && isbtag){
    BtagSFFS = fBTagSFnom->GetFastSimBtagSF(j.flavmc, j.p.Eta(), j.p.Pt(), j.csv);
  }
  return isbtag;
}

void JetSelector::SetSystematics(Jet *j){
  Float_t _csv = j->csv; Int_t _flavmc = j->flavmc; Float_t _pt = j->p.Pt(); Float_t _eta = j->p.Eta();
  if(gIsData) return;
  j->pTJESUp     = rawPt*pt_corrUp;
  j->pTJESDown   = rawPt*pt_corrDown;
  //  j->pTJERUp     = getJetJERpt(*j); // rho
  //  j->pTJERDown   = getJetJERpt(*j); // rho
  j->isBtag_BtagUp      = fBTagSFbUp->IsTagged(_csv, _flavmc, _pt, _eta, (UInt_t)_pt+1);
  j->isBtag_BtagDown    = fBTagSFbDo->IsTagged(_csv, _flavmc, _pt, _eta, (UInt_t)_pt-1);
  j->isBtag_MisTagUp    = fBTagSFlUp->IsTagged(_csv, _flavmc, _pt, _eta, (UInt_t)_pt+3);
  j->isBtag_MisTagDown  = fBTagSFlDo->IsTagged(_csv, _flavmc, _pt, _eta, (UInt_t)_pt-3);
}

