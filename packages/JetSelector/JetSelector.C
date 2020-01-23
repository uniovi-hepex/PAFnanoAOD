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
  MeasType   = "comb";
  minDR      = 0;
  jet_MaxEta = 0;
  jet_MinPt  = 0;
  vetoJet_minPt = 0;
  BtagSFFS   = 1;
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
  year        = GetParam<TString>("year").Atoi();
  gIsData     = GetParam<Bool_t>("IsData");
  selection   = GetParam<TString>("selection");
  gSampleName  = GetParam<TString>("sampleName");
  gOptions     = GetParam<TString>("_options");
  gDoSys       = gOptions.Contains("doSyst")? true : false;
  gSelection   = GetSelection(selection);

  gIsFSRUp = false; gIsFSRDown = false;
  if     (gSampleName.Contains("TTbar_Powheg") && gSampleName.Contains("fsrUp"))   gIsFSRUp = true;
  else if(gSampleName.Contains("TTbar_Powheg") && gSampleName.Contains("fsrDown")) gIsFSRDown = true;

  //---- Select your wp for b-tagging and pt, eta for the jets
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  taggerName="DeepFlav";//"DeepFlav";
  stringWP = "Medium";
  jet_MaxEta = 2.4;
  jet_MinPt  = 30;
  vetoJet_minPt = 30;
  vetoJet_maxEta = 5.0;
  minDR = 0.4;
  if      (gSelection == itt) {
    taggerName="DeepFlav";//"DeepFlav";
    stringWP = "Medium";
    jet_MaxEta = 2.4;
    jet_MinPt  = 30;
    vetoJet_minPt = 30;
    vetoJet_maxEta = 5.0;
    minDR = 0.4;
  }
  else if (gSelection == itWtt) {
    taggerName      = "DeepFlav";
    stringWP        = "Medium";
    jet_MaxEta      = 2.4;
    jet_MinPt       = 30;
    vetoJet_minPt   = 20.;
    vetoJet_maxEta  = 4.7;
    minDR           = 0.4;
  }
  else if (gSelection == itW) {
    taggerName      = "DeepFlav";
    stringWP        = "Medium";
    jet_MaxEta      = 2.4;
    jet_MinPt       = 30;
    vetoJet_minPt   = 20.;
    vetoJet_maxEta  = 4.7;
    minDR           = 0.4;
  }

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  MeasType = "mujets";
  TString pwd  = GetParam<TString>("WorkingDir");
  TString BTagSFPath = Form("%s/packages/BTagSFUtil", pwd.Data());
  
  TString FastSimDataset = "";
  fBTagSFnom = new BTagSFUtil(MeasType, BTagSFPath, taggerName.Data(), stringWP,  0, year);
  if(!gIsData){
    fBTagSFbUp = new BTagSFUtil(MeasType, BTagSFPath, taggerName.Data(), stringWP,  1, year);
    fBTagSFbDo = new BTagSFUtil(MeasType, BTagSFPath, taggerName.Data(), stringWP, -1, year);
    fBTagSFlUp = new BTagSFUtil(MeasType, BTagSFPath, taggerName.Data(), stringWP,  3, year);
    fBTagSFlDo = new BTagSFUtil(MeasType, BTagSFPath, taggerName.Data(), stringWP, -3, year);
  }

  Leptons  = std::vector<Lepton>();
  selJets  = std::vector<Jet>();
  mcJets   = std::vector<Jet>();
  genJets  = std::vector<Jet>();
  vetoJets = std::vector<Jet>();
  selJetsJecCorUp = std::vector<Jet>();
  selJetsJecCorDown = std::vector<Jet>();
  selJetsJecUnCorUp = std::vector<Jet>();
  selJetsJecUnCorDown = std::vector<Jet>();
  selJetsJecUp = std::vector<Jet>();
  selJetsJecDown = std::vector<Jet>();
  selJetsJERUp   = std::vector<Jet>();
  selJetsJERDown   = std::vector<Jet>();
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
  flav = -999999; if(!gIsData) flav = Get<Int_t>("Jet_hadronFlavour", i);
}

void JetSelector::GetGenJetVariables(Int_t i){
  tpJ.SetPtEtaPhiM(Get<Float_t>("GenJet_pt",i), Get<Float_t>("GenJet_eta",i), Get<Float_t>("GenJet_phi", i), Get<Float_t>("GenJet_mass",i));
  eta = Get<Float_t>("GenJet_eta",i);
  pt =  Get<Float_t>("GenJet_pt",i);
  flavmc = Get<Float_t>("GenJet_hadronFlavour", i);
}

void JetSelector::InsideLoop(){
  // Clear vectors...
  selJets.clear();
  selJetsJecCorUp.clear();
  selJetsJecCorDown.clear();
  selJetsJecUnCorUp.clear();
  selJetsJecUnCorDown.clear();
  selJetsJecUp.clear();
  selJetsJecDown.clear();
  selJetsJERUp.clear();
  selJetsJERDown.clear();
  mcJets.clear();
  genJets.clear();
  vetoJets.clear();
  Jets15.clear();
  nBtagJets = 0;
  nBtagJetsJECUp = 0;
  nBtagJetsJECDown = 0;
  Leptons.clear();
  VetoLeptons.clear();
  
  BtagSFFS         = 1.;
  BtagSF           = 1.;
  BtagSFBtagUp     = 1.;
  BtagSFBtagDown   = 1.;
  BtagSFMistagUp   = 1.;
  BtagSFMistagDown = 1.;
  
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
  Float_t errHup    = 0.;
  Float_t errHdn    = 0.;
  Float_t errLup    = 0.;
  Float_t errLdn    = 0.;
  TLorentzVector pJESUp; TLorentzVector pJESDo; TLorentzVector pJERUp; TLorentzVector pJERDo;
  TLorentzVector pJESCorUp; TLorentzVector pJESCorDo;
  TLorentzVector pJESUnCorUp; TLorentzVector pJESUnCorDo;
  Jet tJESUp; Jet tJESDo; Jet tJERUp; Jet tJERDo; Jet tJESCorUp; Jet tJESUnCorUp; Jet tJESCorDo; Jet tJESUnCorDo;

  // Loop over all the jets
  for (Int_t i = 0; i < nJet; i++) {
    // Get all jet variables
    GetJetVariables(i);
    tJ = Jet(tpJ, csv, jetId, flav);
    tJ.SetDeepCSVB(deepcsv);
    tJ.SetDeepCSVC(deepcsvC);
    tJ.SetDeepFlav(deepflav);
    //tJ.isBtag = IsBtag(tJ);

    // Check and clean
    if(tJ.id > 1 && Cleaning(tJ, Leptons, minDR)){ // Jet id == 2, 6

      // MC info
      //if(!gIsData){
      //  GetJetVariables(i);
      //}


      // Fill the selected jets
      tJ.SetIsBtag(fBTagSFnom->IsTagged(deepflav, flav, tJ.p.Pt(), tJ.p.Eta(), (UInt_t)tJ.p.Pt()),  0);
      if(!gIsData){
        tJ.SetIsBtag(fBTagSFbUp->IsTagged(deepflav, flav, tJ.p.Pt(), tJ.p.Eta(), (UInt_t)tJ.p.Pt()),  1);
        tJ.SetIsBtag(fBTagSFbDo->IsTagged(deepflav, flav, tJ.p.Pt(), tJ.p.Eta(), (UInt_t)tJ.p.Pt()), -1);
        tJ.SetIsBtag(fBTagSFlUp->IsTagged(deepflav, flav, tJ.p.Pt(), tJ.p.Eta(), (UInt_t)tJ.p.Pt()),  2);
        tJ.SetIsBtag(fBTagSFlDo->IsTagged(deepflav, flav, tJ.p.Pt(), tJ.p.Eta(), (UInt_t)tJ.p.Pt()), -2);
      }
      if (TMath::Abs(tJ.p.Eta()) < jet_MaxEta){
        if(tJ.p.Pt() > 15 || tJ.pTJESUp > 15 || tJ.pTJESDown > 15 || tJ.pTJERUp > 15 ) Jets15.push_back(tJ);
        if(tJ.p.Pt() > jet_MinPt){
          selJets.push_back(tJ);
          if(tJ.isBtag) nBtagJets++;
        }
      }

      // Get systematics
      //if(gDoSys) SetSystematics(&tJ);
      if(!gIsData){ // JER, JES...
        pJESUp.SetPtEtaPhiM(Get<Float_t>("Jet_pt_jesTotalUp",   i), tJ.p.Eta(), tJ.p.Phi(), Get<Float_t>("Jet_mass_jesTotalUp",   i));
        pJESDo.SetPtEtaPhiM(Get<Float_t>("Jet_pt_jesTotalDown", i), tJ.p.Eta(), tJ.p.Phi(), Get<Float_t>("Jet_mass_jesTotalDown", i));
        pJERUp.SetPtEtaPhiM(Get<Float_t>("Jet_pt_jerUp",   i), tJ.p.Eta(), tJ.p.Phi(), Get<Float_t>("Jet_mass_jerUp",   i));
        pJERDo.SetPtEtaPhiM(Get<Float_t>("Jet_pt_jerDown", i), tJ.p.Eta(), tJ.p.Phi(), Get<Float_t>("Jet_mass_jerDown", i));
        pJESUnCorUp.SetPtEtaPhiM(Get<Float_t>("Jet_pt_jesTotalUnCorrUp",   i), tJ.p.Eta(), tJ.p.Phi(), Get<Float_t>("Jet_mass_jesTotalUnCorrUp",   i));
        pJESUnCorDo.SetPtEtaPhiM(Get<Float_t>("Jet_pt_jesTotalUnCorrDown", i), tJ.p.Eta(), tJ.p.Phi(), Get<Float_t>("Jet_mass_jesTotalUnCorrDown", i));
        pJESCorUp.SetPtEtaPhiM(Get<Float_t>("Jet_pt_jesTotalCorrUp",   i), tJ.p.Eta(), tJ.p.Phi(), Get<Float_t>("Jet_mass_jesTotalCorrUp",   i));
        pJESCorDo.SetPtEtaPhiM(Get<Float_t>("Jet_pt_jesTotalCorrDown", i), tJ.p.Eta(), tJ.p.Phi(), Get<Float_t>("Jet_mass_jesTotalCorrDown", i));
        tJESUp = Jet(pJESUp, csv, jetId, flav);
        tJESDo = Jet(pJESDo, csv, jetId, flav);
        tJERUp = Jet(pJERUp, csv, jetId, flav);
        tJERDo = Jet(pJERDo, csv, jetId, flav);
        tJESCorUp = Jet(pJESCorUp, csv, jetId, flav);
        tJESCorDo = Jet(pJESCorDo, csv, jetId, flav);
        tJESUnCorUp = Jet(pJESUnCorUp, csv, jetId, flav);
        tJESUnCorDo = Jet(pJESUnCorDo, csv, jetId, flav);
        tJESUp.SetIsBtag( fBTagSFnom->IsTagged(deepflav, flav, Get<Float_t>("Jet_pt_jesTotalUp",   i), tJ.p.Eta(), (UInt_t)Get<Float_t>("Jet_pt_jesTotalUp",   i)) );
        tJESDo.SetIsBtag( fBTagSFnom->IsTagged(deepflav, flav, Get<Float_t>("Jet_pt_jesTotalDown", i), tJ.p.Eta(), (UInt_t)Get<Float_t>("Jet_pt_jesTotalDown", i)) );
        tJERUp.SetIsBtag( fBTagSFnom->IsTagged(deepflav, flav, Get<Float_t>("Jet_pt_jerUp",        i), tJ.p.Eta(), (UInt_t)Get<Float_t>("Jet_pt_jerUp",        i)) );
        tJERDo.SetIsBtag( fBTagSFnom->IsTagged(deepflav, flav, Get<Float_t>("Jet_pt_jerDown",      i), tJ.p.Eta(), (UInt_t)Get<Float_t>("Jet_pt_jerDown",      i)) );
        tJESCorUp.SetIsBtag( fBTagSFnom->IsTagged(deepflav, flav, Get<Float_t>("Jet_pt_jesTotalCorrUp",   i), tJ.p.Eta(), (UInt_t)Get<Float_t>("Jet_pt_jesTotalCorrUp",   i)) );
        tJESCorDo.SetIsBtag( fBTagSFnom->IsTagged(deepflav, flav, Get<Float_t>("Jet_pt_jesTotalCorrDown", i), tJ.p.Eta(), (UInt_t)Get<Float_t>("Jet_pt_jesTotalCorrDown", i)) );
        tJESUnCorUp.SetIsBtag( fBTagSFnom->IsTagged(deepflav, flav, Get<Float_t>("Jet_pt_jesTotalUnCorrUp",   i), tJ.p.Eta(), (UInt_t)Get<Float_t>("Jet_pt_jesTotalUnCorrUp",   i)) );
        tJESUnCorDo.SetIsBtag( fBTagSFnom->IsTagged(deepflav, flav, Get<Float_t>("Jet_pt_jesTotalUnCorrDown", i), tJ.p.Eta(), (UInt_t)Get<Float_t>("Jet_pt_jesTotalUnCorrDown", i)) );
        if(TMath::Abs(tJESUp.Eta()) < jet_MaxEta && tJESUp.Pt() > jet_MinPt) selJetsJecUp  .push_back(tJESUp);
        if(TMath::Abs(tJESDo.Eta()) < jet_MaxEta && tJESDo.Pt() > jet_MinPt) selJetsJecDown.push_back(tJESDo);
        if(TMath::Abs(tJERUp.Eta()) < jet_MaxEta && tJERUp.Pt() > jet_MinPt) selJetsJERUp  .push_back(tJERUp);
        if(TMath::Abs(tJERDo.Eta()) < jet_MaxEta && tJERDo.Pt() > jet_MinPt) selJetsJERDown.push_back(tJERDo);
        if(TMath::Abs(tJESCorUp.Eta()) < jet_MaxEta && tJESCorUp.Pt() > jet_MinPt) selJetsJecCorUp  .push_back(tJESCorUp);
        if(TMath::Abs(tJESCorDo.Eta()) < jet_MaxEta && tJESCorDo.Pt() > jet_MinPt) selJetsJecCorDown.push_back(tJESCorDo);
        if(TMath::Abs(tJESUnCorUp.Eta()) < jet_MaxEta && tJESUnCorUp.Pt() > jet_MinPt) selJetsJecUnCorUp  .push_back(tJESUnCorUp);
        if(TMath::Abs(tJESUnCorDo.Eta()) < jet_MaxEta && tJESUnCorDo.Pt() > jet_MinPt) selJetsJecUnCorDown.push_back(tJESUnCorDo);
      }

      // Fill the veto collections
      if (tJ.p.Pt() > vetoJet_minPt && TMath::Abs(tJ.p.Eta()) < vetoJet_maxEta) {
        if (gSelection == itt) {
          if(TMath::Abs(tJ.p.Eta()) > 2.4)  vetoJets.push_back(tJ);
        }
        else if (gSelection == itWtt || gSelection == itW) {
          vetoJets.push_back(tJ);
          
          if (!gIsData) {
            if (tJ.p.Pt() <= 20.001 || TMath::Abs(tJ.p.Eta()) >= 2.4) continue;
            Float_t cutval = -99;

            if      (taggerName == "CSVv2")    cutval = tJ.csv;
            else if (taggerName == "DeepCSV")  cutval = tJ.GetDeepCSVB();
            else if (taggerName == "DeepFlav") cutval = tJ.GetDeepFlav();

            Float_t eff    = fBTagSFnom->JetTagEfficiency(tJ.flavmc, tJ.p.Pt(), tJ.p.Eta());
            Float_t sf     = fBTagSFnom->GetJetSF(cutval, tJ.flavmc, tJ.p.Pt(), tJ.p.Eta());
            Float_t sfHUp  = fBTagSFbUp->GetJetSF(cutval, tJ.flavmc, tJ.p.Pt(), tJ.p.Eta());
            Float_t sfHDn  = fBTagSFbDo->GetJetSF(cutval, tJ.flavmc, tJ.p.Pt(), tJ.p.Eta());
            Float_t sfLUp  = fBTagSFlUp->GetJetSF(cutval, tJ.flavmc, tJ.p.Pt(), tJ.p.Eta());
            Float_t sfLDn  = fBTagSFlDo->GetJetSF(cutval, tJ.flavmc, tJ.p.Pt(), tJ.p.Eta());

            if (tJ.isBtag) {
              mcTag   *= eff;
              dataTag *= eff * sf;

              if (tJ.flavmc == 5 || tJ.flavmc == 4) {
                errHup += (sfHUp - sf ) / sf;
                errHdn += (sf - sfHDn ) / sf;
              }
              else {
                errLup += (sfLUp - sf ) / sf;
                errLdn += (sf - sfLDn ) / sf;
              }
            }
            else {
              mcNoTag   *= ( 1 - eff    );
              dataNoTag *= ( 1 - eff*sf );
              if (tJ.flavmc == 5 || tJ.flavmc == 4) {
                errHup -= eff*(sfHUp - sf ) / (1 - eff*sf);
                errHdn -= eff*(sf - sfHDn ) / (1 - eff*sf);
              }
              else {
                errLup -= eff*(sfLUp - sf ) / (1 - eff*sf);
                errLdn -= eff*(sf - sfLDn ) / (1 - eff*sf);
              }
            }
          }
        }
        else  vetoJets.push_back(tJ);
      }
    }
  }

  BtagSF           *= (dataNoTag * dataTag) / (mcNoTag * mcTag);
  BtagSFBtagUp     *= BtagSF * ( 1 + errHup );
  BtagSFBtagDown   *= BtagSF * ( 1 - errHdn );
  BtagSFMistagUp   *= BtagSF * ( 1 + errLup );
  BtagSFMistagDown *= BtagSF * ( 1 - errLdn );
  
  // Loop over Gen and MC jets...
  if (!gIsData) {
    ngenJet = Get<Int_t>("nGenJet");
    for(Int_t i = 0; i < ngenJet; i++){
      GetGenJetVariables(i);
      tJ = Jet(tpJ, 0, 1, flavmc);
      if (gSelection == itWtt || gSelection == itW) {
        if (tJ.p.Pt() > 20) {
          genJets.push_back(tJ);
        }
      }
      else genJets.push_back(tJ);
    }
  }

  selJets        = SortJetsByPt(selJets);
  selJetsJecCorUp   = SortJetsByPt(selJetsJecCorUp);
  selJetsJecCorDown = SortJetsByPt(selJetsJecCorDown);
  selJetsJecUnCorUp   = SortJetsByPt(selJetsJecUnCorUp);
  selJetsJecUnCorDown   = SortJetsByPt(selJetsJecUnCorDown);
  selJetsJecUp   = SortJetsByPt(selJetsJecUp);
  selJetsJecDown = SortJetsByPt(selJetsJecDown);
  selJetsJERUp   = SortJetsByPt(selJetsJERUp);
  selJetsJERDown = SortJetsByPt(selJetsJERDown);
  Jets15         = SortJetsByPt(Jets15);
  mcJets         = SortJetsByPt(mcJets);
  vetoJets       = SortJetsByPt(vetoJets);

  nSelJets  = selJets.size();
  nJets15   = Jets15.size();
  nVetoJets = vetoJets.size();
  nGenJets  = genJets.size();

  // Set params...
  SetParam("selJets",         selJets);
  SetParam("selJetsJecCorUp",    selJetsJecCorUp);
  SetParam("selJetsJecCorDown",    selJetsJecCorDown);
  SetParam("selJetsJecUnCorUp",    selJetsJecUnCorUp);
  SetParam("selJetsJecUnCorDown",    selJetsJecUnCorDown);
  SetParam("selJetsJecUp",    selJetsJecUp);
  SetParam("selJetsJecDown",  selJetsJecDown);
  SetParam("selJetsJERUp" ,   selJetsJERUp);
  SetParam("selJetsJERDown",  selJetsJERDown);
  SetParam("Jets15",          Jets15);
  SetParam("vetoJets",        vetoJets);
  SetParam("genJets",         genJets);
  SetParam("mcJets",          mcJets);
  SetParam("nSelJets",        nSelJets);
  SetParam("nJets15",         nJets15);
  SetParam("nVetoJets",       nVetoJets);
  SetParam("nGenJets",        nGenJets);
  SetParam("nSelBJets",       nBtagJets);
  SetParam("BtagSF",          BtagSF);
  SetParam("BtagSFBtagUp",    BtagSFBtagUp);
  SetParam("BtagSFBtagDown",  BtagSFBtagDown);
  SetParam("BtagSFMistagUp",  BtagSFMistagUp);
  SetParam("BtagSFMistagDown",BtagSFMistagDown);
  SetParam("BtagSFFS",        BtagSFFS);
}

Bool_t JetSelector::IsBtag(Jet j){
  if(j.Pt() < 20) return false;
  Bool_t isbtag = false;
  // using "weights" as scale factors in the tW analysis :)
  if (gSelection == itWtt || gSelection == itW) isbtag = fBTagSFnom->IsTagged(j.csv, -999999, j.p.Pt(), j.p.Eta(), (UInt_t)j.p.Pt());
  else   isbtag = fBTagSFnom->IsTagged(j.csv,j.flavmc, j.p.Pt(), j.p.Eta(), (UInt_t)j.p.Pt());
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

