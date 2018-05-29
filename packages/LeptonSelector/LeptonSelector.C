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


#include "LeptonSelector.h"

ClassImp(LeptonSelector);
LeptonSelector::LeptonSelector() : PAFChainItemSelector() {}
void LeptonSelector::Summary(){}

void LeptonSelector::Initialise(){
  // Initialise LeptonSelector
  gIsData        = GetParam<Bool_t>("IsData");
  gIsFastSim     = GetParam<Bool_t>("IsFastSim");
  gSelection     = GetParam<Int_t>("iSelection");
  gOptions       = GetParam<TString>("_options");
  gIs2017 = false;
  if(gOptions.Contains("2017")) gIs2017 = true;
  localPath      = GetParam<TString>("WorkingDir");
  LepSF = new LeptonSF(localPath + "/InputFiles/", gOptions);
  ElecScale = new ElecScaleClass(localPath + "/InputFiles/ElecScale.dat");

  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  // Select SFs for each analysis !!!!!!
  if(gSelection == i4tSelec || gSelection == iWZSelec){
    // Los tenemos en funciones, no en histogramas, entonces nos la pela de cargar histogramas
  }

  else if(gSelection == iTopSelec || gSelection == iTWSelec || gSelection == iWWSelec){
    LepSF->loadHisto(iTrigDoubleMuon);
    LepSF->loadHisto(iTrigDoubleElec);
    LepSF->loadHisto(iTrigElMu);
    LepSF->loadHisto(iMuonReco);
    LepSF->loadHisto(iElecReco);
    LepSF->loadHisto(iMuonId,   iTight);
    LepSF->loadHisto(iMuonIsoTightId,   iTight);
    LepSF->loadHisto(iElecId,   iTight);
  }
  else if(gSelection == iStopTopSelec){
    LepSF->loadHisto(iTrigDoubleMuon);
    LepSF->loadHisto(iTrigDoubleElec);
    LepSF->loadHisto(iTrigElMu);
    LepSF->loadHisto(iMuonReco);
    LepSF->loadHisto(iElecReco);
    LepSF->loadHisto(iMuonId,   iTight);
    LepSF->loadHisto(iMuonIsoTightId,   iTight);
    LepSF->loadHisto(iElecId,   iTight);
    if(gIsFastSim){
      LepSF->loadHisto(iElecFastSim);
      LepSF->loadHisto(iMuonFastSim);
    }
  }
  else std::cout << ">>>>>>>>>>>> WRONG SELECTION!!!!!!!!" << std::endl;
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  selLeptons   = std::vector<Lepton>();
  vetoLeptons   = std::vector<Lepton>();
  looseLeptons   = std::vector<Lepton>();
}

//################################################################
//## Definition of wps...
//################################################################
Bool_t LeptonSelector::getSIPcut(Float_t cut){
  if(sip > cut) return false;
  return true;
}

Bool_t LeptonSelector::getGoodVertex(Int_t wp){
  if (gSelection == i4tSelec || gSelection == iWZSelec) {
    if (wp == iTight || wp == iMedium || wp == iLoose){
      if (dxy >= 0.05 || dz >= 0.1) return false;
    }
  } else {
    if(type == 1){ //electrons
      if(wp == iTight){
        if(etaSC <= 1.479 && ((dxy >= 0.05) || (dz  >= 0.10))) return false;
        if(etaSC >  1.479 && ((dxy >= 0.10) || (dz  >= 0.20))) return false;
      }
    }
    else{ // muons
      if(wp == iMedium && (dxy > 0.2  || dz > 0.5)) return false;
      if(wp == iTight  && (dxy > 0.05 || dz > 0.1)) return false;
    }
  }
  return true;
}

Bool_t LeptonSelector::getRelIso03POG(Int_t wp){
  if(type == 1){ // electrons
    if(wp == iWPforStop && RelIso03 > 0.12) return false;
    if(etaSC <= 1.479){
      if(wp == iVeto   && RelIso03 > 0.1750) return false;
      if(wp == iLoose  && RelIso03 > 0.0994) return false;
      if(wp == iMedium && RelIso03 > 0.0695) return false;
      if(wp == iTight  && RelIso03 > 0.0588) return false;
    }
    else if(etaSC > 1.479){
      if(wp == iVeto   && RelIso03 > 0.1590) return false;
      if(wp == iLoose  && RelIso03 > 0.1070) return false;
      if(wp == iMedium && RelIso03 > 0.0821) return false;
      if(wp == iTight  && RelIso03 > 0.0571) return false;
    }
  }
  else{ // muons
    if(wp == iLoose  && RelIso03 > 0.10) return false; 
    if(wp == iTight  && RelIso03 > 0.05) return false; 
    if(wp == iWPforStop && RelIso03 > 0.12) return false; 
    if(wp == iLooseWPforStop && RelIso03 > 0.4) return false; 
  }
  return true;
}

Bool_t LeptonSelector::getRelIso04POG(Int_t wp){ // wps for muons
  if(type == 1) return false; // electrons
  if(wp == iLoose  && RelIso04 > 0.25) return false;
  if(wp == iTight  && RelIso04 > 0.15) return false;
  return true;
}

Bool_t LeptonSelector::getminiRelIso(Int_t wp) {
  if (wp == iTight || wp == iMedium || wp == iLoose) {
    if (miniIso > 0.4) return false;
  }
  return true;
}

Bool_t LeptonSelector::getMuonId(Int_t wp){
  if(wp == iTight   && !tightVar)     return false;
  if(wp == iMedium  && !mediumMuonId) return false;
  return true;
}

Bool_t LeptonSelector::getElecMVA(Int_t wp){
  Float_t point = 0; Float_t aeta = TMath::Abs(etaSC);
  if(wp == iVeryLoose){
    if(aeta < 0.8){
      if     (pt > 10 && pt < 15) point = -0.86;
      else if(pt > 15 && pt < 25) point = -0.86 + (-0.96+0.86)/10*(pt-15);
      else if(pt > 25           ) point = -0.96;
    }
    else if(aeta < 1.479){
      if     (pt > 10 && pt < 15) point = -0.85;
      else if(pt > 15 && pt < 25) point = -0.85 + (-0.96+0.85)/10*(pt-15);
      else if(pt > 25           ) point = -0.96;
    }
    else if(aeta < 2.5){
      if     (pt > 10 && pt < 15) point = -0.81;
      else if(pt > 15 && pt < 25) point = -0.81 + (-0.95+0.81)/10*(pt-15);
      else if(pt > 25           ) point = -0.95;
    }
  }
  if(wp == iLoose){
    if(aeta < 0.8){
      if     (pt > 10 && pt < 15) point = -0.48;
      else if(pt > 15 && pt < 25) point = -0.48 + (-0.85+0.48)/10*(pt-15);
      else if(pt > 25           ) point = -0.85;
    }
    else if(aeta < 1.479){
      if     (pt > 10 && pt < 15) point = -0.67;
      else if(pt > 15 && pt < 25) point = -0.67 + (-0.91+0.67)/10*(pt-15);
      else if(pt > 25           ) point = -0.91;
    }
    else if(aeta < 2.5){
      if     (pt > 10 && pt < 15) point = -0.49;
      else if(pt > 15 && pt < 25) point = -0.49 + (-0.83+0.49)/10*(pt-15);
      else if(pt > 25           ) point = -0.83;
    }
  }
 
  else if(wp == iTight){
    if(aeta < 0.8){
      if     (pt > 10 && pt < 15) point = 0.77;
      else if(pt > 15 && pt < 25) point = 0.77 + (0.52-0.77)/10*(pt-15);
      else if(pt > 25           ) point = 0.52;
    }
    else if(aeta < 1.479){
      if     (pt > 10 && pt < 15) point = 0.56;
      else if(pt > 15 && pt < 25) point = 0.56 + (0.11-0.56)/10*(pt-15);
      else if(pt > 25           ) point = 0.11;
    }
    else if(aeta < 2.5){
      if     (pt > 10 && pt < 15) point = 0.48;
      else if(pt > 15 && pt < 25) point = 0.48 + (-0.01-0.48)/10*(pt-15);
      else if(pt > 25           ) point = -0.01;
    }
  }
  return ( MVAID > point); 
}

Bool_t LeptonSelector::getElecMVAId(Int_t wp, Lepton lep) {
  if (wp == iTight) {
    Float_t   A = -0.86+(-0.85 + 0.86)*(abs(eta) > 0.8)+(-0.81 + 0.86)*(abs(eta) > 1.479);
    Float_t  B = -0.96+(-0.96 + 0.96)*(abs(eta) > 0.8)+(-0.95 + 0.96)*(abs(eta) > 1.479);
    if (pt > 10) {
      if (!(MVAID > min( A , max( B , A+(B-A)/10*(pt-15))) )) return false;
    }
    if (abs(eta) < 0.8) {
      if (sigmaIEtaIEta > 0.011)             return false;
      if (HoE > 0.10)                       return false;
      if (dEtaSC > 0.01)                     return false;
//      if (dPhiSC > 0.04)                     return false;
      if (eImpI < -0.05 || eImpI > 0.0010)   return false;
    }
    else if ((abs(eta) < 1.479) && (abs(eta) >= 0.8)){
      if (sigmaIEtaIEta > 0.011)             return false;
      if (HoE > 0.10)                       return false;
      if (dEtaSC > 0.01)                     return false;
//      if (dPhiSC > 0.04)                     return false;
      if (eImpI < -0.05 || eImpI > 0.0010)   return false;
    }
    else if (abs(eta) >= 1.479) {
      if (sigmaIEtaIEta > 0.030)             return false;
      if (HoE > 0.07)                        return false;
      if (dEtaSC > 0.008)                   return false;
  //    if (dPhiSC > 0.07)                     return false;
      if (eImpI < -0.05 || eImpI > 0.005)   return false;
    }
    if (convVeto == 0)     return false;
    if (lostHits != 0)     return false;
  }

  if (wp == iMedium) {
    if (isGoodLepton(lep)) {
      Float_t   A = -0.86+(-0.85 + 0.86)*(abs(eta) > 0.8)+(-0.81 + 0.86)*(abs(eta) > 1.479);
      Float_t  B = -0.96+(-0.96 + 0.96)*(abs(eta) > 0.8)+(-0.95 + 0.96)*(abs(eta) > 1.479);
      if (pt > 10) {
        if (!(MVAID > min( A , max( B , A+(B-A)/10*(pt-15))) )) return false;
      }
    }
    if (abs(eta) < 0.8) {
      if (!isGoodLepton(lep)) {
        if (MVAID < 0)                          return false;
      }
      if (sigmaIEtaIEta > 0.011)              return false;
      if (HoE > 0.10)                        return false;
      if (dEtaSC > 0.01)                      return false;
//      if (dPhiSC > 0.04)                      return false;
      if (eImpI < -0.05 || eImpI > 0.0010)    return false;
    }
    else if ((abs(eta) < 1.479) && (abs(eta) >= 0.8)) {
      if (!isGoodLepton(lep)) {
        if (MVAID < 0)                          return false;
      }
      if (sigmaIEtaIEta > 0.011)              return false;
      if (HoE > 0.10)                        return false;
      if (dEtaSC > 0.01)                      return false;
//      if (dPhiSC > 0.04)                      return false;
      if (eImpI < -0.05 || eImpI > 0.0010)    return false;
    }
    else if (abs(eta) >= 1.479) {
      if (!isGoodLepton(lep)) {
        if (MVAID < 0.7)                        return false;
      }
        if (sigmaIEtaIEta > 0.030)            return false;
        if (HoE > 0.07)                       return false;
        if (dEtaSC > 0.008)                  return false;
//        if (dPhiSC > 0.07)                    return false;
        if (eImpI < -0.05 || eImpI > 0.005)  return false;
    }
    if (lostHits != 0)                     return false;
  }

  if (wp == iLoose) {
    Float_t   A = -0.86+(-0.85 + 0.86)*(abs(eta) > 0.8)+(-0.81 + 0.86)*(abs(eta) > 1.479);
    Float_t  B = -0.96+(-0.96 + 0.96)*(abs(eta) > 0.8)+(-0.95 + 0.96)*(abs(eta) > 1.479);
    if (pt > 10) {
      if (!(MVAID > min( A , max( B , A+(B-A)/10*(pt-15))) )) return false;
    }
    if (lostHits > 2)     return false;
  }
  return true;
}

Bool_t LeptonSelector::getElecCutBasedId(Int_t wp){
    if(wp == iTight   && tightVar < 4)     return false;
    if(wp == iMedium  && tightVar < 3)     return false;
    if(wp == iLoose   && tightVar < 2)     return false;
    if(wp == iVeto    && tightVar < 1)     return false;
  return true;
}


Bool_t LeptonSelector::getMultiIso(Int_t wp){
  // Needs ptratio and ptrel...
/*  if(wp == iVeryTight) return (miniIso < 0.09 && (ptRatio > 0.84 || ptRel > 7.2));
  if(wp == iTight    ){ 
    DumpVar(evt, "miniIso", miniIso, miniIso < 0.12);
    DumpVar(evt, "(ptRatio > 0.80 || ptRel > 7.2)", ptRatio, (ptRatio > 0.80 || ptRel > 7.2));
    return (miniIso < 0.12 && (ptRatio > 0.80 || ptRel > 7.2));
  }
  if(wp == iMedium   ){ 
    DumpVar(evt, "miniIso", miniIso, miniIso < 0.16);
    DumpVar(evt, "(ptRatio > 0.76 || ptRel > 7.2)", ptRatio, (ptRatio > 0.76 || ptRel > 7.2));
    return (miniIso < 0.16 && (ptRatio > 0.76 || ptRel > 7.2));
  }
  if(wp == iLoose    ) return miniIso < 0.40;
  //if(wp == iLoose    ) return (miniIso < 0.20 && (ptRatio > 0.69 || ptRel > 6.0));
  if(wp == iVeto     ){
    if(type == 1)      return (miniIso < 0.40 && (ptRatio > 0.80 || ptRel > 7.2));
    else               return (miniIso < 0.40 && (ptRatio > 0.76 || ptRel > 7.2));
  }*/    
  return true;
}


//################################################################
//## Lepton definitions for each analysis
//################################################################
//
// Use the functions above to define your objects
//

//============================================================================================
//=============================================================== SELECTED LEPTONS
//============================================================================================
Bool_t LeptonSelector::isGoodLepton(Lepton lep){
  Bool_t passId; Bool_t passIso;
  if(gSelection == iTopSelec || gSelection == iStopTopSelec || gSelection == iTWSelec){
    // Tight cut-based electrons, pT > 20, |eta| < 2.4, RelIso POG, tightIP2D, SIP3D > 4
    // Tight Muon ID, RelIso POG, tightIP2D, SIP3D > 4
    if(lep.isMuon){
      passId  = getMuonId(iTight);
      passIso = getRelIso04POG(iTight);
    }
    if(lep.isElec){
      passId = getElecCutBasedId(iTight) && lostHits <= 1;
      passIso = getRelIso03POG(iTight);
      if(TMath::Abs(etaSC) > 1.4442 && TMath::Abs(etaSC) < 1.566) return false;
    }
    if(lep.p.Pt() < 18 || TMath::Abs(lep.p.Eta()) > 2.4) return false;
    if(passId && passIso && ( (lep.isElec && getGoodVertex(iTight)) || (lep.isMuon && getGoodVertex(iMedium) ))) return true;
    else return false;
  }
  else if(gSelection == iWWSelec){
    if(lep.isMuon){
      passId  = getMuonId(iTight);
      passIso = getRelIso04POG(iTight);
    }
    if(lep.isElec){
      passId = getElecCutBasedId(iTight);
      passIso = getRelIso03POG(iTight);
      if(TMath::Abs(etaSC) > 1.4442 && TMath::Abs(etaSC) < 1.566) return false;
    }
    if(lep.p.Pt() < 20 || TMath::Abs(lep.p.Eta()) > 2.4) return false;
    if(passId && passIso && getGoodVertex(iTight) && getSIPcut(4)) return true;
    else return false;
    return false;
  }
  else if(gSelection == i4tSelec){
    if(lep.isMuon){
      if(lep.p.Pt() < 20) return false;
      if(TMath::Abs(lep.p.Eta()) > 2.4) return false;
      if(!getMuonId(iMedium)) return false;
      if(!getMultiIso(iMedium)) return false;
    }
    if(lep.isElec){
      if(lep.p.Pt() < 20) return false;
      if(TMath::Abs(lep.p.Eta()) > 2.5) return false;
      if(!getElecCutBasedId(iLoose)) return false; // Electron selection criteria used for the emulation
      if(!convVeto) return false;
      if(lostHits != 0) return false;
      if(!getElecMVA(iTight)) return false;
      if(!getMultiIso(iTight)) return false;
    }
    if(!getGoodVertex(iTight)) return false;
    if(!getSIPcut(4)) return false;
    if(TightCharge != 2) return false;
    return true;
  }
  else if(gSelection == iWZSelec){ // Fakeable Objects for the WZ analysis. MVA WP (and related selection criteria) are defined inside the proper analysis
    Bool_t isMVALepton = true;    
    if(lep.isMuon){
      if(!isVetoLepton(lep)) isMVALepton = false;
      if(TightCharge != 2) isMVALepton = false;
      // Added MVA selection in the analysis
    }
    if(lep.isElec){
      if(!isVetoLepton(lep)) isMVALepton = false;
      if(TightCharge != 2) isMVALepton = false;
      if(!convVeto) isMVALepton = false;
      // Added MVA selection in the analysis
    }

    Bool_t isTopLepton = true;
    //Top/Stop ID    
    if(lep.isMuon){
      passId  = getMuonId(iTight);
      passIso = getRelIso04POG(iTight);
    }
    if(lep.isElec){
      passId = getElecCutBasedId(iTight) && lostHits <= 1;
      passIso = getRelIso03POG(iTight);
      if(TMath::Abs(etaSC) > 1.4442 && TMath::Abs(etaSC) < 1.566) isTopLepton = false;
    }
    if(lep.p.Pt() < 20 || TMath::Abs(lep.p.Eta()) > 2.4) isTopLepton = false;
    if(passId && passIso && ( (lep.isElec && getGoodVertex(iTight)) || (lep.isMuon && getGoodVertex(iMedium) ))){
      isTopLepton = true;
    }
    else isTopLepton = false;
    if (isTopLepton && isMVALepton){
      // To be updated when SUSYMVAId available
      //  lepMVASUSYId = getSUSYMVAId(lep, 2) + 10; //This is a terrible but quick way to do it
      return true;
    }
    else if (isTopLepton && !(isMVALepton)){
      lepMVASUSYId = 10; //This is a terrible but quick way to do it
      return true;  
    }
    else if ( (!isTopLepton) && isMVALepton){
      // To be updated when SUSYMVAId available
      //  lepMVASUSYId = getSUSYMVAId(lep, 2); //This is a terrible but quick way to do it
      return true;    
    }
    else return false;
  }
  return true;
}

//============================================================================================
//============================================== VETO LEPTONS
//============================================================================================
Bool_t LeptonSelector::isVetoLepton(Lepton lep){
  Bool_t passId; Bool_t passIso;
  if(gSelection == iStopTopSelec){
    if(lep.isMuon){
      passId = true;
      passIso = getRelIso03POG(iLooseWPforStop);
    }
    else{
      passId = getElecCutBasedId(iVeto);
      passIso = getRelIso03POG(iLooseWPforStop);
    } 
    return passId && passIso && getGoodVertex(iTight) && getSIPcut(4);
  }
  else if(gSelection == iTopSelec || gSelection == iTWSelec || gSelection == iWWSelec){
    return true;
  }
  else if(gSelection == i4tSelec){
    if(lep.isMuon){
      if(lep.p.Pt() < 20) return false;
      if(TMath::Abs(lep.p.Eta()) > 2.4) return false;
      if(!getMuonId(iMedium)) return false;
      if(!getMultiIso(iLoose)) return false;
    }
    if(lep.isElec){
      if(lep.p.Pt() < 20) return false;
      if(TMath::Abs(lep.p.Eta()) > 2.5) return false;
      if(!getElecCutBasedId(iLoose)) return false; // Electron selection criteria used for the emulation
      if(!convVeto) return false;
      if(lostHits != 0) return false;
        if(!getElecMVA(iLoose)) return false;
      if(!getMultiIso(iLoose)) return false;
    }
    //if(!getminiRelIso(iLoose)) return false;
    if(!getGoodVertex(iTight)) return false;
    if(!getSIPcut(4)) return false;
    if(TightCharge != 2) return false;
    return true;
  }
  else if(gSelection == iWZSelec){ // Fakeable Objects for the WZ analysis
    if(lep.isMuon){
      if(lep.p.Pt() < 10) return false;
      if(!isLooseLepton(lep)) return false;
      if(!mediumMuonId) return false;

    }
    if(lep.isElec){
      if(lep.p.Pt() < 10) return false;
      if(!isLooseLepton(lep)) return false;
      if(lostHits > 0) return false;
      if((sigmaIEtaIEta > 0.011 )*(lep.p.Eta() < 0.8) || (sigmaIEtaIEta > 0.011)*(lep.p.Eta() < 1.479)*(lep.p.Eta() > 0.8) || (sigmaIEtaIEta > 0.030)*(lep.p.Eta() > 1.479)) return false;
      if((HoE > 0.1)*(lep.p.Eta() < 0.8) || (HoE > 0.1)*(lep.p.Eta() < 1.479)*(lep.p.Eta() > 0.8) || (HoE > 0.07)*(lep.p.Eta() > 1.479)) return false;
      if((dEtaSC > 0.01)*(lep.p.Eta() < 0.8) || (dEtaSC > 0.01)*(lep.p.Eta() < 1.479)*(lep.p.Eta() > 0.8) || (dEtaSC > 0.008)*(lep.p.Eta() > 1.479)) return false;
      //if((dPhiSC > 0.04)*(lep.p.Eta() < 0.8) || (dPhiSC > 0.04)*(lep.p.Eta() < 1.479)*(lep.p.Eta() > 0.8) || (dPhiSC > 0.07)*(lep.p.Eta() > 1.479)) return false;
      if((eImpI<-0.05) || (eImpI > 0.01)*(lep.p.Eta() < 0.8) || (eImpI > 0.01)*(lep.p.Eta() < 1.479)*(lep.p.Eta() > 0.8) || (eImpI > 0.005)*(lep.p.Eta() > 1.479)) return false;

    }
    // To be updated when SUSYMVAId available
    //lepMVASUSYId = getSUSYMVAId(lep, 1);
    //if (lepMVASUSYId < 0) return false;
    return true;
  }
  return false;
}

//============================================================================================
//============================================== Loose leptons (or other)
//============================================================================================
Bool_t LeptonSelector::isLooseLepton(Lepton lep){
  Bool_t passId; Bool_t passIso;
  if(gSelection == i4tSelec){
    if(lep.isMuon){
      if(lep.p.Pt() < 15) return false;
      if(TMath::Abs(lep.p.Eta()) > 2.4) return false;
      // Smallest muon Id: Loose Id!! 
      if(!getMultiIso(iLoose)) return false;
    }
    if(lep.isElec){
      if(lep.p.Pt() < 15) return false;
      if(TMath::Abs(lep.p.Eta()) > 2.5) return false;
      if(!getElecCutBasedId(iLoose)) return false;
      if(!convVeto) return false;
      if(lostHits > 2) return false;
      if(!getElecMVA(iLoose)) return false;
      if(!getMultiIso(iLoose)) return false;
    }
    //if(!getminiRelIso(iLoose)) return false;
    if(!getGoodVertex(iTight)) return false;
    return true;
  }
  else if(gSelection == iWZSelec){ // Loose leptons are common to all tight WP
    if(lep.isMuon){
      if(lep.p.Pt() < 5) return false;
      if(TMath::Abs(lep.p.Eta()) > 2.4) return false;
      if(!getGoodVertex(iLoose)) return false;
      if(!getSIPcut(8)) return false;
      if(!getminiRelIso(iLoose)) return false;
    }
    if(lep.isElec){
      if(lep.p.Pt() < 7) return false;
      if(TMath::Abs(lep.p.Eta()) > 2.5) return false;
      if(!getGoodVertex(iLoose)) return false;    
      if(!getSIPcut(8)) return false;
      if(!getminiRelIso(iLoose)) return false;
      if(!getElecMVAId(iLoose,lep)) return false;  
      if((MVAID < -0.70)*(lep.p.Eta() < 0.8) || (MVAID < -0.83)*(lep.p.Eta() < 1.479)*(lep.p.Eta() > 0.8) || (MVAID < -0.92)*(lep.p.Eta() > 1.479)) return false;
      if(lostHits > 2) return false;
    }
    return true;
  }
  else if(gSelection == iTopSelec || gSelection == iStopTopSelec || gSelection == iTWSelec){
    // Same as good lepton but no looser cut on pT
    if(lep.isMuon){
      passId  = getMuonId(iTight);
      passIso = getRelIso04POG(iTight);
    }
    if(lep.isElec){
      passId = getElecCutBasedId(iTight) && lostHits <= 1;
      passIso = getRelIso03POG(iTight);
      if(TMath::Abs(etaSC) > 1.4442 && TMath::Abs(etaSC) < 1.566) return false;
    }
    if(lep.p.Pt() < 18 || TMath::Abs(lep.p.Eta()) > 2.4) return false;
    if(passId && passIso && ( (lep.isElec && getGoodVertex(iTight)) || (lep.isMuon && getGoodVertex(iMedium) ))) return true;
    else return false;
  }
    return true;
}


///////////////////////////////////////////////////////////////////////////
// You do not want to change anything below this point
///////////////////////////////////////////////////////////////////////////

void LeptonSelector::InsideLoop(){
  evt = Get<ULong64_t>("event");

  // Clear vectors...
  selLeptons.clear();
  looseLeptons.clear();
  genLeptons.clear();
  vetoLeptons.clear();
  vGenBquarks.clear();

  // Loop over the leptons and select
  nElec     = Get<Int_t>("nElectron");
  nMuon     = Get<Int_t>("nMuon");

  // Loop over the gen leptons and get gen info...
  if(!gIsData){
    ngenLep         = Get<Int_t>("nGenDressedLepton");
    for(Int_t i = 0; i < ngenLep; i++){
      GetGenLeptonVariables(i);
      tL = Lepton(tP, charge, type);
      //if(tL.p.Pt() > 20 && TMath::Abs(tL.p.Eta()) < 2.4) genLeptons.push_back(tL);
      genLeptons.push_back(tL);
    }
  }

  // Loop over reco leptons
  int index = 0; int i; int LepType = 11;
  for(i = 0; i < nElec+nMuon; i++){
    if(i < nElec){
      LepType = 11; 
      index = i;
    }
    else{
      LepType = 13;
      index = i- nElec;
    }
    GetLeptonVariables(index, LepType);
    tL = Lepton(tP, charge, type);
    if(tL.isMuon){
      tL.SetIso(RelIso04);
      tL.SetEnergyUnc(GetMuonEnergyScale());
    }
    else{
      tL.SetIso(RelIso03);
      tL.SetR9(R9);
      tL.SetEnergyUnc(ElecScale->GetUnc(tL.Pt(), tL.Eta(), tL.GetR9()));
    }
/*
    // Set status: Good, conv, flip, fake
    if(!gIsData){
      tL.SetGenMatch(kLGMgood);
      if(abs(mcMatchPDGID) != abs(pdgid))            tL.SetGenMatch(kLGMother);
      if(mcPromptGamma == 1)                         tL.SetGenMatch(kLGMconv);
      if( (tL.GetGenMatch() == kLGMtoGenLep || tL.GetGenMatch() == kLGMgood) && charge*mcMatchPDGID > 0 && abs(mcMatchPDGID) == abs(pdgid))  tL.SetGenMatch(kLGMflip); // flip and no fake
      if(mcMatchID == 0 || mcMatchID == -99)         tL.SetGenMatch(kLGMfake);
//      if(!mcPrompt)   tL.SetGenMatch(kLGMfake);
    }
*/
    if(isGoodLepton(tL)){
      //if(1){
      tL.SetSF(   LepSF->GetLeptonSF(     pt, eta, tL.type) ); // Set SF and error
      tL.SetSFerr(LepSF->GetLeptonSFerror(pt, eta, tL.type) );
      tL.idMVA = lepMVASUSYId;
      selLeptons.push_back(tL);
    }
    if(isVetoLepton(tL)){ // If you need to veto extra leptons...
      //tL.SetSF(1); tL.SetSFerr(1); // To be updated if ever needed
      tL.idMVA = lepMVASUSYId;
      vetoLeptons.push_back(tL);
    }
    if(isLooseLepton(tL)){ // A loose category... used in ttH, for example
      if(gSelection == iStopTopSelec){
        tL.SetSF(   LepSF->GetLeptonSF(     pt, eta, tL.type) ); // Set SF and error
        tL.SetSFerr(LepSF->GetLeptonSFerror(pt, eta, tL.type) );
        looseLeptons.push_back(tL);
      }
      else looseLeptons.push_back(tL);
    }
  }

  nSelLeptons   = selLeptons.size();
  nVetoLeptons   = vetoLeptons.size();
  nLooseLeptons = looseLeptons.size();
  nGenLeptons    = genLeptons.size();

  //=== Trigger SF
  TriggerSF = 1; TriggerSFerr = 0;
  if(gSelection == iTopSelec || gSelection == iStopTopSelec || gSelection == iTWSelec){
    if(nSelLeptons >= 2){
      if     (selLeptons.at(0).isMuon && selLeptons.at(1).isMuon){
        TriggerSF = LepSF->GetTrigDoubleMuSF(    selLeptons.at(0).p.Eta(), selLeptons.at(1).p.Eta());
        TriggerSFerr = LepSF->GetTrigDoubleMuSF_err(selLeptons.at(0).p.Eta(), selLeptons.at(1).p.Eta());
      }
      else if(selLeptons.at(0).isElec && selLeptons.at(1).isElec){
        TriggerSF = LepSF->GetTrigDoubleElSF(    selLeptons.at(0).p.Eta(), selLeptons.at(1).p.Eta());
        TriggerSFerr = LepSF->GetTrigDoubleElSF_err(selLeptons.at(0).p.Eta(), selLeptons.at(1).p.Eta());
      }
      else{
        TriggerSF = LepSF->GetTrigElMuSF(        selLeptons.at(0).p.Eta(), selLeptons.at(1).p.Eta());
        TriggerSFerr = LepSF->GetTrigElMuSF_err(    selLeptons.at(0).p.Eta(), selLeptons.at(1).p.Eta());
      }
    }
  }

  //=== FullSim/FastSim SF
  FSSF = 1; FSSFerr = 0; Float_t id;
  if(gIsFastSim && nSelLeptons >= 1){
    id = 11; if(selLeptons.at(0).isMuon) id = 13;
    FSSF    = LepSF->GetFSSF(selLeptons.at(0).p.Pt(), selLeptons.at(0).p.Eta(), id);
    FSSFerr = LepSF->GetFSSFerr(selLeptons.at(0).p.Pt(), selLeptons.at(0).p.Eta(), id);
  }


 selLeptons   = SortLeptonsByPt(selLeptons);
 vetoLeptons  = SortLeptonsByPt(vetoLeptons);
 looseLeptons = SortLeptonsByPt(looseLeptons);
 genLeptons   = SortLeptonsByPt(genLeptons);

  // Set params for the next selectors
  SetParam("selLeptons",  selLeptons );
  SetParam("vetoLeptons", vetoLeptons);
  SetParam("looseLeptons", looseLeptons);
  SetParam("genLeptons",  genLeptons );
  SetParam("nLeptonsFromTau", nLeptonsFromTau);
  SetParam("nGenLeptons", nGenLeptons);
  SetParam("nSelLeptons", nSelLeptons);
  SetParam("nVetoLeptons", nVetoLeptons);
  SetParam("nLooseLeptons", nVetoLeptons);

  SetParam("TriggerSF",    TriggerSF);
  SetParam("TriggerSFerr", TriggerSFerr);
  SetParam("FSSF",    FSSF);
  SetParam("FSSFerr", FSSFerr);
}

//################################################################
//## Get important variables
//################################################################
void LeptonSelector::GetLeptonVariables(Int_t i, int LepType){ // Once per muon, get all the info
  //tP.SetPxPyPzE(Get<Float_t>("LepGood_px", i), Get<Float_t>("LepGood_py", i), Get<Float_t>("LepGood_pz", i), Get<Float_t>("LepGood_energy", i));
  TString LeptonName = abs(LepType) == 11? "Electron" : "Muon";
  tP.SetPtEtaPhiM(Get<Float_t>(LeptonName + "_pt", i), Get<Float_t>(LeptonName + "_eta", i), Get<Float_t>(LeptonName + "_phi", i), Get<Float_t>(LeptonName + "_mass", i));

  pt = tP.Pt(); eta = tP.Eta(); energy = tP.Energy();
  charge        = Get<Int_t>(LeptonName + "_charge", i);
  type          = TMath::Abs(LepType) == 11 ? 1 : 0;
  pdgid         = LepType;
  dxy           = TMath::Abs(Get<Float_t>(LeptonName + "_dxy", i));
  dz            = TMath::Abs(Get<Float_t>(LeptonName + "_dz", i));
  miniIso       = Get<Float_t>(LeptonName + "_miniPFRelIso_all",i);
  RelIso03      = Get<Float_t>(LeptonName + "_pfRelIso03_all",i);
  sip           = Get<Float_t>(LeptonName + "_sip3d",i);
  MVATTH        = Get<Float_t>(LeptonName + "_mvaTTH",i);       //*
  TightCharge   = Get<Int_t>(LeptonName + "_tightCharge",i);      //*
  etaSC = tP.Eta(); // To be modified

  jetindex      = Get<Int_t>(LeptonName + "_jetIdx", i);  //index of the associated jet (-1 if none)
  genPartIndex  = Get<Int_t>(LeptonName + "_genPartIdx", i);
  genPartFlav   = Get<Int_t>(LeptonName + "_genPartFlav", i);

  RelIso04 = -1; mediumMuonId = -1; SegComp = -1; dEtaSC = -1; HoE = -1; eImpI = -1; 
  lostHits = -1; convVeto = -1; sigmaIEtaIEta = -1; MVAID = -1; R9 = -1;
  SF = 1;
  if(LepType == 13){  // ONLY MUONS
    RelIso04       = Get<Float_t>("Muon_pfRelIso04_all",i);
    mediumMuonId   = Get<Int_t>("Muon_mediumId",i);
    SegComp        = Get<Float_t>("Muon_segmentComp", i);
    tightVar      = Get<Int_t>("Muon_tightId", i);
  }
  else{              // ONLY ELECTRON
    tightVar      = Get<Int_t>("Electron_cutBased", i);
    dEtaSC = Get<Float_t>("Electron_deltaEtaSC", i);
    HoE            = Get<Float_t>("Electron_hoe", i);
    eImpI          = Get<Float_t>("Electron_eInvMinusPInv", i);
    lostHits       = Get<Int_t>("Electron_lostHits", i);
    convVeto       = Get<Int_t>("Electron_convVeto", i);
    sigmaIEtaIEta  = Get<Float_t>("Electron_sieie", i); // Es esta???
    MVAID          = Get<Float_t>("Electron_mvaSpring16GP_WP90",i);  //Electron_mvaSpring16GP_WP80, Electron_mvaSpring16GP_WP90, Electron_mvaSpring16HZZ_WPL
    R9             = Get<Float_t>("Electron_r9",i);
  }

}

void LeptonSelector::GetGenLeptonVariables(Int_t i){
  tP.SetPtEtaPhiM(Get<Float_t>("GenDressedLepton_pt", i), Get<Float_t>("GenDressedLepton_eta", i), Get<Float_t>("GenDressedLepton_phi", i), Get<Float_t>("GenDressedLepton_mass", i));
  pdgid = Get<Int_t>("GenDressedLepton_pdgId",i);
  charge = pdgid < 0 ? 1 : -1;
  type = TMath::Abs(pdgid) == 11 ? 1 : 0;
}


//####################################################
  // Need to compute
//  ptRel         = Get<Float_t>(LeptonName + "_jetPtRelv2",i);
//  ptRatio       = Get<Float_t>(LeptonName + "_jetPtRatiov2",i);
//  jetBTagCSV    = Get<Float_t>(LeptonName + "_jetBTagCSV",i);   //*

  // NO EXISTEN
//  etaSC         = TMath::Abs(Get<Float_t>(LeptonName + "_etaSc",i));
//  dPhiSC         = Get<Float_t>(LeptonName + "_dPhiScTrkIn", i);
//  MVASUSY        = Get<Float_t>(LeptonName + "_mvaSUSY",i);       //*
//  isGlobalMuon = Get<Int_t>(LeptonName + "_isGlobalMuon",i);
//  isTrackerMuon = Get<Int_t>(LeptonName + "_isTrackerMuon",i);

  /*
  // THESE VARIABLES NEED TO BE COMPUTED
  if(!gIsData){
    mcPrompt       = Get<Int_t>(LeptonName + "_mcPrompt", i);
    mcMatchID       = Get<Int_t>(LeptonName + "_mcMatchId", i);
    mcPromptGamma   = Get<Int_t>(LeptonName + "_mcPromptGamma", i);
    mcMatchPDGID    = Get<Int_t>(LeptonName + "_mcMatchPdgId", i);
  }*/


