#pragma once

#include "PAFChainItemSelector.h"
#include "Functions.h"
#include "BTagSFUtil.h"
#include <iostream>
#include <vector>
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"

const Int_t nWeights = 248;


class TWAnalysis : public PAFChainItemSelector{
  public:
    // Main methods
    TWAnalysis();
    virtual ~TWAnalysis(){}
    virtual void InsideLoop();
    virtual void Initialise();
    virtual void Summary();


    // Other methods
    void SetTWVariables();
    void ResetTWVariables();
    void GetLeptonVariables();
    void GetGenLepVariables();
    void GetJetVariables();
    void GetGenJetVariables();
    void GetMETandGenMET();

    void CalculateDressTWVariables();
    void CalculateTWVariables();

    void CalculateSFAndWeights();
    void SetMinimaAndMaxima();

    void DoesItReallyPassDress();
    void DoesItReallyPassReco();

    TLorentzVector GetSysVector(const TString& sys = "");
    Float_t GetTopPtRW();
    Float_t GetHTtot(const TString& sys = "", TString reg = "1j1b");
    Float_t GetLep1Jet1_DR(const TString& sys = "");
    Float_t GetLep12Jet12_DR(const TString& sys = "");
    Float_t GetLep12Jet12MET_DR(const TString& sys = "");

    void    SetTWBDT();

    // Attributes
    std::vector<Lepton> DressLeptons;
    std::vector<Lepton> selLeptons;
    std::vector<Jet>    selJets;
    std::vector<Jet>    selJetsJecUp;
    std::vector<Jet>    selJetsJecDown;
    std::vector<Jet>    selJetsJER;

    std::vector<Jet>    genJets;
    std::vector<Jet>    vetoJets;

    std::vector<Jet>    DressJets;
    std::vector<Jet>    DressLooseCentralJets;
    std::vector<Jet>    DressLooseFwdJets;

    std::vector<Jet> LooseCentralJets, LooseFwdJets;
    std::vector<Jet> LooseCentralJetsJESUp, LooseFwdJetsJESUp;
    std::vector<Jet> LooseCentralJetsJESDown, LooseFwdJetsJESDown;
    std::vector<Jet> LooseCentralJetsJERUp, LooseFwdJetsJERUp;

    TTree* fMiniTree;
    TH1F* fhDummy;
    Float_t   TLHEWeight[254];
    ULong64_t TEvent;

    Char_t TChannel, GenChannel;
    Bool_t TIsSS, TDressIsSS;
    Bool_t TPassMinDressSel, TPassMinRecoSel;
    Bool_t TPassReco, TPassRecoJESUp, TPassRecoJESDown, TPassRecoJERUp;
    Bool_t TPassDress, TPassPart;

    UChar_t TNBJets, TNBJetsJESUp, TNBJetsJESDown, TNBJetsJERUp;
    UChar_t NLooseCentral, NLooseCentralJESUp, NLooseCentralJESDown, NLooseCentralJERUp;
    Float_t FNLooseCentral, FNLooseCentralJESUp, FNLooseCentralJESDown, FNLooseCentralJERUp;
    UChar_t NLooseFwd, NLooseFwdJESUp, NLooseFwdJESDown, NLooseFwdJERUp;
    UChar_t NBLooseCentral, NBLooseCentralJESUp, NBLooseCentralJESDown, NBLooseCentralJERUp;
    Float_t FNBLooseCentral, FNBLooseCentralJESUp, FNBLooseCentralJESDown, FNBLooseCentralJERUp;
    UChar_t DressNJets, DressNLooseCentral, DressNLooseFwd, DressNBJets, DressNBLooseCentral;

    Float_t TLep1_Pt, TLep1_PtJESUp, TLep1_PtJESDown, TLep1_PtJERUp;
    Float_t TLep1_E, TLep1_EJESUp, TLep1_EJESDown, TLep1_EJERUp;
    Float_t TLep1_Phi, TLep1_PhiJESUp, TLep1_PhiJESDown, TLep1_PhiJERUp;
    Float_t TLep1_Eta, TLep1_EtaJESUp, TLep1_EtaJESDown, TLep1_EtaJERUp;
    Float_t TLep2_Pt, TLep2_PtJESUp, TLep2_PtJESDown, TLep2_PtJERUp;
    Float_t TLep2_E, TLep2_EJESUp, TLep2_EJESDown, TLep2_EJERUp;
    Float_t TLep2_Phi, TLep2_PhiJESUp, TLep2_PhiJESDown, TLep2_PhiJERUp;
    Float_t TLep2_Eta, TLep2_EtaJESUp, TLep2_EtaJESDown, TLep2_EtaJERUp;
    Float_t TJet1_Pt,  TJet1_PtJESUp,  TJet1_PtJESDown,  TJet1_PtJERUp;
    Float_t TJet1_E,   TJet1_EJESUp,   TJet1_EJESDown,   TJet1_EJERUp;
    Float_t TJet1_Eta, TJet1_EtaJESUp, TJet1_EtaJESDown, TJet1_EtaJERUp;
    Float_t TJet1_Phi, TJet1_PhiJESUp, TJet1_PhiJESDown, TJet1_PhiJERUp;
    Float_t TJet2_Pt,  TJet2_PtJESUp,  TJet2_PtJESDown,  TJet2_PtJERUp;
    Float_t TJet2_E,   TJet2_EJESUp,   TJet2_EJESDown,   TJet2_EJERUp;
    Float_t TJet2_Eta, TJet2_EtaJESUp, TJet2_EtaJESDown, TJet2_EtaJERUp;
    Float_t TLooseCentral1_Pt,  TLooseCentral1_PtJESUp,  TLooseCentral1_PtJESDown,  TLooseCentral1_PtJERUp;
    Float_t TLooseCentral1_Eta, TLooseCentral1_EtaJESUp, TLooseCentral1_EtaJESDown, TLooseCentral1_EtaJERUp;
    Float_t TLepMuon_Pt,  TLepMuon_PtJESUp,  TLepMuon_PtJESDown,  TLepMuon_PtJERUp;
    Float_t TLepMuon_E,   TLepMuon_EJESUp,   TLepMuon_EJESDown,   TLepMuon_EJERUp;
    Float_t TLepMuon_Phi, TLepMuon_PhiJESUp, TLepMuon_PhiJESDown, TLepMuon_PhiJERUp;
    Float_t TLepMuon_Eta, TLepMuon_EtaJESUp, TLepMuon_EtaJESDown, TLepMuon_EtaJERUp;
    Float_t TLepElec_Pt,  TLepElec_PtJESUp,  TLepElec_PtJESDown,  TLepElec_PtJERUp;
    Float_t TLepElec_E,   TLepElec_EJESUp,   TLepElec_EJESDown,   TLepElec_EJERUp;
    Float_t TLepElec_Phi, TLepElec_PhiJESUp, TLepElec_PhiJESDown, TLepElec_PhiJERUp;
    Float_t TLepElec_Eta, TLepElec_EtaJESUp, TLepElec_EtaJESDown, TLepElec_EtaJERUp;
    Float_t TLep1Lep2_Pt, TLep1Lep2_PtJESUp, TLep1Lep2_PtJESDown, TLep1Lep2_PtJERUp;
    Float_t TLep1Lep2_M, TLep1Lep2_MJESUp, TLep1Lep2_MJESDown, TLep1Lep2_MJERUp;
    Float_t TLep1Lep2_DPhi, TLep1Lep2_DPhiJESUp, TLep1Lep2_DPhiJESDown, TLep1Lep2_DPhiJERUp;
    Float_t TLep1Jet1_M, TLep1Jet1_MJESUp, TLep1Jet1_MJESDown, TLep1Jet1_MJERUp;
    Float_t TLep1Jet1_DPhi, TLep1Jet1_DPhiJESUp, TLep1Jet1_DPhiJESDown, TLep1Jet1_DPhiJERUp;
    Float_t TLep1Jet1_DR, TLep1Jet1_DRJESUp, TLep1Jet1_DRJESDown, TLep1Jet1_DRJERUp;
    Float_t TLep1Jet2_M, TLep1Jet2_MJESUp, TLep1Jet2_MJESDown, TLep1Jet2_MJERUp;
    Float_t TLep1Jet2_DPhi, TLep1Jet2_DPhiJESUp, TLep1Jet2_DPhiJESDown, TLep1Jet2_DPhiJERUp;
    Float_t TLep2Jet1_M, TLep2Jet1_MJESUp, TLep2Jet1_MJESDown, TLep2Jet1_MJERUp;
    Float_t TLep2Jet1_DPhi, TLep2Jet1_DPhiJESUp, TLep2Jet1_DPhiJESDown, TLep2Jet1_DPhiJERUp;
    Float_t TLep2Jet2_M, TLep2Jet2_MJESUp, TLep2Jet2_MJESDown, TLep2Jet2_MJERUp;
    Float_t TLep2Jet2_DPhi, TLep2Jet2_DPhiJESUp, TLep2Jet2_DPhiJESDown, TLep2Jet2_DPhiJERUp;
    Float_t TLep1Lep2Jet1_Pt, TLep1Lep2Jet1_PtJESUp, TLep1Lep2Jet1_PtJESDown, TLep1Lep2Jet1_PtJERUp;
    Float_t TLep12Jet12_DR, TLep12Jet12_DRJESUp, TLep12Jet12_DRJESDown, TLep12Jet12_DRJERUp;
    Float_t TLep12Jet12MET_DR, TLep12Jet12MET_DRJESUp, TLep12Jet12MET_DRJESDown, TLep12Jet12MET_DRJERUp;
    Float_t Sys_Pt, Sys_PtJESUp, Sys_PtJESDown, Sys_PtJERUp;
    Float_t Sys_E, Sys_EJESUp, Sys_EJESDown, Sys_EJERUp;
    Float_t Sys_Eta, Sys_EtaJESUp, Sys_EtaJESDown, Sys_EtaJERUp;
    Float_t Sys_M, Sys_MJESUp, Sys_MJESDown, Sys_MJERUp;
    Float_t Sys_Pz, Sys_PzJESUp, Sys_PzJESDown, Sys_PzJERUp;
    Float_t TSys_PtOverHTtot, TSys_PtOverHTtotJESUp, TSys_PtOverHTtotJESDown, TSys_PtOverHTtotJERUp;
    Float_t TLep1Lep2_HTOverHTtot, TLep1Lep2_HTOverHTtotJESUp, TLep1Lep2_HTOverHTtotJESDown, TLep1Lep2_HTOverHTtotJERUp;
    Float_t THT, THTJESUp, THTJESDown, THTJERUp;
    Float_t THTtot, THTtotJESUp, THTtotJESDown, THTtotJERUp;
    Float_t TMET, TMETJESUp, TMETJESDown, TMETJERUp;
    Float_t TMET_Phi, TMET_PhiJESUp, TMET_PhiJESDown, TMET_PhiJERUp;
    Float_t TBDT, TBDTJESUp, TBDTJESDown, TBDTJERUp;
    Float_t TBDT2j1b, TBDT2j1bJESUp, TBDT2j1bJESDown, TBDT2j1bJERUp;
    Float_t TC_jll, TC_jllJESUp, TC_jllJESDown, TC_jllJERUp;

    Float_t TDressLep1_Pt;
    Float_t TDressLep1_E;
    Float_t TDressLep1_Phi;
    Float_t TDressLep1_Eta;
    Float_t TDressLep2_Pt;
    Float_t TDressLep2_E;
    Float_t TDressLep2_Phi;
    Float_t TDressLep2_Eta;
    Float_t TDressLepMuon_Pt;
    Float_t TDressLepMuon_E;
    Float_t TDressLepMuon_Phi;
    Float_t TDressLepMuon_Eta;
    Float_t TDressLepElec_Pt;
    Float_t TDressLepElec_E;
    Float_t TDressLepElec_Phi;
    Float_t TDressLepElec_Eta;
    Float_t TDressJet1_Pt;
    Float_t TDressJet1_E;
    Float_t TDressJet1_Eta;
    Float_t TDressJet1_Phi;
    Float_t TDressJet2_Pt;
    Float_t TDressJet2_E;
    Float_t TDressJet2_Eta;
    Float_t TDressLep1Lep2_Pt;
    Float_t TDressLep1Lep2_M;
    Float_t TDressLep1Lep2_DPhi;
    Float_t TDressLep1Jet1_M;
    Float_t TDressLep1Jet1_DPhi;
    Float_t TDressLep1Jet2_M;
    Float_t TDressLep1Jet2_DPhi;
    Float_t TDressLep2Jet1_M;
    Float_t TDressLep2Jet1_DPhi;
    Float_t TDressLep2Jet2_M;
    Float_t TDressLep2Jet2_DPhi;
    Float_t DressSys_Pt;
    Float_t DressSys_E;
    Float_t DressSys_Eta;
    Float_t DressSys_M;
    Float_t DressSys_Pz;
    Float_t TDressHT;
    Float_t TDressMET, TDressMET_Phi;

    Float_t TPartMET, TPartMET_Phi;


    UChar_t TNSelLeps, DressNLeps;

    UChar_t TNJets, TNJetsJESUp, TNJetsJESDown, TNJetsJERUp;

    Double_t ElecSF, MuonSF, ElecSFUp, ElecSFDo, MuonSFUp, MuonSFDo, lepSF;
    Double_t TrigSF, TrigSFerr;
    Double_t PUSF, PUSF_Up, PUSF_Down;
    Double_t BtagSF;
    Double_t BtagSFBtagUp, BtagSFBtagDown;
    Double_t BtagSFMistagUp, BtagSFMistagDown;
    Double_t TWeight, TWeight_normal;
    Double_t TWeight_LepEffUp, TWeight_LepEffDown;
    Double_t TWeight_ElecEffUp, TWeight_ElecEffDown;
    Double_t TWeight_MuonEffUp, TWeight_MuonEffDown;
    Double_t TWeight_TrigUp, TWeight_TrigDown;
    Double_t TWeight_FSUp, TWeight_FSDown;
    Double_t TWeight_PUDown, TWeight_PUUp;
    Double_t TWeight_MistagUp, TWeight_MistagDown;
    Double_t TWeight_BtagUp, TWeight_BtagDown;

    TMVA::Reader* BDT1j1b;
    TMVA::Reader* BDT1j1bJESUp;
    TMVA::Reader* BDT1j1bJESDown;
    TMVA::Reader* BDT1j1bJER;
    TMVA::Reader* BDT2j1b;
    TMVA::Reader* BDT2j1bJESUp;
    TMVA::Reader* BDT2j1bJESDown;
    TMVA::Reader* BDT2j1bJER;


  protected:
    // Parameters
    Bool_t   gIsData;
    TString  gSampleName;
    Bool_t   gIsTTbar;
    Bool_t   gIsLHE;
    TString  gOptions;
    Bool_t   gPUWeight;
    Bool_t   passMETfilters;
    Bool_t   passTrigger;
    UShort_t year;

    ClassDef(TWAnalysis, 0);
};
