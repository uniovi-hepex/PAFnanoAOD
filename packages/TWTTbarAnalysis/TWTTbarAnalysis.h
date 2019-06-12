#pragma once

#include "PAFChainItemSelector.h"
#include "Functions.h"
#include "BTagSFUtil.h"
#include <iostream>
#include <vector>

const Int_t nChannels = 3;
enum eSysts   {inorm, nSysts};
const Int_t nWeights = 248;
const TString gChanLabel[nChannels] = {"ElMu", "Muon","Elec"};


class TWTTbarAnalysis : public PAFChainItemSelector{
  public:
    TWTTbarAnalysis();
    virtual ~TWTTbarAnalysis(){}
    virtual void InsideLoop();
    virtual void Initialise();
    virtual void Summary();
    
    std::vector<Lepton> genLeptons  ;
    std::vector<Lepton> selLeptons  ;
    std::vector<Lepton> vetoLeptons ;
    std::vector<Jet>    selJets ;
    std::vector<Jet>    selJetsJecUp   ;
    std::vector<Jet>    selJetsJecDown ;
    std::vector<Jet>    selJetsJER     ;

    std::vector<Jet>    Jets15  ;
    std::vector<Jet>    genJets  ;
    std::vector<Jet>    mcJets  ;
    std::vector<Jet>    vetoJets;

    std::vector<Jet>    GenJets;
    std::vector<Jet>    GenLooseCentralJets;
    std::vector<Jet>    GenLooseFwdJets;
    std::vector<Lepton> GenLeps;
    
    TTree* fMini1j1t;
    
    Bool_t   hasTW;
    Double_t  TLHEWeight[254];

    void GetLeptonVariables(std::vector<Lepton> selLeptons, std::vector<Lepton> VetoLeptons);
    void GetJetVariables(std::vector<Jet> selJets, std::vector<Jet> cleanedJets15, Double_t ptCut = 30);
    void GetGenJetVariables(std::vector<Jet> genJets, std::vector<Jet> mcJets);
    void GetGenLepVariables();
    Double_t getTopPtRW();
    void GetMET();
    void GetGenMET();
    Int_t nFiduJets; Int_t nFidubJets;
    
    void CalculateSFAndWeights();
    void DoesItReallyPassDress();
    void DoesItReallyPassReco();
    void SetMinimaAndMaxima();
    
    Double_t  TrigSF;
    Double_t  TrigSFerr;
    Double_t  PUSF;
    Double_t  PUSF_Up;
    Double_t  PUSF_Down;
    Int_t    gChannel;
    TString  gPar;
    Bool_t   passMETfilters;
    Bool_t   passTrigger;
    Bool_t   isSS;
    Double_t  NormWeight;
    Double_t  BtagSF          ;
    Double_t  BtagSFBtagUp    ;
    Double_t  BtagSFBtagDown  ;
    Double_t  BtagSFMistagUp  ;
    Double_t  BtagSFMistagDown;

    Double_t  TVetoJet1_Pt;
    Double_t  TVetoJet1_Eta;
    Double_t  TVetoJet2_Pt;
    Double_t  TVetoJet2_Eta;
    Double_t  TVetoJet3_Pt;
    Double_t  TVetoJet3_Eta;
  
    void     CalculateTWTTbarVariables();
    void     get20Jets();
    void     ResetTWTTbarVariables();
    void     SetTWTTbarVariables();
    Double_t getDilepMETJetPt(Int_t sys = 0);
    Double_t getDilepPt();
    Double_t getDilepJetPt(const TString& sys = "Norm");
    Double_t getLep1METJetPt(const TString& sys = "Norm");
    Double_t getPtSys(TLorentzVector*, Int_t);
    Double_t getDilepMETJet1Pz(const TString& sys = "Norm");
    Double_t getPzSys(TLorentzVector*, Int_t);
    Double_t getDPtDilep_JetMET(const TString& sys = "Norm");
    Double_t getDPtDilep_MET(const TString& sys = "Norm");
    Double_t getDPtLep1_MET(const TString& sys = "Norm");
    Double_t getDPhiLep1_Lep2();
    Double_t getDPhiLep1_Jet1(const TString& sys = "Norm");
    Double_t getDPhiLep2_Jet1(const TString& sys = "Norm");
    Double_t getDPtL1_L2(const TString& sys = "Norm");
    Double_t getDeltaPt(vector<TLorentzVector>, vector<TLorentzVector>);
    Double_t getDeltaPhi(vector<TLorentzVector>, vector<TLorentzVector>);
    Double_t getSysM(const TString& sys = "Norm");
    Double_t getM(vector<TLorentzVector>);
    Double_t getHTtot(const TString& sys = "Norm");
    Double_t getHTtot2j(const TString& sys = "Norm");

    Double_t getDeltaRDilep_METJets12(Int_t sys = 0);
    Double_t getDeltaRDilep_Jets12(Int_t sys = 0);
    Double_t getDeltaRLep1_Jet1(Int_t sys = 0);

    Double_t getDeltaR(vector<TLorentzVector>, vector<TLorentzVector>);

    Double_t getETSys(Int_t sys = 0);
    Double_t getET_LLJetMET(Int_t sys = 0);

    Double_t getDilepMETPt(Int_t sys= 0);


    //Variables
    
    Double_t ElecSF, MuonSF, ElecSFUp, ElecSFDo, MuonSFUp, MuonSFDo, lepSF;
    Double_t TWeight;
    Double_t TMll;
    Double_t TMllJESUp;
    Double_t TMllJESDown;
    Double_t TMllJERUp;
    Double_t TMET;
    Double_t TGenMET;
    Double_t TMET_Phi;
    Double_t TMET_PhiJESUp;
    Double_t TMET_PhiJESDown;
    Double_t TMET_PhiJERUp;
    Double_t TLeadingJetPt    ;
    Double_t TLeadingJetEta   ;
    Double_t TLeadingJetCSV   ;
    Double_t TSubLeadingJetPt ;
    Double_t TSubLeadingJetEta;
    Double_t TSubLeadingJetCSV;

    Int_t   TNVetoLeps;
    Int_t   TNSelLeps;
    Int_t   TChannel;
    Bool_t   TIsSS;
    Bool_t   TIsFid;
    Double_t TLep_Pt[10];
    Double_t TLep_Eta[10];
    Double_t TLep_Phi[10];
    Double_t TLep_E[10];
    Double_t TLep_Charge[10];

    Int_t   TNJets;
    Int_t   TNVetoJets;
    Int_t   TNBtags;
    Double_t TJet_Pt[20];
    Double_t TJet_Eta[20];
    Double_t TJet_Phi[20];
    Double_t TJet_E[20];
    Int_t   TJet_isBJet[20];
    Double_t THT;       // HT
    Double_t THTtot;
    Double_t THTtotJESUp;
    Double_t THTtotJESDown;
    Double_t THTtotJERUp;
    Double_t THTtot2j;
    Double_t THTtot2jJESUp;
    Double_t THTtot2jJESDown;
    Double_t THTtot2jJERUp;
    Double_t THtRejJ2;

    Double_t TDilepMETPt     ;
    Double_t TETSys          ;
    Double_t TET_LLJetMET    ;
    Double_t TDPtL1_L2       ;
    Double_t TDPtJ2_L2       ;
    Double_t TDR_L1_J1       ;
    Double_t TDR_L1L2_J1J2   ;
    Double_t TDR_L1L2_J1J2MET;

    Double_t TM_LeadingB       ;
    Double_t TM_SubLeadingB    ;
    Double_t TM_LeadingBj2        ;
    Double_t TM_SubLeadingBj2     ;
    Double_t TM_bjetlepton_minmax ;
    Double_t TE_LLB            ;
    Double_t TMT_LLMETB        ;
    Double_t TM_LLB            ;
    Double_t TLeadingJetE       ;
    Double_t TLeadingJetPhi     ;
    Double_t TLeadingLepPhi     ;
    Double_t TLeadingLepE       ;
    Double_t TSubLeadingLepE    ;
    Double_t TSubLeadingLepPhi  ;
    Double_t TLLMETBEta      ;
    Double_t TDPhiLL      ;
    Double_t TDPhiLeadJet ;
    Double_t TDPhiSubLeadJet ;
    
    Double_t TM_LeadingBJESUp       ;
    Double_t TM_SubLeadingBJESUp    ;
    Double_t TE_LLBJESUp            ;
    Double_t TMT_LLMETBJESUp        ;
    Double_t TM_LLBJESUp            ;
    Double_t TLeadingJetPtJESUp    ;
    Double_t TLeadingJetEJESUp       ;
    Double_t TLeadingJetPhiJESUp     ;
    Double_t TLeadingJetEtaJESUp;
    Double_t TLeadingLepPtJESUp;
    Double_t TLeadingLepEJESUp;
    Double_t TLeadingLepPhiJESUp;
    Double_t TLeadingLepEtaJESUp;
    Double_t TSubLeadingLepPtJESUp;
    Double_t TSubLeadingLepEJESUp;
    Double_t TSubLeadingLepPhiJESUp;
    Double_t TSubLeadingLepEtaJESUp;
    Double_t TLLMETBEtaJESUp      ;
    Double_t TDPhiLLJESUp      ;
    Double_t TDPhiLeadJetJESUp ;
    Double_t TDPhiSubLeadJetJESUp ;
    
    Double_t TM_LeadingBJESDown       ;
    Double_t TM_SubLeadingBJESDown    ;
    Double_t TE_LLBJESDown            ;
    Double_t TMT_LLMETBJESDown        ;
    Double_t TM_LLBJESDown            ;
    Double_t TLeadingJetPtJESDown    ;
    Double_t TLeadingJetEJESDown       ;
    Double_t TLeadingJetPhiJESDown     ;
    Double_t TLeadingJetEtaJESDown   ;
    Double_t TLeadingLepPtJESDown;
    Double_t TLeadingLepEJESDown;
    Double_t TLeadingLepPhiJESDown;
    Double_t TLeadingLepEtaJESDown;
    Double_t TSubLeadingLepPtJESDown;
    Double_t TSubLeadingLepEJESDown;
    Double_t TSubLeadingLepPhiJESDown;
    Double_t TSubLeadingLepEtaJESDown;
    Double_t TLLMETBEtaJESDown      ;
    Double_t TDPhiLLJESDown      ;
    Double_t TDPhiLeadJetJESDown ;
    Double_t TDPhiSubLeadJetJESDown ;
    
    Double_t TM_LeadingBJERUp       ;
    Double_t TM_SubLeadingBJERUp    ;
    Double_t TE_LLBJERUp            ;
    Double_t TMT_LLMETBJERUp        ;
    Double_t TM_LLBJERUp            ;
    Double_t TLeadingJetPtJERUp    ;
    Double_t TLeadingJetEJERUp       ;
    Double_t TLeadingJetPhiJERUp     ;
    Double_t TLeadingJetEtaJERUp   ;
    Double_t TLeadingLepPtJERUp;
    Double_t TLeadingLepEJERUp;
    Double_t TLeadingLepPhiJERUp;
    Double_t TLeadingLepEtaJERUp;
    Double_t TSubLeadingLepPtJERUp;
    Double_t TSubLeadingLepEJERUp;
    Double_t TSubLeadingLepPhiJERUp;
    Double_t TSubLeadingLepEtaJERUp;
    Double_t TLLMETBEtaJERUp      ;
    Double_t TDPhiLLJERUp      ;
    Double_t TDPhiLeadJetJERUp ;
    Double_t TDPhiSubLeadJetJERUp ;
    
    Double_t TGenLeadingJetPt   ;
    Double_t TGenLeadingJetE    ;
    Double_t TGenLeadingJetPhi  ;
    Double_t TGenLeadingJetEta  ;
    Double_t TGenLeadingLooseFwdJetPt   ;
    Double_t TGenLeadingLooseFwdJetE    ;
    Double_t TGenLeadingLooseFwdJetPhi  ;
    Double_t TGenLeadingLooseFwdJetEta  ;
    Double_t TGenSubLeadingLooseCentralJetPt   ;
    Double_t TGenSubLeadingLooseCentralJetE    ;
    Double_t TGenSubLeadingLooseCentralJetPhi  ;
    Double_t TGenSubLeadingLooseCentralJetEta  ;
    Double_t TGenLeadingLepPt   ;
    Double_t TGenLeadingLepE    ;
    Double_t TGenLeadingLepPhi  ;
    Double_t TGenLeadingLepEta  ;
    Double_t TGenSubLeadingLepPt ;
    Double_t TGenSubLeadingLepE  ;
    Double_t TGenSubLeadingLepPhi;
    Double_t TGenSubLeadingLepEta;
    Double_t TGenM_LeadingB    ;
    Double_t TGenM_SubLeadingB ;
    Double_t TGenE_LLB         ;
    Double_t TGenM_LeadingBj2  ;
    Double_t TGenM_SubLeadingBj2;
    Double_t TGenM_bjetlepton_minmax;
    Double_t TGenMT_LLMETB     ;
    Double_t TGenM_LLB         ;
    Double_t TGenDilepPt        ;
    Double_t TGenDilepJetPt     ;
    Double_t TGenDilepMETJetPt  ;
    Double_t TGenHTtot          ;
    Double_t TGenMET_Phi        ;
    Double_t TGenDilepMETJet1Pz ;
    Double_t TGenLLMETBEta      ;
    Double_t TGenMSys           ;
    Double_t TGenMll           ;
    Double_t TGenDPhiLL      ;
    Double_t TGenDPhiLeadJet ;
    Double_t TGenDPhiSubLeadJet ;
    Int_t GenChannel           ;
    Bool_t TGenIsSS;
    Bool_t Tpassgen            ;
    Bool_t Tpassreco           ;
    Bool_t TpassrecoJESUp      ;
    Bool_t TpassrecoJESDown    ;
    Bool_t TpassrecoJERUp      ;
    Double_t TWeight_normal     ;
    Int_t nGenJets;
    Int_t nGenbJets;
    Int_t nGenLooseCentralJets;
    Int_t nGenLooseCentralbJets;
    Int_t nGenLooseFwdJets;
    Int_t nGenLeps;
    Int_t nGenJets;
    Int_t nGenLeps;
    Int_t nGenMET;
    Jet tJ;
    TLorentzVector tpJ;
    Lepton tL;
    TLorentzVector tpL;
    TLorentzVector GenMET;
    TLorentzVector tMET;

    // For systematics...
    Int_t   TNJetsJESUp;
    Int_t   TNJetsJESDown;
    Int_t   TNJetsJERUp;
    Int_t   TNBtagsJESUp;
    Int_t   TNBtagsJESDown;
    Int_t   TNBtagsJERUp;
    Double_t TJetJESUp_Pt[20];
    Double_t TJetJESDown_Pt[20];
    Double_t TJetJER_Pt[20];
    Double_t THTJESUp;
    Double_t THTJESDown;
    Double_t THTJERUp;

    Int_t   TNISRJets;
    Double_t TMETJESUp;
    Double_t TMETJESDown;
    Double_t TMETJERUp;
    Double_t TMT2llJESUp;
    Double_t TMT2llJESDown;

    Double_t  TWeight_LepEffUp;
    Double_t  TWeight_LepEffDown;
    Double_t  TWeight_ElecEffUp;
    Double_t  TWeight_ElecEffDown;
    Double_t  TWeight_MuonEffUp;
    Double_t  TWeight_MuonEffDown;
    Double_t  TWeight_TrigUp;
    Double_t  TWeight_TrigDown;
    Double_t  TWeight_FSUp;
    Double_t  TWeight_PUDown;
    Double_t  TWeight_PUUp;
    Double_t  TWeight_FSDown;

    
    Double_t  TWeight_MistagUp  ;
    Double_t  TWeight_MistagDown;
    Double_t  TWeight_BtagUp    ;
    Double_t  TWeight_BtagDown  ;
    
    Double_t  LeadingLeptPt_   ;
    Double_t  LeadingLeptEta_  ;
    Double_t  jetPtSubLeading_ ;
    Double_t  jetEtaSubLeading_;

    Double_t  TLeadingLepPt    ;
    Double_t  TLeadingLepEta   ;
    Double_t  TDilepPt         ;
    Double_t  TSubLeadingLepPt ;
    Double_t  TSubLeadingLepEta;


    Double_t  DilepMETJetPt  , DilepMETJetPtJESUp  , DilepMETJetPtJESDown  , DilepMETJetPtJERUp  ;
    Double_t  Lep1METJetPt   , Lep1METJetPtJESUp   , Lep1METJetPtJESDown   , Lep1METJetPtJERUp   ;
    Double_t  DPtDilep_JetMET, DPtDilep_JetMETJESUp, DPtDilep_JetMETJESDown, DPtDilep_JetMETJERUp;
    Double_t  DPtDilep_MET   , DPtDilep_METJESUp   , DPtDilep_METJESDown   , DPtDilep_METJERUp   ;
    Double_t  DPtLep1_MET    , DPtLep1_METJESUp    , DPtLep1_METJESDown    , DPtLep1_METJERUp    ;
    Double_t  DilepMETJet1Pz , DilepMETJet1PzJESUp , DilepMETJet1PzJESDown , DilepMETJet1PzJERUp ;
    Int_t    nLooseCentral  , nLooseCentralJESUp  , nLooseCentralJESDown  , nLooseCentralJERUp  ;
    Int_t    nLooseFwd      , nLooseFwdJESUp      , nLooseFwdJESDown      , nLooseFwdJERUp      ;
    Int_t    nBLooseCentral , nBLooseCentralJESUp , nBLooseCentralJESDown , nBLooseCentralJERUp ;
    Int_t    nBLooseFwd     , nBLooseFwdJESUp     , nBLooseFwdJESDown     , nBLooseFwdJERUp     ;
    Double_t  TJet2csv       , TJet2csvJESUp       , TJet2csvJESDown       , TJet2csvJERUp       ;
    Double_t  MSys           , MSysJESUp           , MSysJESDown           , MSysJERUp           ;
    Double_t  TJetLoosept    , TJetLooseptJESUp    , TJetLooseptJESDown    , TJetLooseptJERUp    ;
    Double_t  TJetLooseCentralpt, TJetLooseCentralptJESUp    , TJetLooseCentralptJESDown    , TJetLooseCentralptJERUp    ;
    Double_t  TJetLooseFwdpt, TJetLooseFwdptJESUp    , TJetLooseFwdptJESDown    , TJetLooseFwdptJERUp    ;
    Double_t  C_jll          , C_jllJESUp          , C_jllJESDown          , C_jllJERUp          ;
    Double_t  DilepJetPt     , DilepJetPtJESUp     , DilepJetPtJESDown     , DilepJetPtJERUp     ;

    Int_t    nBTotal          , nBTotalJESUp          , nBTotalJESDown          , nBTotalJERUp          ; 
    Double_t  DilepmetjetOverHT, DilepmetjetOverHTJESUp, DilepmetjetOverHTJESDown, DilepmetjetOverHTJERUp; 
    Double_t  HTLepOverHT      , HTLepOverHTJESUp      , HTLepOverHTJESDown      , HTLepOverHTJERUp      ; 
    Double_t  TJet1_pt         , TJet1_ptJESUp         , TJet1_ptJESDown         , TJet1_ptJERUp         ;

    Double_t  TDilepPtJESUp          , TDilepPtJESDown         , TDilepPtJERUp         ; 
    Double_t  TDilepMETPtJESUp       , TDilepMETPtJESDown      , TDilepMETPtJERUp      ;
    Double_t  TETSysJESUp            , TETSysJESDown           , TETSysJERUp           ;
    Double_t  TET_LLJetMETJESUp      , TET_LLJetMETJESDown     , TET_LLJetMETJERUp     ;
    Double_t  THtRejJ2JESUp          , THtRejJ2JESDown         , THtRejJ2JERUp         ;
    Double_t  TDPtL1_L2JESUp         , TDPtL1_L2JESDown        , TDPtL1_L2JERUp        ;
    Double_t  TDPtJ2_L2JESUp         , TDPtJ2_L2JESDown        , TDPtJ2_L2JERUp        ;
    Double_t  TDR_L1_J1JESUp         , TDR_L1_J1JESDown        , TDR_L1_J1JERUp        ;
    Double_t  TDR_L1L2_J1J2JESUp     , TDR_L1L2_J1J2JESDown    , TDR_L1L2_J1J2JERUp    ;
    Double_t  TDR_L1L2_J1J2METJESUp  , TDR_L1L2_J1J2METJESDown , TDR_L1L2_J1J2METJERUp ;
    Double_t  LeadingLeptPt_JESUp    , LeadingLeptPt_JESDown   , LeadingLeptPt_JERUp   ;
    Double_t  LeadingLeptEta_JESUp   , LeadingLeptEta_JESDown  , LeadingLeptEta_JERUp  ;
    Double_t  jetPtSubLeading_JESUp  , jetPtSubLeading_JESDown , jetPtSubLeading_JERUp ;
    Double_t  jetEtaSubLeading_JESUp , jetEtaSubLeading_JESDown, jetEtaSubLeading_JERUp;
    Double_t  TJet2_Pt, TJet2_PtJESUp,TJet2_PtJESDown, TJet2_PtJERUp;


// Histograms
//=====================================================0
    TH1F* fhDummy;
    
  protected:
    Bool_t  gIsData;
    Bool_t  gDoSyst;
    Int_t   gSelection;
    TString gSampleName;
    Bool_t  gIsTTbar;
    Bool_t  gIsTW;
    Bool_t  gIsLHE;
    Bool_t  gIsFSRUp;
    Bool_t  gIsFSRDown;
    
    ClassDef(TWTTbarAnalysis, 0);
};
