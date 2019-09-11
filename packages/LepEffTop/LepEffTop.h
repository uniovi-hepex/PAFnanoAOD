#pragma once

#include "PAFChainItemSelector.h"
#include "Lepton.h"
#include "LeptonSF.h"
#include "Functions.h"
#include <iostream>
#include <vector>
#include "LinkDef.h"

enum ePass{kPass, kFail};
enum eSSOS{kOS, kSS};
enum ePrompt{kPrompt, kNonprompt};

class LepEffTop: public PAFChainItemSelector{

  public:

    LepEffTop(); 
    virtual ~LepEffTop() {}
    virtual void InsideLoop();
    virtual void Initialise();
    virtual void Summary();

    std::vector<Lepton> leptons;
    std::vector<Lepton> genLeptons;

    TH1F* hMuonPt[2][2][2];
    TH1F* hMuonEta[2][2][2];
    TH1F* hMuonPhi[2][2][2];
    TH1F* hMuonIso[2][2][2];
    TH1F* hElecPt[2][2][2];
    TH1F* hElecEta[2][2][2];
    TH1F* hElecPhi[2][2][2];
    TH1F* hElecIso[2][2][2];

    TString GetLabel(Int_t i, Int_t j, Int_t k);
    void InitHistos();
    void FillAll(Int_t i, Int_t j, Int_t k);
    void FillHistos();

  protected:

    Int_t year;
    TString selection;
    Bool_t gIsData;
    Bool_t gIsFastSim;
    Int_t  gSelection;
    Bool_t gDoLepGood; 
    Bool_t gPUWeigth;
    Bool_t gPrefire;
    Bool_t gIs2016;
    Bool_t gIs2017;
    Bool_t gIs2018;
    Int_t  gChannel;
    TString gOptions;
    TString localPath;

    // Trigger     
    Float_t TriggerSF;
    Bool_t  PassMETFilters;
    Bool_t  PassTrigger;
    Bool_t  isSS;
    Float_t NormWeight; 
    Float_t PUSF; Float_t PrefWeight;
    Lepton Elec; Lepton Muon;
    Int_t njets; Int_t nbtags;

    // Lep variables
    Lepton tL;
    Int_t nLep;
    Int_t nElec;
    Int_t nMuon;
    Int_t evt;
    TLorentzVector tP; 
    Float_t pt;
    Float_t eta;
    Float_t energy;
    Int_t   charge; 
    Int_t   type;
    Int_t   pdgid;
    Int_t   tightVar;
    Int_t   mediumMuonId;
    Float_t etaSC;
    Float_t RelIso03;
    Float_t RelIso04;
    Float_t ptRel;
    Float_t ptRatio;
    Float_t miniIso;
    Float_t sigmaIEtaIEta;
    Float_t dEtaSC;
    Float_t dPhiSC;
    Float_t HoE;
    Float_t eImpI;
    Float_t rho;
    Int_t   lostHits;
    Int_t   convVeto;
    Float_t dxy;
    Float_t dz; 
    Float_t sip;
    Float_t SF;
    Float_t MVATTH;
    Float_t MVASUSY;
    Int_t  TightCharge;
    Float_t MVAID;
    Float_t jetBTagCSV;
    Float_t SegComp;
    Int_t isGlobalMuon;
    Int_t isTrackerMuon;
    Int_t lepMVASUSYId;
    Float_t R9;
    Int_t mcMatchID;
    Int_t mcPrompt;
    Int_t mcPromptGamma;
    Int_t mcMatchPDGID;
    Int_t jetindex;
    Int_t genPartIndex;
    Int_t genPartFlav;
    
    // genLeptons
    Int_t ngenLep;
    Int_t gpdgMId;

    Int_t nSelLeptons;
    Int_t nLooseLeptons;

    void GetLeptonVariables(Int_t i, Int_t lepType);
    void GetGenLeptonVariables(Int_t i);

  ClassDef(LepEffTop, 0);
};
