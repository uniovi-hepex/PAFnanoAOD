#pragma once

#include "PAFChainItemSelector.h"
#include "Functions.h"
#include "LeptonSF.h"
//#include "PUWeight.h"


class EventBuilder : public PAFChainItemSelector{

  public:

    std::vector<Float_t> CountLHE;
    std::vector<Float_t> LHEWeights;

    EventBuilder();
    virtual ~EventBuilder();
    virtual void InsideLoop();
    virtual void Initialise();
    virtual void Summary();

    //Bool_t PassHLT_Elec; Bool_t PassHLT_Muon; Bool_t PassHLT_ElMu;
    Bool_t METfilters;
    Bool_t passTrigger;
    Bool_t passTrigger2;
    Bool_t isSS;
    Bool_t gIsFastSim;
    Bool_t makeeffhistos;

    //Float_t  TrigSFElec;      Float_t  TrigSFMuon;      Float_t  TrigSFElMu;
    //Float_t  TrigSFElec_Up;   Float_t  TrigSFElMu_Up;   Float_t  TrigSFMuon_Up;
    //Float_t  TrigSFElec_Down; Float_t  TrigSFMuon_Down; Float_t  TrigSFElMu_Down;

    Float_t TriggerSF; Float_t TriggerSF_Up; Float_t TriggerSF_Down; Float_t TriggerSF_err;
    Float_t PUSF;      Float_t PUSF_Up;      Float_t PUSF_Down;

    TH2F* ElecTrigEffNum;
    TH2F* ElecTrigEffDen;
    TH2F* MuonTrigEffNum;
    TH2F* MuonTrigEffDen;
    TH2F* ElElTrigEffNum;
    TH2F* MuMuTrigEffNum;
    TH2F* ElMuTrigEffNum;
    TH2F* ElElTrigEffDen;
    TH2F* MuMuTrigEffDen;
    TH2F* ElMuTrigEffDen;

    Float_t NormWeight; // Nominal
    Float_t Weight;  // CrossSection/NumberOfGenEvents
    Float_t genWeight;   // For aMCatNLO samples

    Float_t nTrueInt;
    Int_t   gChannel;
    std::vector<Lepton> selLeptons;
    std::vector<Lepton> vetoLeptons;

  protected:

    LeptonSF *TriggSF;
    //PUWeight *fPUWeight;
    //PUWeight *fPUWeightUp;
    //PUWeight *fPUWeightDown;

    Bool_t PassesMETfilters();

    Bool_t TrigElMu();
    Bool_t TrigElEl();
    Bool_t TrigMuMu();
    Bool_t Trig3l4l();

    Bool_t PassesDoubleElecTrigger();
    Bool_t PassesDoubleMuonTrigger();
    Bool_t PassesElMuTrigger();
    Bool_t PassesSingleElecTrigger();
    Bool_t PassesSingleMuonTrigger();
    Bool_t PassesMETtrigger();

    Bool_t gIsSingleMuon;
    Bool_t gIsSingleElec;
    Bool_t gIsDoubleMuon;
    Bool_t gIsDoubleElec;
    Bool_t gIsMuonEG;
    Bool_t gIsMET;
    Bool_t gIsEGamma;
    Bool_t  gIsData;
    Int_t run;
    Int_t year;
    Int_t   gSelection;
    Bool_t gPUWeigth;
    TString selection;
    TString gSampleName;
    TString gPathToHeppyTrees;
    Float_t  gXSec;
    Float_t  gCount;
    Bool_t   gIsMCatNLO;
    Int_t gNEntries;
    Float_t  gSumOfWeights;
    Int_t nEntries;
    Long64_t Count;
    Float_t xsec;
    Int_t nProcessedEvents;

    TString gOptions;
    Bool_t gIs2016;
    Bool_t gIs2017;
    Bool_t gIs2018;
  

    void SetCountLHE();

    ClassDef(EventBuilder, 0);
};
