#pragma once

#include "PAFChainItemSelector.h"
#include "Functions.h"
#include <iostream>
#include <vector>

class CalcEfficiencies : public PAFChainItemSelector{
  public:

    CalcEfficiencies();
    virtual ~CalcEfficiencies(){}
    virtual void InsideLoop();
    virtual void Initialise();
    virtual void Summary();
    
    std::vector<Jet> genBJets;
    std::vector<Jet> taggedJets;

    TH1F* NGenBJets;
    TH1F* NGenBTags;

    void SetJetVariables();
    void SetGenJetVariables();
    Float_t GetNBJets();
    Float_t GetNBTags();

    Int_t gChanel;
    Bool_t passMETfilters;
    Bool_t passTrigger;
    Bool_t isSS;
    Float_t NormWeight;
  
    // Tree Variables
    Float_t TWeight;   // Total nominal weight
    Float_t TMll;      // Invariant mass
    Float_t TMET;      // MET
    Float_t TMET_Phi;  // MET phi

    Int_t   TNVetoLeps;
    Int_t   TNSelLeps;
    Int_t TChannel;
    Float_t TLep_Pt[10];    
    Float_t TLep_Eta[10];
    Float_t TLep_Phi[10];
    Float_t TLep_E[10];
    Float_t TLep_Charge[10];

    Int_t TNJets;            // Jets...
    Int_t TNBtags;
    Float_t TJet_Pt[20];
    Float_t TJet_Eta[20];
    Float_t TJet_Phi[20];
    Float_t TJet_E[20];
    Int_t TJet_isBJet[20];
    Float_t THT;       // HT

    // For systematics...
    Int_t   TNJetsJESUp;
    Int_t   TNJetsJESDown;
    Int_t   TNJetsJER;
    Int_t   TNBtagsUp;
    Int_t   TNBtagsDown;
    Int_t   TNBtagsMisTagUp;
    Int_t   TNBtagsMisTagDown;
    Float_t TJetJESUp_Pt[20];
    Float_t TJetJESDown_Pt[20];
    Float_t TJetJER_Pt[20];
    Float_t THTJESUp;
    Float_t THTJESDown;

    Int_t   TNISRJets;
    Float_t TMETJESUp;
    Float_t TMETJESDown;
    Float_t TMT2llJESUp;
    Float_t TMT2llJESDown;

    Float_t  TWeight_LepEffUp;
    Float_t  TWeight_LepEffDown;
    Float_t  TWeight_TrigUp;
    Float_t  TWeight_TrigDown;
    Float_t  TWeight_FSUp;
    Float_t  TWeight_PUDown;
    Float_t  TWeight_PUUp;
    Float_t  TWeight_FSDown;

  protected:
    Bool_t  gIsData;
    Bool_t  gDoSyst;
    Int_t   gSelection;
    TString gSampleName;

    ClassDef(CalcEfficiencies, 0);
};
