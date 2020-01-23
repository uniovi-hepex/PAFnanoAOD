#pragma once

#include "PAFChainItemSelector.h"
#include "Lepton.h"
#include "Jet.h"
#include "BTagSFUtil.h"
#include "Functions.h"
#include <iostream>
#include <vector>
#include "LinkDef.h"

class JetSelector : public PAFChainItemSelector{

  public:
    TString taggerName;
    TString stringWP;
    JetSelector();
    virtual ~JetSelector();
    virtual void InsideLoop();
    void Initialise();
    void Summary();

    std::vector<Lepton> Leptons;
    std::vector<Lepton> VetoLeptons;
    std::vector<Jet> selJets;
    std::vector<Jet> selJetsJecCorUp;
    std::vector<Jet> selJetsJecCorDown;
    std::vector<Jet> selJetsJecUnCorUp;
    std::vector<Jet> selJetsJecUnCorDown;
    std::vector<Jet> selJetsJecUp;
    std::vector<Jet> selJetsJecDown;
    std::vector<Jet> selJetsJERUp;
    std::vector<Jet> selJetsJERDown;
    std::vector<Jet> Jets15;
    std::vector<Jet> genJets;
    std::vector<Jet> mcJets;
    std::vector<Jet> vetoJets;

    Int_t nSelJets;
    Int_t nGenJets;
    Int_t nVetoJets;
    Int_t nJets15;
    Int_t nBtagJets;
    Int_t nBtagJetsJECUp;
    Int_t nBtagJetsJECDown;
    Int_t nBtagJetsJER;

    Float_t minDR;
    Float_t jet_MaxEta;
    Float_t jet_MinPt;
    Float_t vetoJet_minPt;
    Float_t vetoJet_maxEta;

    Float_t MET_JERUp;
    Float_t MET_JESUp;
    Float_t MET_JESDown;

    Float_t BtagSF;
    Float_t BtagSFBtagUp;
    Float_t BtagSFBtagDown;
    Float_t BtagSFMistagUp;
    Float_t BtagSFMistagDown;

  protected:

    BTagSFUtil *fBTagSFnom ;
    BTagSFUtil *fBTagSFbUp ;
    BTagSFUtil *fBTagSFbDo ;
    BTagSFUtil *fBTagSFlUp ;
    BTagSFUtil *fBTagSFlDo ;

    string MeasType;
    UInt_t evt;

    TString selection;
    Bool_t gIsData;
    Bool_t gDoSys;
    Bool_t gIsFastSim;
    Int_t gSelection;
    TString gSampleName;
    Bool_t gIsFSRUp;
    Bool_t gIsFSRDown;
    Bool_t gIs2017eccop;
    TString gOptions;

    // Jets
    Int_t nJet;
    Jet tJ;
    TLorentzVector tpJ; 
    TLorentzVector tmcJ; 
    Int_t jetId;
    Int_t flavmc;
    Int_t flav;
    Float_t csv;
    Float_t deepcsv;
    Float_t deepcsvC;
    Float_t deepflav;
    Float_t pt; Float_t eta; Float_t rawPt; 
    Float_t pt_corrUp; Float_t pt_corrDown; 
    // genJet
    Int_t ngenJet;
    Float_t BtagSFFS;

    Float_t rho;
    Bool_t gIs2016;
    Bool_t gIs2017;
    Bool_t gIs2018;
    Int_t  year;

    void GetJetVariables(Int_t i);
    void GetGenJetVariables(Int_t i);
    Bool_t IsBtag(Jet j); 
    void SetSystematics(Jet *j);
//    Bool_t Cleaning(Jet j, vector<Lepton> vLep, Float_t minDR = 0.4);

    ClassDef(JetSelector, 0);
};

