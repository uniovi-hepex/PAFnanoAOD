#pragma once

#include "PAFChainItemSelector.h"
#include "Functions.h"
#include "mt2.h"
#include "BTagSFUtil.h"
#include <iostream>
#include <vector>
#include "SUSYnorm.h"
//enum eChannels{iUnkChan, iElMu, iMuon, iElec, nChannels};
const Int_t nChannels = 3;
enum eLevels  {idilepton, iZVeto, iMETcut, i2jets, i1btag, nLevels};
const TString gChanLabel[nChannels] = {"ElMu", "Muon","Elec"};
const TString sCut[nLevels] = {"dilepton", "ZVeto", "MET", "2jets", "1btag"};
const Int_t nPtBins = 14;
const Float_t ptBins[nPtBins+1] = {30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 200, 300, 400, 600};

enum eSysts                   {kNorm, kMuonEffUp,  kMuonEffDown,  kElecEffUp,  kElecEffDown,  kMuonEnergyUp,  kMuonEnergyDown,  kElecEnergyUp,  kElecEnergyDown,  kJESUp,  kJESDown,  kJERUp,  kJERDown,  kPUUp,  kPUDown, kTrigUp, kTrigDown,   kUnclMETUp,   kUnclMETDown, kBtagUp,  kBtagDown,  kMistagUp,  kMistagDown, kISRUp, kISRDown, kFSRUp, kFSRDown, kPrefireUp, kPrefireDown, kSUSYISRUp, kSUSYISRDown, kTopPt, nSysts};
const TString gSyst[nSysts] = {"",    "MuonEffUp", "MuonEffDown", "ElecEffUp", "ElecEffDown", "MuonEnergyUp", "MuonEnergyDown", "ElecEnergyUp", "ElecEnergyDown", "JESUp", "JESDown", "JERUp", "JERDown", "PUUp", "PUDown", "TrigUp", "TrigDown", "UnclMETUp", "UnclMETDown", "BtagUp", "BtagDown", "MistagUp", "MistagDown", "ISRUp", "ISRDown", "FSRUp", "FSRDown", "PrefireUp", "PrefireDown"};


class TopAnalysis : public PAFChainItemSelector{
  public:
    TopAnalysis();
    virtual ~TopAnalysis(){}
    virtual void InsideLoop();
    virtual void Initialise();
    virtual void Summary();
    std::vector<Lepton> genLeptons  ;
    std::vector<Lepton> selLeptons  ;
    std::vector<Lepton> vetoLeptons ;
    std::vector<Jet> selJets ;
    std::vector<Jet> selJetsJecUp   ;
    std::vector<Jet> selJetsJecDown ;
    std::vector<Jet> selJetsJERUp   ;
    std::vector<Jet> selJetsJERDown ;
    std::vector<Jet> selJetsJecCorUp   ;
    std::vector<Jet> selJetsJecCorDown ;
    std::vector<Jet> selJetsJecUnCorUp   ;
    std::vector<Jet> selJetsJecUnCorDown ;
    Int_t JERindex;

    std::vector<Jet> Jets15  ;
    std::vector<Jet> genJets  ;
    std::vector<Jet> mcJets  ;
    std::vector<Jet> vetoJets;

    std::vector<Double_t> SumOfPDFweights;
    std::vector<Double_t> SumOfMEweights;

    BTagSFUtil *fBTagSFnom;
    BTagSFUtil *fBTagSFbUp;
    BTagSFUtil *fBTagSFbDo;
    BTagSFUtil *fBTagSFlUp;
    BTagSFUtil *fBTagSFlDo;
    Int_t era;
    Int_t nPDFweights;
    Int_t nMEweights = 9;
    Int_t   TNVert;
    Float_t   TNVert_pu;
    TTree* fTree;
    TString GetSuffix(int iCh, int iCut, int iSyst = 0);
    void SetLeptonVariables();
    void SetJetVariables();
    void SetEventVariables();
    void SetVariables(int sys = 0);

    Bool_t makeHistos = false;
    Bool_t makeTree   = false;
    Bool_t miniTree   = false;
    vector<int> useSyst = vector<int>();
    Int_t nSyst;

    void GetLeptonVariables(std::vector<Lepton> selLeptons, std::vector<Lepton> VetoLeptons);
    void GetJetVariables(std::vector<Jet> selJets, std::vector<Jet> cleanedJets15, Float_t ptCut = 30);
    void GetGenJetVariables(std::vector<Jet> genJets, std::vector<Jet> mcJets);
    void GetMET();
    void GetWeights();
    Int_t nFiduJets; Int_t nFidubJets; 
    
    Float_t TMuonSF;
    Float_t TElecSF;
    Float_t TrigSF;
    Float_t TrigSFerr;
    Float_t PUSF;
    Float_t PUSF_Up;
    Float_t PUSF_Down;
    Float_t PrefWeight;
    Float_t PrefWeightUp;
    Float_t PrefWeightDo;
    Int_t   gChannel;
    Bool_t  TPassMETFilters;
    Bool_t  TPassTrigger;
    Int_t  TIsHEM;
    Bool_t  isSS;
    Double_t NormWeight;

    void InitHistos();
    void FillDYHistos(Int_t ch);
    void FillHistos(Int_t ch, Int_t cut, Int_t sys);
    void FillCorrHistos();
  
    void get20Jets();
    Float_t GetTopPtWeight(Float_t Pt1, Float_t Pt2);

    Double_t getDilepMETJetPt(const TString& sys = "Norm");
    Double_t getDilepJetPt(const TString& sys = "Norm");
    Double_t getLep1METJetPt(const TString& sys = "Norm");
    Double_t getPtSys(TLorentzVector*, int);
    Double_t getDilepMETJet1Pz(const TString& sys = "Norm");
    Double_t getPzSys(TLorentzVector*, int);
    Double_t getDPtDilep_JetMET(const TString& sys = "Norm");
    Double_t getDPtDilep_MET(const TString& sys = "Norm");
    Double_t getDPtLep1_MET(const TString& sys = "Norm");
    Double_t getDeltaPt(vector<TLorentzVector>, vector<TLorentzVector>);
    Double_t getSysM(const TString& sys = "Norm");
    Double_t getM(vector<TLorentzVector>);

    Int_t CountISRjets();
    Float_t GetISRweight2016(Int_t n);
    Float_t GetISRweight1718(Int_t n);

    Float_t metcut = 50.;
    Float_t mt2cut = 80.;

    //Variables
    ULong64_t event;
    UInt_t    lumiblock;
    Double_t TWeight;   // Total nominal weight
    Double_t TWeight_noPU;
    Float_t TMll;      // Invariant mass
    Float_t TMllMuonUp; Float_t TMllMuonDo; Float_t TMllElecUp; Float_t TMllElecDo;
    Float_t TDilep_Pt;
    Float_t TDeltaEta;
    Float_t TDeltaPhi;
    Float_t TMuonPt;
    Float_t TMuonEta;
    Float_t TMuonPhi;
    Float_t TMuonM;
    Float_t TElecPt;
    Float_t TElecEta;
    Float_t TElecPhi;
    Float_t TElecM;
    Float_t TMuonPtUp;
    Float_t TMuonPtDo;
    Float_t TElecPtUp;
    Float_t TElecPtDo;
    UInt_t  TRun;
    //Int_t   TNVert;
    Float_t TMET;      // MET
    Float_t TMT2;      // MET
    Float_t TMT2lblb; 
    Float_t TMT2orig;
    Float_t TMETorig;
    Float_t TMETpuppi;      // MET
    Float_t TMT2puppi;      // MET
    Float_t TMETpuppi_Phi;  // MET phi

    Float_t TMETPhi;  // MET phi
	Float_t TMETsig; //MET significance
    Float_t TgenTop1Pt = 0;
    Float_t TgenTop2Pt = 0;
    Int_t  TIsOnZ;
    Int_t  TPassDilep;
    Int_t  TPassDilepMuonESUp;
    Int_t  TPassDilepMuonESDo;
    Int_t  TPassDilepElecESUp;
    Int_t  TPassDilepElecESDo;
    Bool_t TPassDilepAny;
    Bool_t TPassJetsAny;
    Bool_t TPassBtagAny;
    Bool_t TPassMETAny;
    Bool_t TPassMT2Any;

    Int_t   TNVetoLeps;
    Int_t   TNSelLeps;
    Int_t   TChannel;
    Int_t   TStatus;
    Int_t   TLep0IsPrompt;
    Int_t   TLep1IsPrompt;
    Int_t   TLep0IsConversion;
    Int_t   TLep1IsConversion;
    Int_t   TIsSS;
    Float_t TLep0Pt;    
    Float_t TLep0Eta;
    Float_t TLep0Phi;
    Float_t TLep0M;
    Int_t   TLep0Id;
    Float_t TLep0Iso;
    Float_t TLep1Pt;    
    Float_t TLep1Eta;
    Float_t TLep1Phi;
    Float_t TLep1M;
    Int_t   TLep1Id;
    Float_t TLep1Iso;

    Float_t TGenLep0Pt;    
    Float_t TGenLep0Eta;
    Float_t TGenLep0Phi;
    Float_t TGenLep1Pt;    
    Float_t TGenLep1Eta;
    Float_t TGenLep1Phi;
    Float_t TGenMET;
    Float_t TGenMET_phi;
    Float_t TGenMT2;

    Float_t TTop0Pt;
    Float_t TTop0Eta;
    Float_t TTop0Phi;
    Float_t TTop1Pt;
    Float_t TTop1Eta;
    Float_t TTop1Phi;

    Int_t TNJets;            // Jets...
    Int_t TNFwdJets; 
    Int_t TNBtags;
    Float_t TJet_Pt[20];
    Float_t TJet_Eta[20];
    Float_t TJet_Phi[20];
    Float_t TJet_M[20];
    Float_t TJet_Csv[20];
    Int_t TNBtagsLoose;
    Float_t TJet0Pt;
    Float_t TJet0Eta;
    Float_t TJet0Csv;

    Float_t TJet0Phi; 
    Float_t TJet0M; 
    Int_t TJet0IsBTag;
    Float_t  TJet1Pt; 
    Float_t TJet1Eta; 
    Float_t TJet1Phi; 
    Float_t TJet1M; 
    Float_t TJet1Csv; 
    Float_t TJet1IsBTag;

    Float_t THT;       // HT

    // For systematics...
    Int_t   TNJetsJESUp;
    Int_t   TNJetsJESDown;
    Int_t   TNJetsJERUp;
    Int_t   TNJetsJERDown;
    Int_t   TNJetsJESCorUp;
    Int_t   TNJetsJESCorDown;
    Int_t   TNJetsJESUnCorUp;
    Int_t   TNJetsJESUnCorDown;
    Int_t   TNBtagsBtagUp;
    Int_t   TNBtagsBtagDown;
    Int_t   TNBtagsMisTagUp;
    Int_t   TNBtagsMisTagDown;
    Int_t   TNBtagsJESUp;
    Int_t   TNBtagsJESDown;
    Float_t TJetJESUp_Pt[20];
    Float_t TJetJESDown_Pt[20];
    Float_t TJetJER_Pt[20];
    Float_t THTJESUp;
    Float_t THTJESDown;
    Float_t THTJESCorUp;
    Float_t THTJESCorDown;
    Float_t THTJESUnCorUp;
    Float_t THTJESUnCorDown;
    Float_t THTJERUp;
    Float_t THTJERDown;
    Float_t TBtagPt;
    Float_t TJet0PtJESUp;
    Float_t TJet0PtJESDown;
    Float_t TJet0PtJERUp;
    Float_t TJet0PtJERDown;
    Float_t TJet0PtJESCorUp;
    Float_t TJet0PtJESCorDown;
    Float_t TJet0PtJESUnCorUp;
    Float_t TJet0PtJESUnCorDown;

    Int_t   TNISRJets;
    Float_t TMETJESUp; Float_t TMETJESDown; Float_t TMETJERUp; Float_t TMETJERDown; Float_t TMETUnclUp; Float_t TMETUnclDown;
    Float_t TMETMuonESUp; Float_t TMETMuonESDown; Float_t TMETElecESUp; Float_t TMETElecESDown;
    Float_t TMETPhiJESUp; Float_t TMETPhiJESDown; Float_t TMETPhiJERUp; Float_t TMETPhiJERDown; Float_t TMETPhiUnclUp; Float_t TMETPhiUnclDown;
    Float_t TMETPhiMuonESUp; Float_t TMETPhiMuonESDown; Float_t TMETPhiElecESUp; Float_t TMETPhiElecESDown;
    Float_t TMT2JESUp; Float_t TMT2JESDown; Float_t TMT2JERUp; Float_t TMT2JERDown; Float_t TMT2UnclUp; Float_t TMT2UnclDown;
    Float_t TMT2MuonESUp; Float_t TMT2MuonESDown; Float_t TMT2ElecESUp; Float_t TMT2ElecESDown;
    Float_t TMT2JESCorUp; Float_t TMT2JESCorDown; Float_t TMT2JESUnCorUp; Float_t TMT2JESUnCorDown; Float_t TMETJESCorUp; Float_t TMETJESCorDown; Float_t TMETJESUnCorUp; Float_t TMETJESUnCorDown; Float_t TMETPhiJESCorUp; Float_t TMETPhiJESCorDown; Float_t TMETPhiJESUnCorUp; Float_t TMETPhiJESUnCorDown;

    
    Float_t m_stop; Float_t m_LSP; Int_t sumWeights; Int_t ngenPart; Int_t iSt; Int_t iLSP; Double_t xsec; Int_t j;Float_t norm;//quitar
    Float_t  TWeightRaw;
    Float_t  TWeight_LepEffUp;
    Float_t  TWeight_LepEffDown;
    Float_t  TWeight_ElecEffUp;
    Float_t  TWeight_ElecEffDown;
    Float_t  TWeight_MuonEffUp;
    Float_t  TWeight_MuonEffDown;
    Float_t  TWeight_TrigUp;
    Float_t  TWeight_TrigDown;
    Float_t  TWeight_PUDown;
    Float_t  TWeight_PUUp;
    Float_t  TWeight_ISRUp;
    Float_t  TWeight_ISRDown;
    Float_t  TWeight_FSRUp;
    Float_t  TWeight_FSRDown;
    Float_t  TWeight_PrefUp;
    Float_t  TWeight_PrefDown;
    Float_t  TWeight_TopPtUp;
    Float_t  TWeight_TopPtDown;
    Float_t  TWeight_ME[9];
    Float_t  TWeight_PDF[110];
    Float_t  TWeight_SUSYISRUp;
    Float_t  TWeight_SUSYISRDown;
    Float_t  TWeight_SUSYISR;
    
    std::vector<Jet> jets;
    Double_t weight;
    Float_t nvert_pu;
    Float_t met, ht, nvert,mt2, invmass; 
    Int_t   nleps, njets, nbtags;
    Float_t lep0pt, lep1pt, lep0eta, lep1eta, lep0iso, lep1iso;
    Float_t dileppt, deltaphi, deltaeta, deltaR;

    TString metvar;
    TString metvarpt;
    TString metvarphi;
    TString JetPt;
    Bool_t gPUWeigth;
    Bool_t gPrefire;
    Bool_t gDoJECunc;
    Bool_t gDoMuonES;
    Bool_t gDoElecES;
    Bool_t gDoPDFunc;
    Bool_t gDoPSunc;
    Bool_t gDoScaleUnc;

// Histograms
//=====================================================0
  TH1F* fHPDFweights[nChannels][nLevels];
  TH1F* fHPSweights[nChannels][nLevels];
  TH1F* fHScaleWeights[nChannels][nLevels];
  TH1F* fHMET[nChannels][nLevels][nSysts];
  TH1F* fHMT2[nChannels][nLevels][nSysts];
  TH1D* fHLep0Eta[nChannels][nLevels][nSysts];
  TH1F* fHLep1Eta[nChannels][nLevels][nSysts];
  TH1F* fHMuonEta[nChannels][nLevels][nSysts];
  TH1F* fHElecEta[nChannels][nLevels][nSysts];
  TH1F* fHDelLepPhi[nChannels][nLevels][nSysts];
  TH1F* fHDelLepEta[nChannels][nLevels][nSysts];
  TH1F* fHDelLepR[nChannels][nLevels][nSysts];
  TH1F* fHR9l0[nChannels][nLevels][nSysts];
  TH1F* fHR9l1[nChannels][nLevels][nSysts];
  TH1F* fHHT[nChannels][nLevels][nSysts];
  TH1F* fHJet0Eta[nChannels][nLevels][nSysts];
  TH1F* fHJet1Eta[nChannels][nLevels][nSysts];
  TH1F* fHJetEta[nChannels][nLevels][nSysts];
  TH1F* fHJetBtagEta[nChannels][nLevels][nSysts];
  TH1F* fHJetBtag0Eta[nChannels][nLevels][nSysts];

  TH1F* fHDiLepPt[nChannels][nLevels][nSysts];
  TH1F* fHLep0Pt[nChannels][nLevels][nSysts];
  TH1F* fHLep1Pt[nChannels][nLevels][nSysts];
  TH1F* fHLep0Iso[nChannels][nLevels][nSysts];
  TH1F* fHLep1Iso[nChannels][nLevels][nSysts];
  TH1F* fHMuonPt[nChannels][nLevels][nSysts];
  TH1F* fHElecPt[nChannels][nLevels][nSysts];
  TH1F* fHMuonIso[nChannels][nLevels][nSysts];
  TH1F* fHElecIso[nChannels][nLevels][nSysts];
  TH1F* fHJetBtag0Pt[nChannels][nLevels][nSysts];
  TH1F* fHJetBtagPt[nChannels][nLevels][nSysts];
  TH1F* fHJetPt[nChannels][nLevels][nSysts];
  TH1F* fHJet0Pt[nChannels][nLevels][nSysts];
  TH1F* fHJet1Pt[nChannels][nLevels][nSysts];
  TH1F* fHNJets[nChannels][nLevels][nSysts];
  TH1F* fHNBtagJets[nChannels][nLevels][nSysts];

  TH1F* fHDYInvMass[nChannels][nLevels][nSysts];
  TH1F* fHDYInvMassSF[nChannels][nLevels][nSysts];
  TH1F* fHInvMass[nChannels][nLevels][nSysts];
  TH1F* fHInvMass2[nChannels][nLevels][nSysts];
  TH1F* fHInvMass2BB[nChannels][nLevels][nSysts];
  TH1F* fHInvMass2EE[nChannels][nLevels][nSysts];
  TH1F* fHInvMass2BE[nChannels][nLevels][nSysts];
  TH1F* fHInvMass2EB[nChannels][nLevels][nSysts];
  TH1F* fHNBtagsNJets[nChannels][nLevels][nSysts];
  TH1F* fHJetCSV[nChannels][nLevels][nSysts];
  TH1F* fHJet0CSV[nChannels][nLevels][nSysts];
  TH1F* fHJet1CSV[nChannels][nLevels][nSysts];
  TH1F* fHJetDeepCSV[nChannels][nLevels][nSysts];
  TH1F* fHJet0DeepCSV[nChannels][nLevels][nSysts];
  TH1F* fHJet1DeepCSV[nChannels][nLevels][nSysts];
  TH1F* fHJetDeepFlav[nChannels][nLevels][nSysts];
  TH1F* fHJet0DeepFlav[nChannels][nLevels][nSysts];
  TH1F* fHJet1DeepFlav[nChannels][nLevels][nSysts];
  TH1F* fHvertices[nChannels][nLevels][nSysts]; 
  TH1F* fHvertices_pu[nChannels][nLevels][nSysts]; 

  TH1F* fhDummy;
  TH1F* fhDummyCh[nChannels];
  TH1F* fhDummy_2leps[nChannels];
  TH1F* fhDummy_trigger[nChannels];
  TH1F* fhDummy_metfilters[nChannels];
  TH1F* fhDummy_OS[nChannels];
  TH1F* fhDummy_minv[nChannels];
  TH1F* fhDummy_lep0pt[nChannels];
  TH1F* fhDummy_mz[nChannels];
  TH1F* fhDummy_met[nChannels];
  TH1F* fhDummy_njets[nChannels];
  TH1F* fhDummy_nbtags[nChannels];
  TH1F*  fHWeightsFidu;
  TH1D*  fHsumwisr;
  TH1D*  fHyields[nChannels][nSysts];
  TH1F*  fHFiduYields[nChannels][nSysts];
  TH1F*  fHSSyields[nChannels][nSysts];

  TH1F* hJetPtReco;
  TH1F* hJetPtGen;
  TH1F* hJetPtRecoB;
  TH1F* hJetPtGenB;
  TH1F* hJetGenRecoPtRatio[nPtBins];
  TH1F* hJetGenRecoPtRatio2[nPtBins];
  Bool_t gIs2016;
  Bool_t gIs2017;
  Bool_t gIs2018;
  protected:

    Int_t year;
    TString selection;
    Bool_t  gIsData;
    Bool_t  gDoSyst;
    Int_t   gSelection;
    TString gSampleName;
    TString sampString;
    
    TString gOptions;
    Bool_t  gIsTTbar;
    Bool_t  gIsSignal; //quitar
    SUSYnorm* snorm;
    Bool_t  gIsTTany;
    Bool_t  gIsLHE;

    ClassDef(TopAnalysis, 0);
};
