#ifndef JET 
#define JET 1

#include <iostream>
#include "TLorentzVector.h"
#include "Lepton.h"

class Jet{
  public:
    Jet(){};
    Jet(TLorentzVector vec, Float_t btag_csv, Int_t Id = -1, Int_t mcFlavour = 0, Float_t btag_deepcsv = -1){
      p = vec;
      csv = btag_csv;
      deepcsv = btag_deepcsv;
      id = Id;
      flavmc = mcFlavour;
      //InitSyst();
    }
    ~Jet(){};

    Bool_t isBtag;
    TLorentzVector p;
    TLorentzVector mcp;
    Int_t id;
    Int_t flavmc;
    Float_t csv; 
    Float_t deepcsv; 
    Float_t deepcsvC;
    Float_t deepflav;
    Float_t genPt;

    // For systematics
    Float_t pTJESUp;
    Float_t pTJESDown;
    Float_t pTJERUp;
    Float_t pTJERDown;
    Bool_t  isBtag_BtagUp;
    Bool_t  isBtag_BtagDown;
    Bool_t  isBtag_MisTagUp;
    Bool_t  isBtag_MisTagDown;

    void InitSyst();
    Float_t Pt(){ return p.Pt();}
    Float_t Eta(){ return p.Eta();}
    Float_t Phi(){ return p.Phi();}
    Float_t E(){ return p.E();}
    void SetMCjet(TLorentzVector p){ mcp = p;}
    void SetGenPt(Float_t pt){ genPt = pt;}
    Float_t GetGenPt(){ return genPt;}
    void SetDeepCSVB(float v);
    void SetDeepCSVC(float v);
    void SetDeepFlav(float v);
    Float_t GetDeepCSVB();
    Float_t GetDeepCSVC();
    Float_t GetDeepFlav();
    void SetIsBtag(Bool_t is = true, Int_t dir = 0){
      if     (dir ==  0) isBtag            = is;
      else if(dir ==  1) isBtag_BtagUp     = is;
      else if(dir == -1) isBtag_BtagDown   = is;
      else if(dir ==  2) isBtag_MisTagUp   = is;
      else if(dir == -2) isBtag_MisTagDown = is;
    }
    Bool_t IsBtag(Float_t dir = 0){
      if     (dir ==  0) return isBtag;
      else if(dir ==  1) return isBtag_BtagUp;
      else if(dir == -1) return isBtag_BtagDown;
      else if(dir ==  2) return isBtag_MisTagUp;
      else if(dir == -2) return isBtag_MisTagDown;
      return isBtag;
    }
};

Int_t GetNJets(std::vector<Jet> jets);
Int_t GetBtags(std::vector<Jet> jets, Int_t dir = 0);
Float_t GetHT(std::vector<Jet> jets);
Float_t GetDeltaHTmiss(TLorentzVector jet, std::vector<Jet> jets, std::vector<Lepton> leptons);

#endif
