#include "Jet.h"

void Jet::InitSyst(){
  pTJESUp   = p.Pt();
  pTJESDown = p.Pt();
  pTJERUp   = p.Pt();
  pTJERDown = p.Pt();
  isBtag_BtagUp     = isBtag;
  isBtag_BtagDown   = isBtag;
  isBtag_MisTagUp   = isBtag;
  isBtag_MisTagDown = isBtag;
}

void Jet::SetDeepCSVB(float v){ deepcsv  = v;}
void Jet::SetDeepCSVC(float v){ deepcsvC = v;}
void Jet::SetDeepFlav(float v){ deepflav = v;}
Float_t Jet::GetDeepCSVB(){ return deepcsv;}
Float_t Jet::GetDeepCSVC(){ return deepcsvC;}
Float_t Jet::GetDeepFlav(){ return deepflav;}


Int_t GetNJets(std::vector<Jet> jets){
  return (int) jets.size();
}

Int_t GetBtags(std::vector<Jet> jets, Int_t dir){
  Int_t nb = 0;
  Int_t njets = GetNJets(jets);
  for(Int_t i = 0; i < njets; i++){
    if(jets.at(i).IsBtag(dir)) nb++;
  }
  return nb;
}

Float_t GetHT(std::vector<Jet> jets){
  Float_t ht = 0;
  Int_t njets = GetNJets(jets);
  for(Int_t i = 0; i < njets; i++) ht += jets.at(i).Pt(); 
  return ht;
}

