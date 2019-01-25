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

