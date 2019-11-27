#ifndef LEPTONSF_H
#define LEPTONSF_H 1

#include <iostream>
#include "Lepton.h"
#include "Functions.h"
#include "TFile.h"
#include "TH2F.h"
#include "TMath.h"
#include "TGraphAsymmErrors.h"
#include "TSystem.h"
#include "TString.h"

//const TString path_to_SF_histos = gSystem->WorkingDirectory() + TString("/InputFiles/");

const Float_t lumiBCDEF = 19.706;
const Float_t lumiGH    = 16.1454;

class SUSYnorm {
 public:
  SUSYnorm(TString path, TString filename);
  ~SUSYnorm() {}
  void loadHisto();
  //Double_t GetStopXSec(Float_t mStop);
  Double_t GetStopXSec(Int_t mStop);
  Double_t GetSUSYnorm(Float_t mStop, Float_t mLSP);
  
  vector<TString> GetAllFiles(TString path, TString  filename = "SMS_T2tt_3J_xqcut_20_top_corridor_2Lfilter_TuneCP5_MLM_p", Bool_t verbose =1);
  TH1* GetHistoFromFiles(vector<TString> Files, TString histoName);
  TH2* GetHistoFromFiles2(vector<TString> Files, TString histoName);
  
 private:
  // Muon SFs
  TH2D*  fhSMS;
  TString path;
  TString filename;
  TString hname;
};
#endif
