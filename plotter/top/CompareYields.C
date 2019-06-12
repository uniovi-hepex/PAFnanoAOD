R__LOAD_LIBRARY(Histo.C+)
R__LOAD_LIBRARY(Looper.C+)
R__LOAD_LIBRARY(Plot.C+)
R__LOAD_LIBRARY(TopHistoReader.C+)

#include "Histo.h"
#include "Looper.h"
#include "Plot.h"

Float_t Lumi = 41.2;
TString chan = "ElMu";
TString level = "1btag";
TString var = "Lep0Eta";

TString sbkg = "TTWJetsToLNu_madspin,TTWJetsToQQ_madspin, TTZToLL_M_1to10, TTZToLLNuNu_M_10, TTZToQQ, WW, ZZ, WZ, DYJetsToLL_M_50, tbarW_noFullHad, tW_noFullHad, TTTo2L2Nu, TTToSemiLeptonic";
TString sdata = "MuonEG, SingleMuon, SingleElectron, DoubleEG, DoubleMuon";
TString sera  = "B,E,F";
TString schannels = "ElMu, Elec, Muon";
vector<TString> bkg  = TStringToVector(sbkg);   Int_t nBkg = bkg.size();
vector<TString> data = TStringToVector(sdata); Int_t nData = data.size();
vector<TString> era  = TStringToVector(sera); Int_t nEras = era.size();
vector<TString> channels = TStringToVector(schannels); Int_t nChan = channels.size();

TString path =  "/pool/ciencias/userstorage/juanr/top/2017/sep03/";
TopHistoReader *t = new TopHistoReader(path);

Float_t GetLumiForEra(TString era);
Float_t GetBkgYield(TString era, TString ch = "");
Float_t GetDataYield(TString era, TString ch = "");
TString FixSpaces(TString t, Int_t n = 10);
TString PrintCompareYields(TString era, TString chan = "");
void PrintDummy(TString era, TString dataset, TString ch = "");

Float_t GetBkgYield(TString era, TString ch){
  t->SetPath(path);
  t->SetLevel(level);
  t->SetVar(var); t->SetRebin(1);
  if(ch == "") ch = chan; t->SetChan(ch);

  t->SetIsData(false);
  t->SetLumi(GetLumiForEra(era)*1000);
  
  Float_t y = 0;
  TString pr = "";
  for(Int_t i = 0; i < nBkg; i++)  y += t->GetYield(bkg.at(i), "", "", "", var);
  return y;
}

Float_t GetDataYield(TString era, TString ch){
  //t->SetPath("/nfs/fanae/user/juanr/nanoAOD/Top_temp/");//path+"/data/");
  t->SetPath(path);//+"/perEra/");
  t->SetLevel(level);
  t->SetVar(var); t->SetRebin(1);
  if(ch == "") ch = chan; t->SetChan(ch);

  t->SetIsData(true);
  t->SetLumi(GetLumiForEra(era)*1000);
  
  Float_t y = 0;
  TString pr = "";
  for(Int_t i = 0; i < nData; i++)  y += t->GetYield(data.at(i)+"_Run2017"+era, "", "", "",var);
  return y;
}


Float_t GetLumiForEra(TString era){
  if(era == 'B') return 4.823;
  if(era == 'C') return 9.664;
  if(era == 'D') return 4.252;
  if(era == 'E') return 9.278;
  if(era == 'F') return 13.54;
  return 0;
}

TString FixSpaces(TString t, Int_t n){
  while(t.Sizeof()-1 < n){
    t += " ";
  }
  return t;
}

TString PrintCompareYields(TString era, TString chan){
  float ybkg = GetBkgYield(era, chan);
  float ydat = GetDataYield(era, chan);
  float r = fabs(ybkg-ydat)/ybkg;
  float ratio = ydat/ybkg;
  TString out = Form("\033[1;32mChannel \033[0;32m%s\033[0m, \033[1;34mEra \033[0;34m%s\033[0m, [\033[1;35mdata\033[0m, \033[1;36mbkg]\033[0m = [\033[1;35m%1.1f\033[0m, \033[1;36m%1.1f\033[0m]", t->GetChan().Data(), era.Data(), ydat, ybkg);
  out += Form(", data/MC = %1.2f\n",ratio);
  TString a = Form("%1.2f, %1.2f", ydat, ybkg);

  //if     (r < 0.05) out += Form(" (%1.2f)\n",ratio);
  //else if(r < 0.10) out += " \033[1;33mWARNING (" + TString(Form("%1.2f",ratio)) + ")\n\033[0m";
  //else                      out += " \033[1;31mWARNING (" + TString(Form("%1.2f",ratio)) + ")\n\033[0m";
  cout << out;
  return a;
}

void PrintDummy(TString era, TString dataset, TString ch){
//  cout << "#### Dataset: " << dataset + " " + era << ", channel: " << ch << endl;
  t->SetPath(path+"/data/");
  t->SetLevel(level);
  t->SetVar(var); t->SetRebin(1);
  if(ch == "") ch = chan; t->SetChan(ch);
  t->SetIsData(true);
  TString sdlevels = ",2leps,trigger,metfilters,OS,minv,lep0pt,mz,met,njets,nbtags"; 
  vector<TString> dlevels = TStringToVector(sdlevels); Int_t nDLev = dlevels.size();
  TString lev; TString name;
  TString datatag = dataset+" "+era+" "+ch; datatag = datatag.ReplaceAll("Electron", "Elec").ReplaceAll("Double", "Dou").ReplaceAll("Single", "Sin");
  TString outnum = FixSpaces(datatag,14) + " |";
  TString outlev = FixSpaces("",14) + " |";
  Float_t norm = t->GetNamedHisto("fhDummy_"+ ch, dataset + "_Run2017" + era)->GetEntries();
  for(int i = 1; i < nDLev; i++){
    lev = dlevels.at(i);
    name = "fhDummy_" + lev + ch;
    outlev += FixSpaces(lev) + " | ";
    outnum += FixSpaces(Form("%1.2f", t->GetNamedHisto(name, dataset + "_Run2017" + era)->GetEntries()/norm*100)) + " | ";
  }
  cout << outlev << endl;
  cout << outnum << endl;
}

void CompareYields(){
  TString iCh = ""; TString iEr = "";
  for(Int_t j = 0; j < nChan; j++){
    iCh = channels.at(j); 
    /*for(Int_t k = 0; k < nData; k++){
      TString dataset = data.at(k);
      if(iCh == "ElMu" && (dataset.BeginsWith("Double"))) continue;
      if(iCh == "Elec" && (!dataset.Contains("Elec") && !(dataset == "DoubleEG"))) continue;
      if(iCh == "Muon" && (!dataset.Contains("Muon") || dataset == "MuonEG")) continue;
      //PrintDummy("E", dataset, iCh);
    }*/
    Float_t totData = 0; Float_t totMC = 0;
    TString v;
    vector<Float_t> f;
    for(Int_t i = 0; i < nEras; i++){
      iEr = era.at(i);
      path = Form("/pool/ciencias/userstorage/juanr/top/2017/nov22/run%s/", iEr.Data());
      v = PrintCompareYields(iEr, iCh);
      f = TStringToFloat(v);
      totData += f[0];
      totMC   += f[1];
    }
    cout << Form("Channel = %s, [data, bkg] = [%1.2f, %1.2f], data/MC = %1.2f\n", iCh.Data(), totData, totMC, totData/totMC);
  }
}

