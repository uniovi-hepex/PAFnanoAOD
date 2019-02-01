R__LOAD_LIBRARY(Histo.C+)
R__LOAD_LIBRARY(Looper.C+)
R__LOAD_LIBRARY(Plot.C+)
R__LOAD_LIBRARY(TopHistoReader.C+)

#include "Histo.h"
#include "Looper.h"
#include "Plot.h"

bool doOnlyDoubleLepDatasets = false;
TString Era = "F";
Float_t Lumi =  0.308545;
TString path =  "/pool/ciencias/userstorage/juanr/top/5TeV/oct10/";
TString plotDir = "/nfs/fanae/user/juanr/www/plots5TeV/temp2/";


//TString chan = "Elec";
//TString level = "dilepton";
TString sChannels    = "ElMu, Elec, Muon";
TString sLevels      = "dilepton, ZVeto, MET, 2jets, 1btag";
TString plotsSel     = "Jets,Leptons,Global";
vector<TString> vChan= TStringToVector(sChannels);
vector<TString> vLev = TStringToVector(sLevels  );
Int_t   nChan        = vChan.size();
Int_t   nLev         = vLev .size();

void DrawStack(TopHistoReader *t, TString var, int rebin = 1, TString xtit = "", TString ytit = "Events", TString subfolder = "");
void draw(TString level = "1btag", TString chan = "ElMu");


void p5TeVplot(){
  int ichan; int ilev; TString chan; TString lev;
  for(ichan = 0; ichan < nChan; ichan++){
    chan = vChan.at(ichan);
    for(ilev = 0; ilev < nLev; ilev++){
      lev = vLev.at(ilev);
      draw(lev, chan);
    }
  }
}

void draw(TString level, TString chan){
  TopHistoReader *t = new TopHistoReader(path);
  t->SetLumi(Lumi*1000);
  TString f = "";

  t->SetLevel(level);
  t->SetChan(chan);
  // Leptons
  if(plotsSel.Contains("Leptons")){
  f = "Leptons";
  DrawStack(t, "Lep0Pt",  100, "Leading lepton p_{T} (GeV)", "Events", f);
  DrawStack(t, "Lep1Pt",  100, "Subleading lepton p_{T} (GeV)", "Events", f);
  DrawStack(t, "Lep0Eta",  5, "Leading lepton #eta", "Events", f);
  DrawStack(t, "Lep1Eta",  5, "Subleading lepton #eta", "Events", f);
  DrawStack(t, "DiLepPt", 160, "p_{T}^{ll} (GeV)", "Events", f);
  DrawStack(t, "DelLepPhi", 10, "#Delta#phi(l1,l2) (rad/#pi)", "Events", f);
  //DrawStack(t, "Lep0Iso",  2, "Leading lepton RelIso", "Events", f);
  //DrawStack(t, "Lep1Iso",  2, "Subleading lepton RelIso", "Events", f);
  //DrawStack(t, "DelLepEta", 4, "#Delta#eta(l1,l2) ", "Events", f);
  if(chan == "ElMu"){
    DrawStack(t, "ElecPt",  100, "Electron p_{T} (GeV)", "Events", f);
    DrawStack(t, "MuonPt",  100, "Muon p_{T} (GeV)", "Events", f);
    DrawStack(t, "ElecEta",  5, "Electron #eta", "Events", f);
    DrawStack(t, "MuonEta",  5, "Muon #eta", "Events", f);
    //DrawStack(t, "ElecIso",  2, "Electron RelIso03", "Events", f);
    //DrawStack(t, "MuonIso",  2, "Muon RelIso04", "Events", f);
  }
  DrawStack(t, "InvMass", 20, "m_{ll} (GeV)", "Events", f);
  if(chan != "ElMu") DrawStack(t, "InvMass2", 2, "m_{ll} (GeV)", "Events", f);
  }

  // Global
  if(plotsSel.Contains("Global")){
  f = "Global";
  DrawStack(t, "NJets",  1, "Jet multiplicity", "Events", f);
  //DrawStack(t, "NBtagsNJets",  1, "NBtagNJets", "Events", f);
  DrawStack(t, "NBtagsNJets2", 1, "#Jets/#Btags", "Events", f);
  DrawStack(t, "NBtagJets",  1, "Btag multiplicity", "Events", f);
  DrawStack(t, "MET", 100, "MET (GeV)", "Events", f);
  DrawStack(t, "HT",  100, "HT (GeV)", "Events", f);
  DrawStack(t, "Vtx", 1, "nVtx", "Events", f);
  }

  // Jets
  if(plotsSel.Contains("Jets")){
  f = "Jets";
  DrawStack(t, "Jet0Pt", 270, "Leading Jet p_{T} (GeV)", "Events", f);
  DrawStack(t, "Jet1Pt", 220, "Subleading Jet p_{T} (GeV)", "Events", f);
  DrawStack(t, "Jet0CSV", 10, "Leading Jet CSV", "Events", f);
  DrawStack(t, "Jet1CSV", 10, "Subleading CSV", "Events", f);
  DrawStack(t, "Jet0DeepCSV", 10, "Leading Jet Deep CSV B", "Events", f);
  DrawStack(t, "Jet1DeepCSV", 10, "Subleading Deep CSV B", "Events", f);
  DrawStack(t, "Jet0Eta", 5, "Leading Jet #eta", "Events", f);
  DrawStack(t, "Jet1Eta", 5, "Subleading Jet #eta", "Events", f);
  //DrawStack(t, "Jet0DeepFlav", 2, "Leading Jet Deep Flavour B", "Events", f);
  //DrawStack(t, "Jet1DeepFlav", 2, "Subleading Deep Flavour B", "Events", f);
  
  //DrawStack(t, "JetAllPt", 100, "All Jet p_{T} (GeV)", "Events", f);
  //DrawStack(t, "JetAllCSV", 2, "All Jet CSV", "Events", f);
  //DrawStack(t, "JetAllDeepCSV", 2, "All Jet Deep CSV B", "Events", f);
  //DrawStack(t, "JetAllDeepFlav", 2, "All Jet Deep Flavour B", "Events", f);
  //DrawStack(t, "JetAllEta", 2, "All Jet #eta", "Events", f);
  }
}

void DrawStack(TopHistoReader *t, TString var, int rebin, TString xtit, TString ytit, TString subfolder){
  t->SetVar(var); t->SetRebin(rebin);
  t->SetIsData(false);

  TH1F* hprobe = t->GetHisto("TT");
  Int_t nbins = hprobe->GetNbinsX();
  Float_t b0 = hprobe->GetBinLowEdge(1);
  Float_t bN = hprobe->GetBinLowEdge(nbins+1);

  Plot* p = new Plot("plot", "", "All", nbins, b0, bN, "");


  // VV
  p->PrepareHisto(t->GetHisto("WWTo2L2Nu"), "WWTo2L2Nu", "VV", itBkg, kGreen-8);
  p->PrepareHisto(t->GetHisto("WZTo3LNU"), "WZTo3LNU", "VV", itBkg);
  //p->PrepareHisto(t->GetHisto("ZZ"), "ZZ", "ttV + VV", itBkg);

  // DY
  p->PrepareHisto(t->GetHisto("DYJetsToLL_M_10to50"), "DYJetsToLL_M_10to50", "DY", itBkg, kAzure+2);
  p->PrepareHisto(t->GetHisto("DYJetsToLL_MLL-50"), "DYJetsToLL_MLL-50", "DY", itBkg, kAzure+2);

  // NonW/Z
  p->PrepareHisto(t->GetHisto("TT_semilep"), "TT_Semilep", "NonW/Z leptons", itBkg, kGray+1);
  p->PrepareHisto(t->GetHisto("WJetsToLNu"), "WJetsToLNu", "NonW/Z leptons", itBkg, kGray+1);

  // tW
  p->PrepareHisto(t->GetHisto("tbarW_noFullHad"), "tbarW_noFullHad", "tW", itBkg, kViolet+2);
  p->PrepareHisto(t->GetHisto("tW_noFullHad"), "tW_noFullHad", "tW", itBkg);

  // ttbar
  p->PrepareHisto(t->GetHisto("TT"), "TT", "ttbar", itBkg, kRed);

  // Uncertainties
  /*
  vector<TString> sa = TStringToVector("TTWJetsToLNu_madspin, TTWJetsToQQ_madspin, TTZToLL_M_1to10, TTZToLLNuNu_M_10, TTZToQQ, WW, WZ, ZZ, DYJetsToLL_M_50, tbarW_noFullHad, tW_noFullHad, TTTo2L2Nu");
  vector<TString> pr = TStringToVector("ttV + VV, ttV + VV, ttV + VV, ttV + VV, ttV + VV, ttV + VV, ttV + VV, ttV + VV, DY, tW, tW, ttbar");
  vector<TString> unc = TStringToVector("MuonEffUp, MuonEffDown, ElecEffUp, ElecEffDown, TrigUp, TrigDown, JESUp, JESDown, JERUp, JERDown, PUUp, PUDown, BtagUp, BtagDown, MistagUp, MistagDown");

  for(int iu = 0; iu < (int) unc.size(); iu++){
    for(int is = 0; is < (int) sa.size(); is++){
      TString sample = sa.at(is); TString process = pr.at(is); TString sys = unc.at(iu);
      t->SetSyst(sys);
      p->PrepareHisto(t->GetHisto(sample), sample+"_"+sys, process, itSys, 1, sys);
    }
  }
  */

  // Data
  t->SetSyst("");
  t->SetIsData(true);
  p->PrepareHisto(t->GetHisto("HighEGJet"), "HighEGJet", "Data", itData);
  p->PrepareHisto(t->GetHisto("DoubleMuon"), "DoubleMuon", "Data", itData);
  p->PrepareHisto(t->GetHisto("SingleMuon"), "SingleMuon", "Data", itData);
  t->SetPath(path);

  if(var == "NBtagsNJets2"){
    p->SetBinLabels("0/0, 1/0, #splitline{1/1}{(blinded)}, 2/0, 2/1, 2/2, #geq 3/#geq 0");
    p->SetData();
    p->hData->SetBinContent(3, 0);
  }
  
  // Drawing options and labels...
  p->SetRatioErrorColor(kViolet-2);
  p->SetRatioErrorStyle(3244);
  p->SetStackErrorStyle(3244);
  p->doYieldsInLeg = false;
  
  p->SetTitleX(xtit, 0.18);
  p->SetTitleY(ytit, 0.075);
  p->SetYaxisOffset(0.65, 0.06);
  p->SetRatioYtitle("Data/MC", 0.17);
  p->SetYratioOffset(0.27, 0.15);
  p->SetPadPlotMargins("0.08, 0.12, 0.02, 0.10");
  p->SetRatioMin(0.5); p->SetRatioMax(1.5);
  p->SetLumi(Lumi*1000);
  p->SetTextForLumi("%3.2f pb^{-1} (5.02 TeV)");
 
  TString levch = "";
  if     (t->GetChan() == "Elec") levch = "ee";
  else if(t->GetChan() == "Muon") levch = "#mu#mu";
  else if(t->GetChan() == "ElMu") levch = "e#mu";
  p->SetChLabel(levch);
  if(t->GetLevel() == "2jets") p->SetChLabel(levch + " + #geq 2jets");
  if(t->GetLevel() == "1btag") p->SetChLabel(levch + " #geq 2jets + #geq 1btag");

  p->SetChLabelPos(0.25, 0.95);
  p->SetLegendTextSize(0.04);
  p->SetLegendNCol(3);
  p->SetLegendPosition(0.42, 0.79, 0.90, 0.90);
  
 
  // Save the stack plot
  TString plotFolder = plotDir+ "/" + subfolder + "/" + t->GetLevel()+"/";
  gSystem->mkdir(plotDir, true);
  gSystem->mkdir(plotDir+ "/" + subfolder + "/", true);
  gSystem->mkdir(plotFolder, true);
  p->SetPlotFolder(plotFolder);
  p->SetVarName(var+"_"+t->GetChan()+"_"+t->GetLevel());
  p->doSetLogy = false;
  p->doSignal = false;
  p->doSys = true;
  p->DrawStack();
  delete p;
  return;
}
