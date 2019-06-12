R__LOAD_LIBRARY(Histo.C+)
R__LOAD_LIBRARY(Looper.C+)
R__LOAD_LIBRARY(Plot.C+)
R__LOAD_LIBRARY(TopHistoReader.C+)

#include "Histo.h"
#include "Looper.h"
#include "Plot.h"

TString Era = "";
Float_t Lumi = 13.54;//41.2;
TString path18 =  "/pool/ciencias/userstorage/juanr/top/2018/sep25/";
TString path17 =  "/pool/ciencias/userstorage/juanr/top/2017/sep05/";
TString path17perEra =  "/pool/ciencias/userstorage/juanr/top/2017/sep05/perEra/";
TString plotDir = "/nfs/fanae/user/juanr/www/plots2018/";


//TString chan = "Elec";
//TString level = "dilepton";
TString sChannels    = "ElMu, Elec, Muon";
TString sLevels      = "dilepton, ZVeto, MET, 2jets, 1btag";
vector<TString> vChan= TStringToVector(sChannels);
vector<TString> vLev = TStringToVector(sLevels  );
Int_t   nChan        = vChan.size();
Int_t   nLev         = vLev .size();

void DrawStack(TopHistoReader *t, TString var, int rebin = 1, TString xtit = "", TString ytit = "Events", TString subfolder = "");
void draw(TString level = "1btag", TString chan = "ElMu");


void CompYears(){
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
  TopHistoReader *t = new TopHistoReader(path18);
  t->SetLumi(Lumi*1000);
  TString f = "";

  t->SetLevel(level);
  t->SetChan(chan);
  // Leptons
  f = "Leptons";
  //DrawStack(t, "DelLepPhi", 4, "#Delta#phi(l1,l2) (rad/#pi)", "Events", f);
  //DrawStack(t, "DelLepEta", 4, "#Delta#eta(l1,l2) ", "Events", f);

  DrawStack(t, "Lep0Pt",  100, "Leading lepton p_{T} (GeV)", "Events", f);
  DrawStack(t, "Lep1Pt",  100, "Subleading lepton p_{T} (GeV)", "Events", f);
  DrawStack(t, "Lep0Eta",  2, "Leading lepton #eta", "Events", f);
  DrawStack(t, "Lep1Eta",  2, "Subleading lepton #eta", "Events", f);
  DrawStack(t, "Lep0Iso",  2, "Leading lepton RelIso", "Events", f);
  DrawStack(t, "Lep1Iso",  2, "Subleading lepton RelIso", "Events", f);
  DrawStack(t, "DiLepPt", 80, "p_{T}^{ll} (GeV)", "Events", f);
  DrawStack(t, "DelLepPhi", 4, "#Delta#phi(l1,l2) (rad/#pi)", "Events", f);
  DrawStack(t, "DelLepEta", 4, "#Delta#eta(l1,l2) ", "Events", f);
  if(chan == "ElMu"){
    DrawStack(t, "ElecPt",  100, "Electron p_{T} (GeV)", "Events", f);
    DrawStack(t, "MuonPt",  100, "Muon p_{T} (GeV)", "Events", f);
    DrawStack(t, "ElecEta",  2, "Electron #eta", "Events", f);
    DrawStack(t, "MuonEta",  2, "Muon #eta", "Events", f);
    DrawStack(t, "ElecIso",  2, "Electron RelIso03", "Events", f);
    DrawStack(t, "MuonIso",  2, "Muon RelIso04", "Events", f);
  }
  DrawStack(t, "InvMass", 20, "m_{ll} (GeV)", "Events", f);
  if(chan != "ElMu") DrawStack(t, "InvMass2", 2, "m_{ll} (GeV)", "Events", f);

  // Global
  f = "Global";
  DrawStack(t, "NJets",  1, "Jet multiplicity", "Events", f);
  DrawStack(t, "NBtagsNJets",  1, "NBtagNJets", "Events", f);
  DrawStack(t, "NBtagJets",  1, "Btag multiplicity", "Events", f);
  DrawStack(t, "MET", 100, "MET (GeV)", "Events", f);
  DrawStack(t, "MT2", 60, "M_{T2} (GeV)", "Events", f);
  DrawStack(t, "HT",  100, "HT (GeV)", "Events", f);
  DrawStack(t, "Vtx", 1, "nVtx", "Events", f);

  // Jets
  f = "Jets";
  DrawStack(t, "Jet0Pt", 100, "Leading Jet p_{T} (GeV)", "Events", f);
  DrawStack(t, "Jet1Pt", 100, "Subleading Jet p_{T} (GeV)", "Events", f);
  DrawStack(t, "Jet0CSV", 2, "Leading Jet CSV", "Events", f);
  DrawStack(t, "Jet1CSV", 2, "Subleading CSV", "Events", f);
  DrawStack(t, "Jet0DeepCSV", 2, "Leading Jet Deep CSV B", "Events", f);
  DrawStack(t, "Jet1DeepCSV", 2, "Subleading Deep CSV B", "Events", f);
  DrawStack(t, "Jet0DeepFlav", 2, "Leading Jet Deep Flavour B", "Events", f);
  DrawStack(t, "Jet1DeepFlav", 2, "Subleading Deep Flavour B", "Events", f);
  DrawStack(t, "Jet0Eta", 2, "Leading Jet #eta", "Events", f);
  DrawStack(t, "Jet1Eta", 2, "Subleading Jet #eta", "Events", f);
  
  //DrawStack(t, "JetAllPt", 100, "All Jet p_{T} (GeV)", "Events", f);
  //DrawStack(t, "JetAllCSV", 2, "All Jet CSV", "Events", f);
  //DrawStack(t, "JetAllDeepCSV", 2, "All Jet Deep CSV B", "Events", f);
  //DrawStack(t, "JetAllDeepFlav", 2, "All Jet Deep Flavour B", "Events", f);
  //DrawStack(t, "JetAllEta", 2, "All Jet #eta", "Events", f);
}

void DrawStack(TopHistoReader *t, TString var, int rebin, TString xtit, TString ytit, TString subfolder){
  t->SetVar(var); t->SetRebin(rebin);
  t->SetIsData(false);

  t->SetPath(path18);
  TH1F* hprobe = t->GetHisto("DoubleMuon_Run2018");
  Int_t nbins = hprobe->GetNbinsX();
  Float_t b0 = hprobe->GetBinLowEdge(1);
  Float_t bN = hprobe->GetBinLowEdge(nbins+1);

  Plot* p = new Plot("plot", "", "All", nbins, b0, bN, "");

  // Data
  t->SetIsData(true);
/*
  p->PrepareHisto(t->GetHisto("MuonEG"), "MuonEG", "Data", itSignal);
  p->PrepareHisto(t->GetHisto("DoubleMuon"), "DoubleMuon", "Data", itSignal);
  p->PrepareHisto(t->GetHisto("DoubleEG"), "DoubleEG", "Data", itSignal);
  if(!doOnlyDoubleLepDatasets){
    p->PrepareHisto(t->GetHisto("SingleMuon"), "SingleMuon", "Data", itSignal);
    p->PrepareHisto(t->GetHisto("SingleElectron"), "SingleElectron", "Data", itSignal);
  }
*/
 
  // 2017 BCDE
  t->SetPath(path17perEra);
  vector<TString> BCDE = TStringToVector("B, C, D, E");
  for(int i = 0; i < (int) BCDE.size(); i++){
    Era = BCDE.at(i);
    p->PrepareHisto(t->GetHisto("MuonEG_Run2017"+Era), "MuonEG Run 2017 BCDE", "Run 2017 BCDE", itSignal, 1);
    p->PrepareHisto(t->GetHisto("DoubleMuon_Run2017"+Era), "DoubleMuon Run 2017 BCDE", "Run 2017 BCDE", itSignal);
    p->PrepareHisto(t->GetHisto("DoubleEG_Run2017"+Era), "DoubleEG Run 2017 BCDE", "Run 2017 BCDE", itSignal);
    p->PrepareHisto(t->GetHisto("SingleMuon_Run2017"+Era), "SingleMuon Run 2017 BCDE", "Run 2017 BCDE", itSignal);
    p->PrepareHisto(t->GetHisto("SingleElectron_Run2017"+Era), "SingleElectron Run 2017 BCDE", "Run 2017 BCDE", itSignal);
  }

  // 2017 F
  t->SetPath(path17perEra);
  Era = "F";
  p->PrepareHisto(t->GetHisto("MuonEG_Run2017"+Era), "MuonEG Run 2017 "+Era, "Run 2017 "+Era, itSignal, kRed+1);
  p->PrepareHisto(t->GetHisto("DoubleMuon_Run2017"+Era), "DoubleMuon Run 2017", "Run 2017 "+Era, itSignal);
  p->PrepareHisto(t->GetHisto("DoubleEG_Run2017"+Era), "DoubleEG Run 2017"+Era, "Run 2017 "+Era, itSignal);
  p->PrepareHisto(t->GetHisto("SingleMuon_Run2017"+Era), "SingleMuon Run 2017"+Era, "Run 2017 "+Era, itSignal);
  p->PrepareHisto(t->GetHisto("SingleElectron_Run2017"+Era), "SingleElectron Run 2017"+Era, "Run 2017 "+Era, itSignal);

  // 2018 A, B, C
  vector<TString> eras = TStringToVector("A, B, C");
  t->SetPath(path18);
  Int_t Colors2018[3] = {kTeal+2, kAzure-2, kViolet+2};
  Int_t color = 0;
  for(int i = 0; i < (int) eras.size(); i++){
    Era = eras.at(i);
    color = Colors2018[i];
    p->PrepareHisto(t->GetHisto("MuonEG_Run2018"+Era), "MuonEG Run 2018"+Era, "Run 2018 " + Era, itSignal, color);
    p->PrepareHisto(t->GetHisto("DoubleMuon_Run2018"+Era), "DoubleMuon Run 2018"+Era, "Run 2018 " + Era, itSignal, color);
    p->PrepareHisto(t->GetHisto("DoubleEG_Run2018"+Era), "DoubleEG Run 2018"+Era, "Run 2018 " + Era, itSignal, color);
    p->PrepareHisto(t->GetHisto("SingleMuon_Run2018"+Era), "SingleMuon Run 2018"+Era, "Run 2018 " + Era, itSignal, color);
    p->PrepareHisto(t->GetHisto("SingleElectron_Run2018"+Era), "SingleElectron Run 2018"+Era, "Run 2018 " + Era, itSignal, color);
  }

  // Drawing options and labels...
  p->SetRatioErrorColor(kViolet-2);
  p->SetRatioErrorStyle(3244);
  p->SetStackErrorStyle(3244);
  p->doYieldsInLeg = false;
  p->doLegend = true;

  p->SetPadPlotLimits( "0.02, 0.45, 0.98, 0.95");
  p->SetPadRatioLimits("0.02, 0.1, 0.98, 0.45");
  p->SetPadPlotMargins( "0.00, 0.08, 0, 0.08");
  p->SetPadRatioMargins("0.00, 0.05, 0, 0.08");

  p->SetTitleX(xtit, 0.12);
  p->SetXaxisOffset(1, 0.11);

  p->SetTitleY(ytit, 0.12);
  p->SetYaxisOffset(0.4, 0.09);

  p->SetRatioYtitle("Data / 2017 BCDE", 0.12); // 0.14
  p->SetYratioOffset(0.32, 0.11);
  p->SetRatioMin(0.5); p->SetRatioMax(1.5);

  p->SetLumi(Lumi);
 
  TString levch = "";
  if     (t->GetChan() == "Elec") levch = "ee";
  else if(t->GetChan() == "Muon") levch = "#mu#mu";
  else if(t->GetChan() == "ElMu") levch = "e#mu";
  p->SetChLabel(levch);
  if(t->GetLevel() == "2jets") p->SetChLabel(levch + " + #geq 2jets");
  if(t->GetLevel() == "1btag") p->SetChLabel(levch + " #geq 2jets + #geq 1btag");

  //p->SetCMSlabel(TString t);
  p->SetCMSmodeLabel("Preliminary");
  p->SetCMSLabelPos(0.1, 1.05, 0.08);
  p->SetCMSmodeLabelPos(0.17, 1.045, 0.09);

  p->SetTextForLumi("13 TeV", 0.85, 1.05, 0.08);
  p->SetChLabelPos(0.35, 0.95, 0.08);
  
  p->SetLegendTextSize(0.06);
  p->SetLegendNCol(1);
  p->SetLegendPosition(0.72, 0.70, 0.95, 0.96);
  
 
  // Save the stack plot
  TString plotFolder = plotDir+ "/" + subfolder + "/" + t->GetLevel()+"/";
  gSystem->mkdir(plotDir, true);
  gSystem->mkdir(plotDir+ "/" + subfolder + "/", true);
  gSystem->mkdir(plotFolder, true);
  p->SetPlotFolder(plotFolder);
  p->SetVarName(var+"_"+t->GetChan()+"_"+t->GetLevel());
  p->doSetLogy = false;
  p->doSignal = true;
  p->doSys = false;
  p->DrawComp("", 1, "hist");
  delete p;
  return;
}
