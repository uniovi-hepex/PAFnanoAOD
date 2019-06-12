R__LOAD_LIBRARY(Histo.C+)
R__LOAD_LIBRARY(Looper.C+)
R__LOAD_LIBRARY(Plot.C+)
R__LOAD_LIBRARY(TResultsTable.C+)
R__LOAD_LIBRARY(CrossSection.C+)
R__LOAD_LIBRARY(TopHistoReader.C+)

#include "Histo.h"
#include "Looper.h"
#include "Plot.h"

Float_t Lumi = 41.2;
TString path =  "/pool/ciencias/userstorage/juanr/top/2017/sep05/";
TString plotDir = "/nfs/fanae/user/juanr/www/plots2017/temp/";


//TString chan = "Elec";
//TString level = "dilepton";
TString sChannels    = "ElMu, Elec, Muon";
TString sLevels      = "dilepton, ZVeto, MET, 2jets, 1btag";
vector<TString> vChan= TStringToVector(sChannels);
vector<TString> vLev = TStringToVector(sLevels  );
Int_t   nChan        = vChan.size();
Int_t   nLev         = vLev .size();

void draw(TString level = "1btag", TString chan = "ElMu");


void xsec(){
  cout << ">> ElMu" << endl;
  draw("1btag", "ElMu");
  cout << ">> Elec" << endl;
  draw("1btag", "Elec");
  cout << ">> Muon" << endl;
  draw("1btag", "Muon");
  return;
}

void draw(TString level, TString chan){
  TopHistoReader *t = new TopHistoReader(path);
  t->SetLumi(Lumi*1000);
  t->SetLevel(level);
  t->SetChan(chan);
  t->SetRebin(1);
  t->SetVar("Lep0Eta");
  t->SetIsData(false);

  TH1F* hprobe = t->GetHisto("WW");
  Int_t nbins = hprobe->GetNbinsX();
  Float_t b0 = hprobe->GetBinLowEdge(1);
  Float_t bN = hprobe->GetBinLowEdge(nbins+1);
  Plot* p = new Plot("plot", "", "All", nbins, b0, bN, "");

  p->PrepareHisto(t->GetHisto("TTWJetsToLNu_madspin"), "TTWJetsToLNu_madspin", "ttV", itBkg, kOrange+1);
  p->PrepareHisto(t->GetHisto("TTWJetsToQQ_madspin"), "TTWJetsToQQ_madspin", "ttV", itBkg);
  p->PrepareHisto(t->GetHisto("TTZToLL_M_1to10"), "TTZToLL_M_1to10", "ttV", itBkg);
  p->PrepareHisto(t->GetHisto("TTZToLLNuNu_M_10"), "TTZToLLNuNu_M_10", "ttV", itBkg);
  p->PrepareHisto(t->GetHisto("TTZToQQ"), "TTZToQQ", "ttV", itBkg);
  p->PrepareHisto(t->GetHisto("WW"), "WW", "VV", itBkg, kGreen-8);
  p->PrepareHisto(t->GetHisto("WZ"), "WZ", "VV", itBkg);
  p->PrepareHisto(t->GetHisto("ZZ"), "ZZ", "VV", itBkg);

  p->PrepareHisto(t->GetHisto("DYJetsToLL_M_50"), "DYJetsToLL_M_50", "DY", itBkg, kAzure+2);
  p->PrepareHisto(t->GetHisto("TTToSemiLeptonic"), "TTToSemiLeptonic", "NonW/Z leptons", itBkg, kGray+1);
  p->PrepareHisto(t->GetHisto("WJetsToLNu_MLM"), "WJetsToLNu_MLM", "NonW/Z leptons", itBkg, kGray+1);
  p->PrepareHisto(t->GetHisto("tbarW_noFullHad"), "tbarW_noFullHad", "tW", itBkg, kViolet+2);
  p->PrepareHisto(t->GetHisto("tW_noFullHad"), "tW_noFullHad", "tW", itBkg);
  p->PrepareHisto(t->GetHisto("TTTo2L2Nu"), "TTTo2L2Nu", "ttbar", itBkg, kRed);

  // Uncertainties
  vector<TString> sa = TStringToVector("TTWJetsToLNu_madspin, TTWJetsToQQ_madspin, TTZToLL_M_1to10, TTZToLLNuNu_M_10, TTZToQQ, WW, WZ, ZZ, DYJetsToLL_M_50, tbarW_noFullHad, tW_noFullHad, TTTo2L2Nu");
  vector<TString> pr = TStringToVector("ttV, ttV, ttV, ttV, ttV, VV, VV, VV, DY, tW, tW, ttbar");
  vector<TString> unc = TStringToVector("MuonEffUp, MuonEffDown, ElecEffUp, ElecEffDown, TrigUp, TrigDown, JESUp, JESDown, JERUp, JERDown, PUUp, PUDown, BtagUp, BtagDown, MistagUp, MistagDown");

  p->PrepareHisto(t->GetHisto("TTTo2L2Nu_TuneCP5up"),   "TTTo2L2Nu_TuneCP5Up",   "ttbar", itSys, 1, "UEUp");
  p->PrepareHisto(t->GetHisto("TTTo2L2Nu_TuneCP5down"), "TTTo2L2Nu_TuneCP5Down", "ttbar", itSys, 1, "UEDown");
  p->PrepareHisto(t->GetHisto("TTTo2L2Nu_hdampUP"),   "TTTo2L2Nu_hdampUP",   "ttbar", itSys, 1, "hdampUp");
  p->PrepareHisto(t->GetHisto("TTTo2L2Nu_hdampDOWN"),   "TTTo2L2Nu_hdampDOWN",   "ttbar", itSys, 1, "hdampDown");

  for(int iu = 0; iu < (int) unc.size(); iu++){
    for(int is = 0; is < (int) sa.size(); is++){
      TString sample = sa.at(is); TString process = pr.at(is); TString sys = unc.at(iu);
      t->SetSyst(sys);
      p->PrepareHisto(t->GetHisto(sample), sample+"_"+sys, process, itSys, 1, sys);
    }
  }

  p->AddSystematic("stat"); 
  p->AddToSystematicLabels("MuonEff, ElecEff, Trig, PU, Btag, Mistag, stat, UE, hdamp");

  t->SetSyst("");
  t->SetIsData(true);
  p->PrepareHisto(t->GetHisto("MuonEG_Run2017"), "MuonEG", "Data", itData);
  p->PrepareHisto(t->GetHisto("DoubleMuon_Run2017"), "DoubleMuon", "Data", itData);
  p->PrepareHisto(t->GetHisto("DoubleEG_Run2017"), "DoubleEG", "Data", itData);
  p->PrepareHisto(t->GetHisto("SingleMuon_Run2017"), "SingleMuon", "Data", itData);
  p->PrepareHisto(t->GetHisto("SingleElectron_Run2017"), "SingleElectron", "Data", itData);
  p->SetLumi(Lumi);
  p->SetLumiUnc(0.023);
 
  // Save the stack plot
  TString plotFolder = plotDir;
  gSystem->mkdir(plotFolder, true);
  p->SetPlotFolder(plotFolder);
  p->SetVarName("yields_"+t->GetChan()+"_"+t->GetLevel());
  p->doSignal = false;
  p->PrintSystYields();
  p->PrintYields("", "", "", "tex, txt, html");


  // Cross section
  TString xch  = t->GetChan();
  TString xlev = t->GetLevel();
  CrossSection *x = new CrossSection(p, "ttbar");
  x->SetTheoXsec(831.8);
  x->SetChannelTag(xch);
  x->SetLevelTag(xlev);
  x->SetBR(0.03263);
  //x->SetNFiducialEvents(1.15245e+06);
  //x->SetNSimulatedEvents(77229341);

  x->SetEfficiencySyst("JES, Btag, Mistag, ElecEff, MuonEff, LepEff, Trig, PU, JER");
  x->SetAcceptanceSyst("stat, UE, hdamp, scale, pdf, isr, fsr, q2, ME, CR");

  x->SwitchLabel("VV", "Dibosons");
  x->SwitchLabel("DY", "Drell-Yan");
  x->SwitchLabel("UE", "Underlying Event");
  x->SwitchLabel("fsr", "FSR scale");
  x->SwitchLabel("isr", "ISR scale");
  x->SwitchLabel("hdamp", "ME/PS matching (hdamp)");
  x->SwitchLabel("JES", "Jet Energy Scale");
  x->SwitchLabel("JER", "Jet Energy Resolution");
  x->SwitchLabel("Btag", "b-tagging efficiency");
  x->SwitchLabel("Mistag", "Mistagging efficiency");
  x->SwitchLabel("ElecEff", "Electron efficiencies");
  x->SwitchLabel("MuonEff", "Muon efficiencies");
  x->SwitchLabel("Trig", "Trigger efficiencies");
  x->SwitchLabel("PU", "Pile-Up");
  x->SwitchLabel("Scale", "ME scales");
  x->SwitchLabel("pdf", "PDF");
  x->SwitchLabel("stat", "MC stat");
  x->SwitchLabel("CR", "Color reconnection");

  x->SetUnc("Drell-Yan", 0.30);
  x->SetUnc("NonW/Z leptons",     0.50);
  x->SetUnc("Dibosons",  0.30);
  x->SetUnc("tW",        0.12);
  x->SetUnc("ttV",       0.30);

  // Adding extrapolation unc for muon eff.
  Float_t extElUnc2(0), extMuUnc2(0);
  if     (xch == "ElMu") extMuUnc2 = 0.005*0.005;
  else if(xch == "Muon") extMuUnc2 = 0.005*0.005*2;
  else if(xch == "Elec") extElUnc2 = 0;

  Float_t y = x->GetYield("ttbar");
  Float_t stat =  TMath::Abs(x->GetUnc("Muon efficiencies") - y)/y;
  Float_t MuonUnc = TMath::Sqrt(stat*stat + extMuUnc2);
  MuonUnc = MuonUnc*y + y;
  x->SetUnc("Muon efficiencies",  MuonUnc);

  x->SetOutputFolder(plotDir + "/xsec/");

  // Cross section from MC.............................................
  x->SetTableName("xsec_unc_MC");
  x->SetYieldsTableName("yields_MC");
  x->SetXsecTableName("xsec_MC");
  x->PrintSystematicTable("html,tex,txt");
  x->PrintYieldsTable("html,tex,txt");
  x->PrintCrossSection("html");

  delete x;
  delete p;
  return;
}
