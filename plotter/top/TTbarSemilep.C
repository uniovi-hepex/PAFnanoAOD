R__LOAD_LIBRARY(Histo.C+)
R__LOAD_LIBRARY(Looper.C+)
R__LOAD_LIBRARY(Plot.C+)
R__LOAD_LIBRARY(TopHistoReader.C+)

#include "Histo.h"
#include "Looper.h"
#include "Plot.h"

TString outpath = "/nfs/fanae/user/juanr/www/stop/CombTTbarSemilep/jun28/";
TString level = "Lep";

void DrawStack(TopHistoReader *t, TString var, int rebin = 1, TString xtit = "", TString ytit = "Events");
void NormPlot(Plot* p);

void TTbarSemilep(){
  TString path = "/nfs/fanae/user/juanr/nanoAOD/Top_temp/jun28/";//CombTTbarSemilep/temp/";
  TopHistoReader *t = new TopHistoReader(path);
  t->SetChan("All");
  t->SetLumi(35850);

  // level = "Lep";  DrawStack(t, "MET", 10, "MET (GeV)");

  TString ch[2] = {"Elec", "Muon"};
  TString lev[5] = {"Lep", "g1btag", "exact1btag", "fwdJet", "g2btag"};

  //t->SetChan(ch[1]); level = "g2btag"; 
  //t->SetLevel(level);
  //DrawStack(t, "METPhi", -1, "MET #phi (rad)");

  for(int j = 0; j < 2; j++){
    t->SetChan(ch[j]);
    for(int k = 4; k < 5; k++){
      level = lev[k];
      t->SetLevel(level);
      DrawStack(t, "MET", 20, "MET (GeV)");
      DrawStack(t, "NJets", 1, "Number of jets");
      DrawStack(t, "JetPt", 10, "Jet p_{T} (GeV)", "Jets");
      DrawStack(t, "HT", 50, "HT (GeV)");
      //DrawStack(t, "MT", 20, "m_{T}(lep,MET) (GeV)");
      if(ch[j] == "Muon") DrawStack(t, "MuonPt", 8, "Muon p_{T} (GeV)");
      if(ch[j] == "Elec") DrawStack(t, "ElecPt", 20, "Election p_{T} (GeV)");
      //DrawStack(t, "LepPhi", 1, "Lepton #phi (rad/#pi)");
      //DrawStack(t, "METPhi", 2, "MET #phi (rad/#pi)");
      DrawStack(t, "DeltaPhi", 2, "#Delta #phi (rad/#pi)");
      DrawStack(t, "InvMass", 2, "#Dilepton invariant mass (GeV)");
      //DrawStack(t, "DPhiLepMET", 1, "#Delta#phi(e/#mu, MET) (rad/#pi)");
      DrawStack(t, "DPhiLepMET2", 1, "#Delta#phi(e/#mu, MET) (rad/#pi)");
      //DrawStack(t, "NPV", 1, "Number of primary verteces");
      //DrawStack(t, "MT2", 1, "M_{T2} (GeV)");
      DrawStack(t, "METComb", 10, "MET (GeV)");
      //DrawStack(t, "METCombPhi", 2, "MET #phi (rad/#pi)");
      DrawStack(t, "METCombDPhi", 1, "MET1-MET2 #Delta #phi (rad/#pi)");
      if(k==4) DrawStack(t, "BestDijetWMass", 20, "Best W m_{jj} (GeV)");
      //DrawStack(t, "NBtags", 1, "Number of b-tagged jets");
      //DrawStack(t, "Jet0Pt", 10, "Leading jet p_{T} (GeV)");
      //DrawStack(t, "TrueInt", 1, "nTrueInt");
      //DrawStack(t, "Jet0PtNoBtag", 10, "Leading not b-tagged jet p_{T} (GeV)");
      //DrawStack(t, "FHT", 50, "HT (including fwd jets) (GeV)");
      //DrawStack(t, "FwdJet0Pt", 20, "Leading fwd jet p_{T} (GeV)");
      //DrawStack(t, "FwdJetPt", 20, "Forward Jet p_{T} (GeV)", "Jets");
      //DrawStack(t, "NJetsNBtag", 1, "[NJets, NBtags]");
      //DrawStack(t, "NFwdJet", 1, "Number of forward jets");
    }
  }
}

void DrawStack(TopHistoReader *t, TString var, int rebin, TString xtit, TString ytit){
  t->SetVar(var); t->SetRebin(rebin);
  t->SetIsData(false);

  TH1F* hprobe = t->GetHisto("TT");
  Int_t nbins = hprobe->GetNbinsX();
  double b0 = hprobe->GetBinLowEdge(1);
  double bN = hprobe->GetBinLowEdge(nbins+1);

  Plot* p = new Plot("plot", "", "All", nbins, b0, bN, "");

  // DY
  p->PrepareHisto(t->GetHisto("DYJetsToLL_M_50"), "DYJetsToLL_M_50", "DY", itBkg, kAzure+2);
  
  // WJets
  p->PrepareHisto(t->GetHisto("WJetsToLNu"), "WJetsToLNu", "W+Jets", itBkg, kGreen-7);

  // s-channel
  p->PrepareHisto(t->GetHisto("ST_s_channel"), "ST_s_channel", "single top", itBkg, kMagenta+2);
 
  // t-channel
  p->PrepareHisto(t->GetHisto("ST_t_channel_antitop"), "ST_t_channel_antitop", "single top", itBkg, kOrange+1);
  p->PrepareHisto(t->GetHisto("ST_t_channel_top"),     "ST_t_channel_top",     "single top", itBkg, kOrange+1);

  // tW
  p->PrepareHisto(t->GetHisto("tbarW"), "TbarW", "single top", itBkg, kViolet+2);
  p->PrepareHisto(t->GetHisto("tW"),    "TW", "single top", itBkg);

 // ttbar
  p->PrepareHisto(t->GetHisto("TT"), "TTTo2L2Nu", "ttbar", itBkg, kRed);

  // Data
  t->SetIsData(true);
  p->PrepareHisto(t->GetHisto("SingleMuon"), "SingleMuon", "Data", itData);
  p->PrepareHisto(t->GetHisto("SingleElectron"), "SingleElectron", "Data", itData);

  p->ScaleProcess("ttbar", p->GetHisto("Data")->GetYield()/p->GetHisto("ttbar")->GetYield());
  
  // Drawing options and labels...
  p->SetRatioErrorColor(kGray);
  p->SetRatioErrorStyle(3244);
  p->SetStackErrorStyle(3244);
  p->doYieldsInLeg = true;
  
  p->SetTitleX(xtit, 0.18);
  p->SetTitleY(ytit, 0.08);
  p->SetYaxisOffset(0.45, 0.06);
  p->SetYratioOffset(0.22, 0.17);
  p->SetPadPlotMargins("0.08, 0.12, 0.03, 0.10");
  p->SetRatioMin(0.6); p->SetRatioMax(1.4);

  p->SetLegendTextSize(0.05);
  p->SetLegendPosition(0.60, 0.30, 0.90, 0.90);

  if(var.Contains("Phi") || var.Contains("phi")){
    p->SetLegendNCol(3);
    p->SetLegendPosition(0.33, 0.78, 0.93, 0.9);
    p->SetLegendTextSize(0.04);
  }
  
 
  // Save the stack plot
  p->doSignal = false;
  p->SetVarName(var+"_"+t->GetChan());
  p->SetPlotFolder(outpath + "/" + level + "/");// +t->GetChan() + "/");
  p->doSetLogy = false;
  p->DrawStack();

  p->SetVarName(var+"_"+t->GetChan()+"_log");
  p->SetPlotFolder(outpath + "/" + level + "/log/");// +t->GetChan() + "/");
  p->doSetLogy = true;
  p->SetPlotMinimum(10);
  p->DrawStack();

  delete p;
  return;
}

void NormPlot(Plot* p){
  TString name;
  Int_t nhistos = p->VBkgs.size();
  for(int i = 0; i < nhistos; i++){
    name = p->VBkgs.at(i)->GetProcess();
    cout << " >> " << name << endl;
  }
}
