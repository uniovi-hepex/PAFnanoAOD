#include "Looper.h"
#include "TString.h"


class TopHistoReader{
  //const TString DefinedChannels[] = {"ElMu", "Elec", "Muon"}; Int_t nDefinedChannels = 3;

  bool doStackOverflow = true;
  bool isData = false;
  bool verbose = false;

  public:
  TopHistoReader(TString _path = "./", TString _process = "", TString _var = "", TString _chan = "", TString _level = ""){
    path = _path;
    process = _process;
    var = _var;
    chan = _chan;
    level = _level;
    if(!path.EndsWith("/")) path += "/";
    lumi = 39000;
    fname = "";
    rebin = 1;
  }

  TString GetHistoName(){ 
    TString name = "H_" + var + "_" + chan + "_" + level;
    if(syst != "") name += "_" + syst;
    return name;
  }

  Histo* GetNamedHisto(TString name, TString pr = "", Int_t rb = -1){
    // Load sum of histograms...
    if(pr.Contains(',')){
      vector<TString> vpr = TStringToVector(pr);
      Histo* h = GetNamedHisto(name, vpr.at(0), rb);
      for(Int_t j = 1; j < (Int_t) vpr.size(); j++) h->Add(GetNamedHisto(name, vpr.at(j), rb));
      h->SetDirectory(0);
      return h;
    }
    if(name.Contains(',')){
      vector<TString> vname = TStringToVector(name);
      Histo* h = GetNamedHisto(vname.at(0), pr, rb);
      for(Int_t j = 1; j < (Int_t) vname.size(); j++) h->Add(GetNamedHisto(vname.at(j), pr, rb));
      h->SetDirectory(0);
      return h;
    }

    if(pr != "") process = pr; if(rb != -1) rebin = rb;
    // Check if file exists...
    TString filename = path + "Tree_" + process + ".root";
    if(gSystem->AccessPathName(path + "Tree_" + process + ".root")){
      if(!process.EndsWith(".root")) process += ".root";
      filename = path + process;
    }

    if(verbose) cout << "Opening file: " << filename << endl;
    TFile* file = TFile::Open(filename);
    TH1F* h;
    if(verbose) cout << "Looking for histo: " << name << endl;
    file->GetObject(name,h);
    if(h == 0) cout << "ERROR: histo \"" << name << "\" not found in file: " << filename << endl;

    int nb = h->GetNbinsX(); 
    if(doStackOverflow){
      h->SetBinContent(nb, h->GetBinContent(nb+2));
      h->SetBinContent(nb+2, 0);
    }
    h->Rebin(rebin);
    if(!isData) h->Scale(lumi);
    h->SetLineColor(1); h->SetFillStyle(0);
    h->SetLineWidth(2); h->SetStats(0);
    h->SetTitle("");
    h->GetXaxis()->SetTitle("");
    h->GetYaxis()->SetTitle("");
    h->SetDirectory(0);
    delete file;
    Histo* histo = new Histo(*h);
    histo->SetStyle();
    return histo;
  }

  Int_t GetBinNumberForLevel(TString lev = ""){
    if(lev != "") level = lev;
    TString sLevels = "dilepton, ZVeto, MET, 2jets, 1btag";
    vector<TString> levels = TStringToVector(sLevels);
    Int_t nLevels = levels.size();
    for(int i = 0; i < nLevels; i++) if(level == levels.at(i)) return i;
    cout << Form("WARNING: level name \"%s\" unknown!!\n", level.Data());
    return 0;
  }

  Histo* GetYieldHisto(TString ch = "", TString pr = ""){
    if(ch != "") chan = ch;
    return GetNamedHisto("H_Yields_" + chan, pr);
  }

  Histo* GetYieldSSHisto(TString ch = "", TString pr = ""){
    if(ch != "") chan = ch;
    return GetNamedHisto("H_SSYields_" + chan, pr);
  }

  Histo* GetHisto(TString pr = "", TString v = "", Int_t rb = -1, TString lev = "", TString ch = "", TString s = ""){ 
    if(pr != "") process = pr; if(rb != -1) rebin = rb;
    if(v != "") var = v; if(lev != "") level = lev; if(ch != "") chan = ch; if(s != "") syst = s;
    TString hname = GetHistoName();
    return GetNamedHisto(hname, pr, rb);
  }

  Float_t GetYield(TString pr = "", TString lev = "", TString ch = "", TString s = "", TString v = "", bool isSS = false){ 
    Histo* h; Float_t y = 0;
    if(v != ""){
      h = GetHisto(pr, v, 1, lev, ch, s);
      y = h->Integral(); 
    }
    else{
      if(isSS) h = GetYieldHisto(ch, pr);
      else     h = GetYieldSSHisto(ch, pr);
      Int_t bin = GetBinNumberForLevel(lev);
      y = h->GetBinContent(bin);
    }
    delete h;
    return y;
  }

  Float_t GetYieldSS(TString pr = "", TString lev = "", TString ch = "", TString s = "", TString v = ""){ 
    return GetYield(pr, lev, ch, s, v, true);
  }

  Histo* GetSystHisto(TString syst = ""){ 
    return GetHisto("", "", -1,"", "", syst); 
  }

  vector<TString> vprocess;

  public:
  TString GetVar(){return var;}
  TString GetChan(){return chan;}
  TString GetLevel(){return level;}
  TString GetSyst(){return syst;}
  TString GetProcess(){return process;}
  TString GetPath(){return path;}
  float GetLumi(){return lumi;}

  void SetVar(TString t){var = t;}
  void SetChan(TString t){chan = t;}
  void SetLevel(TString t){level = t;}
  void SetSyst(TString t){syst = t;}
  void SetProcess(TString t){process = t;}
  void SetPath(TString t){path = t;}
  void SetLumi(float f){lumi = f;}
  void SetIsData(bool val = true){isData = val;}
  void SetRebin(int v){ rebin = v;}
  void SetVerbose(bool v){verbose = v;}

  void SetSamples(TString samples){vprocess = TStringToVector(samples);}
 

  protected:
  TString var;
  TString chan;
  TString level;
  TString syst;
  TString process;
  TString path;
  TFile* file;
  Int_t rebin;
  Float_t lumi;
  TString fname;

};
