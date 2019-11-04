#include "TString.h"
#include <iostream>
#include <fstream>
#include "TProof.h"

std::vector<TString> TStringToVector(TString t, char separator = ','){
  std::vector<TString> v;
  t.ReplaceAll(" ", "");
  Int_t n = t.CountChar(separator);
  TString element;
  for(Int_t i = 0; i < n; i++){
    element = t(0, t.First(separator));
    t = t(t.First(separator)+1, t.Sizeof());
    v.push_back(element);
  }
  v.push_back(t);
  return v;
}

void run(TString samp, TString selection, Double_t xsec, Double_t sumofweights, Int_t year, TString outname, Int_t nSlots = 1, TString outpath = "", TString options = "", Bool_t isamcatnlo = false, Bool_t isData = false, Long64_t nEvents = 0, Long64_t FirstEvent = 0, TString workingdir = "", TString path = "", Bool_t debug = false) {



  if (debug) {
    gProofDebugMask = TProofDebug::kAll;
    gProofDebugLevel = 5;
  }

  PAFProject* myProject = 0;
  vector<TString> samples = TStringToVector(samp);
  if(path != ""){
    if(!path.EndsWith("/")) path += "/";
    Int_t nFiles = samples.size(); Int_t i;
    for(i = 0; i < nFiles; i++) samples.at(i) = path + samples.at(i);
  }

  // PAF mode selection (based on number of slots)
  PAFIExecutionEnvironment* pafmode = 0;
  if   (nSlots <=1 ) pafmode = new PAFSequentialEnvironment();
  else               pafmode = new PAFPROOFLiteEnvironment(nSlots);


  myProject = new PAFProject(pafmode);
  myProject->AddDataFiles(samples);
  myProject->SetDefaultTreeName("/Events");
  
  // Deal with first and last event
  if (nEvents > 0   ) myProject->SetNEvents(nEvents);
  if (FirstEvent > 0) myProject->SetFirstEvent(FirstEvent);
  
  // Set output file
  if(outpath == "") outpath = "./"+selection+"_temp";
  gSystem->mkdir(outpath, 1);
  TString outfilename = Form("%s", outname.Data());
  if(!outfilename.EndsWith(".root")) outfilename += ".root";
  myProject->SetOutputFile(outpath + "/" + outfilename);
  
  Double_t tmpw = xsec/sumofweights;

  // Parameters for the analysis
  if(workingdir == "") workingdir = gSystem->pwd();
  myProject->SetInputParam("sampleName", outname);
  myProject->SetInputParam("IsData",     isData);
  myProject->SetInputParam("weight",     tmpw);
  myProject->SetInputParam("IsMCatNLO",  isamcatnlo);
  myProject->SetInputParam("selection",  selection);
  myProject->SetInputParam("WorkingDir", workingdir);
  myProject->SetInputParam("xsec",       xsec);
  myProject->SetInputParam("path", path); //quitar
  myProject->SetInputParam("_options",   options);
  myProject->SetInputParam("year",       TString(Form("%i",year)));
  
  if (selection == "TWAnalysis") {
    TString tmvalibpath = gSystem->Getenv("ROOTSYS");
    tmvalibpath += "/lib/libTMVA.so";
    myProject->AddLibrary(tmvalibpath);
  }

  // Adding packages
  myProject->AddSelectorPackage("LeptonSelector");
  myProject->AddSelectorPackage("JetSelector");
  myProject->AddSelectorPackage("EventBuilder");
  
  // Analysis selector
  if (selection != "")
    myProject->AddSelectorPackage(selection);
  if (selection == "TopAnalysis")
    myProject->AddSelectorPackage("LepEffTop");
  
  // Additional packages
  myProject->AddPackage("Lepton");
  myProject->AddPackage("Jet");
  myProject->AddPackage("Functions");
  if (selection == "TopAnalysis") myProject->AddPackage("mt2");
  myProject->AddPackage("LeptonSF");
  myProject->AddPackage("BTagSFUtil");
  
  if (selection == "TopAnalysis")
    myProject->AddPackage("SUSYnorm");

  myProject->Run();
}
