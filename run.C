#include "TString.h"
#include <iostream>
#include <fstream>

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

void run(TString samp, TString selection, Float_t xsec, Float_t sumofweights, Int_t year, TString outname, Int_t nSlots = 1, TString outpath = "", TString options = "", Bool_t isamcatnlo = false, Bool_t isData = false, Int_t nEvents = 0, Int_t FirstEvent = 0){

  PAFProject* myProject = 0;
  vector<TString> samples = TStringToVector(samp);

  // PAF mode selection (based on number of slots)
  PAFIExecutionEnvironment* pafmode = 0;
  if      (nSlots <=1 )  pafmode = new PAFSequentialEnvironment();
  else                   pafmode = new PAFPROOFLiteEnvironment(nSlots);

  myProject = new PAFProject(pafmode);
  myProject->AddDataFiles(samples);
  myProject->SetDefaultTreeName("Events");
  
  // Deal with first and last event
  if (nEvents > 0   ) myProject->SetNEvents(nEvents);
  if (FirstEvent > 0) myProject->SetFirstEvent(FirstEvent);
  
  // Set output file
  if(outpath == "") outpath = "./"+selection+"_temp";
  gSystem->mkdir(outpath, 1);
  if(!outname.EndsWith(".root")) outname += ".root";
  myProject->SetOutputFile(outpath + "/" + outname);
  
  // Parameters for the analysis
  myProject->SetInputParam("sampleName",        outname);
  myProject->SetInputParam("IsData",            isData    );
  myProject->SetInputParam("weight",            xsec/sumofweights);
  myProject->SetInputParam("IsMCatNLO",         isamcatnlo);
  myProject->SetInputParam("selection",         selection);
  myProject->SetInputParam("WorkingDir",        gSystem->pwd());
  myProject->SetInputParam("xsec",              xsec);
  myProject->SetInputParam("_options",          options);
  myProject->SetInputParam("year",              TString(Form("%i",year)));
  
  // Adding packages
  myProject->AddSelectorPackage("LeptonSelector");
  myProject->AddSelectorPackage("JetSelector");
  myProject->AddSelectorPackage("EventBuilder");
  
  // Analysis selector
  if(selection != "")
    myProject->AddSelectorPackage(selection);
  
  // Additional packages
  myProject->AddPackage("Lepton");
  myProject->AddPackage("Jet");
  myProject->AddPackage("mt2");
  myProject->AddPackage("Functions");
  myProject->AddPackage("LeptonSF");
  myProject->AddPackage("BTagSFUtil");
  
  myProject->Run();
}
