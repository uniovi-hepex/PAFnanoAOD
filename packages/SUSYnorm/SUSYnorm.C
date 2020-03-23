#include "SUSYnorm.h"

#include "PAF.h"

#include <iostream>

SUSYnorm::SUSYnorm(TString _path, TString _filename): fhSMS(0) {
  path  = _path;
  filename = _filename;

  hname = "hSMS";
  loadHisto();
};

void SUSYnorm::loadHisto(){
  vector<TString> files = GetAllFiles(path, filename);
  for(int i=0; i<(int)files.size(); i++){
    cout << Form("[%i] ",i) << files[i] << endl;
  }
  fhSMS = (TH2D*)GetHistoFromFiles2(files, hname);
}

vector<TString> SUSYnorm::GetAllFiles(TString path, TString  filename, Bool_t verbose) {
  if(verbose) cout << Form("[GetAllFiles]: Obtaining files of form %s in folder %s\n", filename.Data(), path.Data());
  vector<TString> theFiles;
  if(!path.EndsWith("/")) path += "/";

  if(filename.Contains(",")){
    theFiles = TStringToVector(filename);
    for(int i=0; i<(int)theFiles.size(); i++){
      theFiles[i] = path+theFiles[i];
    }
    return theFiles;
  }

  TString command("ls ");
  if(filename != "")
    command += 
      path + "/" + filename + " " +
      path + "/" + filename + ".root " +
      path + "/" + filename + "_[0-9].root " +
      path + "/" + filename + "_[0-9][0-9].root " +
      path + "/Tree_" + filename + ".root " +
      path + "/Tree_" + filename + "_[0-9].root " +
      path + "/Tree_" + filename + "_[0-9][0-9].root";
  else command += path;

  command += " 2> /dev/null";
  if(verbose) cout << "[GetAllFiles] " << Form("Executing command: %s\n", command.Data());

  //We cannot use GetFromPipe because it is too verbose, so we implement
  //the full code
  //    TString result=gSystem->GetFromPipe(command);
  TString result;
  FILE *pipe = gSystem->OpenPipe(command, "r");
  if (!pipe) cout << Form("ERROR [GetAllFiles] Cannot run command \"%s\"\n", command.Data());
  else {
    TString line;
    while (line.Gets(pipe)) {
      //if(line.Contains("13")) continue;
      if (result != "")	result += "\n";
      result += line;
    }
    gSystem->ClosePipe(pipe);
  }
  
  if (result != "" ) {
    TObjArray* filesfound = result.Tokenize(TString('\n'));
    if (!filesfound) cout << ("ERROR [GetAllFiles]: Could not parse output while finding files\n");
    else {
      for (int i = 0; i < filesfound->GetEntries(); i++) theFiles.push_back(filesfound->At(i)->GetName());
      filesfound->Clear();
      delete filesfound;
    }
  }

  if (theFiles.size() == 0) cout << ("ERROR [GetAllFiles]: Could not find data!\n");
  // Removing duplicated files
  for(int i = 1; i < (int) theFiles.size(); i++){
    if(theFiles.at(i) == theFiles.at(i-1)){
      theFiles.erase(theFiles.begin()+i);
      i--;
    }
  }
  return theFiles;
}

TH1* SUSYnorm::GetHistoFromFiles(vector<TString> Files, TString histoName){
  Int_t nFiles = Files.size(); TFile *f;
  TH1 *h = nullptr; TH1 *htemp = nullptr;
  for(Int_t i = 0; i < nFiles; i++){
    f = TFile::Open(Files.at(i));
    if(!f->FindKey(histoName)){
      cout << "WARNING: not found histogram \"" << histoName << "\"" << endl;
      return h;
    }
    if(i == 0) {
      f->GetObject(histoName, h); 
    }
    else{
      f->GetObject(histoName, htemp); 
      h->Add(htemp);
    }
  }
  h->SetDirectory(0);
  f->Close();
  return h;
}

TH2* SUSYnorm::GetHistoFromFiles2(vector<TString> Files, TString histoName){
  Int_t nFiles = Files.size(); TFile *f;
  TH2 *h = nullptr; TH2 *htemp = nullptr;
  for(Int_t i = 0; i < nFiles; i++){
    f = TFile::Open(Files.at(i));
    if(!f->FindKey(histoName)){
      cout << "WARNING: not found histogram \"" << histoName << "\"" << endl;
      return h;
    }
    if(i == 0) {
      f->GetObject(histoName, h); 
    }
    else{
      f->GetObject(histoName, htemp); 
      h->Add(htemp);
    }
  }
  h->SetDirectory(0);
  f->Close();
  return h;
}


Double_t SUSYnorm::GetSUSYnorm(Float_t mStop, Float_t mLSP){
  Int_t bin = fhSMS->FindBin(mLSP, mStop, 1);
  return fhSMS->GetBinContent(bin);
}



Double_t SUSYnorm::GetStopXSec(Int_t StopMass){
  if      (StopMass == 100) return 1770;
  else if (StopMass == 105) return 1450;
  else if (StopMass == 110) return 1200;
  else if (StopMass == 115) return 998;
  else if (StopMass == 120) return 832;
  else if (StopMass == 125) return 697;
  else if (StopMass == 130) return 586;
  else if (StopMass == 135) return 495;
  else if (StopMass == 140) return 419;
  else if (StopMass == 145) return 357;
  else if (StopMass == 150) return 304;
  else if (StopMass == 155) return 261;
  else if (StopMass == 160) return 224;
  else if (StopMass == 165) return 194;
  else if (StopMass == 170) return 168;
  else if (StopMass == 175) return 146;
  else if (StopMass == 180) return 127;
  else if (StopMass == 185) return 111;
  else if (StopMass == 190) return 97.3;
  else if (StopMass == 195) return 85.6;
  else if (StopMass == 200) return 75.5;
  else if (StopMass == 205) return 66.8;
  else if (StopMass == 210) return 59.3;
  else if (StopMass == 215) return 52.7;
  else if (StopMass == 220) return 47.0;
  else if (StopMass == 225) return 42.0;
  else if (StopMass == 230) return 37.7;
  else if (StopMass == 235) return 33.8;
  else if (StopMass == 240) return 30.5;
  else if (StopMass == 245) return 27.5;
  else if (StopMass == 250) return 24.8;
  else if (StopMass == 255) return 22.5;
  else if (StopMass == 260) return 20.4;
  else if (StopMass == 265) return 18.6;
  else if (StopMass == 270) return 16.9;
  else if (StopMass == 275) return 15.5;
  else if (StopMass == 280) return 14.1;
  else if (StopMass == 285) return 12.9;
  else if (StopMass == 290) return 11.9;
  else if (StopMass == 295) return 10.9;
  else if (StopMass == 300) return 10.0;
  else if (StopMass == 305) return 9.18;
  else if (StopMass == 310) return 8.43;
  else if (StopMass == 315) return 7.75;
  else if (StopMass == 320) return 7.13;
  else if (StopMass == 325) return 6.57;
  else if (StopMass == 330) return 6.06;
  else if (StopMass == 335) return 5.59;
  else if (StopMass == 340) return 5.17;
  else if (StopMass == 345) return 4.78;
  else if (StopMass == 350) return 4.43;
  else if (StopMass == 355) return 4.10;
  else if (StopMass == 360) return 3.81;
  else if (StopMass == 365) return 3.54;
  else if (StopMass == 370) return 3.29;
  else if (StopMass == 375) return 3.06;
  else if (StopMass == 380) return 2.85;
  else if (StopMass == 385) return 2.65;
  else if (StopMass == 390) return 2.47;
  else if (StopMass == 395) return 2.31;
  else if (StopMass == 400) return 2.15;
  else if (StopMass == 405) return 2.01;
  else if (StopMass == 410) return 1.88;
  else if (StopMass == 415) return 1.76;
  else if (StopMass == 420) return 1.64;
  else if (StopMass == 425) return 1.54;
  else if (StopMass == 430) return 1.44;
  else if (StopMass == 435) return 1.35;
  else if (StopMass == 440) return 1.26;
  else if (StopMass == 445) return 1.19;
  else if (StopMass == 450) return 1.11;
  else if (StopMass == 455) return 1.05;
  else if (StopMass == 460) return 0.983;
  else if (StopMass == 465) return 0.925;
  else if (StopMass == 470) return 0.870;
  else if (StopMass == 475) return 0.819;
  else if (StopMass == 480) return 0.771;
  else if (StopMass == 485) return 0.727;
  else if (StopMass == 490) return 0.685;
  else if (StopMass == 495) return 0.646;
  else if (StopMass == 500) return 0.609;
  else if (StopMass == 505) return 0.575;
  else if (StopMass == 510) return 0.543;
  else if (StopMass == 515) return 0.513;
  else if (StopMass == 520) return 0.484;
  else if (StopMass == 525) return 0.458;
  else if (StopMass == 530) return 0.433;
  else if (StopMass == 535) return 0.409;
  else if (StopMass == 540) return 0.387;
  else if (StopMass == 545) return 0.367;
  else if (StopMass == 550) return 0.347;
  else if (StopMass == 555) return 0.329;
  else if (StopMass == 560) return 0.312;
  else if (StopMass == 565) return 0.296;
  else if (StopMass == 570) return 0.280;
  else if (StopMass == 575) return 0.266;
  else if (StopMass == 580) return 0.252;
  else if (StopMass == 585) return 0.240;
  else if (StopMass == 590) return 0.228;
  else if (StopMass == 595) return 0.216;
  else if (StopMass == 600) return 0.205;
  else if (StopMass == 605) return 0.195;
  else if (StopMass == 610) return 0.186;
  else if (StopMass == 615) return 0.177;
  else if (StopMass == 620) return 0.168;
  else if (StopMass == 625) return 0.160;
  else if (StopMass == 630) return 0.152;
  else if (StopMass == 635) return 0.145;
  else if (StopMass == 640) return 0.138;
  else if (StopMass == 645) return 0.131;
  else if (StopMass == 650) return 0.125;
  else if (StopMass == 655) return 0.119;
  else if (StopMass == 660) return 0.114;
  else if (StopMass == 665) return 0.108;
  else if (StopMass == 670) return 0.103;
  else if (StopMass == 675) return 0.987E-01;
  else if (StopMass == 680) return 0.942E-01;
  else if (StopMass == 685) return 0.899E-01;
  else if (StopMass == 690) return 0.858E-01;
  else if (StopMass == 695) return 0.820E-01;
  else if (StopMass == 700) return 0.783E-01;
  else if (StopMass == 705) return 0.748E-01;
  else if (StopMass == 710) return 0.715E-01;
  else if (StopMass == 715) return 0.683E-01;
  else if (StopMass == 720) return 0.653E-01;
  else if (StopMass == 725) return 0.624E-01;
  else if (StopMass == 730) return 0.597E-01;
  else if (StopMass == 735) return 0.571E-01;
  else if (StopMass == 740) return 0.546E-01;
  else if (StopMass == 745) return 0.523E-01;
  else if (StopMass == 750) return 0.500E-01;
  else if (StopMass == 755) return 0.479E-01;
  else if (StopMass == 760) return 0.459E-01;
  else if (StopMass == 765) return 0.439E-01;
  else if (StopMass == 770) return 0.421E-01;
  else if (StopMass == 775) return 0.403E-01;
  else if (StopMass == 780) return 0.386E-01;
  else if (StopMass == 785) return 0.370E-01;
  else if (StopMass == 790) return 0.355E-01;
  else if (StopMass == 795) return 0.340E-01;
  else if (StopMass == 800) return 0.326E-01;
  else if (StopMass == 805) return 0.313E-01;
  else if (StopMass == 810) return 0.300E-01;
  else if (StopMass == 815) return 0.288E-01;
  else if (StopMass == 820) return 0.276E-01;
  else if (StopMass == 825) return 0.265E-01;
  else if (StopMass == 830) return 0.254E-01;
  else if (StopMass == 835) return 0.244E-01;
  else if (StopMass == 840) return 0.234E-01;
  else if (StopMass == 845) return 0.225E-01;
  else if (StopMass == 850) return 0.216E-01;
  else if (StopMass == 855) return 0.208E-01;
  else if (StopMass == 860) return 0.199E-01;
  else if (StopMass == 865) return 0.192E-01;
  else if (StopMass == 870) return 0.184E-01;
  else if (StopMass == 875) return 0.177E-01;
  else if (StopMass == 880) return 0.170E-01;
  else if (StopMass == 885) return 0.164E-01;
  else if (StopMass == 890) return 0.157E-01;
  else if (StopMass == 895) return 0.151E-01;
  else if (StopMass == 900) return 0.145E-01;
  else if (StopMass == 905) return 0.140E-01;
  else if (StopMass == 910) return 0.135E-01;
  else if (StopMass == 915) return 0.129E-01;
  else if (StopMass == 920) return 0.125E-01;
  else if (StopMass == 925) return 0.120E-01;
  else if (StopMass == 930) return 0.115E-01;
  else if (StopMass == 935) return 0.111E-01;
  else if (StopMass == 940) return 0.107E-01;
  else if (StopMass == 945) return 0.103E-01;
  else if (StopMass == 950) return 0.991E-02;
  else if (StopMass == 955) return 0.954E-02;
  else if (StopMass == 960) return 0.919E-02;
  else if (StopMass == 965) return 0.885E-02;
  else if (StopMass == 970) return 0.853E-02;
  else if (StopMass == 975) return 0.822E-02;
  else if (StopMass == 980) return 0.792E-02;
  else if (StopMass == 985) return 0.763E-02;
  else if (StopMass == 990) return 0.735E-02;
  else if (StopMass == 995) return 0.709E-02;
  else if (StopMass == 1000) return 0.683E-02;
  else{ 
    cout << Form("Warning [GetStopXSec]: No Cross Section for that mass!! (mStop = %i) Extrapolating...\n", StopMass);
    Float_t v0; Float_t vf; Int_t pass;
    Float_t pmass; Float_t nmass;
    pass = 5; if(StopMass > 400) pass = 25;

    pmass = StopMass - StopMass%pass;
    nmass = StopMass - StopMass%pass + pass;

    v0 = GetStopXSec(pmass);
    vf = GetStopXSec(nmass);

    Float_t x  = float(StopMass%pass)/pass;
    Float_t newXsec = v0 + (vf-v0)*x;

    cout << Form("xsec(%g) = %g; xsec(%g) = %g --> xsec(%g) = %g\n", pmass, v0, nmass, vf, (Float_t) StopMass, newXsec);

    return newXsec;
  }
}

std::vector<TString> SUSYnorm::TStringToVector(TString t, char separator){
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

