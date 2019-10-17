#include "TopAnalysis.h"

ClassImp(TopAnalysis);

TopAnalysis::TopAnalysis() : PAFChainItemSelector() {
  fTree           = 0;
  fHWeightsFidu   = 0;
  TPassMETFilters = 0;
  TPassTrigger    = 0;
  isSS            = 0;
  event           = 0;
  lumiblock       = 0;
  //sumWeights      = 0; //quitar



  for(Int_t ch = 0; ch < nChannels; ch++){
    for(Int_t cut = 0; cut < nLevels; cut++){
      if(!gIsData){
        fHPSweights[ch][cut]  = 0;
        fHPDFweights[ch][cut]  = 0;
        fHScaleWeights[ch][cut]  = 0;
      }
      for(Int_t sys = 0; sys < nSyst; sys++){
        if(gIsData && sys!=0) break;
        fHMET[ch][cut][sys]         = 0;
        fHMT2[ch][cut][sys]         = 0;
        fHLep0Eta[ch][cut][sys]     = 0;
        fHLep1Eta[ch][cut][sys]     = 0;
        fHMuonEta[ch][cut][sys]     = 0;
        fHElecEta[ch][cut][sys]     = 0;
        fHDelLepPhi[ch][cut][sys]   = 0;
        fHDelLepEta[ch][cut][sys]   = 0;
        fHHT[ch][cut][sys]          = 0;
        fHJetEta[ch][cut][sys]     = 0;
        fHJetBtagEta[ch][cut][sys]     = 0;
        fHJetBtag0Eta[ch][cut][sys]     = 0;
        fHJet0Eta[ch][cut][sys]     = 0;
        fHJet1Eta[ch][cut][sys]     = 0;
        fHDYInvMass[ch][cut][sys]       = 0;
        fHDYInvMassSF[ch][cut][sys]     = 0;
        fHInvMass[ch][cut][sys]       = 0;
        fHInvMass2[ch][cut][sys]      = 0;
        fHNBtagsNJets[ch][cut][sys]   = 0;
        fHNJets[ch][cut][sys]        = 0;
        fHNBtagJets[ch][cut][sys]    = 0;
        fHJet0Pt[ch][cut][sys]       = 0;
        fHJet1Pt[ch][cut][sys]       = 0;
        fHJetPt[ch][cut][sys]       = 0;
        fHJetBtagPt[ch][cut][sys]       = 0;
        fHJetBtag0Pt[ch][cut][sys]       = 0;
        fHDiLepPt[ch][cut][sys]      = 0;
        fHLep0Iso[ch][cut][sys]       = 0;
        fHLep1Iso[ch][cut][sys]       = 0;
        fHLep0Pt[ch][cut][sys]       = 0;
        fHLep1Pt[ch][cut][sys]       = 0;
        fHMuonIso[ch][cut][sys]       = 0;
        fHElecIso[ch][cut][sys]       = 0;
        fHMuonPt[ch][cut][sys]       = 0;
        fHElecPt[ch][cut][sys]       = 0;
        fHJetCSV[ch][cut][sys]       = 0;
        fHJet0CSV[ch][cut][sys]       = 0;
        fHJet1CSV[ch][cut][sys]       = 0;
        fHJetDeepCSV[ch][cut][sys]       = 0;
        fHJet0DeepCSV[ch][cut][sys]       = 0;
        fHJet1DeepCSV[ch][cut][sys]       = 0;
        fHJetDeepFlav[ch][cut][sys]       = 0;
        fHJet0DeepFlav[ch][cut][sys]       = 0;
        fHJet1DeepFlav[ch][cut][sys]       = 0;
        fHvertices[ch][cut][sys]     = 0;
        if(cut == 0){
          fHyields[ch][sys]     = 0;
          fHFiduYields[ch][sys]     = 0;
          fHSSyields[ch][sys]   = 0;
        }
      }
    }
  }
}
void TopAnalysis::Summary(){}

void TopAnalysis::Initialise(){
  gIsData      = GetParam<Bool_t>("IsData");
  selection    = GetParam<TString>("selection");
  gSampleName  = GetParam<TString>("sampleName");
  gOptions     = GetParam<TString>("_options");
  gDoSyst      = true;// gOptions.Contains("doSyst")? true : false;
  year         = GetParam<TString>("year").Atoi();
  gIsTTbar     = false;
  gIsTTany     = false;
  gSelection   = GetSelection(selection);
  gDoJECunc    = gOptions.Contains("JECunc")? true : false;
  gDoPDFunc    = gOptions.Contains("PDF")? true : false;
  gDoPSunc     = gOptions.Contains("PS")? true : false;
  gDoScaleUnc  = gOptions.Contains("Scale")? true : false;
  gPUWeigth    = gOptions.Contains("PUweight")? true : false;
  gPrefire     = gOptions.Contains("prefire")? true : false;
  JetPt        = gOptions.Contains("JetPtNom")? "Jet_pt_nom" : "Jet_pt";
  if ((gSampleName == "TT" || gSampleName.BeginsWith("TT_")) && year == 2016) gIsTTbar = true;
  if (gIsTTbar || gSampleName.BeginsWith("TTTo2L2Nu")) gIsTTany = true;


  nPDFweights = gIsTTbar ? 100 : 33;

  makeTree   = true;
  miniTree = true;
  makeHistos = false;
  gIsSignal= false; //quitar
  if(makeTree){
	if ((gSampleName.BeginsWith("stop"))) gIsSignal = true; //quitar
    fTree   = CreateTree("MiniTree","Created with PAF");
    SetLeptonVariables();
    SetJetVariables();
    SetEventVariables();

  }
  
  
  if(gIsSignal){
	  TString path = GetParam<TString>("path");
	  //snorm = new SUSYnorm(path, gSampleName);}
	  cout<<path<<endl;
	  snorm = new SUSYnorm(path , "SMS_T2tt_3J_xqcut_20_top_corridor_2Lfilter_TuneCUETP8M2T4_madgra");}
     
  // Uncertainties
  useSyst.push_back(kNorm);
  if(gDoSyst && !gIsData){
    useSyst.push_back(kMuonEffUp);
    useSyst.push_back(kMuonEffDown);
    useSyst.push_back(kElecEffUp);
    useSyst.push_back(kElecEffDown);
    useSyst.push_back(kBtagUp);
    useSyst.push_back(kBtagDown);
    useSyst.push_back(kMistagUp);
    useSyst.push_back(kMistagDown);
    useSyst.push_back(kTrigUp);
    useSyst.push_back(kTrigDown);

    if(gDoJECunc){
      useSyst.push_back(kJESUp);
      useSyst.push_back(kJESDown);
      useSyst.push_back(kJERUp);
      useSyst.push_back(kJERDown);
    }
    if(gPUWeigth){
      useSyst.push_back(kPUUp);
      useSyst.push_back(kPUDown);
    }
    if(gPrefire){
      useSyst.push_back(kPrefireUp);
      useSyst.push_back(kPrefireDown);
    }
    if(gDoPSunc){
      useSyst.push_back(kISRUp);
      useSyst.push_back(kISRDown);
      useSyst.push_back(kFSRUp);
      useSyst.push_back(kFSRDown);
    }
  }
  nSyst = useSyst.size();
  InitHistos();
  metvar   = year == 2017? "METFixEE2017" : "MET";
  metvarpt = year == 2017? "METFixEE2017_pt" : "MET_pt";
  if(year != 2017 and gOptions.Contains("JetPtNom")) metvarpt = "MET_pt_nom";

  gIs2017 = false; gIs2016 = false; gIs2018 = false;
  if     (year == 2017) gIs2017 = true;
  else if(year == 2018) gIs2018 = true;
  else if(year == 2016) gIs2016 = true;
  // b tagging
  TString pwd  = GetParam<TString>("WorkingDir");
  TString BTagSFPath = Form("%s/packages/BTagSFUtil", pwd.Data());
  TString taggerName="DeepFlav"; //"CSVv2"; //"DeepCSV"; // DeepFlav   //quitar dejar DeepFlav
  TString MeasType = "mujets";
  TString stringWP = "Medium";//"Medium"; //quitar dejar medium
  //if(taggerName == "DeepFlav" && year == 2017) MeasType = "comb";
  fBTagSFnom = new BTagSFUtil(MeasType.Data(), BTagSFPath, taggerName.Data(), stringWP,  0, year);
  if(!gIsData){
    fBTagSFbUp = new BTagSFUtil(MeasType.Data(), BTagSFPath, taggerName.Data(), stringWP,  1, year);
    fBTagSFbDo = new BTagSFUtil(MeasType.Data(), BTagSFPath, taggerName.Data(), stringWP, -1, year);
    fBTagSFlUp = new BTagSFUtil(MeasType.Data(), BTagSFPath, taggerName.Data(), stringWP,  3, year);
    fBTagSFlDo = new BTagSFUtil(MeasType.Data(), BTagSFPath, taggerName.Data(), stringWP, -3, year);
  }
}

void TopAnalysis::InsideLoop(){
  event     = Get<ULong64_t>("event");
  lumiblock = Get<UInt_t>("luminosityBlock");
//aqui quitar

  if (gIsSignal){
	  ngenPart=Get<Int_t>("nGenPart");

  int i;
  iSt = 0;
  iLSP = 0;
  for(Int_t i = 0; i < ngenPart; i++){
    Int_t genpdgId = Get<Int_t>("GenPart_pdgId", i);

    if(genpdgId == 1000006)  iSt=i; 
    if(genpdgId == 1000022)  iLSP=i; }
  m_stop = Get<Float_t>("GenPart_mass", iSt);
  m_LSP = Get<Float_t>("GenPart_mass", iLSP);

}
  
  // Vectors with the objects
  genLeptons     = GetParam<vector<Lepton>>("genLeptons");
  selLeptons     = GetParam<vector<Lepton>>("selLeptons");
  vetoLeptons    = GetParam<vector<Lepton>>("vetoLeptons");
  selJets        = GetParam<vector<Jet>>("selJets");
  //selJetsJecUp   = GetParam<vector<Jet>>("selJetsJecUp");
  //selJetsJecDown = GetParam<vector<Jet>>("selJetsJecDown");
  //Jets15         = GetParam<vector<Jet>>("Jets15");
  //vetoJets       = GetParam<vector<Jet>>("vetoJets");
  //genJets        = GetParam<vector<Jet>>("genJets");
  //mcJets         = GetParam<vector<Jet>>("mcJets");

  // Weights and SFs
  NormWeight     = GetParam<Double_t>("NormWeight");
  TrigSF         = GetParam<Float_t>("TriggerSF");
  TrigSFerr      = GetParam<Float_t>("TriggerSFerr");

  if(gIsSignal){
	  Float_t norm = snorm->GetSUSYnorm(m_stop,m_LSP);
	  xsec = snorm->GetStopXSec(m_stop);
	  NormWeight = norm != 0 ? xsec/norm : 0;
	  //cout << Form("m_stop:%f y m_LSP=%f  xsec: %f, %f\n", m_stop, m_LSP,  xsec,NormWeight);
  }
  if(!gIsData && gPUWeigth){
    PUSF         = Get<Float_t>("puWeight");
    PUSF_Up      = Get<Float_t>("puWeightUp");
    PUSF_Down    = Get<Float_t>("puWeightDown");
  }
  else{PUSF = 1; PUSF_Up = 1; PUSF_Down = 1;}
  if(!gIsData && gPrefire){
    PrefWeight   = Get<Float_t>("PrefireWeight");
    PrefWeightUp = Get<Float_t>("PrefireWeight_Up");
    PrefWeightDo = Get<Float_t>("PrefireWeight_Down");}
  else{PrefWeight = 1; PrefWeightUp = 1; PrefWeightDo = 1;}

  // Event variables
  gChannel       = GetParam<Int_t>("gChannel");
  TPassMETFilters = GetParam<Bool_t>("METfilters");
  TPassTrigger    = GetParam<Bool_t>("passTrigger");
  isSS           = GetParam<Bool_t>("isSS");

  // Leptons and Jets
  GetLeptonVariables(selLeptons, vetoLeptons);
  GetGenJetVariables(genJets, mcJets);
  GetMET();
  GetWeights();
  GetJetVariables(selJets, Jets15); //quitar: esto estaba comentado (lo descomente pa minitrees de stop)

  fhDummy->Fill(1);
  if(gIsTTbar) FillCorrHistos(); 

  if(gIsTTbar && genLeptons.size() < 2) return; // Dilepton selection for ttbar!!!
  //if(gIsTTbar && genLeptons.size() >= 2) return; // Dilepton selection for ttbar!!! //quitar y dejar la de arriba

  // Number of events in fiducial region
  if(!gIsData && makeHistos) {
    if(genLeptons.size() >= 2){ // MIND THE POSSIBLE SKIM (on reco leptons) IN THE SAMPLE!!
      Int_t GenChannel = -1;
      if(genLeptons.at(0).isElec && genLeptons.at(1).isMuon) GenChannel = iElMu;
      if(genLeptons.at(0).isMuon && genLeptons.at(1).isElec) GenChannel = iElMu;
      if(genLeptons.at(0).isMuon && genLeptons.at(1).isMuon) GenChannel = iMuon;
      if(genLeptons.at(0).isElec && genLeptons.at(1).isElec) GenChannel = iElec;
      if( ( (genLeptons.at(0).p.Pt() > 25 && genLeptons.at(1).p.Pt() > 20) || (genLeptons.at(0).p.Pt() > 20 && genLeptons.at(1).p.Pt() > 25) )
          && (TMath::Abs(genLeptons.at(0).p.Eta()) < 2.4 && TMath::Abs(genLeptons.at(1).p.Eta()) < 2.4) 
          && ( (genLeptons.at(0).p + genLeptons.at(1).p).M() > 20 ) ){
        fHFiduYields[GenChannel-1][0] -> Fill(idilepton);
        if(GenChannel == iElMu || (TMath::Abs((genLeptons.at(0).p + genLeptons.at(1).p).M() - 91) > 15) ){
          fHFiduYields[GenChannel-1][0] -> Fill(iZVeto);
          if(GenChannel == iElMu || TGenMET > 40){   // MET > 40 in ee, µµ
            fHFiduYields[GenChannel-1][0] -> Fill(iMETcut);
            if(nFiduJets >= 2){ //At least 2 jets
              fHFiduYields[GenChannel-1][0] -> Fill(i2jets);
              if(nFidubJets >= 1){ // At least 1 b-tag
                fHFiduYields[GenChannel-1][0] -> Fill(i1btag);

              }
            }
          }
        }
      }
    }
  }

  // Event Selection
  // ===================================================================================================================
  if(gIsTTbar){
    SetVariables();
    FillCorrHistos();
  }
  if (TNSelLeps >= 2 && TPassTrigger && TPassMETFilters) {
    for(Int_t sys = 0; sys < nSyst; sys++){
      if(!gDoSyst && sys > 0) break;
      if(gIsData  && sys > 0) break;
     
      // Get values or the corresponding variation
      SetVariables(useSyst.at(sys));
      if (invmass > 20 && lep0pt > 25 && lep1pt > 20) {
        if(isSS) fHSSyields[gChannel][sys] -> Fill(idilepton, weight);
        else {
          fHyields[gChannel][sys] -> Fill(idilepton, weight);
          FillHistos(gChannel, idilepton, sys);
          if(sys == 0) FillDYHistos(gChannel); // Only once
        }

        if (TChannel == iElMu || (TMath::Abs(invmass - 91) > 15)  ){ //  Z Veto in ee, µµ
          if (isSS) fHSSyields[gChannel][sys] -> Fill(iZVeto, weight);
          else {      fHyields[gChannel][sys] -> Fill(iZVeto, weight);
            FillHistos(gChannel, iZVeto, sys);}

          if(TChannel == iElMu || met > 40){   // MET > 40 in ee, µµ
            if(isSS) fHSSyields[gChannel][sys] -> Fill(iMETcut, weight);
            else{      fHyields[gChannel][sys] -> Fill(iMETcut, weight);
              FillHistos(gChannel, iMETcut, sys);}

            if(njets > 1){ //At least 2 jets
              if(isSS) fHSSyields[gChannel][sys] -> Fill(i2jets, weight);
              else{      fHyields[gChannel][sys] -> Fill(i2jets, weight);
                FillHistos(gChannel, i2jets, sys); }

              if(nbtags > 0){ // At least 1 b-tag
                if(isSS) fHSSyields[gChannel][sys] -> Fill(i1btag, weight);
                else{   
                  fHyields[gChannel][sys] -> Fill(i1btag, weight);
                  FillHistos(gChannel, i1btag, sys);
                }
                if (!isSS && sys == 0 && makeTree && TChannel == iElMu) fTree->Fill();
              }
            }
          }
        }
      }
    }
  }
  SetParam("NJets",  njets);
  SetParam("NBtags", nbtags);
}



//#####################################################################
// Functions
//------------------------------------------------------------------
void TopAnalysis::GetLeptonVariables(std::vector<Lepton> selLeptons, std::vector<Lepton> VetoLeptons){
  TNSelLeps = selLeptons.size();
  Int_t nVetoLeptons = VetoLeptons.size();
  TNVetoLeps = (nVetoLeptons == 0) ? TNSelLeps : nVetoLeptons;
  if(TNSelLeps < 2) gChannel = -1;
  else if(selLeptons.at(0).isMuon && selLeptons.at(1).isElec) gChannel = iElMu;
  else if(selLeptons.at(0).isElec && selLeptons.at(1).isMuon) gChannel = iElMu;
  else if(selLeptons.at(0).isMuon && selLeptons.at(1).isMuon) gChannel = iMuon;
  else if(selLeptons.at(0).isElec && selLeptons.at(1).isElec) gChannel = iElec;
  if(TNSelLeps > 1){
    TMll = (selLeptons.at(0).p + selLeptons.at(1).p).M();
    TDilep_Pt = (selLeptons.at(0).p + selLeptons.at(1).p).Pt();
    TDeltaPhi = selLeptons.at(0).p.DeltaPhi(selLeptons.at(1).p);
    TDeltaEta = selLeptons.at(0).p.Eta() - selLeptons.at(1).p.Eta();
  }
  TChannel = gChannel;
  TIsSS = isSS;
  gChannel = gChannel -1; // gchannel used for chan index of histograms

  Int_t id0 = 0; Int_t id1 = 0;
  TLep0Pt = 0; TLep0Eta = 0; TLep0Phi = 0; TLep0M = 0; TLep0Id = 0;
  TLep1Pt = 0; TLep1Eta = 0; TLep1Phi = 0; TLep1M = 0; TLep1Id = 0;
  if(TNSelLeps >= 2){
    id0 = selLeptons.at(0).isMuon? 13 : 11;
    id1 = selLeptons.at(1).isMuon? 13 : 11;
    TLep0Pt  = selLeptons.at(0).Pt();
    TLep0Eta = selLeptons.at(0).Eta();
    TLep0Phi = selLeptons.at(0).Phi();
    TLep0M   = selLeptons.at(0).p.M();
    TLep0Id  = id0*selLeptons.at(0).charge;
    TLep1Pt  = selLeptons.at(1).Pt();
    TLep1Eta = selLeptons.at(1).Eta();
    TLep1Phi = selLeptons.at(1).Phi();
    TLep1M   = selLeptons.at(1).p.M();
    TLep1Id  = id1*selLeptons.at(1).charge;
  }
}

void TopAnalysis::GetJetVariables(std::vector<Jet> selJets, std::vector<Jet> cleanedJets15, Float_t ptCut){
  TNJets = selJets.size(); THT = 0;
  TNFwdJets   = vetoJets.size() - selJets.size();
  TNBtags = 0; TNBtagsBtagUp = 0; TNBtagsBtagDown = 0;
  TNBtagsMisTagUp = 0;  TNBtagsMisTagDown = 0;
  THTJESUp = 0; THTJESDown = 0;
  TBtagPt = 0;

  for(Int_t i = 0; i < TNJets; i++){
    TJet_Pt[i]     = selJets.at(i).Pt();
    TJet_Eta[i]    = selJets.at(i).Eta();
    TJet_Phi[i]    = selJets.at(i).Phi();
    TJet_M[i]      = selJets.at(i).p.M();
    TJet_Csv[i]    = selJets.at(i).csv;
    THT += selJets.at(i).p.Pt();
    if(selJets.at(i).isBtag){
      TNBtags++;
      if(TBtagPt == 0) TBtagPt = selJets.at(i).Pt();
    }
  }
  SetParam("THT",THT);
  TJet0Pt = 0; TJet0Eta = 0; TJet0Phi = 0; TJet0M = 0; TJet0Csv = -1; TJet0IsBTag = -1;
  TJet1Pt = 0; TJet1Eta = 0; TJet1Phi = 0; TJet1M = 0; TJet1Csv = -1; TJet1IsBTag = -1;
  if(TNJets >= 1){
    TJet0Pt  = selJets.at(0).Pt();
    TJet0Eta = selJets.at(0).Eta();
    TJet0Phi = selJets.at(0).Phi();
    TJet0M   = selJets.at(0).p.M();
    TJet0Csv = selJets.at(0).csv;
    TJet0IsBTag = selJets.at(0).isBtag;
  }
  if(TNJets >= 2){
    TJet1Pt  = selJets.at(1).Pt();
    TJet1Eta = selJets.at(1).Eta();
    TJet1Phi = selJets.at(1).Phi();
    TJet1M   = selJets.at(1).p.M();
    TJet1Csv = selJets.at(1).csv;
    TJet1IsBTag = selJets.at(1).isBtag;
  }

  if(gIsData) return;  // For systematics...
  for(Int_t i = 0; i < TNJets; i++){
    if(selJets.at(i).isBtag_BtagUp    ) TNBtagsBtagUp++;
    if(selJets.at(i).isBtag_BtagDown  ) TNBtagsBtagDown++;
    if(selJets.at(i).isBtag_MisTagUp  ) TNBtagsMisTagUp++;
    if(selJets.at(i).isBtag_MisTagDown) TNBtagsMisTagDown++;
  }
  TNJetsJESUp    = 0;
  TNJetsJESDown  = 0;
  TNBtagsJESUp    = 0;
  TNBtagsJESDown  = 0;
  TNJetsJERUp      = 0;  
  for(Int_t i = 0; i < (Int_t) cleanedJets15.size(); i++){
    if(cleanedJets15.at(i).pTJESUp > ptCut){
      THTJESUp += cleanedJets15.at(i).pTJESUp;
      TNJetsJESUp++;
      if(cleanedJets15.at(i).isBtag) TNBtagsJESUp++;
      TJetJESUp_Pt[i] = cleanedJets15.at(i).pTJESUp;
    }
    if(cleanedJets15.at(i).pTJESDown > ptCut){
      THTJESDown += cleanedJets15.at(i).pTJESDown;
      TNJetsJESDown++;
      if(cleanedJets15.at(i).isBtag) TNBtagsJESDown++;
      TJetJESDown_Pt[i] = cleanedJets15.at(i).pTJESDown;
    }
    if(cleanedJets15.at(i).pTJERUp > ptCut){
      TNJetsJERUp++;
      TJetJER_Pt[i] = cleanedJets15.at(i).pTJERUp;
    }
  }

}

void TopAnalysis::GetGenJetVariables(std::vector<Jet> genJets, std::vector<Jet> mcJets){
  if(gIsData) return;
  nFiduJets = 0; nFidubJets = 0; 
  Int_t nGenJets = genJets.size();
  Int_t nmcJets = mcJets.size();
  for(Int_t i = 0; i < nGenJets; i++) if(genJets.at(i).p.Pt() > 30 && TMath::Abs(genJets.at(i).p.Eta()) < 2.4)                         nFiduJets++;
  for(Int_t i = 0; i <  nmcJets; i++) if(mcJets.at(i).p.Pt()  > 30 && TMath::Abs(mcJets.at(i).Eta())    < 2.4 && mcJets.at(i).isBtag)  nFidubJets++;
}

void TopAnalysis::GetMET(){
    TRun        = gIsData ? Get<UInt_t>("run") : 1;
    TMET        = Get<Float_t>(metvarpt); // MET_pt
    TMET_Phi    = Get<Float_t>(metvar+"_phi");  // MET phi
    if((Int_t) selLeptons.size() >= 2) TMT2 = getMT2ll(selLeptons.at(0), selLeptons.at(1), TMET, TMET_Phi);
    TMETJESUp = 0; TMETJESDown = 0; TGenMET = 0; TgenTop1Pt = 0; TgenTop2Pt = 0;
    TMETJERUp = 0; TMETJERDown = 0; TMETUnclUp = 0; TMETUnclDown = 0;
    TMT2JESUp = 0; TMT2JESDown = 0; TMT2JERUp = 0; TMT2JERDown = 0; TMT2UnclUp = 0; TMT2UnclDown = 0;
    TMT2MESUp = 0; TMT2MESDown = 0; TMT2EESUp = 0; TMT2EESDown = 0;
    if(gIsData) TNVert = Get<Int_t>("PV_npvs");
    if(gIsData) return;
    TNVert      = Get<Int_t>("Pileup_nPU");
    TGenMET     = Get<Float_t>("GenMET_pt");
    //TMETJESUp    = Get<Float_t>("met_jecUp_pt"  );
    //TMETJESDown  = Get<Float_t>("met_jecDown_pt");
    //TMETJERUp    = Get<Float_t>("met_jecUp_pt"  );
    //TMETJERDown  = Get<Float_t>("met_jecUp_pt"  );
    //TMETUnclUp   = Get<Float_t>("met_jecUp_pt"  );
    //TMETUnclDown = Get<Float_t>("met_jecUp_pt"  );
}

void TopAnalysis::GetWeights(){
  TWeight_ElecEffUp = 1; TWeight_ElecEffDown = 1; TWeight_MuonEffUp = 1; TWeight_MuonEffDown = 1;
  TWeight_TrigUp = 1; TWeight_TrigDown = 1; TWeight_PUDown = 1; TWeight_PUUp = 1; TWeight = 1;
  TWeight_ISRUp   = 1; TWeight_ISRDown = 1; TWeight_FSRUp   = 1; TWeight_FSRDown = 1; TWeight_PrefUp = 1; TWeight_PrefDown = 1;
  if(gIsData) return;
  if(TNSelLeps < 2) return;
  Float_t lepSF   = selLeptons.at(0).GetSF( 0)*selLeptons.at(1).GetSF( 0);
  Float_t ElecSF = 1; Float_t MuonSF = 1;
  Float_t ElecSFUp = 1; Float_t ElecSFDo = 1; Float_t MuonSFUp = 1; Float_t MuonSFDo = 1;
  Float_t stat = 0; 
  Float_t fsrUp = 1; Float_t fsrDo = 1;  Float_t isrUp = 1; Float_t isrDo = 1; 

  if(gDoPSunc){
    // [0] is ISR=0.5 FSR=1; [1] is ISR=1 FSR=0.5; [2] is ISR=2 FSR=1; [3] is ISR=1 FSR=2
    isrDo = Get<Float_t>("PSWeight", 0);
    fsrDo = Get<Float_t>("PSWeight", 1);
    isrUp = Get<Float_t>("PSWeight", 2);
    fsrUp = Get<Float_t>("PSWeight", 3);
  }
  

  //For muons
  //https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonReferenceEffsRun2
  // OLD --> Additional 1% for ID + 0.5% for Isolation + 0.5% single muon triggers
  // 0.5% for iso (AN from Sergio)

  if(TChannel == iElec){
    ElecSF   = selLeptons.at(0).GetSF( 0)*selLeptons.at(1).GetSF( 0);
    ElecSFUp = selLeptons.at(0).GetSF( 1)*selLeptons.at(1).GetSF( 1);
    ElecSFDo = selLeptons.at(0).GetSF(-1)*selLeptons.at(1).GetSF(-1);
    MuonSFUp = 1; MuonSFDo = 1; MuonSF = 1;
  }
  else if(TChannel == iMuon){
    MuonSF   = selLeptons.at(0).GetSF( 0)*selLeptons.at(1).GetSF( 0);
    MuonSFUp = selLeptons.at(0).GetSF( 1)*selLeptons.at(1).GetSF( 1);
    MuonSFDo = selLeptons.at(0).GetSF(-1)*selLeptons.at(1).GetSF(-1);
    ElecSFUp = 1; ElecSFDo = 1; ElecSF = 1;
  }
  else{
    if(selLeptons.at(0).isMuon){
      MuonSF   *= selLeptons.at(0).GetSF( 0);
      MuonSFUp *= selLeptons.at(0).GetSF( 1)+MuonSF*0.005;
      MuonSFDo *= selLeptons.at(0).GetSF(-1)-MuonSF*0.005;
    }
    else{
      ElecSF   *= selLeptons.at(0).GetSF( 0);
      ElecSFUp *= selLeptons.at(0).GetSF( 1);
      ElecSFDo *= selLeptons.at(0).GetSF(-1);
    }
    if(selLeptons.at(1).isMuon){
      MuonSF   *= selLeptons.at(1).GetSF( 0);
      MuonSFUp *= selLeptons.at(1).GetSF( 1)+MuonSF*0.005;
      MuonSFDo *= selLeptons.at(1).GetSF(-1)-MuonSF*0.005;
    }
    else{
      ElecSF   *= selLeptons.at(1).GetSF( 0);
      ElecSFUp *= selLeptons.at(1).GetSF( 1);
      ElecSFDo *= selLeptons.at(1).GetSF(-1);
    }
  }
  TWeight             = NormWeight*ElecSF*MuonSF*TrigSF*PUSF*PrefWeight;
  TWeight_ElecEffUp   = NormWeight*ElecSFUp*MuonSF*TrigSF*PUSF*PrefWeight;
  TWeight_ElecEffDown = NormWeight*ElecSFDo*MuonSF*TrigSF*PUSF*PrefWeight;
  TWeight_MuonEffUp   = NormWeight*ElecSF*MuonSFUp*TrigSF*PUSF*PrefWeight;
  TWeight_MuonEffDown = NormWeight*ElecSF*MuonSFDo*TrigSF*PUSF*PrefWeight;
  TWeight_TrigUp     = NormWeight*lepSF*(TrigSF+TrigSFerr)*PUSF*PrefWeight;
  TWeight_TrigDown   = NormWeight*lepSF*(TrigSF-TrigSFerr)*PUSF*PrefWeight;
  TWeight_PUDown     = NormWeight*lepSF*TrigSF*PUSF_Up*PrefWeight;
  TWeight_PUUp       = NormWeight*lepSF*TrigSF*PUSF_Down*PrefWeight;
  TWeight_ISRUp      = NormWeight*ElecSF*MuonSF*TrigSF*PUSF*PrefWeight*isrUp;
  TWeight_ISRDown    = NormWeight*ElecSF*MuonSF*TrigSF*PUSF*PrefWeight*isrDo;
  TWeight_FSRUp      = NormWeight*ElecSF*MuonSF*TrigSF*PUSF*PrefWeight*fsrUp;
  TWeight_FSRDown    = NormWeight*ElecSF*MuonSF*TrigSF*PUSF*PrefWeight*fsrDo;
  TWeight_PrefUp     = NormWeight*ElecSF*MuonSF*TrigSF*PUSF*PrefWeightUp;
  TWeight_PrefDown   = NormWeight*ElecSF*MuonSF*TrigSF*PUSF*PrefWeightDo;
}

void TopAnalysis::InitHistos(){
  fhDummy = CreateH1F("fhDummy","fhDummy", 1, 0, 2);
  for(Int_t ch = 0; ch < nChannels; ch++){
    fhDummyCh[ch] = CreateH1F("fhDummy_"+gChanLabel[ch],"fhDummy_"+gChanLabel[ch], 1, 0, 2);
    fhDummy_2leps[ch] = CreateH1F("fhDummy_2leps"+gChanLabel[ch],"fhDummy_2leps"+gChanLabel[ch], 1, 0, 2);
    fhDummy_trigger[ch] = CreateH1F("fhDummy_trigger"+gChanLabel[ch],"fhDummy_trigger"+gChanLabel[ch], 1, 0, 2);
    fhDummy_metfilters[ch] = CreateH1F("fhDummy_metfilters"+gChanLabel[ch],"fhDummy_metfilters"+gChanLabel[ch], 1, 0, 2);
    fhDummy_OS[ch] = CreateH1F("fhDummy_OS"+gChanLabel[ch],"fhDummy_OS"+gChanLabel[ch], 1, 0, 2);
    fhDummy_minv[ch] = CreateH1F("fhDummy_minv"+gChanLabel[ch],"fhDummy_minv"+gChanLabel[ch], 1, 0, 2);
    fhDummy_lep0pt[ch] = CreateH1F("fhDummy_lep0pt"+gChanLabel[ch],"fhDummy_lep0pt"+gChanLabel[ch], 1, 0, 2);
    fhDummy_mz[ch] = CreateH1F("fhDummy_mz"+gChanLabel[ch],"fhDummy_mz"+gChanLabel[ch], 1, 0, 2);
    fhDummy_met[ch] = CreateH1F("fhDummy_met"+gChanLabel[ch],"fhDummy_met"+gChanLabel[ch], 1, 0, 2);
    fhDummy_njets[ch] = CreateH1F("fhDummy_njets"+gChanLabel[ch],"fhDummy_njets"+gChanLabel[ch], 1, 0, 2);
    fhDummy_nbtags[ch] = CreateH1F("fhDummy_nbtags"+gChanLabel[ch],"fhDummy_nbtags"+gChanLabel[ch], 1, 0, 2);

    TString suf; Int_t indexSyst;
    for(Int_t sys = 0; sys < nSyst; sys++){
      suf = gChanLabel[ch];
      indexSyst = useSyst.at(sys);
      if(sys > 0) suf += "_" + gSyst[indexSyst];
      fHyields[ch][sys]     = CreateH1F("H_Yields_"+suf,"", nLevels, -0.5, nLevels-0.5);
      fHFiduYields[ch][sys]     = CreateH1F("H_FiduYields_"+suf,"", nLevels, -0.5, nLevels-0.5);
      fHSSyields[ch][sys]   = CreateH1F("H_SSYields_"+suf,"", nLevels, -0.5, nLevels-0.5);
    }
  }
  if(!makeHistos) return;
  hJetPtReco  = CreateH1F("H_JetPtReco",  "H_JetPtReco",  nPtBins, (Float_t*) ptBins);
  hJetPtGen   = CreateH1F("H_JetPtGen",   "H_JetPtGen",   nPtBins, (Float_t*)ptBins);
  hJetPtRecoB = CreateH1F("H_JetPtRecoB", "H_JetPtRecoB", nPtBins, (Float_t*)ptBins);
  hJetPtGenB  = CreateH1F("H_JetPtGenB",  "H_JetPtGenB",  nPtBins, (Float_t*)ptBins);
  for(Int_t i = 0; i < nPtBins; i++){
    hJetGenRecoPtRatio[i]  = CreateH1F(Form("H_JetGenRecoPtRatio_%i", i),  Form("H_JetGenRecoPtRatio_%i", i),  30, 0.5, 1.5);
    hJetGenRecoPtRatio2[i] = CreateH1F(Form("H_JetGenRecoPtRatio2_%i", i), Form("H_JetGenRecoPtRatio2_%i", i), 30, 0.5, 1.5);
  }
  
  TString suffix;
  for(Int_t ch = 0; ch < nChannels; ch++){
    for(Int_t cut = 0; cut < nLevels; cut++){
      suffix = GetSuffix(ch, cut, 0);
      if(!gIsData){
        if(gDoPDFunc)   fHPDFweights[ch][cut]       = CreateH1F("H_PDFweights_"  +suffix, "PDFweights", nPDFweights, 0.5, nPDFweights+0.5); // 33 as nominal... 100 for 2016 old tune
        if(gDoScaleUnc) fHScaleWeights[ch][cut]     = CreateH1F("H_ScaleWeights_"+suffix, "ScaleWeights", 9, 0.5, 9.5);
        if(gDoPSunc)    fHPSweights[ch][cut]        = CreateH1F("H_PSweights_"   +suffix, "PSweights", 4, 0.5, 4.5);
      }
      for(Int_t sys = 0; sys < nSyst; sys++){
        if(gIsData && sys > 0) break;
        suffix = GetSuffix(ch, cut, sys);
        fHMET[ch][cut][sys]         = CreateH1F("H_MET_"        +suffix, "MET"       , 3000, 0,300);
        fHMT2[ch][cut][sys]         = CreateH1F("H_MT2_"        +suffix, "MT2"       , 3000, 0,300);
        fHMuonEta[ch][cut][sys]     = CreateH1F("H_MuonEta_"    +suffix, "Lep0Eta"   , 50  ,-2.5 ,2.5);
        fHElecEta[ch][cut][sys]     = CreateH1F("H_ElecEta_"    +suffix, "Lep1Eta"   , 50  ,-2.5 ,2.5);
        fHLep0Eta[ch][cut][sys]     = CreateH1F("H_Lep0Eta_"    +suffix, "Lep0Eta"   , 50  ,-2.5 ,2.5);
        fHLep1Eta[ch][cut][sys]     = CreateH1F("H_Lep1Eta_"    +suffix, "Lep1Eta"   , 50  ,-2.5 ,2.5);
        fHDelLepPhi[ch][cut][sys]   = CreateH1F("H_DelLepPhi_"  +suffix, "DelLepPhi" , 100, 0, 1);
        fHDelLepEta[ch][cut][sys]   = CreateH1F("H_DelLepEta_"  +suffix, "DelLepEta" ,  96,0, 4.8);
        fHHT[ch][cut][sys]          = CreateH1F("H_HT_"         +suffix, "HT"        , 4700,30,500);
        fHJet0Eta[ch][cut][sys]     = CreateH1F("H_Jet0Eta_"  +suffix, "Jet0Eta"   , 50,-2.5,2.5);
        fHJet1Eta[ch][cut][sys]     = CreateH1F("H_Jet1Eta_"  +suffix, "Jet1Eta"   , 50,-2.5,2.5);
        fHDYInvMass[ch][cut][sys]   = CreateH1F("H_DY_InvMass_"    +suffix, "InvMass"   ,  300,  0., 300.);
        fHDYInvMassSF[ch][cut][sys] = CreateH1F("H_DY_SF_InvMass_"    +suffix, "InvMass"   ,  300,  0., 300.);
        fHInvMass[ch][cut][sys]     = CreateH1F("H_InvMass_"    +suffix, "InvMass"   ,  300,  0., 300.);
        fHInvMass2[ch][cut][sys]    = CreateH1F("H_InvMass2_"   +suffix, "InvMass2"  ,  400, 70., 110.);
        fHNBtagsNJets[ch][cut][sys] = CreateH1F("H_NBtagsNJets_"+suffix, "NBtagsNJets"   ,15 , -0.5, 14.5);
        fHNJets[ch][cut][sys]       = CreateH1F("H_NJets_"      +suffix, "NJets"     , 8 ,-0.5, 7.5);
        fHNBtagJets[ch][cut][sys]   = CreateH1F("H_NBtagJets_"  +suffix, "NBtagJets" , 4 ,-0.5, 3.5);
        fHJet0Pt[ch][cut][sys]      = CreateH1F("H_Jet0Pt_"     +suffix, "Jet0Pt"    , 2700,30,300);
        fHJet1Pt[ch][cut][sys]      = CreateH1F("H_Jet1Pt_"     +suffix, "Jet1Pt"    , 2200,30,250);
        fHDiLepPt[ch][cut][sys]     = CreateH1F("H_DiLepPt_"    +suffix, "DiLepPt"   , 1600,20,180);
        fHMuonPt[ch][cut][sys]      = CreateH1F("H_MuonPt_"     +suffix, "Lep0Pt"    , 1800,20,200);
        fHElecPt[ch][cut][sys]      = CreateH1F("H_ElecPt_"     +suffix, "Lep1Pt"    , 1800,20,200);
        fHMuonIso[ch][cut][sys]     = CreateH1F("H_MuonIso_"     +suffix, "RelIso"    , 50,0,0.15);
        fHElecIso[ch][cut][sys]     = CreateH1F("H_ElecIso_"     +suffix, "RelIso"    , 50,0,0.15);
        fHLep0Pt[ch][cut][sys]      = CreateH1F("H_Lep0Pt_"     +suffix, "Lep0Pt"    , 1800,20,200);
        fHLep1Pt[ch][cut][sys]      = CreateH1F("H_Lep1Pt_"     +suffix, "Lep1Pt"    , 1800,20,200);
        fHLep0Iso[ch][cut][sys]     = CreateH1F("H_Lep0Iso_"     +suffix, "RelIso"    , 50,0,0.15);
        fHLep1Iso[ch][cut][sys]     = CreateH1F("H_Lep1Iso_"     +suffix, "RelIso"    , 50,0,0.15);
        fHJetEta[ch][cut][sys]      = CreateH1F("H_JetAllEta_"     +suffix, "JetAllEta"    , 50,-2.5,2.5);
        fHJetPt[ch][cut][sys]       = CreateH1F("H_JetAllPt_"     +suffix, "JetAllPt"    , 2700,30,300);
        fHJetCSV[ch][cut][sys]      = CreateH1F("H_JetAllCSV_" +suffix, "CSV" , 100,0, 1.0);
        fHJet0CSV[ch][cut][sys]     = CreateH1F("H_Jet0CSV_" +suffix, "Jet0CSV" , 100,0, 1.0);
        fHJet1CSV[ch][cut][sys]     = CreateH1F("H_Jet1CSV_" +suffix, "Jet1CSV" , 100,0, 1.0);
        fHJetDeepCSV[ch][cut][sys]  = CreateH1F("H_JetAllDeepCSV_" +suffix, "DeepCSV" , 100,0, 1.0);
        fHJet0DeepCSV[ch][cut][sys] = CreateH1F("H_Jet0DeepCSV_" +suffix, "Jet0DeepCSV" , 100,0, 1.0);
        fHJet1DeepCSV[ch][cut][sys] = CreateH1F("H_Jet1DeepCSV_" +suffix, "Jet1DeepCSV" , 100,0, 1.0);
        fHJetDeepFlav[ch][cut][sys]  = CreateH1F("H_JetAllDeepFlav_" +suffix, "DeepFlav" , 100,0, 1.0);
        fHJet0DeepFlav[ch][cut][sys] = CreateH1F("H_Jet0DeepFlav_" +suffix, "Jet0DeepFlav" , 100,0, 1.0);
        fHJet1DeepFlav[ch][cut][sys] = CreateH1F("H_Jet1DeepFlav_" +suffix, "Jet1DeepFlav" , 100,0, 1.0);
        fHvertices[ch][cut][sys]    = CreateH1F("H_Vtx_"+suffix, "Vtx", 101, -0.5, 100.5); 
      }
    }
  }
}

void TopAnalysis::FillDYHistos(Int_t ch){
  if(!makeHistos) return;
  Int_t sys = 0;
  Int_t cut;
  Float_t EventWeight = weight;
  cut = idilepton;
  fHDYInvMass[ch][cut][sys]       -> Fill(invmass, EventWeight);
  fHDYInvMassSF[ch][cut][sys]     -> Fill(invmass, EventWeight);
  cut = iZVeto;
  fHDYInvMass[ch][cut][sys]       -> Fill(invmass, EventWeight);
  fHDYInvMassSF[ch][cut][sys]     -> Fill(invmass, EventWeight);
  cut = iMETcut;
  if(met > 40) fHDYInvMassSF[ch][cut][sys]     -> Fill(invmass, EventWeight);
  fHDYInvMass[ch][cut][sys]       -> Fill(invmass, EventWeight);
  //  }
  if(njets > 1){
    cut = i2jets;
    if(met > 40) fHDYInvMassSF[ch][cut][sys]     -> Fill(invmass, EventWeight);
    fHDYInvMass[ch][cut][sys]       -> Fill(invmass, EventWeight);
    if(nbtags > 0){
      cut = i1btag;
      if(met > 40) fHDYInvMassSF[ch][cut][sys]     -> Fill(invmass, EventWeight);
      fHDYInvMass[ch][cut][sys]       -> Fill(invmass, EventWeight);
    }
  }
}

void TopAnalysis::FillHistos(Int_t ch, Int_t cut, Int_t sys){
  if(gIsData && sys != 0) return;
  if(!makeHistos) return;

  if(!gIsData && sys == 0){
    Int_t i = 0;
    if(gDoPSunc){
      for(i = 0; i < Get<Int_t>("nPSWeight"); i++)
        fHPSweights[ch][cut]->Fill(i+1, Get<Float_t>("PSWeight",i)*weight);
    }
    if(gDoPDFunc){
      for(i = 0; i < Get<Int_t>("nLHEPdfWeight"); i++)
        fHPDFweights[ch][cut]->Fill(i+1, Get<Float_t>("LHEPdfWeight",i)*weight);
    }
    if(gDoScaleUnc){
      for(i = 0; i < Get<Int_t>("nLHEScaleWeight"); i++)
        fHScaleWeights[ch][cut]->Fill(i+1, Get<Float_t>("LHEScaleWeight",i)*weight);
    }
  } 

  // Global
  fHMET[ch][cut][sys]         -> Fill(met, weight);
  fHMT2[ch][cut][sys]         -> Fill(mt2, weight);
  fHHT[ch][cut][sys]          -> Fill(ht, weight);
  fHNJets[ch][cut][sys]        -> Fill(njets, weight);
  fHNBtagJets[ch][cut][sys]    -> Fill(nbtags, weight);
  fHvertices[ch][cut][sys]      -> Fill(nvert, weight);

  // Leptons
  fHLep0Eta[ch][cut][sys]     -> Fill(lep0eta, weight);
  fHLep1Eta[ch][cut][sys]     -> Fill(lep1eta, weight);
  fHLep0Pt[ch][cut][sys]      -> Fill(lep0pt, weight);
  fHLep1Pt[ch][cut][sys]      -> Fill(lep1pt, weight);
  fHLep0Iso[ch][cut][sys]      -> Fill(lep0iso, weight);
  fHLep1Iso[ch][cut][sys]      -> Fill(lep1iso, weight);
  if(selLeptons.at(0).isMuon){
    fHMuonPt[ch][cut][sys]    -> Fill(lep0pt, weight);
    fHMuonEta[ch][cut][sys]   -> Fill(lep0eta, weight);
    fHMuonIso[ch][cut][sys]   -> Fill(lep0iso, weight);
  }
  else{
    fHElecPt[ch][cut][sys]    -> Fill(lep0pt, weight);
    fHElecEta[ch][cut][sys]   -> Fill(lep0eta, weight);
    fHElecIso[ch][cut][sys]   -> Fill(lep0iso, weight);
  }
  if(selLeptons.at(1).isMuon){
    fHMuonPt[ch][cut][sys]    -> Fill(lep1pt, weight);
    fHMuonEta[ch][cut][sys]   -> Fill(lep1eta, weight);
    fHMuonIso[ch][cut][sys]   -> Fill(lep1iso, weight);
  }
  else{
    fHElecPt[ch][cut][sys]    -> Fill(lep1pt, weight);
    fHElecEta[ch][cut][sys]   -> Fill(lep1eta, weight);
    fHElecIso[ch][cut][sys]   -> Fill(lep1iso, weight);
  }
  fHDiLepPt[ch][cut][sys]     -> Fill(dileppt, weight);
  fHDelLepPhi[ch][cut][sys]   -> Fill(TMath::Abs(deltaphi)/3.141592, weight);
  fHDelLepEta[ch][cut][sys]   -> Fill(TMath::Abs(deltaeta), weight);
  fHInvMass[ch][cut][sys]     -> Fill(invmass, weight);
  fHInvMass2[ch][cut][sys]    -> Fill(invmass, weight);

  // Jets
  if(njets > 0){ 
    fHJet0Eta[ch][cut][sys]     -> Fill(jets.at(0).Eta(), weight);
    fHJet0Pt [ch][cut][sys]     -> Fill(jets.at(0).Pt(), weight);
    fHJet0CSV[ch][cut][sys]     -> Fill(jets.at(0).csv, weight);
    fHJet0DeepCSV[ch][cut][sys] -> Fill(jets.at(0).GetDeepCSVB(), weight);
    fHJet0DeepFlav[ch][cut][sys] -> Fill(jets.at(0).GetDeepFlav(), weight);
  }
  if(njets > 1){
    fHJet1Eta[ch][cut][sys]     -> Fill(jets.at(1).Eta(), weight);
    fHJet1Pt [ch][cut][sys]     -> Fill(jets.at(1).Pt(), weight);
    fHJet1CSV[ch][cut][sys]     -> Fill(jets.at(1).csv, weight);
    fHJet1DeepCSV[ch][cut][sys] -> Fill(jets.at(1).GetDeepCSVB(), weight);
    fHJet1DeepFlav[ch][cut][sys] -> Fill(jets.at(1).GetDeepFlav(), weight);
  }
  for(Int_t i = 0; i < (Int_t) jets.size(); i++){ 
    fHJetPt[ch][cut][sys]  ->Fill(jets.at(i).Pt());
    fHJetEta[ch][cut][sys] ->Fill(jets.at(i).Eta());
    fHJetCSV[ch][cut][sys] -> Fill(jets.at(i).csv, weight);
    fHJetDeepCSV[ch][cut][sys] -> Fill(jets.at(i).GetDeepCSVB(), weight);
    fHJetDeepFlav[ch][cut][sys] -> Fill(jets.at(i).GetDeepFlav(), weight);
  }
  // Btag histos
  /*
  bool firstBtag = false;
  for(Int_t i = 0; i < (Int_t) jets.size(); i++){
    if(jets.at(i).isBtag){
      if(!firstBtag){
        fHJetBtag0Pt[ch][cut][sys]  ->Fill(jets.at(i).Pt(), weight);
        fHJetBtag0Eta[ch][cut][sys] ->Fill(jets.at(i).Eta(), weight);
      }
      firstBtag = true;
      fHJetBtagPt[ch][cut][sys]  ->Fill(jets.at(i).Pt(), weight);
      fHJetBtagEta[ch][cut][sys] ->Fill(jets.at(i).Eta(), weight);
    }
  }*/

  // NBtagNJets
  if(njets == 0) fHNBtagsNJets[ch][cut][sys]   -> Fill(nbtags,    weight);
  if(njets == 1) fHNBtagsNJets[ch][cut][sys]   -> Fill(nbtags+1,  weight);
  if(njets == 2) fHNBtagsNJets[ch][cut][sys]   -> Fill(nbtags+3,  weight);
  if(njets == 3) fHNBtagsNJets[ch][cut][sys]   -> Fill(nbtags+6,  weight);
  if(njets == 4) fHNBtagsNJets[ch][cut][sys]   -> Fill(nbtags+10, weight);
}

void TopAnalysis::FillCorrHistos(){
  if(!makeHistos) return;
  Float_t pTreco; Float_t pTgen; Float_t dR = 0.3; Float_t dPtoPt; Bool_t isBtag; Bool_t isBJet;
  Float_t bin0; Float_t bin1;
  Int_t njets = jets.size();
  for(Int_t i = 0; i < njets; i++){
    pTreco = jets.at(i).Pt(); pTgen = jets.at(i).GetGenPt();
    dPtoPt = fabs(pTreco - pTgen)/pTreco;
    isBtag = jets.at(i).isBtag; isBJet = jets.at(i).flavmc == 5;

    hJetPtReco ->Fill(pTreco);
    hJetPtGen  ->Fill(pTgen);
    if(isBtag) hJetPtRecoB->Fill(pTreco);
    if(isBJet) hJetPtGenB ->Fill(pTreco);
    for(Int_t b = 0; b < nPtBins; b++){
      bin0 = ptBins[b]; bin1 = ptBins[b+1];
      if(pTreco >= bin0 && pTreco < bin1){
        hJetGenRecoPtRatio[b]->Fill(pTreco/pTgen); 
        if(dR < 0.2) hJetGenRecoPtRatio2[b]->Fill(pTreco/pTgen);
      }
    }
  }
}

void TopAnalysis::SetLeptonVariables(){
  //fTree->Branch("TNVetoLeps",   &TNVetoLeps,   "TNVetoLeps/I");
  if(!miniTree){
    fTree->Branch("TChannel",     &TChannel,     "TChannel/I");
    fTree->Branch("TIsSS",        &TIsSS,        "TIsSS/B");
    fTree->Branch("TLep0Phi",     &TLep0Phi,     "TLep0Phi/F");
    fTree->Branch("TLep0M",       &TLep0M,       "TLep0M/F");
    fTree->Branch("TLep0Id",      &TLep0Id,      "TLep0Id/I");
    fTree->Branch("TLep1Phi",     &TLep1Phi,     "TLep1Phi/F");
    fTree->Branch("TLep1M",       &TLep1M,       "TLep1M/F");
    fTree->Branch("TLep1Id",      &TLep1Id,      "TLep1Id/I");
  }
    fTree->Branch("TMll",         &TMll,         "TMll/F");
    fTree->Branch("TDilep_Pt",    &TDilep_Pt,    "TDilep_Pt/F");
    fTree->Branch("TLep0Pt",      &TLep0Pt,      "TLep0Pt/F");
    fTree->Branch("TLep0Eta",     &TLep0Eta,     "TLep0Eta/F");
    fTree->Branch("TNSelLeps",    &TNSelLeps,    "TNSelLeps/I");
    fTree->Branch("TDeltaPhi",    &TDeltaPhi,    "TDeltaPhi/F");
    fTree->Branch("TDeltaEta",    &TDeltaEta,    "TDeltaEta/F");
    fTree->Branch("TLep1Pt",      &TLep1Pt,      "TLep1Pt/F");
    fTree->Branch("TLep1Eta",     &TLep1Eta,     "TLep1Eta/F");
}

void TopAnalysis::SetJetVariables(){
  if(!miniTree){
    //fTree->Branch("TNFwdJets",     &TNFwdJets,   "TNFwdJets/I");

    /*fTree->Branch("TJet_Csv",      TJet_Csv,     "TJet_Csv[TNJets]/F");
      fTree->Branch("TJet_Pt",       TJet_Pt,      "TJet_Pt[TNJets]/F");
      fTree->Branch("TJet_Eta",      TJet_Eta,     "TJet_Eta[TNJets]/F");
      fTree->Branch("TJet_Phi",      TJet_Phi,     "TJet_Phi[TNJets]/F");
      fTree->Branch("TJet_M",        TJet_M,       "TJet_M[TNJets]/F");*/
    fTree->Branch("TJet0Phi",      &TJet0Phi,    "TJet0Phi/F");
    fTree->Branch("TJet0M",        &TJet0M,      "TJet0M/F");
    fTree->Branch("TJet0Csv",      &TJet0Csv,    "TJet0Csv/F");
    fTree->Branch("TJet0IsBTag",   &TJet0IsBTag, "TJet0IsBTag/I");
    fTree->Branch("TJet1Phi",      &TJet1Phi,    "TJet1Phi/F");
    fTree->Branch("TJet1M",        &TJet1M,      "TJet1M/F");
    fTree->Branch("TJet1Csv",      &TJet1Csv,    "TJet1Csv/F");
    fTree->Branch("TJet1IsBTag",   &TJet1IsBTag, "TJet1IsBTag/I");
    fTree->Branch("TBtagPt",       &TBtagPt,     "TBtagPt/F");

    fTree->Branch("TNJetsJESUp",           &TNJetsJESUp,         "TNJetsJESUp/I");
    fTree->Branch("TNJetsJESDown",           &TNJetsJESDown,         "TNJetsJESDown/I");
    fTree->Branch("TNJetsJERUp",           &TNJetsJERUp,         "TNJetsJERUp/I");

    fTree->Branch("TNBtagsBtagUp",     &TNBtagsBtagUp,   "TNBtagsBtagUp/I");
    fTree->Branch("TNBtagsBtagDown",   &TNBtagsBtagDown, "TNBtagsBtagDown/I");
    fTree->Branch("TNBtagsMisTagUp",     &TNBtagsMisTagUp,   "TNBtagsMisTagUp/I");
    fTree->Branch("TNBtagsMisTagDown",   &TNBtagsMisTagDown, "TNBtagsMisTagDown/I");

    fTree->Branch("TNBtagsJESUp",   &TNBtagsJESUp, "TNBtagsJESUp/I");
    fTree->Branch("TNBtagsJESDown",  &TNBtagsJESDown, "TNBtagsJESDown/I");
    /*fTree->Branch("TJetJESUp_Pt",      TJetJESUp_Pt,      "TJetJESUp_Pt[TNJetsJESUp]/F");
      fTree->Branch("TJetJESDown_Pt",    TJetJESDown_Pt,    "TJetJESDown_Pt[TNJetsJESDown]/F");
      fTree->Branch("TJetJER_Pt",        TJetJER_Pt,        "TJetJER_Pt[TNJetsJERUp]/F");*/
    fTree->Branch("THTJESUp",     &THTJESUp,     "THTJESUp/F");
    fTree->Branch("THTJESDown",   &THTJESDown,   "THTJESDown/F");
  }

    fTree->Branch("TJet0Pt",       &TJet0Pt,     "TJet0Pt/F");
    fTree->Branch("TJet0Eta",      &TJet0Eta,    "TJet0Eta/F");
    fTree->Branch("TJet1Pt",       &TJet1Pt,     "TJet1Pt/F");
    fTree->Branch("TJet1Eta",      &TJet1Eta,    "TJet1Eta/F");
    fTree->Branch("THT",          &THT,          "THT/F");
    fTree->Branch("TNJets",        &njets,      "TNJets/I");
    fTree->Branch("TNBtags",       &nbtags,     "TNBtags/I");
}

void TopAnalysis::SetEventVariables(){
  if(!miniTree){
    fTree->Branch("TWeight_LepEffUp",      &TWeight_LepEffUp,      "TWeight_LepEffUp/F");
    fTree->Branch("TWeight_LepEffDown",    &TWeight_LepEffDown,    "TWeight_LepEffDown/F");
    fTree->Branch("TWeight_ElecEffUp",      &TWeight_ElecEffUp,      "TWeight_ElecEffUp/F");
    fTree->Branch("TWeight_ElecEffDown",    &TWeight_ElecEffDown,    "TWeight_ElecEffDown/F");
    fTree->Branch("TWeight_MuonEffUp",      &TWeight_MuonEffUp,      "TWeight_MuonEffUp/F");
    fTree->Branch("TWeight_MuonEffDown",    &TWeight_MuonEffDown,    "TWeight_MuonEffDown/F");
    fTree->Branch("TWeight_TrigUp",        &TWeight_TrigUp,        "TWeight_TrigUp/F");
    fTree->Branch("TWeight_TrigDown",      &TWeight_TrigDown,      "TWeight_TrigDown/F");
    fTree->Branch("TWeight_PUUp",        &TWeight_PUUp,        "TWeight_PUUp/F");
    fTree->Branch("TWeight_PUDown",        &TWeight_PUDown,        "TWeight_PUDown/F");

    fTree->Branch("TEvent",          &event,           "TEvent/l");
    fTree->Branch("TLuminosityBlock",&lumiblock,       "TLuminosityBlock/i");
    fTree->Branch("TPassMETFilters", &TPassMETFilters, "TPassMETFilters/B");
    fTree->Branch("TPassTrigger",    &TPassTrigger,    "TPassTrigger/B");
    fTree->Branch("TRun",            &TRun,            "TRun/i");
    fTree->Branch("TNVert",          &TNVert,          "TNVert/I");
    //fTree->Branch("TGenMET",         &TGenMET,         "TGenMET/F");
    //fTree->Branch("TgenTop1Pt",   &TgenTop1Pt,   "TgenTop1Pt/F");
    //fTree->Branch("TgenTop2Pt",   &TgenTop2Pt,   "TgenTop2Pt/F");
    fTree->Branch("TMET_Phi",     &TMET_Phi,     "TMET_Phi/F");
    fTree->Branch("TMETJESUp",    &TMETJESUp,    "TMETJESUp/F");
    fTree->Branch("TMETJESDown",  &TMETJESDown,  "TMETJESDown/F");
  }
  fTree->Branch("TWeight",      &TWeight,      "TWeight/F");
  fTree->Branch("TMET",            &TMET,            "TMET/F");
  fTree->Branch("TMT2",            &TMT2,            "TMT2/F");

  if (gIsSignal){
    fTree->Branch("Tm_stop", &m_stop, "Tm_stop/F");  //quitar
    fTree->Branch("Tm_LSP", &m_LSP, "Tm_LSP/F");  //quitar
  }
}

TString TopAnalysis::GetSuffix(int iCh, int iCut, int iSyst){
   Int_t indexSyst = useSyst.at(iSyst);
   TString t = iSyst > 0 ? gChanLabel[iCh]+"_"+sCut[iCut]+"_"+gSyst[indexSyst] : gChanLabel[iCh]+"_"+sCut[iCut];
   return t;
}

void TopAnalysis::SetVariables(int sys){
  // Global variables to fill the histos where
  // the value changes according to the variation
  // that you choose
  // Global
  nleps = selLeptons.size(); weight = TWeight;
  met = TMET; mt2 = TMT2; ht = THT; nvert = TNVert; invmass = TMll;
  // Leptons
  if(nleps >= 2){
    lep0pt = TLep0Pt; lep1pt = TLep1Pt; lep0eta = TLep0Eta; lep1eta = TLep1Eta; 
    lep0iso = selLeptons.at(0).GetIso(); lep1iso = selLeptons.at(1).GetIso();
    dileppt = (selLeptons.at(0).p + selLeptons.at(1).p).Pt();
    deltaphi = selLeptons.at(0).p.DeltaPhi(selLeptons.at(1).p);
    deltaeta = selLeptons.at(0).p.Eta() - selLeptons.at(1).p.Eta();
  }

  // Jets
  njets = 0; nbtags = 0;
  jets.clear(); Jet jet; TLorentzVector t; 
  float pt, eta, phi, m; bool isbtag; float csv, deepcsv, deepflav;
  Int_t nJets = Get<Int_t>("nJet");
  Int_t jetid, flav;
  for(int i = 0; i < nJets; i++){
    pt = Get<Float_t>(JetPt,i); eta = Get<Float_t>("Jet_eta",i); phi = Get<Float_t>("Jet_phi", i); m = Get<Float_t>("Jet_mass",i);
    csv = Get<Float_t>("Jet_btagCSVV2", i); deepcsv = Get<Float_t>("Jet_btagDeepB", i); deepflav = Get<Float_t>("Jet_btagDeepFlavB", i);
    jetid = Get<Int_t>("Jet_jetId",i);
    Int_t nGenJets = Get<Int_t>("nGenJet");
    Int_t genJetIdx = Get<Int_t>("Jet_genJetIdx", i);
    Float_t genPt = -999;
    if(genJetIdx > 0 and genJetIdx < nJets) genPt = Get<Float_t>("GenJet_pt", genJetIdx);
    flav = -999999; if(!gIsData) flav = Get<Int_t>("Jet_hadronFlavour", i);

    // JES and JER ----> Update this to propagate to MET, MT2, etc!!
    // MET_pt_jerUp, MET_phi_jerUp, MET_pt_jesTotalUp, MET_phi_jesTotalUp 
    if     (sys == kJESUp){
      pt = Get<Float_t>("Jet_pt_jesTotalUp", i);
      m  = Get<Float_t>("Jet_mass_jesTotalUp", i);
    }
    else if(sys == kJESDown){
      pt = Get<Float_t>("Jet_pt_jesTotalDown", i);
      m  = Get<Float_t>("Jet_mass_jesTotalDown", i);
    }
    else if(sys == kJERUp){
      pt = Get<Float_t>("Jet_pt_jerUp", i);
      m  = Get<Float_t>("Jet_mass_jerUp", i);
    }
    else if(sys == kJERDown){
      pt = Get<Float_t>("Jet_pt_jerDown", i);
      m  = Get<Float_t>("Jet_mass_jerDown", i);
    }
    Float_t alg = deepflav; 
    if     (sys == kBtagUp)      isbtag = fBTagSFbUp->IsTagged(alg, flav, pt, eta);
    else if(sys == kBtagDown)    isbtag = fBTagSFbDo->IsTagged(alg, flav, pt, eta);
    else if(sys == kMistagUp)    isbtag = fBTagSFlUp->IsTagged(alg, flav, pt, eta);
    else if(sys == kMistagDown)  isbtag = fBTagSFlDo->IsTagged(alg, flav, pt, eta);
    else                         isbtag = fBTagSFnom->IsTagged(alg, flav, pt, eta);
    if(pt > 30 && TMath::Abs(eta) <= 2.4 && jetid > 1){ // pt > 30, |eta| < 2.4, id > 1 (2, 6)
      t.SetPtEtaPhiM(pt, eta, phi, m);
      jet = Jet(t, csv);
      jet.isBtag = isbtag;
      jet.SetDeepCSVB(deepcsv);
      jet.SetDeepFlav(deepflav);
      jet.SetGenPt(genPt);

      if(Cleaning(jet, selLeptons, 0.4)){
        njets++;
        if(jet.isBtag) nbtags++; 
        jets.push_back(jet);
      }
    }
  }
  TNJets = njets; TNBtags = nbtags;

  if     (sys == kMuonEffUp  ) weight = TWeight_MuonEffUp;
  else if(sys == kMuonEffDown) weight = TWeight_MuonEffDown;
  else if(sys == kElecEffUp  ) weight = TWeight_ElecEffUp;
  else if(sys == kElecEffDown) weight = TWeight_ElecEffDown;
  else if(sys == kTrigUp     ) weight = TWeight_TrigUp;
  else if(sys == kTrigDown   ) weight = TWeight_TrigDown;
  else if(sys == kPUUp       ) weight = TWeight_PUUp;   
  else if(sys == kPUDown     ) weight = TWeight_PUDown;
  else if(sys == kISRUp      ) weight = TWeight_ISRUp;
  else if(sys == kISRDown    ) weight = TWeight_ISRDown;
  else if(sys == kFSRUp      ) weight = TWeight_FSRUp;
  else if(sys == kFSRDown    ) weight = TWeight_FSRDown;
  else if(sys == kPrefireUp  ) weight = TWeight_PrefUp;
  else if(sys == kPrefireDown) weight = TWeight_PrefDown;
}

  
