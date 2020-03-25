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
        fHDelLepR[ch][cut][sys]   = 0;
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
        fHInvMass2EE[ch][cut][sys]      = 0;
        fHInvMass2BB[ch][cut][sys]      = 0;
        fHInvMass2EB[ch][cut][sys]      = 0;
        fHInvMass2BE[ch][cut][sys]      = 0;
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
        fHvertices_pu[ch][cut][sys]     = 0;
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
  sampString   = GetParam<TString>("sampString");
  gSampleName  = GetParam<TString>("sampleName");
  gOptions     = GetParam<TString>("_options");
  gDoSyst      = true;// gOptions.Contains("doSyst")? true : false;
  year         = GetParam<TString>("year").Atoi();
  gIsTTbar     = false;
  gIsTTany     = false;
  gSelection   = GetSelection(selection);
  gDoElecES    = false;
  gDoMuonES    = gOptions.Contains("MuonES")? true : false;
  gDoJECunc    = gOptions.Contains("JECunc")? true : false;
  gDoPDFunc    = gOptions.Contains("PDF")? true : false;
  gDoPSunc     = gOptions.Contains("PS")? true : false;
  gDoScaleUnc  = gOptions.Contains("Scale")? true : false;
  gPUWeigth    = gOptions.Contains("PUweight")? true : false;
  gPrefire     = gOptions.Contains("prefire")? true : false;
  JetPt        = gOptions.Contains("JetPtNom")? "Jet_pt" : "Jet_pt";
  if ((gSampleName == "TT" || gSampleName.BeginsWith("TT_")) && year == 2016) gIsTTbar = true;
  if (gIsTTbar || gSampleName.BeginsWith("TTTo2L2Nu")) gIsTTany = true;

  nPDFweights = gIsTTbar ? 100 : 33;

  makeTree   = true;
  miniTree = true;
  makeHistos = true;

  gIsSignal= false; //quitar
	if ((gSampleName.BeginsWith("stop"))) gIsSignal = true; //quitar
  if(makeTree){
    fTree   = CreateTree("MiniTree","Created with PAF");
    SetLeptonVariables();
    SetJetVariables();
    SetEventVariables();
  }
  
  
  if(gIsSignal){
	  TString path = GetParam<TString>("path");
    snorm = new SUSYnorm(path ,  sampString);// "SMS_T2tt_3J_xqcut_20_top_corridor_2Lfilter_TuneCP5_MLM_p, SMS_T2tt_3J_xqcut_20_top_corridor_2Lfilter_TuneCUETP8M2T4_madgra");
  }
     
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
    if(gIsTTany){
      useSyst.push_back(kTopPt);
    }
  }
  nSyst = useSyst.size();
  InitHistos();
  metvar    = year == 2017? "METFixEE2017"     : "MET";
  metvarphi = year == 2017? "METFixEE2017_phi" : "MET_phi_nom";
  metvarpt  = year == 2017? "METFixEE2017_pt_nom"  : "MET_pt_nom";
  if(gIsData){
    metvarphi = "MET_phi";
    metvarpt  = "MET_pt";
  }
  //if(year != 2017 and gOptions.Contains("JetPtNom")){
  //  metvarpt   = "MET_pt_nom";
  //  metvarphi  = "MET_phi_nom";
  //}
  //metvar="MET";
  //metvarpt="MET_pt";
  
  
  
  gIs2017 = false; gIs2016 = false; gIs2018 = false;
  if     (year == 2017) gIs2017 = true;
  else if(year == 2018) gIs2018 = true;
  else if(year == 2016) gIs2016 = true;
  // b tagging
  TString pwd  = GetParam<TString>("WorkingDir");
  TString BTagSFPath = Form("%s/packages/BTagSFUtil", pwd.Data());
  TString taggerName="DeepFlav"; //"CSVv2"; //"DeepCSV"; // DeepFlav   //quitar dejar DeepFlav
  TString MeasType = "mujets"; //quitar: dejar mujets
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
  TTop0Pt = 0; TTop0Eta = 0; TTop0Phi = 0;
  TTop1Pt = 0; TTop1Eta = 0; TTop1Phi = 0;
  if (gIsSignal || gIsTTany){
    ngenPart=Get<Int_t>("nGenPart");
    int i; iSt = 0; iLSP = 0;
    for(Int_t i = 0; i < ngenPart; i++){
      Int_t genpdgId = Get<Int_t>("GenPart_pdgId", i);
      if     (genpdgId == -6 and Get<Int_t>("GenPart_status", i) == 22){
        TTop0Pt  = Get<Float_t>("GenPart_pt",  i);
        TTop0Eta = Get<Float_t>("GenPart_eta", i);
        TTop0Phi = Get<Float_t>("GenPart_phi", i);
      }
      else if(genpdgId ==  6 and Get<Int_t>("GenPart_status", i) == 22){
        TTop1Pt  = Get<Float_t>("GenPart_pt",  i);
        TTop1Eta = Get<Float_t>("GenPart_eta", i);
        TTop1Phi = Get<Float_t>("GenPart_phi", i);
      }
      if(genpdgId == 1000006)  iSt=i; 
      if(genpdgId == 1000022)  iLSP=i; 
    }
    if(gIsSignal){
      m_stop = Get<Float_t>("GenPart_mass", iSt);
      m_LSP  = Get<Float_t>("GenPart_mass", iLSP);
    }
  }

  // Vectors with the objects
  genLeptons     = GetParam<vector<Lepton>>("genLeptons");
  selLeptons     = GetParam<vector<Lepton>>("selLeptons");
  vetoLeptons    = GetParam<vector<Lepton>>("vetoLeptons");
  selJets        = GetParam<vector<Jet>>("selJets");
  selJetsJecUp   = GetParam<vector<Jet>>("selJetsJecUp");
  selJetsJecDown = GetParam<vector<Jet>>("selJetsJecDown");
  selJetsJERUp   = GetParam<vector<Jet>>("selJetsJERUp");
  selJetsJERDown = GetParam<vector<Jet>>("selJetsJERDown");
  selJetsJecCorUp   = GetParam<vector<Jet>>("selJetsJecCorUp");
  selJetsJecCorDown = GetParam<vector<Jet>>("selJetsJecCorDown");
  selJetsJecUnCorUp   = GetParam<vector<Jet>>("selJetsJecUnCorUp");
  selJetsJecUnCorDown = GetParam<vector<Jet>>("selJetsJecUnCorDown");
  //Jets15         = GetParam<vector<Jet>>("Jets15");
  //vetoJets       = GetParam<vector<Jet>>("vetoJets");
  //genJets        = GetParam<vector<Jet>>("genJets");
  //mcJets         = GetParam<vector<Jet>>("mcJets");
  fhDummy->Fill(1);

  if( (int)selLeptons.size() < 2) return;
  
  // Weights and SFs
  NormWeight     = GetParam<Double_t>("NormWeight");
  TrigSF         = GetParam<Float_t>("TriggerSF");
  TrigSFerr      = GetParam<Float_t>("TriggerSFerr");

  if(gIsSignal){
	  Float_t norm = snorm->GetSUSYnorm(m_stop,m_LSP);
	  xsec = snorm->GetStopXSec(m_stop)*(3*0.108)*(3*0.108);
	  NormWeight = norm != 0 ? xsec/norm*Get<Float_t>("genWeight") : 0;
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
  //PUSF = 1; PUSF_Up = 1; PUSF_Down = 1;

  // Event variables
  gChannel       = GetParam<Int_t>("gChannel");
  TPassMETFilters = GetParam<Bool_t>("METfilters");
  TPassTrigger    = GetParam<Bool_t>("passTrigger");
  isSS           = GetParam<Bool_t>("isSS");
  TIsHEM          = GetParam<Bool_t>("isHEM");

  // Leptons and Jets
  GetLeptonVariables(selLeptons, vetoLeptons);
  GetGenJetVariables(genJets, mcJets);
  GetMET();
  GetWeights();
  GetJetVariables(selJets, Jets15); //quitar: esto estaba comentado (lo descomente pa minitrees de stop)

  if(gIsTTbar && makeHistos) FillCorrHistos(); 

  if(gIsTTbar && genLeptons.size() < 2) return; // Dilepton selection for ttbar!!! 
  //if(gIsTTbar && genLeptons.size() >= 2) return; // Semilep selection for ttbar!!! //quitar y dejar la de arriba

  //if(gIsSignal){
  //  if(m_stop != 205) return;
  //  if(m_LSP  !=  20) return;
  //}

  // GenLeptons
  TGenLep0Pt = 0; TGenLep0Eta = 0; TGenLep0Phi = 0;
  TGenLep1Pt = 0; TGenLep1Eta = 0; TGenLep1Phi = 0;

  TLorentzVector gLep0 = TLorentzVector();
  TLorentzVector gLep1 = TLorentzVector();

  if(genLeptons.size() >= 1){
    TGenLep0Pt  = genLeptons.at(0).Pt(); 
    TGenLep0Eta = genLeptons.at(0).Eta(); 
    TGenLep0Phi = genLeptons.at(0).Phi(); 
  }
  if(genLeptons.size() >= 2){
    TGenLep1Pt  = genLeptons.at(1).Pt(); 
    TGenLep1Eta = genLeptons.at(1).Eta(); 
    TGenLep1Phi = genLeptons.at(1).Phi(); 
  }
  if(!gIsData){
    TGenMET = Get<Float_t>("GenMET_pt"); TGenMET_phi = Get<Float_t>("GenMET_phi");
    TGenMT2 = genLeptons.size() >= 2 ? getMT2ll(genLeptons.at(0), genLeptons.at(1), TGenMET, TGenMET_phi) : -1;    
}

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
    if(makeHistos) FillCorrHistos();
  }

  if (TNSelLeps >= 2 && TPassTrigger && TPassMETFilters) { 
  //if (TPassTrigger && TPassMETFilters) {
    for(Int_t sys = 0; sys < nSyst; sys++){
      if(!gDoSyst && sys > 0) break;
      if(gIsData  && sys > 0) break;
      // Get values or the corresponding variation

      SetVariables(useSyst.at(sys));
      //if (sys == 0 && makeTree && TChannel == iElMu && TPassDilepAny && TPassJetsAny && TPassBtagAny && TPassMETAny && TPassMT2Any) fTree->Fill();
      if (sys == 0 && makeTree && TChannel == iElMu && TPassDilepAny && TPassJetsAny && TPassBtagAny) fTree->Fill(); //quitar:sincro, dejar esta o la otra  
      //if (sys == 0 && makeTree && TChannel == iElMu && TPassDilepAny) fTree->Fill();
      
      // ee, mm
      //if (sys == 0 && makeTree && TChannel == iElec && !TIsOnZ && TPassDilepAny && TPassJetsAny && TPassBtagAny) fTree->Fill(); //quitar:sincro, dejar esta o la otra  
      //if (sys == 0 && makeTree && TChannel == iMuon && !TIsOnZ && TPassDilepAny && TPassJetsAny && TPassBtagAny) fTree->Fill(); //quitar:sincro, dejar esta o la otra  
      
      //if (!isSS) fHyields[gChannel][sys] -> Fill(iZVeto, weight); //quitar: sincro
      if (invmass > 20 && lep0pt > 25 && lep1pt > 20) {

        if(isSS) fHSSyields[gChannel][sys] -> Fill(idilepton, weight);
        else {
          fHyields[gChannel][sys] -> Fill(idilepton, weight);
          FillHistos(gChannel, idilepton, sys);
          if(sys == 0) FillDYHistos(gChannel); // Only once
        }

        if (TChannel == iElMu || (TMath::Abs(invmass - 91) > 15)  ){ //  Z Veto in ee, µµ

          if (isSS) fHSSyields[gChannel][sys] -> Fill(iZVeto, weight);
          else {      fHyields[gChannel][sys] -> Fill(iZVeto, weight); //quitar: sincro, descomentar
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

                //if (!isSS && sys == 0 && makeTree && TChannel == iElMu) fTree->Fill();
                //if (!isSS && sys == 0 && makeTree && TChannel == iElMu && met >= 50 && mt2 >=80) fTree->Fill();
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
  TIsOnZ = false;
  gChannel = gChannel -1; // gchannel used for chan index of histograms

  Int_t id0 = 0; Int_t id1 = 0; Int_t indexMuon = 0; Int_t indexElec = 1;
  TLep0Pt = 0; TLep0Eta = 0; TLep0Phi = 0; TLep0M = 0; TLep0Id = 0; TLep0Phi = 0;
  TLep1Pt = 0; TLep1Eta = 0; TLep1Phi = 0; TLep1M = 0; TLep1Id = 0; TLep1Phi = 0;
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
    indexMuon = id0 == 13? 0 : 1;
    indexElec = id0 == 13? 1 : 0;
    TMuonPt  = selLeptons.at(indexMuon).Pt();
    TMuonEta = selLeptons.at(indexMuon).Eta();
    TMuonPhi = selLeptons.at(indexMuon).Phi();
    TMuonM   = selLeptons.at(indexMuon).M();
    TElecPt  = selLeptons.at(indexElec).Pt();
    TElecEta = selLeptons.at(indexElec).Eta();
    TElecPhi = selLeptons.at(indexElec).Phi();
    TElecM   = selLeptons.at(indexElec).M();

    TStatus = -1;
    if      ( selLeptons.at(0).IsPrompt() && selLeptons.at(1).IsPrompt() ) TStatus = 1;
    else if ( selLeptons.at(0).IsConversion() || selLeptons.at(1).IsConversion() ) TStatus = 22;
    else if ( selLeptons.at(0).IsFromB() || selLeptons.at(1).IsFromB() ) TStatus = 5;
    else if ( selLeptons.at(0).IsFromC() || selLeptons.at(1).IsFromC() ) TStatus = 4;
    else if ( selLeptons.at(0).IsFromL() || selLeptons.at(1).IsFromL() ) TStatus = 3;
    else TStatus = 0;

    TLorentzVector muon = TLorentzVector();
    TLorentzVector elec = TLorentzVector();
    TMllMuonUp = TMll; TMllMuonDo = TMll; TMllElecUp = TMll; TMllElecDo = TMll;
    TMuonPtUp = TMuonPt; TMuonPtDo = TMuonPt; TElecPtUp = TElecPt; TElecPtDo = TElecPt;
    if(gDoMuonES){
      TMuonPtUp = selLeptons.at(indexMuon).GetPtUp();
      TMuonPtDo = selLeptons.at(indexMuon).GetPtDo();
      elec.SetPtEtaPhiM(TElecPt,   TElecEta, TElecPhi, TElecM);
      muon.SetPtEtaPhiM(TMuonPtUp, TMuonEta, TMuonPhi, TMuonM);
      TMllMuonUp = (muon+elec).M();
      muon.SetPtEtaPhiM(TMuonPtDo, TMuonEta, TMuonPhi, TMuonM);
      TMllMuonDo = (muon+elec).M();
    }
    if(gDoElecES){
      //GetEnergyUnc
      TElecPtUp = selLeptons.at(indexMuon).Pt();//.PtUp();
      TElecPtDo = selLeptons.at(indexMuon).Pt();//.PtDo();
      muon.SetPtEtaPhiM(TMuonPt,   TMuonEta, TMuonPhi, TMuonM);
      elec.SetPtEtaPhiM(TElecPtUp, TElecEta, TElecPhi, TElecM);
      TMllElecUp = (muon+elec).M();
      elec.SetPtEtaPhiM(TElecPtDo, TElecEta, TElecPhi, TElecM);
      TMllElecDo = (muon+elec).M();
    }
    TIsOnZ             = abs( (selLeptons.at(0).p+selLeptons.at(1).p).M() - 91 ) < 15;
  }
  TPassDilep         = ( (TMuonPt   > 25 && TElecPt   > 20) || (TMuonPt   > 20 && TElecPt   > 25) ) && TMll       > 20;
  if(gIsData){
    TPassDilepAny = TPassDilep;
    TPassMETAny = TMET > metcut;
    TPassMT2Any = TMT2 > mt2cut;
  }
  else{
    TPassDilepMuonESUp = ( (TMuonPtUp > 25 && TElecPt   > 20) || (TMuonPtUp > 20 && TElecPt   > 25) ) && TMllMuonUp > 20;
    TPassDilepMuonESDo = ( (TMuonPtDo > 25 && TElecPt   > 20) || (TMuonPtDo > 20 && TElecPt   > 25) ) && TMllMuonDo > 20; 
    TPassDilepElecESUp = ( (TMuonPt   > 25 && TElecPtUp > 20) || (TMuonPt   > 20 && TElecPtUp > 25) ) && TMllElecUp > 20; 
    TPassDilepElecESDo = ( (TMuonPt   > 25 && TElecPtDo > 20) || (TMuonPt   > 20 && TElecPtDo > 25) ) && TMllElecDo > 20; 
    TPassDilepAny = TPassDilep || TPassDilepMuonESUp || TPassDilepMuonESDo || TPassDilepElecESUp || TPassDilepElecESDo;
  }
}

void TopAnalysis::GetJetVariables(std::vector<Jet> selJets, std::vector<Jet> cleanedJets15, Float_t ptCut){
  TNJets = selJets.size(); THT = 0;
  TNFwdJets   = vetoJets.size() - selJets.size();
  TNBtags = 0; TNBtagsBtagUp = 0; TNBtagsBtagDown = 0;
  TNBtagsMisTagUp = 0;  TNBtagsMisTagDown = 0;
  THTJESUp = 0; THTJESDown = 0; THTJESCorUp = 0; THTJESCorDown = 0; THTJESUnCorUp = 0; THTJESUnCorDown = 0;
  THTJERUp = 0; THTJERDown = 0;
  TBtagPt = 0;

  TNBtags = GetBtags(selJets);
  // BtagPt
  if(TNBtags >= 1){
    for(Int_t i = 0; i < TNJets; i++){
      if(selJets.at(i).isBtag){
        TBtagPt = selJets.at(i).Pt();
        break;
      }
    }
  }
  THT = GetHT(selJets);

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

  if(gIsData){
    TPassJetsAny = (TNJets >= 2);
    TPassBtagAny = (TNBtags >= 1);
    return;  // For systematics...
  }
  TNBtagsBtagUp     = GetBtags(selJets,  1);
  TNBtagsBtagDown   = GetBtags(selJets, -1);
  TNBtagsMisTagUp   = GetBtags(selJets,  3);
  TNBtagsMisTagDown = GetBtags(selJets, -3);

  // NJets
  TNJetsJESUp       = GetNJets(selJetsJecUp);
  TNJetsJESDown     = GetNJets(selJetsJecDown);
  TNJetsJESCorUp    = GetNJets(selJetsJecCorUp);
  TNJetsJESCorDown  = GetNJets(selJetsJecCorDown);
  TNJetsJESUnCorUp  = GetNJets(selJetsJecUnCorUp);
  TNJetsJESUnCorDown= GetNJets(selJetsJecUnCorDown);
  TNJetsJERUp       = GetNJets(selJetsJERUp);
  TNJetsJERDown     = GetNJets(selJetsJERDown);

  // HT
  THTJESUp          = GetHT(selJetsJecUp);
  THTJESDown        = GetHT(selJetsJecDown);
  THTJERUp          = GetHT(selJetsJERUp);
  THTJERDown        = GetHT(selJetsJERDown);
  THTJESCorUp       = GetHT(selJetsJecCorUp);
  THTJESCorDown     = GetHT(selJetsJecCorDown);
  THTJESUnCorUp     = GetHT(selJetsJecUnCorUp);
  THTJESUnCorDown   = GetHT(selJetsJecUnCorDown);

  // Jet0Pt
  TJet0PtJERUp = 0; TJet0PtJERDown = 0;
  TJet0PtJESUp = 0; TJet0PtJESDown = 0;
  TJet0PtJESCorUp = 0; TJet0PtJESCorDown = 0;
  TJet0PtJESUnCorUp = 0; TJet0PtJESUnCorDown = 0;
  if(TNJetsJESUp        > 0) TJet0PtJESUp        = selJetsJecUp.at(0).Pt();
  if(TNJetsJESDown      > 0) TJet0PtJESDown      = selJetsJecDown.at(0).Pt();
  if(TNJetsJERUp        > 0) TJet0PtJERUp        = selJetsJERUp.at(0).Pt();
  if(TNJetsJERDown      > 0) TJet0PtJERDown      = selJetsJERDown.at(0).Pt();
  if(TNJetsJESCorUp     > 0) TJet0PtJESCorUp     = selJetsJecCorUp.at(0).Pt();
  if(TNJetsJESCorDown   > 0) TJet0PtJESCorDown   = selJetsJecCorDown.at(0).Pt();
  if(TNJetsJESUnCorUp   > 0) TJet0PtJESUnCorUp   = selJetsJecUnCorUp.at(0).Pt();
  if(TNJetsJESUnCorDown > 0) TJet0PtJESUnCorDown = selJetsJecUnCorDown.at(0).Pt();

  TPassJetsAny = (TNJets >= 2) || (TNJetsJESCorUp >= 2) || (TNJetsJESCorDown >= 2) || (TNJetsJESUnCorUp >= 2) || (TNJetsJESUnCorDown >= 2) || (TNJetsJERUp >= 2) || (TNJetsJERDown >= 2);
  TPassBtagAny = (TNBtags >= 1) || (TNBtagsBtagUp >= 1) || (TNBtagsBtagDown >= 1) || (TNBtagsMisTagUp >= 1) || (TNBtagsMisTagDown >= 1);
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
  TMT2 = 0;
  TRun        = gIsData ? Get<UInt_t>("run") : 1;
  TMET        = Get<Float_t>(metvarpt); // MET_pt
  TMETPhi    = Get<Float_t>(metvarphi);  // MET phi
  TMETorig    = Get<Float_t>("MET_pt"); TMT2orig = 0;
  TMETsig = Get<Float_t>("MET_significance");
  TLorentzVector pmet; pmet.SetPtEtaPhiM(TMET, 0, TMETPhi, 0);
  if((Int_t) selLeptons.size() >= 2){
    TMT2 = getMT2ll(selLeptons.at(0), selLeptons.at(1), TMET, TMETPhi);
    if(gIs2017) TMT2orig = getMT2ll(selLeptons.at(0), selLeptons.at(1), TMETorig, Get<Float_t>("MET_phi"));
    
    TMT2lblb = getMT2lb(selLeptons.at(0).p, selLeptons.at(1).p, pmet, selJets);
  }
  

  if(gIs2018){
    TMETpuppi        = Get<Float_t>("PuppiMET_pt"); // MET_pt
    TMETpuppi_Phi    = Get<Float_t>("PuppiMET_phi");  // MET phi
    if((Int_t) selLeptons.size() >= 2){
      TMT2puppi=getMT2ll(selLeptons.at(0), selLeptons.at(1), TMETpuppi, TMETpuppi_Phi);
    }
  }

  if(gIsData){
    TPassMETAny = TMET > metcut;
    TPassMT2Any = TMT2 > mt2cut;
  }


  TNVert_pu = 0; TgenTop1Pt = 0; TgenTop2Pt = 0; TGenMET = 0; 
  TMETJESUp = 0; TMETJESDown = 0; TMETJERUp = 0; TMETJERDown = 0; TMETUnclUp = 0; TMETUnclDown = 0; TMETMuonESUp = 0; TMETMuonESDown = 0; TMETElecESUp = 0; TMETElecESDown = 0; TMETPhiJESUp = 0; TMETPhiJESDown = 0; TMETPhiJERUp = 0; TMETPhiJERDown = 0; TMETPhiUnclUp = 0; TMETPhiUnclDown = 0; TMETPhiMuonESUp = 0; TMETPhiMuonESDown = 0; TMETPhiElecESUp = 0; TMETPhiElecESDown = 0; TMT2JESUp = 0; TMT2JESDown = 0; TMT2JERUp = 0; TMT2JERDown = 0; TMT2UnclUp = 0; TMT2UnclDown = 0; TMT2MuonESUp = 0; TMT2MuonESDown = 0; TMT2ElecESUp = 0; TMT2ElecESDown = 0;
TMT2JESCorUp = 0; TMT2JESCorDown = 0; TMT2JESUnCorUp = 0; TMT2JESUnCorDown = 0; TMETJESCorUp = 0; TMETJESCorDown = 0; TMETJESUnCorUp = 0; TMETJESUnCorDown = 0; TMETPhiJESCorUp = 0; TMETPhiJESCorDown = 0; TMETPhiJESUnCorUp = 0; TMETPhiJESUnCorDown = 0;

  if(gIsData) TNVert = Get<Int_t>("PV_npvsGood"); //PV_npvs
  else        TNVert = Get<Int_t>("PV_npvsGood"); //Pileup_nPU
  if(gIsData) return;
  TNVert_pu      = Get<Float_t>("Pileup_nTrueInt"); 
  TGenMET     = Get<Float_t>("GenMET_pt");
  if(gIs2017){
	  TMETUnclUp      = Get<Float_t>("METFixEE2017_pt_unclustEnUp");
	  TMETUnclDown    = Get<Float_t>("METFixEE2017_pt_unclustEnDown");
	  TMETPhiUnclUp   = Get<Float_t>("METFixEE2017_phi_unclustEnUp");
	  TMETPhiUnclDown = Get<Float_t>("METFixEE2017_phi_unclustEnDown");}
  else{
	  TMETUnclUp      = Get<Float_t>("MET_pt_unclustEnUp");
	  TMETUnclDown    = Get<Float_t>("MET_pt_unclustEnDown");
	  TMETPhiUnclUp   = Get<Float_t>("MET_phi_unclustEnUp");
	  TMETPhiUnclDown = Get<Float_t>("MET_phi_unclustEnDown");}
  TMT2UnclUp      = getMT2ll(selLeptons.at(0), selLeptons.at(1), TMETUnclUp,   TMETPhiUnclUp);
  TMT2UnclDown    = getMT2ll(selLeptons.at(0), selLeptons.at(1), TMETUnclDown, TMETPhiUnclDown);
  if(gDoJECunc){
    
    if(gIs2017){
		TMETJESUp      = Get<Float_t>("METFixEE2017_pt_jesTotalUp");
		TMETPhiJESUp   = Get<Float_t>("METFixEE2017_phi_jesTotalUp");
		TMETJESDown    = Get<Float_t>("METFixEE2017_pt_jesTotalDown");
		TMETPhiJESDown = Get<Float_t>("METFixEE2017_phi_jesTotalDown");
		TMETJESCorUp        = Get<Float_t>("METFixEE2017_pt_jesTotalCorrUp");
		TMETPhiJESCorUp     = Get<Float_t>("METFixEE2017_phi_jesTotalCorrUp");
		TMETJESCorDown      = Get<Float_t>("METFixEE2017_pt_jesTotalCorrDown");
		TMETPhiJESCorDown   = Get<Float_t>("METFixEE2017_phi_jesTotalCorrDown");
		TMETJESUnCorUp        = Get<Float_t>("METFixEE2017_pt_jesTotalUnCorrUp");
		TMETPhiJESUnCorUp     = Get<Float_t>("METFixEE2017_phi_jesTotalUnCorrUp");
		TMETJESUnCorDown      = Get<Float_t>("METFixEE2017_pt_jesTotalUnCorrDown");
		TMETPhiJESUnCorDown   = Get<Float_t>("METFixEE2017_phi_jesTotalUnCorrDown");
		TMETJERUp      = Get<Float_t>("METFixEE2017_pt_jerUp");
		TMETPhiJERUp   = Get<Float_t>("METFixEE2017_phi_jerUp");
		TMETJERDown    = Get<Float_t>("METFixEE2017_pt_jerDown");
		TMETPhiJERDown = Get<Float_t>("METFixEE2017_phi_jerDown");}
	else{
		TMETJESUp      = Get<Float_t>("MET_pt_jesTotalUp");
		TMETPhiJESUp   = Get<Float_t>("MET_phi_jesTotalUp");
		TMETJESDown    = Get<Float_t>("MET_pt_jesTotalDown");
		TMETPhiJESDown = Get<Float_t>("MET_phi_jesTotalDown");
		TMETJESCorUp        = Get<Float_t>("MET_pt_jesTotalCorrUp");
		TMETPhiJESCorUp     = Get<Float_t>("MET_phi_jesTotalCorrUp");
		TMETJESCorDown      = Get<Float_t>("MET_pt_jesTotalCorrDown");
		TMETPhiJESCorDown   = Get<Float_t>("MET_phi_jesTotalCorrDown");
		TMETJESUnCorUp        = Get<Float_t>("MET_pt_jesTotalUnCorrUp");
		TMETPhiJESUnCorUp     = Get<Float_t>("MET_phi_jesTotalUnCorrUp");
		TMETJESUnCorDown      = Get<Float_t>("MET_pt_jesTotalUnCorrDown");
		TMETPhiJESUnCorDown   = Get<Float_t>("MET_phi_jesTotalUnCorrDown");
		TMETJERUp      = Get<Float_t>("MET_pt_jerUp");
		TMETPhiJERUp   = Get<Float_t>("MET_phi_jerUp");
		TMETJERDown    = Get<Float_t>("MET_pt_jerDown");
		TMETPhiJERDown = Get<Float_t>("MET_phi_jerDown");}
    TMT2JESCorUp      = getMT2ll(selLeptons.at(0), selLeptons.at(1), TMETJESCorUp,   TMETPhiJESCorUp);
    TMT2JESCorDown    = getMT2ll(selLeptons.at(0), selLeptons.at(1), TMETJESCorDown, TMETPhiJESCorDown);
    TMT2JESUnCorUp      = getMT2ll(selLeptons.at(0), selLeptons.at(1), TMETJESUnCorUp,   TMETPhiJESUnCorUp);
    TMT2JESUnCorDown    = getMT2ll(selLeptons.at(0), selLeptons.at(1), TMETJESUnCorDown, TMETPhiJESUnCorDown);
    TMT2JERUp      = getMT2ll(selLeptons.at(0), selLeptons.at(1), TMETJERUp,   TMETPhiJERUp);
    TMT2JERDown    = getMT2ll(selLeptons.at(0), selLeptons.at(1), TMETJERDown, TMETPhiJERDown);
  }
  if(gDoMuonES){
    TLorentzVector muon; TLorentzVector muonup; TLorentzVector muondo; TLorentzVector elec;
    muon  .SetPtEtaPhiM(TMuonPt,   TMuonEta, TMuonPhi, TMuonM);
    muonup.SetPtEtaPhiM(TMuonPtUp, TMuonEta, TMuonPhi, TMuonM);
    muondo.SetPtEtaPhiM(TMuonPtDo, TMuonEta, TMuonPhi, TMuonM);
    TLorentzVector metmuup = LEStoMET(pmet, muon, muonup);
    TLorentzVector metmudo = LEStoMET(pmet, muon, muondo);
    TMETMuonESDown    = metmudo.Pt();
    TMETPhiMuonESDown = metmudo.Phi();
    TMETMuonESUp      = metmuup.Pt();
    TMETPhiMuonESUp   = metmuup.Phi();
    if((Int_t) selLeptons.size() >= 2){
      elec.SetPtEtaPhiM(TElecPt,   TElecEta, TElecPhi, TElecM);
      TMT2MuonESUp   = getMT2ll(muonup, elec, TMETMuonESUp,   TMETPhiMuonESUp);
      TMT2MuonESDown = getMT2ll(muondo, elec, TMETMuonESDown, TMETPhiMuonESDown);
    }
  }
  if(gDoElecES){
    TLorentzVector elec; TLorentzVector elecup; TLorentzVector elecdo; TLorentzVector muon;
    elec  .SetPtEtaPhiM(TElecPt,   TElecEta, TElecPhi, TElecM);
    elecup.SetPtEtaPhiM(TElecPtUp, TElecEta, TElecPhi, TElecM);
    elecdo.SetPtEtaPhiM(TElecPtDo, TElecEta, TElecPhi, TElecM);
    TLorentzVector metelup = LEStoMET(pmet, elec, elecup);
    TLorentzVector meteldo = LEStoMET(pmet, elec, elecdo);
    TMETElecESDown    = meteldo.Pt();
    TMETPhiElecESDown = meteldo.Pt();
    TMETElecESUp      = metelup.Phi();
    TMETPhiElecESUp   = metelup.Phi();
    if((Int_t) selLeptons.size() >= 2){
      muon  .SetPtEtaPhiM(TMuonPt,   TMuonEta, TMuonPhi, TMuonM);
      TMT2ElecESUp   = getMT2ll(muon, elecup, TMETElecESUp,   TMETPhiElecESUp);
      TMT2ElecESDown = getMT2ll(muon, elecdo, TMETElecESDown, TMETPhiElecESDown);
    }

  }
  if(!gIsData){
    TPassMETAny = (TMET > metcut || TMETJESUp > metcut || TMETJESDown > metcut || TMETJERUp > metcut || TMETJERDown > metcut || TMETMuonESUp > metcut || TMETMuonESDown > metcut || TMETElecESUp > metcut || TMETElecESDown > metcut || TMETUnclUp > metcut || TMETUnclDown > metcut || TMETJESCorUp > metcut || TMETJESCorDown > metcut || TMETJESUnCorUp > metcut || TMETJESUnCorDown > metcut);
    TPassMT2Any = (TMT2 > mt2cut || TMT2JESUp > mt2cut || TMT2JESDown > mt2cut || TMT2JERUp > mt2cut || TMT2JERDown > mt2cut || TMT2MuonESUp > mt2cut || TMT2MuonESDown > mt2cut || TMT2ElecESUp > mt2cut || TMT2ElecESDown > mt2cut || TMT2UnclUp > mt2cut || TMT2UnclDown > mt2cut || TMT2JESCorUp > mt2cut || TMT2JESCorDown > mt2cut || TMT2JESUnCorUp > mt2cut || TMT2JESUnCorDown > mt2cut);
  }
}




void TopAnalysis::GetWeights(){
  TWeight_ElecEffUp = 1; TWeight_ElecEffDown = 1; TWeight_MuonEffUp = 1; TWeight_MuonEffDown = 1;
  TWeight_TrigUp = 1; TWeight_TrigDown = 1; TWeight_PUDown = 1; TWeight_PUUp = 1; TWeight = 1; TWeight_noPU = 1; 
  TWeight_ISRUp   = 1; TWeight_ISRDown = 1; TWeight_FSRUp   = 1; TWeight_FSRDown = 1; TWeight_PrefUp = 1; TWeight_PrefDown = 1; TWeight_TopPtUp = 1; TWeight_TopPtDown = 1;
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

  //if(gIsSignal) PUSF = 1;
  TWeight_noPU             = NormWeight*ElecSF*MuonSF*TrigSF*PrefWeight;

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
  TWeight_TopPtUp   = TWeight;
  TWeight_TopPtDown = TWeight;
  if(gIsTTany) TWeight_TopPtUp    = NormWeight*ElecSF*MuonSF*TrigSF*PUSF*PrefWeight*GetTopPtWeight(TTop0Pt, TTop1Pt);
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
      fHyields[ch][sys]     = CreateH1D("H_Yields_"+suf,"", nLevels, -0.5, nLevels-0.5);
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
        fHLep0Eta[ch][cut][sys]     = CreateH1D("H_Lep0Eta_"    +suffix, "Lep0Eta"   , 50  ,-2.5 ,2.5);
        fHLep1Eta[ch][cut][sys]     = CreateH1F("H_Lep1Eta_"    +suffix, "Lep1Eta"   , 50  ,-2.5 ,2.5);
        fHDelLepPhi[ch][cut][sys]   = CreateH1F("H_DelLepPhi_"  +suffix, "DelLepPhi" , 500, 0, 5);
        fHDelLepEta[ch][cut][sys]   = CreateH1F("H_DelLepEta_"  +suffix, "DelLepEta" ,  96,0, 4.8);
        fHDelLepR[ch][cut][sys]   = CreateH1F("H_DelLepR_"  +suffix, "DelLepR" ,  96,0, 4.8);
        fHHT[ch][cut][sys]          = CreateH1F("H_HT_"         +suffix, "HT"        , 4700,30,500);
        fHJet0Eta[ch][cut][sys]     = CreateH1F("H_Jet0Eta_"  +suffix, "Jet0Eta"   , 50,-2.5,2.5);
        fHJet1Eta[ch][cut][sys]     = CreateH1F("H_Jet1Eta_"  +suffix, "Jet1Eta"   , 50,-2.5,2.5);
        fHDYInvMass[ch][cut][sys]   = CreateH1F("H_DY_InvMass_"    +suffix, "InvMass"   ,  300,  0., 300.);
        fHDYInvMassSF[ch][cut][sys] = CreateH1F("H_DY_SF_InvMass_"    +suffix, "InvMass"   ,  300,  0., 300.);
        fHInvMass[ch][cut][sys]     = CreateH1F("H_InvMass_"    +suffix, "InvMass"   ,  300,  0., 300.);
        fHInvMass2[ch][cut][sys]    = CreateH1F("H_InvMass2_"   +suffix, "InvMass2"  ,  400, 70., 110.);
        fHInvMass2BE[ch][cut][sys]    = CreateH1F("H_InvMass2BE_"   +suffix, "InvMass2BE"  ,  400, 70., 110.);
        fHInvMass2EB[ch][cut][sys]    = CreateH1F("H_InvMass2EB_"   +suffix, "InvMass2EB"  ,  400, 70., 110.);
        fHInvMass2BB[ch][cut][sys]    = CreateH1F("H_InvMass2BB_"   +suffix, "InvMass2BB"  ,  400, 70., 110.);
        fHInvMass2EE[ch][cut][sys]    = CreateH1F("H_InvMass2EE_"   +suffix, "InvMass2EE"  ,  400, 70., 110.);
        fHNBtagsNJets[ch][cut][sys] = CreateH1F("H_NBtagsNJets_"+suffix, "NBtagsNJets"   ,15 , -0.5, 14.5);
        fHNJets[ch][cut][sys]       = CreateH1F("H_NJets_"      +suffix, "NJets"     , 16 ,-0.5, 15.5);
        fHNBtagJets[ch][cut][sys]   = CreateH1F("H_NBtagJets_"  +suffix, "NBtagJets" , 8 ,-0.5, 7.5);
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
        fHvertices_pu[ch][cut][sys]    = CreateH1F("H_VtxPu_"+suffix, "Vtx", 101, -0.5, 100.5);

      }
    }
  }
}

void TopAnalysis::FillDYHistos(Int_t ch){
  if(!makeHistos) return;
  Int_t sys = 0;
  Int_t cut;
  Double_t EventWeight = weight;
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
  fHvertices_pu[ch][cut][sys]      -> Fill(nvert_pu, weight);
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
  fHDelLepR[ch][cut][sys]     -> Fill(deltaR, weight);
  fHInvMass2[ch][cut][sys]    -> Fill(invmass, weight);
  
  //Endcaps and barrel
  if(lep0eta <= 1.479 & lep1eta <= 1.479) fHInvMass2BB[ch][cut][sys]     -> Fill(invmass, weight);
  if(lep0eta <= 1.479 & lep1eta > 1.479) fHInvMass2BE[ch][cut][sys]     -> Fill(invmass, weight);
  if(lep0eta > 1.479 & lep1eta <= 1.479) fHInvMass2EB[ch][cut][sys]     -> Fill(invmass, weight);
  if(lep0eta > 1.479 & lep1eta > 1.479) fHInvMass2EE[ch][cut][sys]     -> Fill(invmass, weight);

  // Jets
  if(njets > 0){ 
    fHJet0Eta[ch][cut][sys]     -> Fill(selJets.at(0).Eta(), weight);
    fHJet0Pt [ch][cut][sys]     -> Fill(selJets.at(0).Pt(), weight);
    fHJet0CSV[ch][cut][sys]     -> Fill(selJets.at(0).csv, weight);
    fHJet0DeepCSV[ch][cut][sys] -> Fill(selJets.at(0).GetDeepCSVB(), weight);
    fHJet0DeepFlav[ch][cut][sys] -> Fill(selJets.at(0).GetDeepFlav(), weight);
  }
  if(njets > 1){
    fHJet1Eta[ch][cut][sys]     -> Fill(selJets.at(1).Eta(), weight);
    fHJet1Pt [ch][cut][sys]     -> Fill(selJets.at(1).Pt(), weight);
    fHJet1CSV[ch][cut][sys]     -> Fill(selJets.at(1).csv, weight);
    fHJet1DeepCSV[ch][cut][sys] -> Fill(selJets.at(1).GetDeepCSVB(), weight);
    fHJet1DeepFlav[ch][cut][sys] -> Fill(selJets.at(1).GetDeepFlav(), weight);
  }
  for(Int_t i = 0; i < (Int_t) selJets.size(); i++){ 
    fHJetPt[ch][cut][sys]  ->Fill(selJets.at(i).Pt());
    fHJetEta[ch][cut][sys] ->Fill(selJets.at(i).Eta());
    fHJetCSV[ch][cut][sys] -> Fill(selJets.at(i).csv, weight);
    fHJetDeepCSV[ch][cut][sys] -> Fill(selJets.at(i).GetDeepCSVB(), weight);
    fHJetDeepFlav[ch][cut][sys] -> Fill(selJets.at(i).GetDeepFlav(), weight);
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
  cout << "PASAAAAAAA FillCorr" << endl;
  Float_t pTreco; Float_t pTgen; Float_t dR = 0.3; Float_t dPtoPt; Bool_t isBtag; Bool_t isBJet;
  Float_t bin0; Float_t bin1;
  Int_t njets = selJets.size();
  for(Int_t i = 0; i < njets; i++){
    pTreco = selJets.at(i).Pt(); pTgen = selJets.at(i).GetGenPt();
    dPtoPt = fabs(pTreco - pTgen)/pTreco;
    isBtag = selJets.at(i).isBtag; isBJet = selJets.at(i).flavmc == 5;

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
    fTree->Branch("TLep0M",       &TLep0M,       "TLep0M/F");
    fTree->Branch("TLep1M",       &TLep1M,       "TLep1M/F");
  }
  fTree->Branch("TIsSS",        &TIsSS,        "TIsSS/I");
  fTree->Branch("TStatus",      &TStatus,      "TStatus/I");
  fTree->Branch("TLep0Id",      &TLep0Id,      "TLep0Id/I");
  fTree->Branch("TLep1Id",      &TLep1Id,      "TLep1Id/I");
  fTree->Branch("TMll",         &TMll,         "TMll/F");
  fTree->Branch("TDilep_Pt",    &TDilep_Pt,    "TDilep_Pt/F");
  fTree->Branch("TLep0Pt",      &TLep0Pt,      "TLep0Pt/F");
  fTree->Branch("TLep0Eta",     &TLep0Eta,     "TLep0Eta/F");
  fTree->Branch("TLep0Phi",     &TLep0Phi,     "TLep0Phi/F");
  fTree->Branch("TDeltaPhi",    &TDeltaPhi,    "TDeltaPhi/F");
  fTree->Branch("TDeltaEta",    &TDeltaEta,    "TDeltaEta/F");
  fTree->Branch("TLep1Pt",      &TLep1Pt,      "TLep1Pt/F");
  fTree->Branch("TLep1Eta",     &TLep1Eta,     "TLep1Eta/F");
  fTree->Branch("TLep1Phi",     &TLep1Phi,     "TLep1Phi/F");

  fTree->Branch("TGenLep0Pt",      &TGenLep0Pt,      "TGenLep0Pt/F");
  fTree->Branch("TGenLep0Eta",     &TGenLep0Eta,     "TGenLep0Eta/F");
  fTree->Branch("TGenLep0Phi",     &TGenLep0Phi,     "TGenLep0Phi/F");
  fTree->Branch("TGenLep1Pt",      &TGenLep1Pt,      "TGenLep1Pt/F");
  fTree->Branch("TGenLep1Eta",     &TGenLep1Eta,     "TGenLep1Eta/F");
  fTree->Branch("TGenLep1Phi",     &TGenLep1Phi,     "TGenLep1Phi/F");
  fTree->Branch("TPassDilep",     &TPassDilep,     "TPassDilep/I");
  if(!gIsData){
    fTree->Branch("TPassDilepMuonESUp",     &TPassDilepMuonESUp,     "TPassDilepMuonESUp/I");
    fTree->Branch("TPassDilepMuonESDo",     &TPassDilepMuonESDo,     "TPassDilepMuonESDown/I");
    fTree->Branch("TPassDilepElecESUp",     &TPassDilepElecESUp,     "TPassDilepElecESUp/I");
    fTree->Branch("TPassDilepElecESDo",     &TPassDilepElecESDo,     "TPassDilepElecESDown/I");
  }
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
  }

  if(!gIsData){
    // JER / JES uncertainties
    fTree->Branch("TNJetsJESUp",       &TNJetsJESUp,     "TNJetsJESUp/I");
    fTree->Branch("TNJetsJESDown",     &TNJetsJESDown,   "TNJetsJESDown/I");
    fTree->Branch("TNJetsJESCorUp",       &TNJetsJESCorUp,     "TNJetsJESCorUp/I");
    fTree->Branch("TNJetsJESCorDown",     &TNJetsJESCorDown,   "TNJetsJESCorDown/I");
    fTree->Branch("TNJetsJESUnCorUp",       &TNJetsJESUnCorUp,     "TNJetsJESUnCorUp/I");
    fTree->Branch("TNJetsJESUnCorDown",     &TNJetsJESUnCorDown,   "TNJetsJESUnCorDown/I");
    fTree->Branch("TNJetsJERUp",       &TNJetsJERUp,     "TNJetsJERUp/I");
    fTree->Branch("TNJetsJERDown",     &TNJetsJERDown,   "TNJetsJERDown/I");
    fTree->Branch("THTJESUp",     &THTJESUp,     "THTJESUp/F");
    fTree->Branch("THTJESDown",   &THTJESDown,   "THTJESDown/F");
    fTree->Branch("THTJESCorUp",     &THTJESCorUp,     "THTJESCorUp/F");
    fTree->Branch("THTJESCorDown",   &THTJESCorDown,   "THTJESCorDown/F");
    fTree->Branch("THTJESUnCorUp",     &THTJESUnCorUp,     "THTJESUnCorUp/F");
    fTree->Branch("THTJESUnCorDown",   &THTJESUnCorDown,   "THTJESUnCorDown/F");
    fTree->Branch("THTJERUp",     &THTJERUp,     "THTJERUp/F");
    fTree->Branch("THTJERDown",   &THTJERDown,   "THTJERDown/F");
    fTree->Branch("TJet0PtJESCorUp",       &TJet0PtJESCorUp,     "TJet0PtJESCorUp/F");
    fTree->Branch("TJet0PtJESCorDown",     &TJet0PtJESCorDown,   "TJet0PtJESCorDown/F");
    fTree->Branch("TJet0PtJESUnCorUp",       &TJet0PtJESUnCorUp,     "TJet0PtJESUnCorUp/F");
    fTree->Branch("TJet0PtJESUnCorDown",     &TJet0PtJESUnCorDown,   "TJet0PtJESUnCorDown/F");

    // B-tagging uncertainties
    fTree->Branch("TNBtagsBtagUp",     &TNBtagsBtagUp,   "TNBtagsBtagUp/I");
    fTree->Branch("TNBtagsBtagDown",   &TNBtagsBtagDown, "TNBtagsBtagDown/I");
    fTree->Branch("TNBtagsMisTagUp",     &TNBtagsMisTagUp,   "TNBtagsMisTagUp/I");
    fTree->Branch("TNBtagsMisTagDown",   &TNBtagsMisTagDown, "TNBtagsMisTagDown/I");
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
  fTree->Branch("TWeight",      &TWeight,      "TWeight/D");
  fTree->Branch("TWeight_noPU",      &TWeight_noPU,      "TWeight_noPU/D");

  fTree->Branch("TMET",            &TMET,            "TMET/F");
  fTree->Branch("TMETPhi",     &TMETPhi,     "TMETPhi/F");
  fTree->Branch("TMT2",            &TMT2,            "TMT2/F");
  fTree->Branch("TMT2lblb",       &TMT2lblb,          "TMT2lblb/F");
  fTree->Branch("TMETsig", &TMETsig, "TMETsig/F");
  if(!gIsData){
    fTree->Branch("TWeight_ElecEffUp",      &TWeight_ElecEffUp,     "TWeight_ElecEffUp/F");
    fTree->Branch("TWeight_ElecEffDown",    &TWeight_ElecEffDown,   "TWeight_ElecEffDown/F");
    fTree->Branch("TWeight_MuonEffUp",      &TWeight_MuonEffUp,     "TWeight_MuonEffUp/F");
    fTree->Branch("TWeight_MuonEffDown",    &TWeight_MuonEffDown,   "TWeight_MuonEffDown/F");
    fTree->Branch("TWeight_TrigUp",         &TWeight_TrigUp,        "TWeight_TrigUp/F");
    fTree->Branch("TWeight_TrigDown",       &TWeight_TrigDown,      "TWeight_TrigDown/F");
    fTree->Branch("TWeight_PUUp",           &TWeight_PUUp,          "TWeight_PUUp/F");
    fTree->Branch("TWeight_PUDown",         &TWeight_PUDown,        "TWeight_PUDown/F");
    fTree->Branch("TWeight_ISRUp",          &TWeight_ISRUp,          "TWeight_ISRUp/F");
    fTree->Branch("TWeight_ISRDown",        &TWeight_ISRDown,        "TWeight_ISRDown/F");
    fTree->Branch("TWeight_FSRUp",           &TWeight_FSRUp,          "TWeight_FSRUp/F");
    fTree->Branch("TWeight_FSRDown",         &TWeight_FSRDown,        "TWeight_FSRDown/F");
    fTree->Branch("TWeight_PrefUp",           &TWeight_PrefUp,          "TWeight_PrefUp/F");
    fTree->Branch("TWeight_PrefDown",         &TWeight_PrefDown,        "TWeight_PrefDown/F");
    if(gIsTTany) fTree->Branch("TWeight_TopPtUp",         &TWeight_TopPtUp,        "TWeight_TopPtUp/F");
    if(gIsTTany) fTree->Branch("TWeight_TopPtDown",         &TWeight_TopPtDown,        "TWeight_TopPtDown/F");
    fTree->Branch("TNVert_pu",       &TNVert_pu,          "TNVert_pu/F");

    fTree->Branch("TMETJESUp",      &TMETJESUp,       "TMETJESUp/F");
    fTree->Branch("TMETJESDown",    &TMETJESDown,     "TMETJESDown/F");
    fTree->Branch("TMETJERUp",      &TMETJERUp,       "TMETJERUp/F");
    fTree->Branch("TMETJERDown",    &TMETJERDown,     "TMETJERDown/F");
    fTree->Branch("TMETJESCorUp",      &TMETJESCorUp,       "TMETJESCorUp/F");
    fTree->Branch("TMETJESCorDown",    &TMETJESCorDown,     "TMETJESCorDown/F");
    fTree->Branch("TMETJESUnCorUp",      &TMETJESUnCorUp,       "TMETJESUnCorUp/F");
    fTree->Branch("TMETJESUnCorDown",    &TMETJESUnCorDown,     "TMETJESUnCorDown/F");

    fTree->Branch("TMT2JESUp",      &TMT2JESUp,       "TMT2JESUp/F");
    fTree->Branch("TMT2JESDown",    &TMT2JESDown,     "TMT2JESDown/F");
    fTree->Branch("TMT2JERUp",      &TMT2JERUp,       "TMT2JERUp/F");
    fTree->Branch("TMT2JERDown",    &TMT2JERDown,     "TMT2JERDown/F");
    fTree->Branch("TMT2JESCorUp",      &TMT2JESCorUp,       "TMT2JESCorUp/F");
    fTree->Branch("TMT2JESCorDown",    &TMT2JESCorDown,     "TMT2JESCorDown/F");
    fTree->Branch("TMT2JESUnCorUp",      &TMT2JESUnCorUp,       "TMT2JESUnCorUp/F");
    fTree->Branch("TMT2JESUnCorDown",    &TMT2JESUnCorDown,     "TMT2JESUnCorDown/F");

    if(gDoMuonES){
      fTree->Branch("TMETMuonESUp",   &TMETMuonESUp,    "TMETMuonESUp/F");
      fTree->Branch("TMETMuonESDown", &TMETMuonESDown,  "TMETMuonESDown/F");
      fTree->Branch("TMT2MuonESUp",   &TMT2MuonESUp,    "TMT2MuonESUp/F");
      fTree->Branch("TMT2MuonESDown", &TMT2MuonESDown,  "TMT2MuonESDown/F");
    }
    if(gDoElecES){
      fTree->Branch("TMETElecESUp",   &TMETElecESUp,    "TMETElecESUp/F");
      fTree->Branch("TMETElecESDown", &TMETElecESDown,  "TMETElecESDown/F");
      fTree->Branch("TMT2ElecESUp",   &TMT2ElecESUp,    "TMT2ElecESUp/F");
      fTree->Branch("TMT2ElecESDown", &TMT2ElecESDown,  "TMT2ElecESDown/F");
    }
    fTree->Branch("TMETUnclUp",   &TMETUnclUp,    "TMETUnclUp/F");
    fTree->Branch("TMETUnclDown", &TMETUnclDown,  "TMETUnclDown/F");
    fTree->Branch("TMT2UnclUp",   &TMT2UnclUp,    "TMT2UnclUp/F");
    fTree->Branch("TMT2UnclDown", &TMT2UnclDown,  "TMT2UnclDown/F");
  }

  if(gIs2018) fTree->Branch("TIsHEM",          &TIsHEM,    "TIsHEM/B");
  fTree->Branch("TNVert",          &TNVert,          "TNVert/I");
  if(gIs2017){
    fTree->Branch("TMETorig",            &TMETorig,            "TMETorig/F");
    fTree->Branch("TMT2orig",            &TMT2orig,            "TMT2orig/F");
  }
  if(gIs2018){
    fTree->Branch("TMETpuppi",            &TMETpuppi,            "TMETpuppi/F");
    fTree->Branch("TMT2puppi",            &TMT2puppi,            "TMT2puppi/F");
    fTree->Branch("TMETpuppi_Phi",            &TMETpuppi_Phi,            "TMETpuppi_Phi/F");

  }
  if(!gIsData){
    fTree->Branch("TGenMET",         &TGenMET,         "TGenMET/F");
    fTree->Branch("TGenMET_phi",     &TGenMET_phi,      "TGenMET_phi/F");
    fTree->Branch("TGenMT2",         &TGenMT2,         "TGenMT2/F");
  }
  if(!miniTree){
    fTree->Branch("TEvent",          &event,           "TEvent/l");
    fTree->Branch("TLuminosityBlock",&lumiblock,       "TLuminosityBlock/i");
    fTree->Branch("TPassMETFilters", &TPassMETFilters, "TPassMETFilters/B");
    fTree->Branch("TPassTrigger",    &TPassTrigger,    "TPassTrigger/B");
    fTree->Branch("TRun",            &TRun,            "TRun/i");
    //fTree->Branch("TgenTop1Pt",   &TgenTop1Pt,   "TgenTop1Pt/F");
    //fTree->Branch("TgenTop2Pt",   &TgenTop2Pt,   "TgenTop2Pt/F");
    fTree->Branch("TMETJESUp",    &TMETJESUp,    "TMETJESUp/F");
    fTree->Branch("TMETJESDown",  &TMETJESDown,  "TMETJESDown/F");
  }

  if (gIsSignal || gIsTTany){
    fTree->Branch("TTop0Pt",  &TTop0Pt,  "TTop0Pt/F");
    fTree->Branch("TTop0Eta", &TTop0Eta, "TTop0Eta/F");
    fTree->Branch("TTop0Phi", &TTop0Phi, "TTop0Phi/F");
    fTree->Branch("TTop1Pt",  &TTop1Pt,  "TTop1Pt/F");
    fTree->Branch("TTop1Eta", &TTop1Eta, "TTop1Eta/F");
    fTree->Branch("TTop1Phi", &TTop1Phi, "TTop1Phi/F");
    if(gIsSignal){
      fTree->Branch("Tm_stop", &m_stop, "Tm_stop/F");
      fTree->Branch("Tm_LSP", &m_LSP, "Tm_LSP/F");  
    }
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
  met = TMET; mt2 = TMT2; nvert = TNVert; invmass = TMll; nvert_pu = TNVert_pu;
  // Leptons
  if(nleps >= 2){
    lep0pt = TLep0Pt; lep1pt = TLep1Pt; lep0eta = TLep0Eta; lep1eta = TLep1Eta; 
    lep0iso = selLeptons.at(0).GetIso(); lep1iso = selLeptons.at(1).GetIso();
    dileppt = (selLeptons.at(0).p + selLeptons.at(1).p).Pt();
    deltaphi = selLeptons.at(0).p.DeltaPhi(selLeptons.at(1).p);
    deltaeta = selLeptons.at(0).p.Eta() - selLeptons.at(1).p.Eta();
    deltaR = selLeptons.at(0).p.DeltaR(selLeptons.at(1).p);
  }

  // Jets
  njets = 0; nbtags = 0; ht = 0;
  jets.clear(); Jet jet; TLorentzVector t; 
  float pt, eta, phi, m; bool isbtag; float csv, deepcsv, deepflav;
  Int_t nJets = Get<Int_t>("nJet");
  Int_t jetid, flav;
  for(int i = 0; i < nJets; i++){
    pt = Get<Float_t>(JetPt,i); eta = Get<Float_t>("Jet_eta",i); phi = Get<Float_t>("Jet_phi", i); m = Get<Float_t>("Jet_mass",i);
    csv = Get<Float_t>("Jet_btagCSVV2", i); deepcsv = Get<Float_t>("Jet_btagDeepB", i); deepflav = Get<Float_t>("Jet_btagDeepFlavB", i);
    jetid = Get<Int_t>("Jet_jetId",i);

    Float_t genPt = -999;
    if(!gIsData){
      Int_t nGenJets = Get<Int_t>("nGenJet");
      Int_t genJetIdx = Get<Int_t>("Jet_genJetIdx", i);
      if(genJetIdx > 0 and genJetIdx < nJets) genPt = Get<Float_t>("GenJet_pt", genJetIdx);    
    }
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
    Float_t alg = deepflav; //quitar: dejar deepflav
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
        ht += jet.Pt();
      }
    }
  }
 //TNJets = njets; TNBtags = nbtags; THT = ht;
  njets= TNJets ; nbtags=TNBtags ; ht=THT;
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
  else if(sys == kTopPt      ) weight = TWeight_TopPtUp;
}

Float_t TopAnalysis::GetTopPtWeight(Float_t Pt1, Float_t Pt2){
  Float_t normWeight = 1.00298;
  Float_t a = 0.0615; Float_t b = 0.0005;
  Float_t SF1; Float_t SF2;
  SF1 = TMath::Exp(a - Pt1*b);
  SF2 = TMath::Exp(a - Pt2*b);
  return TMath::Sqrt(SF1*SF2)*normWeight;
}

