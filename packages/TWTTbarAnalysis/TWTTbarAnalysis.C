#include "TWTTbarAnalysis.h"
ClassImp(TWTTbarAnalysis);

bool GreaterThan(Float_t i, Float_t j){ return (i > j);}

//#####################################################################
// Core PAF methods
//---------------------------------------------------------------------

TWTTbarAnalysis::TWTTbarAnalysis() : PAFChainItemSelector() {
  fhDummy        = 0;
  passMETfilters = 0;
  passTrigger    = 0;
  isSS           = 0;
  
  for(Int_t i = 0; i < 254; i++) TLHEWeight[i] = 0;
}



void TWTTbarAnalysis::Initialise() {
  gIsData     = GetParam<Bool_t>("IsData");
  gSelection  = GetParam<Int_t>("iSelection");
  gSampleName = GetParam<TString>("sampleName");
  gDoSyst     = GetParam<Bool_t>("doSyst");
  gPar        = GetParam<TString>("_options");
  if (gPar.Contains("Semi")) {
    cout << "> Running the semileptonic ttbar sample" << endl;
  }
  gIsTTbar    = false;
  gIsLHE      = false;

  if (gSampleName.Contains("TTbar") || gSampleName.Contains("TTJets")) gIsTTbar = true;
  if (gSampleName == "TTbar_Powheg")   gIsLHE = true;
  
  if (gSampleName.Contains("_")) {
    TObjArray *tx = gSampleName.Tokenize("_");
    if (((TObjString*) tx->Last())->GetString().IsDigit()) gIsLHE = true;
  }
  
  fhDummy = CreateH1F("fhDummy", "fhDummy", 1, 0, 2);
  
  fMiniTree = CreateTree("fMiniTree", "MiniTree");
  SetTWTTbarVariables();
  
  genLeptons      = std::vector<Lepton>();
  selLeptons      = std::vector<Lepton>();
  vetoLeptons     = std::vector<Lepton>();
  selJets         = std::vector<Jet>();
  selJetsJecUp    = std::vector<Jet>();
  selJetsJecDown  = std::vector<Jet>();
  selJetsJER      = std::vector<Jet>();
  Jets15          = std::vector<Jet>();
  genJets         = std::vector<Jet>();
  mcJets          = std::vector<Jet>();
  vetoJets        = std::vector<Jet>();
  
  GenJets         = std::vector<Jet>();
  GenLooseCentralJets = std::vector<Jet>();
  GeNLooseFwdJets     = std::vector<Jet>();
  GenLeps             = std::vector<Lepton>();
}



void TWTTbarAnalysis::InsideLoop() {
  ResetTWTTbarVariables();
  // Vectors with the objects
  genLeptons        = GetParam<vector<Lepton>>("genLeptons");
  selLeptons        = GetParam<vector<Lepton>>("selLeptons");
//   vetoLeptons       = GetParam<vector<Lepton>>("vetoLeptons");
  selJets           = GetParam<vector<Jet>>("selJets");
  selJetsJecUp      = GetParam<vector<Jet>>("selJetsJecUp");
  selJetsJecDown    = GetParam<vector<Jet>>("selJetsJecDown");
  selJetsJER        = GetParam<vector<Jet>>("selJetsJER");
//   Jets15            = GetParam<vector<Jet>>("Jets15");
  vetoJets          = GetParam<vector<Jet>>("vetoJets");
  genJets           = GetParam<vector<Jet>>("genJets");
  mcJets            = GetParam<vector<Jet>>("mcJets");
  
  // Weights and SFs
  NormWeight        = GetParam<Float_t>("NormWeight");
  TrigSF            = GetParam<Float_t>("TriggerSF");
  TrigSFerr         = GetParam<Float_t>("TriggerSFerr");
  PUSF              = GetParam<Float_t>("PUSF");
  PUSF_Up           = GetParam<Float_t>("PUSF_Up");
  PUSF_Down         = GetParam<Float_t>("PUSF_Down");
  BtagSF            = GetParam<Float_t>("BtagSF");
  BtagSFBtagUp      = GetParam<Float_t>("BtagSFBtagUp");
  BtagSFBtagDown    = GetParam<Float_t>("BtagSFBtagDown");
  BtagSFMistagUp    = GetParam<Float_t>("BtagSFMistagUp");
  BtagSFMistagDown  = GetParam<Float_t>("BtagSFMistagDown");

  // Event variables
  gChannel          = GetParam<Int_t>("gChannel");
  passMETfilters    = GetParam<Bool_t>("METfilters");
  passTrigger       = GetParam<Bool_t>("passTrigger");
  isSS              = GetParam<Bool_t>("isSS");
  year              = GetParam<Int_t>("Year");
  
  // Leptons and Jets
  GetLeptonVariables();
  GetGenLepVariables();
  GetJetVariables();
  GetGenJetVariables();
  GetMETandGenMET();
  
  
  fhDummy->Fill(1);

  if (gPar.Contains("Semi")) {
    if (gIsTTbar && genLeptons.size() > 1 ) return;
  } else {
    if (gIsTTbar && genLeptons.size() < 2 ) return; // Dilepton selection for ttbar!!!
  }
  
  TWeight_normal = NormWeight;
  
  if ((gPar.Contains("Unfolding")) && (nGenLeps >= 2) && (nGenJets == 2) && (nGenbJets == 2)) { // Checking if we pass the selection with gen things
    if(GenLeps.at(0).isElec && GenLeps.at(1).isMuon) GenChannel = iElMu; // ...but first, let's redefine the GenChannel to get it right
    if(GenLeps.at(0).isMuon && GenLeps.at(1).isElec) GenChannel = iElMu;
    if(GenLeps.at(0).isMuon && GenLeps.at(1).isMuon) GenChannel = iMuon;
    if(GenLeps.at(0).isElec && GenLeps.at(1).isElec) GenChannel = iElec;
    TGenIsSS            = (GenLeps.at(0).charge * GenLeps.at(1).charge) > 0;
    
    TDressLep1Jet1_M     = (GenJets.at(0).p + GenLeps.at(0).p).M();
    TDressLep2Jet1_M   = (GenJets.at(0).p + GenLeps.at(1).p).M();
    TDressLep1Jet2_M    = (GenJets.at(1).p + GenLeps.at(0).p).M();
    TDressLep2Jet2_M = (GenJets.at(1).p + GenLeps.at(1).p).M();
    if(GreaterThan(TDressLep1Jet1_M,TDressLep2Jet2_M)) {
        if(GreaterThan(TDressLep2Jet1_M,TDressLep1Jet2_M)) {
            if(!GreaterThan(TDressLep1Jet1_M,TDressLep2Jet1_M)) {TDressMiniMax = TDressLep1Jet1_M;}
            else                                               {TDressMiniMax = TDressLep2Jet1_M;}
        }
        else {
            if(!GreaterThan(TDressLep1Jet1_M,TDressLep1Jet2_M))  {TDressMiniMax = TDressLep1Jet1_M; }
            else                                               {TDressMiniMax = TDressLep1Jet2_M;}
        }
    }
    else {
        if(GreaterThan(TDressLep2Jet1_M,TDressLep1Jet2_M)) {
            if(!GreaterThan(TDressLep2Jet2_M,TDressLep2Jet1_M)) {TDressMiniMax = TDressLep2Jet2_M;}
            else                                                    {TDressMiniMax = TDressLep2Jet1_M;}
        }
        else {
            if(!GreaterThan(TDressLep2Jet2_M,TDressLep1Jet2_M))  {TDressMiniMax = TDressLep2Jet2_M;}
            else                                                    {TDressMiniMax = TDressLep1Jet2_M;}
        }  
    }
    TDressLep1Lep2Jet1_E           = (GenJets.at(0).p + GenLeps.at(0).p + GenLeps.at(1).p).E();
    TGenLep1Lep2METJet1_Mt       = (GenJets.at(0).p + GenLeps.at(0).p + GenLeps.at(1).p + GenMET).Mt();
    TGenLep1Lep2Jet1_M           = (GenJets.at(0).p + GenLeps.at(0).p + GenLeps.at(1).p).M();
    TGenDilepPt         = (GenLeps.at(0).p + GenLeps.at(1).p).Pt();
    TGenLep1Lep2Jet1_Pt      = (GenLeps.at(0).p + GenLeps.at(1).p + GenJets.at(0).p).Pt();
    TGenLep1Lep2METJet1_Pt   = (GenLeps.at(0).p + GenLeps.at(1).p + GenJets.at(0).p + GenMET).Pt();
    TGenHTtot           = (GenLeps.at(0).Pt() + GenLeps.at(1).Pt() + GenJets.at(0).Pt() + GenMET.Pt());
    TGenLep1Lep2METJet1_Pz  = (GenLeps.at(0).p + GenLeps.at(1).p + GenJets.at(0).p + GenMET).Pz();
    TGenLep1Lep2METJet1_Eta       = (GenLeps.at(0).p + GenLeps.at(1).p + GenJets.at(0).p + GenMET).Eta();
    TGenMSys            = (GenLeps.at(0).p + GenLeps.at(1).p + GenJets.at(0).p + GenMET).M();
    TGenLep1Lep2_M             = (GenLeps.at(0).p + GenLeps.at(1).p).M();
    TGenDPhiLL          = GenLeps.at(0).p.DeltaPhi(GenLeps.at(1).p);
    TGenDPhiLeadJet     = GenLeps.at(0).p.DeltaPhi(GenJets.at(0).p);
    TGenDPhiSubLeadJet  = GenLeps.at(1).p.DeltaPhi(GenJets.at(0).p);
    
    if (((GenLeps.at(0).p + GenLeps.at(1).p).M() > 20) && (GenLeps.at(0).p.Pt() > 25) && 
        (GenChannel == iElMu || GenChannel == iMuon || GenChannel == iElec) && !TGenIsSS && (nGenLooseCentralJets == 2)) {
      if (GenChannel == iMuon || GenChannel == iElec) {
        if ((TGenMET > 20) && (abs(TGenLep1Lep2_M-90.19)>15)) {TPassGen = 1;}   
      }
      else {TPassGen = 1;}
    }
  }
  
  
  
  if ((TNSelLeps >= 2) && passTrigger && passMETfilters && ((selLeptons.at(0).p + selLeptons.at(1).p).M() > 20) &&
      (selLeptons.at(0).p.Pt() > 25) && (TChannel == iElMu || TChannel == iElec || TChannel == iMuon) && (!isSS)) {
    Float_t lepSF     = selLeptons.at(0).GetSF( 0) * selLeptons.at(1).GetSF( 0);
    Float_t ElecSF    = 1; Float_t MuonSF   = 1;
    Float_t ElecSFUp  = 1; Float_t ElecSFDo = 1; Float_t MuonSFUp = 1; Float_t MuonSFDo = 1;
    if (TChannel == iElec) {
      ElecSF   = selLeptons.at(0).GetSF( 0)*selLeptons.at(1).GetSF( 0);
      ElecSFUp = selLeptons.at(0).GetSF( 1)*selLeptons.at(1).GetSF( 1);
      ElecSFDo = selLeptons.at(0).GetSF(-1)*selLeptons.at(1).GetSF(-1);
      MuonSFUp = 1; MuonSFDo = 1; MuonSF = 1;
    }
    else if (TChannel == iMuon) {
      MuonSFUp = selLeptons.at(0).GetSF( 1)*selLeptons.at(1).GetSF( 1);
      MuonSFDo = selLeptons.at(0).GetSF(-1)*selLeptons.at(1).GetSF(-1);
      ElecSFUp = 1; ElecSFDo = 1; ElecSF = 1;
    }
    else {
      if(selLeptons.at(0).isMuon) {
        MuonSF   *= selLeptons.at(0).GetSF( 0);
        MuonSFUp *= selLeptons.at(0).GetSF( 1);
        MuonSFDo *= selLeptons.at(0).GetSF(-1);
      }
      else {
        ElecSF   *= selLeptons.at(0).GetSF( 0);
        ElecSFUp *= selLeptons.at(0).GetSF( 1);
        ElecSFDo *= selLeptons.at(0).GetSF(-1);
      }
      if (selLeptons.at(1).isMuon) {
        MuonSF   *= selLeptons.at(1).GetSF( 0);
        MuonSFUp *= selLeptons.at(1).GetSF( 1);
        MuonSFDo *= selLeptons.at(1).GetSF(-1);
      }
      else {
        ElecSF   *= selLeptons.at(1).GetSF( 0);
        ElecSFUp *= selLeptons.at(1).GetSF( 1);
        ElecSFDo *= selLeptons.at(1).GetSF(-1);
      }
    }
    
    TWeight               = NormWeight * ElecSF   * MuonSF * TrigSF * PUSF * BtagSF;
    TWeight_ElecEffUp     = NormWeight * ElecSFUp * MuonSF * TrigSF * PUSF * BtagSF;
    TWeight_ElecEffDown   = NormWeight * ElecSFDo * MuonSF * TrigSF * PUSF * BtagSF;
    TWeight_MuonEffUp     = NormWeight * ElecSF   * (MuonSF+TMath::Sqrt(TMath::Power(MuonSFUp-MuonSF,2)+TMath::Power(MuonSF*0.0122,2)))*TrigSF*PUSF*BtagSF;
    TWeight_MuonEffDown   = NormWeight * ElecSF   * (MuonSF-TMath::Sqrt(TMath::Power(MuonSFDo-MuonSF,2)+TMath::Power(MuonSF*0.0122,2)))*TrigSF*PUSF*BtagSF;
    TWeight_TrigUp        = NormWeight * lepSF    * (TrigSF+TrigSFerr) * PUSF   * BtagSF;
    TWeight_TrigDown      = NormWeight * lepSF    * (TrigSF-TrigSFerr) * PUSF   * BtagSF;
    TWeight_PUDown        = NormWeight * lepSF    * TrigSF * PUSF_Up   * BtagSF;
    TWeight_PUUp          = NormWeight * lepSF    * TrigSF * PUSF_Down * BtagSF;
    TWeight_BtagUp        = NormWeight * ElecSF   * MuonSF * TrigSF * PUSF * BtagSFBtagUp    ;
    TWeight_BtagDown      = NormWeight * ElecSF   * MuonSF * TrigSF * PUSF * BtagSFBtagDown  ;
    TWeight_MistagUp      = NormWeight * ElecSF   * MuonSF * TrigSF * PUSF * BtagSFMistagUp  ;
    TWeight_MistagDown    = NormWeight * ElecSF   * MuonSF * TrigSF * PUSF * BtagSFMistagDown;

    if(gIsData) TWeight = 1;
    
    CalculateTWTTbarVariables();
    
    if ((TNJets == 2) && (TNBJets == 2) && (nLooseCentral == 2)) {
      if (TChannel == iMuon || TChannel == iElec) {
          if ((TMET > 20) && (abs(TLep1Lep2_M-90.19)>15)) {TPassReco = 1;}
      }
      else {TPassReco = 1;}
    }
    if ((TNJetsJESUp == 1) && (TNBJetsJESUp == 1) && (NLooseCentralJESUp == 1)) {
      TPassRecoJESUp = 1;
    }
    if ((TNJetsJESDown == 1) && (TNBJetsJESDown == 1) && (NLooseCentralJESDown == 1)) {
      TPassRecoJESDown = 1;
    }
    if ((TNJetsJERUp == 1) && (TNBJetsJERUp == 1) && (NLooseCentralJERUp == 1)) {
      TPassRecoJERUp = 1;
    }
  }
  
  //   Setting maximum value of the unfolding candidate variables.
  if (TMET                 >= 200)         TMET                 = 199.999;
  if (TLep1Jet1_M          >= 400)         TLep1Jet1_M          = 399.999;
  if (TLep2Jet1_M       >= 300)         TLep2Jet1_M       = 299.999;
  if (TLep1Jet2_M        >= 400)         TLep1Jet2_M        = 399.999;
  if (TLep2Jet2_M     >= 300)         TLep2Jet2_M     = 299.999;
  if (TMiniMax >= 420)         TMiniMax = 419.999;
  if (TLep1Lep2Jet1_M               >= 400)         TLep1Lep2Jet1_M               = 399.999;
  if (TLep1Lep2METJet1_Mt           >= 500)         TLep1Lep2METJet1_Mt           = 499.999;
  if (TLep1Lep2Jet1_E               >= 700)         TLep1Lep2Jet1_E               = 699.999;
  if (TLep1Lep2_Pt             >= 200)         TLep1Lep2_Pt             = 199.999;
  if (Lep1Lep2Jet1_Pt           >= 200)         Lep1Lep2Jet1_Pt           = 199.999;
  if (Lep1Lep2METJet1_Pt        >= 150)         Lep1Lep2METJet1_Pt        = 149.999;
  if (THTtot               >= 600)         THTtot               = 599.999;
  if (TJet1_Pt        >= 300)         TJet1_Pt        = 299.999;
  if (TJet1_E         >= 400)         TJet1_E         = 399.999;
  if (TJet1_Phi       >= TMath::Pi()) TJet1_Phi       = 3.14;
  if (TJet1_Eta       >= 2.4)         TJet1_Eta       = 2.39999;
  if (TLep1_Pt        >= 300)         TLep1_Pt        = 299.999;
  if (TLep1_E         >= 350)         TLep1_E         = 349.999;
  if (TLep1_Phi       >= TMath::Pi()) TLep1_Phi       = 3.14;
  if (TLep1_Eta       >= 2.4)         TLep1_Eta       = 2.39999;
  if (TLep2_Pt     >= 150)         TLep2_Pt     = 149.999;
  if (TLep2_E      >= 250)         TLep2_E      = 249.999;
  if (TLep2_Phi    >= TMath::Pi()) TLep2_Phi    = 3.14;
  if (TLep2_Eta    >= 2.4)         TLep2_Eta    = 2.39999;
  if (Lep1Lep2METJet1_Pz       >= 450)         Lep1Lep2METJet1_Pz       = 449.999;
  if (TLep1Lep2METJet1_Eta           >= 5)           TLep1Lep2METJet1_Eta           = 4.999;
  if (MSys                 >= 700)         MSys                 = 699.999;
  if (TLep1Lep2_M                 >= 300)         TLep1Lep2_M                 = 299.999;
  if (TLep1Lep2_DPhi              >= TMath::Pi()) TLep1Lep2_DPhi              = 3.14;
  if (TLep1Jet1_DPhi         >= TMath::Pi()) TLep1Jet1_DPhi         = 3.14;
  if (TLep2Jet1_DPhi      >= TMath::Pi()) TLep2Jet1_DPhi      = 3.14;
  
  if (TMET                 < 0)            TMET                 = 0;
  if (TLep1Jet1_M          < 0)            TLep1Jet1_M          = 0;
  if (TLep2Jet1_M       < 0)            TLep2Jet1_M       = 0;
  if (TLep1Jet2_M        < 0)            TLep1Jet2_M        = 0;
  if (TLep2Jet2_M     < 0)            TLep2Jet2_M     = 0;
  if (TMiniMax < 0)            TMiniMax = 0;
  if (TLep1Lep2Jet1_M               < 0)            TLep1Lep2Jet1_M               = 0;
  if (TLep1Lep2METJet1_Mt           < 0)            TLep1Lep2METJet1_Mt           = 0;
  if (TLep1Lep2Jet1_E               < 0)            TLep1Lep2Jet1_E               = 0;
  if (TLep1Lep2_Pt             < 0)            TLep1Lep2_Pt             = 0;
  if (Lep1Lep2Jet1_Pt           < 0)            Lep1Lep2Jet1_Pt           = 0;
  if (Lep1Lep2METJet1_Pt        < 0)            Lep1Lep2METJet1_Pt        = 0;
  if (THTtot               < 0)            THTtot               = 0;
  if (TJet1_Pt        < 0)            TJet1_Pt        = 0;
  if (TJet1_E         < 0)            TJet1_E         = 0;
  if (TJet1_Phi       < -TMath::Pi()) TJet1_Phi       = -3.14;
  if (TJet1_Eta       < -2.4)         TJet1_Eta       = -2.39999;
  if (TLep1_Pt        < 0)            TLep1_Pt        = 0;
  if (TLep1_E         < 0)            TLep1_E         = 0;
  if (TLep1_Phi       < -TMath::Pi()) TLep1_Phi       = -3.14;
  if (TLep1_Eta       < -2.4)         TLep1_Eta       = -2.39999;
  if (TLep2_Pt     < 0)            TLep2_Pt     = 0;
  if (TLep2_E      < 0)            TLep2_E      = 0;
  if (TLep2_Phi    < -TMath::Pi()) TLep2_Phi    = -3.14;
  if (TLep2_Eta    < -2.4)         TLep2_Eta    = -2.39999;
  if (Lep1Lep2METJet1_Pz       < -450)         Lep1Lep2METJet1_Pz       = -449.999;
  if (TLep1Lep2METJet1_Eta           < -5)           TLep1Lep2METJet1_Eta           = -4.999;
  if (MSys                 < 0)            MSys                 = 0;
  if (TLep1Lep2_M                 < 0)            TLep1Lep2_M                 = 0;
  if (TLep1Lep2_DPhi              < -TMath::Pi()) TLep1Lep2_DPhi              = -3.14;
  if (TLep1Jet1_DPhi         < -TMath::Pi()) TLep1Jet1_DPhi         = -3.14;
  if (TLep2Jet1_DPhi      < -TMath::Pi()) TLep2Jet1_DPhi      = -3.14;
  
  
  if (TMETJESUp                 >= 200)         TMETJESUp                 = 199.999;
  if (TLep1Jet1_MJESUp          >= 400)         TLep1Jet1_MJESUp          = 399.999;
  if (TLep2Jet1_MJESUp       >= 300)         TLep2Jet1_MJESUp       = 299.999;
  if (TLep1Lep2Jet1_MJESUp               >= 400)         TLep1Lep2Jet1_MJESUp               = 399.999;
  if (TLep1Lep2METJet1_MtJESUp           >= 500)         TLep1Lep2METJet1_MtJESUp           = 499.999;
  if (TLep1Lep2Jet1_EJESUp               >= 700)         TLep1Lep2Jet1_EJESUp               = 699.999;
  if (TLep1Lep2_PtJESUp             >= 200)         TLep1Lep2_PtJESUp             = 199.999;
  if (Lep1Lep2Jet1_PtJESUp           >= 200)         Lep1Lep2Jet1_PtJESUp           = 199.999;
  if (Lep1Lep2METJet1_PtJESUp        >= 150)         Lep1Lep2METJet1_PtJESUp        = 149.999;
  if (THTtotJESUp               >= 600)         THTtotJESUp               = 599.999;
  if (TJet1_PtJESUp        >= 150)         TJet1_PtJESUp        = 149.999;
  if (TJet1_EJESUp         >= 400)         TJet1_EJESUp         = 399.999;
  if (TJet1_PhiJESUp       >= TMath::Pi()) TJet1_PhiJESUp       = 3.14;
  if (TJet1_EtaJESUp       >= 2.4)         TJet1_EtaJESUp       = 2.39999;
  if (TLep1_PtJESUp        >= 150)         TLep1_PtJESUp        = 149.999;
  if (TLep1_EJESUp         >= 350)         TLep1_EJESUp         = 349.999;
  if (TLep1_PhiJESUp       >= TMath::Pi()) TLep1_PhiJESUp       = 3.14;
  if (TLep1_EtaJESUp       >= 2.4)         TLep1_EtaJESUp       = 2.39999;
  if (TLep2_PtJESUp     >= 150)         TLep2_PtJESUp     = 149.999;
  if (TLep2_EJESUp      >= 250)         TLep2_EJESUp      = 249.999;
  if (TLep2_PhiJESUp    >= TMath::Pi()) TLep2_PhiJESUp    = 3.14;
  if (TLep2_EtaJESUp    >= 2.4)         TLep2_EtaJESUp    = 2.39999;
  if (Lep1Lep2METJet1_PzJESUp       >= 450)         Lep1Lep2METJet1_PzJESUp       = 449.999;
  if (TLep1Lep2METJet1_EtaJESUp           >= 5)           TLep1Lep2METJet1_EtaJESUp           = 4.999;
  if (MSysJESUp                 >= 700)         MSysJESUp                 = 699.999;
  if (TLep1Lep2_MJESUp                 >= 300)         TLep1Lep2_MJESUp                 = 299.999;
  if (TLep1Lep2_DPhiJESUp              >= TMath::Pi()) TLep1Lep2_DPhiJESUp              = 3.14;
  if (TLep1Jet1_DPhiJESUp         >= TMath::Pi()) TLep1Jet1_DPhiJESUp         = 3.14;
  if (TLep2Jet1_DPhiJESUp      >= TMath::Pi()) TLep2Jet1_DPhiJESUp      = 3.14;
    
  if (TMETJESUp                 < 0)            TMETJESUp                 = 0;
  if (TLep1Jet1_MJESUp          < 0)            TLep1Jet1_MJESUp          = 0;
  if (TLep2Jet1_MJESUp       < 0)            TLep2Jet1_MJESUp       = 0;
  if (TLep1Lep2Jet1_MJESUp               < 0)            TLep1Lep2Jet1_MJESUp               = 0;
  if (TLep1Lep2METJet1_MtJESUp           < 0)            TLep1Lep2METJet1_MtJESUp           = 0;
  if (TLep1Lep2Jet1_EJESUp               < 0)            TLep1Lep2Jet1_EJESUp               = 0;
  if (TLep1Lep2_PtJESUp             < 0)            TLep1Lep2_PtJESUp             = 0;
  if (Lep1Lep2Jet1_PtJESUp           < 0)            Lep1Lep2Jet1_PtJESUp           = 0;
  if (Lep1Lep2METJet1_PtJESUp        < 0)            Lep1Lep2METJet1_PtJESUp        = 0;
  if (THTtotJESUp               < 0)            THTtotJESUp               = 0;
  if (TJet1_PtJESUp        < 0)            TJet1_PtJESUp        = 0;
  if (TJet1_EJESUp         < 0)            TJet1_EJESUp         = 0;
  if (TJet1_PhiJESUp       < -TMath::Pi()) TJet1_PhiJESUp       = -3.14;
  if (TJet1_EtaJESUp       < -2.4)         TJet1_EtaJESUp       = -2.39999;
  if (TLep1_PtJESUp        < 0)            TLep1_PtJESUp        = 0;
  if (TLep1_EJESUp         < 0)            TLep1_EJESUp         = 0;
  if (TLep1_PhiJESUp       < -TMath::Pi()) TLep1_PhiJESUp       = -3.14;
  if (TLep1_EtaJESUp       < -2.4)         TLep1_EtaJESUp       = -2.39999;
  if (TLep2_PtJESUp     < 0)            TLep2_PtJESUp     = 0;
  if (TLep2_EJESUp      < 0)            TLep2_EJESUp      = 0;
  if (TLep2_PhiJESUp    < -TMath::Pi()) TLep2_PhiJESUp    = -3.14;
  if (TLep2_EtaJESUp    < -2.4)         TLep2_EtaJESUp    = -2.39999;
  if (Lep1Lep2METJet1_PzJESUp       < -450)         Lep1Lep2METJet1_PzJESUp       = -449.999;
  if (TLep1Lep2METJet1_EtaJESUp           < -5)           TLep1Lep2METJet1_EtaJESUp           = -4.999;
  if (MSysJESUp                 < 0)            MSysJESUp                 = 0;
  if (TLep1Lep2_MJESUp                 < 0)            TLep1Lep2_MJESUp                 = 0;
  if (TLep1Lep2_DPhiJESUp              < -TMath::Pi()) TLep1Lep2_DPhiJESUp              = -3.14;
  if (TLep1Jet1_DPhiJESUp         < -TMath::Pi()) TLep1Jet1_DPhiJESUp         = -3.14;
  if (TLep2Jet1_DPhiJESUp      < -TMath::Pi()) TLep2Jet1_DPhiJESUp      = -3.14;
  
  
  if (TMETJESDown                 >= 200)         TMETJESDown                 = 199.999;
  if (TLep1Jet1_MJESDown          >= 400)         TLep1Jet1_MJESDown          = 399.999;
  if (TLep2Jet1_MJESDown       >= 300)         TLep2Jet1_MJESDown       = 299.999;
  if (TLep1Lep2Jet1_MJESDown               >= 400)         TLep1Lep2Jet1_MJESDown               = 399.999;
  if (TLep1Lep2METJet1_MtJESDown           >= 500)         TLep1Lep2METJet1_MtJESDown           = 499.999;
  if (TLep1Lep2Jet1_EJESDown               >= 700)         TLep1Lep2Jet1_EJESDown               = 699.999;
  if (TLep1Lep2_PtJESDown             >= 200)         TLep1Lep2_PtJESDown             = 199.999;
  if (Lep1Lep2Jet1_PtJESDown           >= 200)         Lep1Lep2Jet1_PtJESDown           = 199.999;
  if (Lep1Lep2METJet1_PtJESDown        >= 150)         Lep1Lep2METJet1_PtJESDown        = 149.999;
  if (THTtotJESDown               >= 600)         THTtotJESDown               = 599.999;
  if (TJet1_PtJESDown        >= 150)         TJet1_PtJESDown        = 149.999;
  if (TJet1_EJESDown         >= 400)         TJet1_EJESDown         = 399.999;
  if (TJet1_PhiJESDown       >= TMath::Pi()) TJet1_PhiJESDown       = 3.14;
  if (TJet1_EtaJESDown       >= 2.4)         TJet1_EtaJESDown       = 2.39999;
  if (TLep1_PtJESDown        >= 150)         TLep1_PtJESDown        = 149.999;
  if (TLep1_EJESDown         >= 350)         TLep1_EJESDown         = 349.999;
  if (TLep1_PhiJESDown       >= TMath::Pi()) TLep1_PhiJESDown       = 3.14;
  if (TLep1_EtaJESDown       >= 2.4)         TLep1_EtaJESDown       = 2.39999;
  if (TLep2_PtJESDown     >= 150)         TLep2_PtJESDown     = 149.999;
  if (TLep2_EJESDown      >= 250)         TLep2_EJESDown      = 249.999;
  if (TLep2_PhiJESDown    >= TMath::Pi()) TLep2_PhiJESDown    = 3.14;
  if (TLep2_EtaJESDown    >= 2.4)         TLep2_EtaJESDown    = 2.39999;
  if (Lep1Lep2METJet1_PzJESDown       >= 450)         Lep1Lep2METJet1_PzJESDown       = 449.999;
  if (TLep1Lep2METJet1_EtaJESDown           >= 5)           TLep1Lep2METJet1_EtaJESDown           = 4.999;
  if (MSysJESDown                 >= 700)         MSysJESDown                 = 699.999;
  if (TLep1Lep2_MJESDown                 >= 300)         TLep1Lep2_MJESDown                 = 299.999;
  if (TLep1Lep2_DPhiJESDown              >= TMath::Pi()) TLep1Lep2_DPhiJESDown              = 3.14;
  if (TLep1Jet1_DPhiJESDown         >= TMath::Pi()) TLep1Jet1_DPhiJESDown         = 3.14;
  if (TLep2Jet1_DPhiJESDown      >= TMath::Pi()) TLep2Jet1_DPhiJESDown      = 3.14;
    
  if (TMETJESDown                 < 0)            TMETJESDown                 = 0;
  if (TLep1Jet1_MJESDown          < 0)            TLep1Jet1_MJESDown          = 0;
  if (TLep2Jet1_MJESDown       < 0)            TLep2Jet1_MJESDown       = 0;
  if (TLep1Lep2Jet1_MJESDown               < 0)            TLep1Lep2Jet1_MJESDown               = 0;
  if (TLep1Lep2METJet1_MtJESDown           < 0)            TLep1Lep2METJet1_MtJESDown           = 0;
  if (TLep1Lep2Jet1_EJESDown               < 0)            TLep1Lep2Jet1_EJESDown               = 0;
  if (TLep1Lep2_PtJESDown             < 0)            TLep1Lep2_PtJESDown             = 0;
  if (Lep1Lep2Jet1_PtJESDown           < 0)            Lep1Lep2Jet1_PtJESDown           = 0;
  if (Lep1Lep2METJet1_PtJESDown        < 0)            Lep1Lep2METJet1_PtJESDown        = 0;
  if (THTtotJESDown               < 0)            THTtotJESDown               = 0;
  if (TJet1_PtJESDown        < 0)            TJet1_PtJESDown        = 0;
  if (TJet1_EJESDown         < 0)            TJet1_EJESDown         = 0;
  if (TJet1_PhiJESDown       < -TMath::Pi()) TJet1_PhiJESDown       = -3.14;
  if (TJet1_EtaJESDown       < -2.4)         TJet1_EtaJESDown       = -2.39999;
  if (TLep1_PtJESDown        < 0)            TLep1_PtJESDown        = 0;
  if (TLep1_EJESDown         < 0)            TLep1_EJESDown         = 0;
  if (TLep1_PhiJESDown       < -TMath::Pi()) TLep1_PhiJESDown       = -3.14;
  if (TLep1_EtaJESDown       < -2.4)         TLep1_EtaJESDown       = -2.39999;
  if (TLep2_PtJESDown     < 0)            TLep2_PtJESDown     = 0;
  if (TLep2_EJESDown      < 0)            TLep2_EJESDown      = 0;
  if (TLep2_PhiJESDown    < -TMath::Pi()) TLep2_PhiJESDown    = -3.14;
  if (TLep2_EtaJESDown    < -2.4)         TLep2_EtaJESDown    = -2.39999;
  if (Lep1Lep2METJet1_PzJESDown       < -450)         Lep1Lep2METJet1_PzJESDown       = -449.999;
  if (TLep1Lep2METJet1_EtaJESDown           < -5)           TLep1Lep2METJet1_EtaJESDown           = -4.999;
  if (MSysJESDown                 < 0)            MSysJESDown                 = 0;
  if (TLep1Lep2_MJESDown                 < 0)            TLep1Lep2_MJESDown                 = 0;
  if (TLep1Lep2_DPhiJESDown              < -TMath::Pi()) TLep1Lep2_DPhiJESDown              = -3.14;
  if (TLep1Jet1_DPhiJESDown         < -TMath::Pi()) TLep1Jet1_DPhiJESDown         = -3.14;
  if (TLep2Jet1_DPhiJESDown      < -TMath::Pi()) TLep2Jet1_DPhiJESDown      = -3.14;
  
  
  if (TMETJERUp                 >= 200)         TMETJERUp                 = 199.999;
  if (TLep1Jet1_MJERUp          >= 400)         TLep1Jet1_MJERUp          = 399.999;
  if (TLep2Jet1_MJERUp       >= 300)         TLep2Jet1_MJERUp       = 299.999;
  if (TLep1Lep2Jet1_MJERUp               >= 400)         TLep1Lep2Jet1_MJERUp               = 399.999;
  if (TLep1Lep2METJet1_MtJERUp           >= 500)         TLep1Lep2METJet1_MtJERUp           = 499.999;
  if (TLep1Lep2Jet1_EJERUp               >= 700)         TLep1Lep2Jet1_EJERUp               = 699.999;
  if (TLep1Lep2_PtJERUp             >= 200)         TLep1Lep2_PtJERUp             = 199.999;
  if (Lep1Lep2Jet1_PtJERUp           >= 200)         Lep1Lep2Jet1_PtJERUp           = 199.999;
  if (Lep1Lep2METJet1_PtJERUp        >= 150)         Lep1Lep2METJet1_PtJERUp        = 149.999;
  if (THTtotJERUp               >= 600)         THTtotJERUp               = 599.999;
  if (TJet1_PtJERUp        >= 150)         TJet1_PtJERUp        = 149.999;
  if (TJet1_EJERUp         >= 400)         TJet1_EJERUp         = 399.999;
  if (TJet1_PhiJERUp       >= TMath::Pi()) TJet1_PhiJERUp       = 3.14;
  if (TJet1_EtaJERUp       >= 2.4)         TJet1_EtaJERUp       = 2.39999;
  if (TLep1_PtJERUp        >= 150)         TLep1_PtJERUp        = 149.999;
  if (TLep1_EJERUp         >= 350)         TLep1_EJERUp         = 349.999;
  if (TLep1_PhiJERUp       >= TMath::Pi()) TLep1_PhiJERUp       = 3.14;
  if (TLep1_EtaJERUp       >= 2.4)         TLep1_EtaJERUp       = 2.39999;
  if (TLep2_PtJERUp     >= 150)         TLep2_PtJERUp     = 149.999;
  if (TLep2_EJERUp      >= 250)         TLep2_EJERUp      = 249.999;
  if (TLep2_PhiJERUp    >= TMath::Pi()) TLep2_PhiJERUp    = 3.14;
  if (TLep2_EtaJERUp    >= 2.4)         TLep2_EtaJERUp    = 2.39999;
  if (Lep1Lep2METJet1_PzJERUp       >= 450)         Lep1Lep2METJet1_PzJERUp       = 449.999;
  if (TLep1Lep2METJet1_EtaJERUp           >= 5)           TLep1Lep2METJet1_EtaJERUp           = 4.999;
  if (MSysJERUp                 >= 700)         MSysJERUp                 = 699.999;
  if (TLep1Lep2_MJERUp                 >= 300)         TLep1Lep2_MJERUp                 = 299.999;
  if (TLep1Lep2_DPhiJERUp              >= TMath::Pi()) TLep1Lep2_DPhiJERUp              = 3.14;
  if (TLep1Jet1_DPhiJERUp         >= TMath::Pi()) TLep1Jet1_DPhiJERUp         = 3.14;
  if (TLep2Jet1_DPhiJERUp      >= TMath::Pi()) TLep2Jet1_DPhiJERUp      = 3.14;
    
  if (TMETJERUp                 < 0)            TMETJERUp                 = 0;
  if (TLep1Jet1_MJERUp          < 0)            TLep1Jet1_MJERUp          = 0;
  if (TLep2Jet1_MJERUp       < 0)            TLep2Jet1_MJERUp       = 0;
  if (TLep1Lep2Jet1_MJERUp               < 0)            TLep1Lep2Jet1_MJERUp               = 0;
  if (TLep1Lep2METJet1_MtJERUp           < 0)            TLep1Lep2METJet1_MtJERUp           = 0;
  if (TLep1Lep2Jet1_EJERUp               < 0)            TLep1Lep2Jet1_EJERUp               = 0;
  if (TLep1Lep2_PtJERUp             < 0)            TLep1Lep2_PtJERUp             = 0;
  if (Lep1Lep2Jet1_PtJERUp           < 0)            Lep1Lep2Jet1_PtJERUp           = 0;
  if (Lep1Lep2METJet1_PtJERUp        < 0)            Lep1Lep2METJet1_PtJERUp        = 0;
  if (THTtotJERUp               < 0)            THTtotJERUp               = 0;
  if (TJet1_PtJERUp        < 0)            TJet1_PtJERUp        = 0;
  if (TJet1_EJERUp         < 0)            TJet1_EJERUp         = 0;
  if (TJet1_PhiJERUp       < -TMath::Pi()) TJet1_PhiJERUp       = -3.14;
  if (TJet1_EtaJERUp       < -2.4)         TJet1_EtaJERUp       = -2.39999;
  if (TLep1_PtJERUp        < 0)            TLep1_PtJERUp        = 0;
  if (TLep1_EJERUp         < 0)            TLep1_EJERUp         = 0;
  if (TLep1_PhiJERUp       < -TMath::Pi()) TLep1_PhiJERUp       = -3.14;
  if (TLep1_EtaJERUp       < -2.4)         TLep1_EtaJERUp       = -2.39999;
  if (TLep2_PtJERUp     < 0)            TLep2_PtJERUp     = 0;
  if (TLep2_EJERUp      < 0)            TLep2_EJERUp      = 0;
  if (TLep2_PhiJERUp    < -TMath::Pi()) TLep2_PhiJERUp    = -3.14;
  if (TLep2_EtaJERUp    < -2.4)         TLep2_EtaJERUp    = -2.39999;
  if (Lep1Lep2METJet1_PzJERUp       < -450)         Lep1Lep2METJet1_PzJERUp       = -449.999;
  if (TLep1Lep2METJet1_EtaJERUp           < -5)           TLep1Lep2METJet1_EtaJERUp           = -4.999;
  if (MSysJERUp                 < 0)            MSysJERUp                 = 0;
  if (TLep1Lep2_MJERUp                 < 0)            TLep1Lep2_MJERUp                 = 0;
  if (TLep1Lep2_DPhiJERUp              < -TMath::Pi()) TLep1Lep2_DPhiJERUp              = -3.14;
  if (TLep1Jet1_DPhiJERUp         < -TMath::Pi()) TLep1Jet1_DPhiJERUp         = -3.14;
  if (TLep2Jet1_DPhiJERUp      < -TMath::Pi()) TLep2Jet1_DPhiJERUp      = -3.14;
  
  
  if (TGenMET              >= 200)         TGenMET              = 199.999;
  if (TDressLep1Jet1_M       >= 400)         TDressLep1Jet1_M       = 399.999;
  if (TDressLep2Jet1_M    >= 300)         TDressLep2Jet1_M    = 299.999;
  if (TGenLep1Lep2Jet1_M            >= 400)         TGenLep1Lep2Jet1_M            = 399.999;
  if (TDressLep1Jet2_M     >= 400)         TDressLep1Jet2_M     = 399.999;
  if (TDressLep2Jet2_M  >= 300)         TDressLep2Jet2_M  = 299.999;
  if (TDressMiniMax>=420)        TDressMiniMax=419.999;
  if (TGenLep1Lep2METJet1_Mt        >= 500)         TGenLep1Lep2METJet1_Mt        = 499.999;
  if (TDressLep1Lep2Jet1_E            >= 700)         TDressLep1Lep2Jet1_E            = 699.999;
  if (TGenDilepPt          >= 200)         TGenDilepPt          = 199.999;
  if (TGenLep1Lep2Jet1_Pt       >= 200)         TGenLep1Lep2Jet1_Pt       = 199.999;
  if (TGenLep1Lep2METJet1_Pt    >= 150)         TGenLep1Lep2METJet1_Pt    = 149.999;
  if (TGenHTtot            >= 600)         TGenHTtot            = 599.999;
  if (TDressJet1_Pt     >= 300)         TDressJet1_Pt     = 299.999;
  if (TDressJet1_E      >= 400)         TDressJet1_E      = 399.999;
  if (TDressJet1_Phi    >= TMath::Pi()) TDressJet1_Phi    = 3.14;
  if (TDressJet1_Eta    >= 2.4)         TDressJet1_Eta    = 2.39999;
  if (TDressLep1_Pt     >= 300)         TDressLep1_Pt     = 299.999;
  if (TDressLep1_E      >= 350)         TDressLep1_E      = 349.999;
  if (TDressLep1_Phi    >= TMath::Pi()) TDressLep1_Phi    = 3.14;
  if (TDressLep1_Eta    >= 2.4)         TDressLep1_Eta    = 2.39999;
  if (TDressLep2_Pt  >= 150)         TDressLep2_Pt  = 149.999;
  if (TDressLep2_E   >= 250)         TDressLep2_E   = 249.999;
  if (TDressLep2_Phi >= TMath::Pi()) TDressLep2_Phi = 3.14;
  if (TDressLep2_Eta >= 2.4)         TDressLep2_Eta = 2.39999;
  if (TGenLep1Lep2METJet1_Pz   >= 450)         TGenLep1Lep2METJet1_Pz   = 449.999;
  if (TGenLep1Lep2METJet1_Eta        >= 5)           TGenLep1Lep2METJet1_Eta        = 4.999;
  if (TGenMSys             >= 700)         TGenMSys             = 699.999;
  if (TGenLep1Lep2_M              >= 300)         TGenLep1Lep2_M              = 299.999;
  if (TGenDPhiLL           >= TMath::Pi()) TGenDPhiLL           = 3.14;
  if (TGenDPhiLeadJet      >= TMath::Pi()) TGenDPhiLeadJet      = 3.14;
  if (TGenDPhiSubLeadJet   >= TMath::Pi()) TGenDPhiSubLeadJet   = 3.14;
  
  if (TGenMET              < 0)            TGenMET              = 0;
  if (TDressLep1Jet1_M       < 0)            TDressLep1Jet1_M       = 0;
  if (TDressLep2Jet1_M    < 0)            TDressLep2Jet1_M    = 0;
  if (TGenLep1Lep2Jet1_M            < 0)            TGenLep1Lep2Jet1_M            = 0;
  if (TDressLep1Jet2_M     < 0)            TDressLep1Jet2_M     = 0;
  if (TDressLep2Jet2_M  < 0)            TDressLep2Jet2_M  = 0;
  if (TDressMiniMax<0)           TDressMiniMax=0;
  if (TGenLep1Lep2METJet1_Mt        < 0)            TGenLep1Lep2METJet1_Mt        = 0;
  if (TDressLep1Lep2Jet1_E            < 0)            TDressLep1Lep2Jet1_E            = 0;
  if (TGenDilepPt          < 0)            TGenDilepPt          = 0;
  if (TGenLep1Lep2Jet1_Pt       < 0)            TGenLep1Lep2Jet1_Pt       = 0;
  if (TGenLep1Lep2METJet1_Pt    < 0)            TGenLep1Lep2METJet1_Pt    = 0;
  if (TGenHTtot            < 0)            TGenHTtot            = 0;
  if (TDressJet1_Pt     < 0)            TDressJet1_Pt     = 0;
  if (TDressJet1_E      < 0)            TDressJet1_E      = 0;
  if (TDressJet1_Phi    < -TMath::Pi()) TDressJet1_Phi    = -3.14;
  if (TDressJet1_Eta    < -2.4)         TDressJet1_Eta    = -2.39999;
  if (TDressLep1_Pt     < 0)            TDressLep1_Pt     = 0;
  if (TDressLep1_E      < 0)            TDressLep1_E      = 0;
  if (TDressLep1_Phi    < -TMath::Pi()) TDressLep1_Phi    = -3.14;
  if (TDressLep1_Eta    < -2.4)         TDressLep1_Eta    = -2.39999;
  if (TDressLep2_Pt  < 0)            TDressLep2_Pt  = 0;
  if (TDressLep2_E   < 0)            TDressLep2_E   = 0;
  if (TDressLep2_Phi < -TMath::Pi()) TDressLep2_Phi = -3.14;
  if (TDressLep2_Eta < -2.4)         TDressLep2_Eta = -2.39999;
  if (TGenLep1Lep2METJet1_Pz   < -450)         TGenLep1Lep2METJet1_Pz   = -449.999;
  if (TGenLep1Lep2METJet1_Eta        < -5)           TGenLep1Lep2METJet1_Eta        = -4.999;
  if (TGenMSys             < 0)            TGenMSys             = 0;
  if (TGenLep1Lep2_M              < 0)            TGenLep1Lep2_M              = 0;
  if (TGenDPhiLL           < -TMath::Pi()) TGenDPhiLL           = -3.14;
  if (TGenDPhiLeadJet      < -TMath::Pi()) TGenDPhiLeadJet      = -3.14;
  if (TGenDPhiSubLeadJet   < -TMath::Pi()) TGenDPhiSubLeadJet   = -3.14;
  
  //   Setting dummy value for gen events that don't pass the reco selection for
  // unfolding procedures.
  if (gPar.Contains("Unfolding")) {
    if (TPassGen && !TPassReco) {
        TMET              = 99999;
        TMET_Phi          = 99999;
        TLep1Jet1_M       = 99999;
        TLep2Jet1_M    = 99999;
        TLep1Lep2Jet1_M            = 99999;
        TLep1Jet2_M     = 99999;
        TLep2Jet2_M  = 99999;
        TMiniMax=99999;
        TLep1Lep2METJet1_Mt        = 99999;
        TLep1Lep2Jet1_E            = 99999;
        TLep1Lep2_Pt          = 99999;
        Lep1Lep2Jet1_Pt        = 99999;
        Lep1Lep2METJet1_Pt     = 99999;
        THTtot            = 99999;
        TJet1_Pt     = 99999;
        TJet1_E      = 99999;
        TJet1_Phi    = 99999;
        TJet1_Eta    = 99999;
        TLep1_Pt     = 99999;
        TLep1_E      = 99999;
        TLep1_Phi    = 99999;
        TLep1_Eta    = 99999;
        TLep2_Pt  = 99999;
        TLep2_E   = 99999;
        TLep2_Phi = 99999;
        TLep2_Eta = 99999;
        Lep1Lep2METJet1_Pz    = 99999;
        TLep1Lep2METJet1_Eta        = 99999;
        MSys              = 99999;
        TLep1Lep2_M              = 99999;
        TLep1Lep2_DPhi           = 99999;
        TLep1Jet1_DPhi      = 99999;
        TLep2Jet1_DPhi   = 99999;
    }
    if (TPassGen && !TPassRecoJESUp) {
        TMETJESUp              = 99999;
        TMET_PhiJESUp          = 99999;
        TLep1Jet1_MJESUp       = 99999;
        TLep2Jet1_MJESUp    = 99999;
        TLep1Lep2Jet1_MJESUp            = 99999;
        TLep1Lep2METJet1_MtJESUp        = 99999;
        TLep1Lep2Jet1_EJESUp            = 99999;
        TLep1Lep2_PtJESUp          = 99999;
        Lep1Lep2Jet1_PtJESUp        = 99999;
        Lep1Lep2METJet1_PtJESUp     = 99999;
        THTtotJESUp            = 99999;
        TJet1_PtJESUp     = 99999;
        TJet1_EJESUp      = 99999;
        TJet1_PhiJESUp    = 99999;
        TJet1_EtaJESUp    = 99999;
        TLep1_PtJESUp     = 99999;
        TLep1_EJESUp      = 99999;
        TLep1_PhiJESUp    = 99999;
        TLep1_EtaJESUp    = 99999;
        TLep2_PtJESUp  = 99999;
        TLep2_EJESUp   = 99999;
        TLep2_PhiJESUp = 99999;
        TLep2_EtaJESUp = 99999;
        Lep1Lep2METJet1_PzJESUp    = 99999;
        TLep1Lep2METJet1_EtaJESUp        = 99999;
        MSysJESUp              = 99999;
        TLep1Lep2_MJESUp              = 99999;
        TLep1Lep2_DPhiJESUp           = 99999;
        TLep1Jet1_DPhiJESUp      = 99999;
        TLep2Jet1_DPhiJESUp   = 99999;
    }
    if (TPassGen && !TPassRecoJESDown) {
        TMETJESDown              = 99999;
        TMET_PhiJESDown          = 99999;
        TLep1Jet1_MJESDown       = 99999;
        TLep2Jet1_MJESDown    = 99999;
        TLep1Lep2Jet1_MJESDown            = 99999;
        TLep1Lep2METJet1_MtJESDown        = 99999;
        TLep1Lep2Jet1_EJESDown            = 99999;
        TLep1Lep2_PtJESDown          = 99999;
        Lep1Lep2Jet1_PtJESDown        = 99999;
        Lep1Lep2METJet1_PtJESDown     = 99999;
        THTtotJESDown            = 99999;
        TJet1_PtJESDown     = 99999;
        TJet1_EJESDown      = 99999;
        TJet1_PhiJESDown    = 99999;
        TJet1_EtaJESDown    = 99999;
        TLep1_PtJESDown     = 99999;
        TLep1_EJESDown      = 99999;
        TLep1_PhiJESDown    = 99999;
        TLep1_EtaJESDown    = 99999;
        TLep2_PtJESDown  = 99999;
        TLep2_EJESDown   = 99999;
        TLep2_PhiJESDown = 99999;
        TLep2_EtaJESDown = 99999;
        Lep1Lep2METJet1_PzJESDown    = 99999;
        TLep1Lep2METJet1_EtaJESDown        = 99999;
        MSysJESDown              = 99999;
        TLep1Lep2_MJESDown              = 99999;
        TLep1Lep2_DPhiJESDown           = 99999;
        TLep1Jet1_DPhiJESDown      = 99999;
        TLep2Jet1_DPhiJESDown   = 99999;
    }
    if (TPassGen && !TPassRecoJERUp) {
        TMETJERUp              = 99999;
        TMET_PhiJERUp          = 99999;
        TLep1Jet1_MJERUp       = 99999;
        TLep2Jet1_MJERUp    = 99999;
        TLep1Lep2Jet1_MJERUp            = 99999;
        TLep1Lep2METJet1_MtJERUp        = 99999;
        TLep1Lep2Jet1_EJERUp            = 99999;
        TLep1Lep2_PtJERUp          = 99999;
        Lep1Lep2Jet1_PtJERUp        = 99999;
        Lep1Lep2METJet1_PtJERUp     = 99999;
        THTtotJERUp            = 99999;
        TJet1_PtJERUp     = 99999;
        TJet1_EJERUp      = 99999;
        TJet1_PhiJERUp    = 99999;
        TJet1_EtaJERUp    = 99999;
        TLep1_PtJERUp     = 99999;
        TLep1_EJERUp      = 99999;
        TLep1_PhiJERUp    = 99999;
        TLep1_EtaJERUp    = 99999;
        TLep2_PtJERUp  = 99999;
        TLep2_EJERUp   = 99999;
        TLep2_PhiJERUp = 99999;
        TLep2_EtaJERUp = 99999;
        Lep1Lep2METJet1_PzJERUp    = 99999;
        TLep1Lep2METJet1_EtaJERUp        = 99999;
        MSysJERUp              = 99999;
        TLep1Lep2_MJERUp              = 99999;
        TLep1Lep2_DPhiJERUp           = 99999;
        TLep1Jet1_DPhiJERUp      = 99999;
        TLep2Jet1_DPhiJERUp   = 99999;
    }
  }
//   if (TPassReco || TPassGen) {
  
  if (TPassGen || TPassReco || TPassRecoJESUp || TPassRecoJESDown || TPassRecoJERUp 
      || (TNJets == 2 && TNBJets == 2) || (TNJetsJESUp == 1 && TNBJetsJESUp == 1) 
      || (TNJetsJESDown == 1 && TNBJetsJESDown == 1) || (TNJetsJERUp == 1 && TNBJetsJERUp == 1))  { // If needed, filling.
    ffMiniTree->Fill();
  }
}



void TWTTbarAnalysis::Summary(){}            //=============== Summary




//#####################################################################
// Functions
//---------------------------------------------------------------------
void TWTTbarAnalysis::GetLeptonVariables(){
  TNSelLeps = selLeptons.size();
  for(Int_t i = 0; i < TNSelLeps; i++){
    TLep_Pt[i]     = selLeptons.at(i).Pt();
    TLep_Eta[i]    = selLeptons.at(i).Eta();
    TLep_Phi[i]    = selLeptons.at(i).Phi();
    TLep_E[i]      = selLeptons.at(i).E();
    TLep_Charge[i] = selLeptons.at(i).charge;
  }
  if(TNSelLeps < 2) gChannel = -1;
  else if(selLeptons.at(0).isMuon && selLeptons.at(1).isElec) gChannel = iElMu;
  else if(selLeptons.at(0).isElec && selLeptons.at(1).isMuon) gChannel = iElMu;
  else if(selLeptons.at(0).isMuon && selLeptons.at(1).isMuon) gChannel = iMuon;
  else if(selLeptons.at(0).isElec && selLeptons.at(1).isElec) gChannel = iElec;
  if(TNSelLeps > 1) TLep1Lep2_M = (selLeptons.at(0).p + selLeptons.at(1).p).M();
  TChannel = gChannel;
  TIsSS = isSS;
  gChannel = gChannel -1; // gchannel used for chan index of histograms
  
  bool TIsOSDilep = false;
  if (TNSelLeps >= 2)
    TIsOSDilep = passTrigger && passMETfilters && (!isSS) && ((selLeptons.at(0).p + selLeptons.at(1).p).M() > 20) && selLeptons.at(0).Pt() > 25;
  else
    TIsOSDilep = false;
  
  if (TNSelLeps >= 1) {
    TLep1_Pt          = selLeptons.at(0).Pt();
    TLep1_E           = selLeptons.at(0).E();
    TLep1_Phi         = selLeptons.at(0).Phi();
    TLep1_Eta         = selLeptons.at(0).Eta();
    TLep1_PtJESUp     = selLeptons.at(0).Pt();
    TLep1_EJESUp      = selLeptons.at(0).E();
    TLep1_PhiJESUp    = selLeptons.at(0).Phi();
    TLep1_EtaJESUp    = selLeptons.at(0).Eta();
    TLep1_PtJESDown   = selLeptons.at(0).Pt();
    TLep1_EJESDown    = selLeptons.at(0).E();
    TLep1_PhiJESDown  = selLeptons.at(0).Phi();
    TLep1_EtaJESDown  = selLeptons.at(0).Eta();
    TLep1_PtJERUp     = selLeptons.at(0).Pt();
    TLep1_EJERUp      = selLeptons.at(0).E();
    TLep1_PhiJERUp    = selLeptons.at(0).Phi();
    TLep1_EtaJERUp    = selLeptons.at(0).Eta();
    if (TNSelLeps >= 2) {
      TLep2_Pt          = selLeptons.at(1).Pt();
      TLep2_E           = selLeptons.at(1).E();
      TLep2_Phi         = selLeptons.at(1).Phi();
      TLep2_Eta         = selLeptons.at(1).Eta();
      TLep2_PtJESUp     = selLeptons.at(1).Pt();
      TLep2_EJESUp      = selLeptons.at(1).E();
      TLep2_PhiJESUp    = selLeptons.at(1).Phi();
      TLep2_EtaJESUp    = selLeptons.at(1).Eta();
      TLep2_PtJESDown   = selLeptons.at(1).Pt();
      TLep2_EJESDown    = selLeptons.at(1).E();
      TLep2_PhiJESDown  = selLeptons.at(1).Phi();
      TLep2_EtaJESDown  = selLeptons.at(1).Eta();
      TLep2_PtJERUp     = selLeptons.at(1).Pt();
      TLep2_EJERUp      = selLeptons.at(1).E();
      TLep2_PhiJERUp    = selLeptons.at(1).Phi();
      TLep2_EtaJERUp    = selLeptons.at(1).Eta();
    }
  }
  
  SetParam("TIsOSDilep", TIsOSDilep);
}


void TWTTbarAnalysis::GetJetVariables() {
  TNJets = selJets.size(); THT = 0;
  TNBJets = 0;
  THTJESUp = 0; THTJESDown = 0; 
  for (Int_t i = 0; i < TNJets; i++) {
    TJet_Pt[i]     = selJets.at(i).Pt();
    TJet_Eta[i]    = selJets.at(i).Eta();
    TJet_Phi[i]    = selJets.at(i).Phi();
    TJet_E[i]      = selJets.at(i).E();
    TJet_isBJet[i] = selJets.at(i).isBtag;
    THT += TJet_Pt[i];
    if(selJets.at(i).isBtag) TNBJets++;
  }
  
  if (TNJets > 0) {
    TJet1_Pt  = selJets.at(0).Pt();
    TJet1_E   = selJets.at(0).E();
    TJet1_Phi = selJets.at(0).Phi();
    TJet1_Eta = selJets.at(0).Eta();
    TJet1_CSV = selJets.at(0).csv;
  }
  if (TNJets > 1) {
    TJet2_Pt  = selJets.at(1).Pt();
    TJet2_Eta = selJets.at(1).Eta();
    TJet2_CSV = selJets.at(1).csv;
  }
  
  TNVetoJets = vetoJets.size();
  if (TNVetoJets > 0) {
    TVetoJet1_Pt     = vetoJets.at(0).Pt();
    TVetoJet1_Eta    = vetoJets.at(0).Eta();
  }
  else{
    TVetoJet1_Pt     = -99.;
    TVetoJet1_Eta    = -99.;
  }

  if (TNVetoJets > 1){
    TVetoJet2_Pt     = vetoJets.at(1).Pt();
    TVetoJet2_Eta    = vetoJets.at(1).Eta();
  }
  else{
    TVetoJet2_Pt     = -99.;
    TVetoJet2_Eta    = -99.;
  }

  if (TNVetoJets > 2){
    TVetoJet3_Pt     = vetoJets.at(2).Pt();
    TVetoJet3_Eta    = vetoJets.at(2).Eta();
  }
  else{
    TVetoJet3_Pt     = -99.;
    TVetoJet3_Eta    = -99.;
  }

  SetParam("THT",THT);

  if (gIsData) return;
  
  TNJetsJESUp     = 0;
  TNJetsJESDown   = 0;
  TNBJetsJESUp    = 0;
  TNBJetsJESDown  = 0;
  TNBJetsJERUp    = 0;
  TNJetsJERUp     = 0;

  TNJetsJESUp   = selJetsJecUp.size();
  TNJetsJESDown = selJetsJecDown.size();
  TNJetsJERUp   = selJetsJER.size();
  
  for (auto& jet : selJetsJecUp) {
    THTJESUp += jet.pTJESUp;
    if (jet.isBtag) TNBJetsJESUp++;
  }
  for (auto& jet : selJetsJecDown) {
    THTJESDown += jet.pTJESDown;
    if (jet.isBtag) TNBJetsJESDown++;
  }
  for (auto& jet : selJetsJER) {
    if (jet.isBtag) TNBJetsJERUp++;
  }
  
  if (!gIsData) {
    if (TNJetsJESUp >= 1) {
      TJet1_PtJESUp  = selJetsJecUp.at(0).Pt();
      TJet1_EJESUp   = selJetsJecUp.at(0).E();
      TJet1_PhiJESUp = selJetsJecUp.at(0).Phi();
      TJet1_EtaJESUp = selJetsJecUp.at(0).Eta();
    }
    if (TNJetsJESDown >= 1) {
      TJet1_PtJESDown  = selJetsJecDown.at(0).Pt();
      TJet1_JESDown    = selJetsJecDown.at(0).E();
      TJet1_PhiJESDown = selJetsJecDown.at(0).Phi();
      TJet1_EtaJESDown = selJetsJecDown.at(0).Eta();
    }
    if (TNJetsJERUp >= 1) {
      TJet1_PtJERUp  = selJetsJER.at(0).Pt();
      TJet1_JERUp    = selJetsJER.at(0).E();
      TJet1_PhiJERUp = selJetsJER.at(0).Phi();
      TJet1_EtaJERUp = selJetsJER.at(0).Eta();
    }
  }
}


void TWTTbarAnalysis::GetGenJetVariables() {  // TERMINAR DE REHACER DESDE JET SELECTOR
  if (gIsData) return;
  nGenJets = genJets.size();
  DressNJets = 0; DressNBJets = 0; DressNLooseCentral = 0; DressNBLooseCentral = 0; DressNLooseFwd = 0;
  
  DressJets.clear();
  DressLooseCentralJets.clear();
  DressLooseFwdJets.clear();
  
  for (UShort_t i = 0; i < (UShort_t)nGenJets; i++) {
    if (TMath::Abs(genJets.at(i).p.Eta()) < 2.4 && Cleaning(genJets.at(i), GenLeps, 0.4)) {
      DressLooseCentralJets.push_back(genJets.at(i))
      if (genJets.at(i).isBtag) nGenLooseCentralbJets++;
      if (genJets.at(i).p.Pt() > 30) {
        DressJets.push_back(genJets.at(i));
        if (TMath::Abs(genJets.at(i).flavmc) == 5) DressNBJets++;
      }
    }
    else if (TMath::Abs(genJets.at(i).p.Eta()) < 4.7) {
      DressLooseFwdJets.push_back(genJets.at(i));
    }
  }
  
  DressNJets          = DressJets.size();
  DressNLooseCentral  = DressLooseCentralJets.size();
  DressNLooseFwd      = DressLooseFwdJets.size();
  
  if (DressNLooseCentral >= 2) {
    TDressLooseCentral2_Pt   = DressLooseCentralJets.at(1).Pt();
    TDressLooseCentral2_E    = DressLooseCentralJets.at(1).E();
    TDressLooseCentral2_Phi  = DressLooseCentralJets.at(1).Phi();
    TDressLooseCentral2_Eta  = DressLooseCentralJets.at(1).Eta();
  }
  if (DressNJets >= 1) {
    TDressJet1_Pt   = DressJets.at(0).Pt();
    TDressJet1_E    = DressJets.at(0).E();
    TDressJet1_Phi  = DressJets.at(0).Phi();
    TDressJet1_Eta  = DressJets.at(0).Eta();
  }
  if (DressNLooseFwd >= 1) {
    TDressLooseFwd1_Pt  = DressLooseFwdJets.at(0).Pt();
    TDressLooseFwd1_E   = DressLooseFwdJets.at(0).E();
    TDressLooseFwd1_Phi = DressLooseFwdJets.at(0).Phi();
    TDressLooseFwd1_Eta = DressLooseFwdJets.at(0).Eta();
  }
}


void TWTTbarAnalysis::GetGenLepVariables() {
  if (gIsData) return;
  nGenLeps = genLeptons.size();
  
  if (nGenLeps >= 1) {
    TDressLep1_Pt   = genLeptons.at(0).Pt();
    TDressLep1_E    = genLeptons.at(0).E();
    TDressLep1_Phi  = genLeptons.at(0).Phi();
    TDressLep1_Eta  = genLeptons.at(0).Eta();
    if (nGenLeps >= 2) {
      TDressLep2_Pt   = genLeptons.at(1).Pt();
      TDressLep2_E    = genLeptons.at(1).E();
      TDressLep2_Phi  = genLeptons.at(1).Phi();
      TDressLep2_Eta  = genLeptons.at(1).Eta();
    }
  }
}


Float_t TWTTbarAnalysis::getTopPtRW() { // REVISAR
  if (!gIsTTbar) return 1;
  
  Float_t SF = 1;

  if (Get<Int_t>("nGenTop") > 1) {
    Float_t pt1 = Get<Float_t>("GenTop_pt"  , 0);
    Float_t pt2 = Get<Float_t>("GenTop_pt"  , 1);
    SF *= TMath::Exp(0.0615 - 0.0005 * pt1);
    SF *= TMath::Exp(0.0615 - 0.0005 * pt2);
    return TMath::Sqrt(SF);
  }
  else {
    cout << "Error Error Error get top pt rw" << endl;
    return -1;
  }
}


void TWTTbarAnalysis::GetMETandGenMET() {
  if      (year == 2017) {
    TMET        = Get<Float_t>("METFixEE2017_pt");
    TMET_Phi    = Get<Float_t>("METFixEE2017_phi");
  }
  else if (year == 2018) { // CAMBIAR PA QUE SEA LISTO Y DETECTE EL NOM CUANDO LO HAYA
    TMET        = Get<Float_t>("MET_pt_nom");
    TMET_Phi    = Get<Float_t>("MET_phi_nom");
  }
  
  if (gIsData) return;
//   TMETJESUp   = Get<Float_t>("met_jecUp_pt"  );
//   TMETJESDown = Get<Float_t>("met_jecDown_pt");
//   TMET_PhiJESUp   = Get<Float_t>("met_jecUp_phi"  );
//   TMET_PhiJESDown = Get<Float_t>("met_jecDown_phi");
//   Float_t  diff_MET_JER_phi = GetParam<Float_t>("diff_MET_JER_phi");
//   Float_t  diff_MET_JER_pt  = GetParam<Float_t>("diff_MET_JER_pt");
// 
//   TLorentzVector diff_MET_JER; diff_MET_JER.SetPtEtaPhiM(diff_MET_JER_pt, 0.,diff_MET_JER_phi, 0.);
//   TLorentzVector vMET; vMET.SetPtEtaPhiM(TMET, 0., TMET_Phi, 0);
// 
//   TMET_PhiJERUp     = (vMET + diff_MET_JER).Phi();
//   TMETJERUp         = (vMET + diff_MET_JER).Pt();

  TPartMET      = Get<Float_t>("GenMET_pt");
  TPartMET_Phi  = Get<Float_t>("GenMET_phi");
  TDressMET     = Get<Float_t>("MET_fiducialGenPt");
  TDressMET_Phi = Get<Float_t>("MET_fiducialGenPhi");
  if (gIsLHE)  for(Int_t i = 0; i < Get<Int_t>("nLHEweight"); i++)   TLHEWeight[i] = Get<Float_t>("LHEweight_wgt", i);
}


void TWTTbarAnalysis::SetTWTTbarVariables() {
  // Detector level variables
  fMiniTree->Branch("TChannel",              &TChannel,              "TChannel/I");
  fMiniTree->Branch("TIsSS",                 &TIsSS,                 "TIsSS/B");
  fMiniTree->Branch("TWeight",               &TWeight,               "TWeight/F");
  fMiniTree->Branch("TWeight_ElecEffUp",     &TWeight_ElecEffUp,     "TWeight_ElecEffUp/F");
  fMiniTree->Branch("TWeight_ElecEffDown",   &TWeight_ElecEffDown,   "TWeight_ElecEffDown/F");
  fMiniTree->Branch("TWeight_MuonEffUp",     &TWeight_MuonEffUp,     "TWeight_MuonEffUp/F");
  fMiniTree->Branch("TWeight_MuonEffDown",   &TWeight_MuonEffDown,   "TWeight_MuonEffDown/F");
  fMiniTree->Branch("TWeight_TrigUp",        &TWeight_TrigUp,        "TWeight_TrigUp/F");
  fMiniTree->Branch("TWeight_TrigDown",      &TWeight_TrigDown,      "TWeight_TrigDown/F");
  fMiniTree->Branch("TWeight_PUUp",          &TWeight_PUUp,          "TWeight_PUUp/F");
  fMiniTree->Branch("TWeight_PUDown",        &TWeight_PUDown,        "TWeight_PUDown/F");
  fMiniTree->Branch("TWeight_BtagUp",        &TWeight_BtagUp,        "TWeight_BtagUp/F");
  fMiniTree->Branch("TWeight_BtagDown",      &TWeight_BtagDown,      "TWeight_BtagDown/F");
  fMiniTree->Branch("TWeight_MistagUp",      &TWeight_MistagUp,      "TWeight_MistagUp/F");
  fMiniTree->Branch("TWeight_MistagDown",    &TWeight_MistagDown,    "TWeight_MistagDown/F");
  fMiniTree->Branch("TWeight_normal",        &TWeight_normal,        "TWeight_normal/F");
  fMiniTree->Branch("TLHEWeight",            TLHEWeight,             "TLHEWeight[254]/F");
  
  
  fMiniTree->Branch("TPassReco",             &TPassReco,             "TPassReco/B");
  fMiniTree->Branch("TPassRecoJESUp",        &TPassRecoJESUp,        "TPassRecoJESUp/B");
  fMiniTree->Branch("TPassRecoJESDown",      &TPassRecoJESDown,      "TPassRecoJESDown/B");
  fMiniTree->Branch("TPassRecoJERUp",        &TPassRecoJERUp,        "TPassRecoJERUp/B");
  fMiniTree->Branch("TNJets"       ,         &TNJets,                "TNJets/I"        );
  fMiniTree->Branch("TNJetsJESUp"  ,         &TNJetsJESUp,           "TNJetsJESUp/I"   );
  fMiniTree->Branch("TNJetsJESDown",         &TNJetsJESDown,         "TNJetsJESDown/I" );
  fMiniTree->Branch("TNJetsJERUp",           &TNJetsJERUp,           "TNJetsJERUp/I" );
  fMiniTree->Branch("TNBJets"       ,        &TNBJets,               "TNBJets/I"        );
  fMiniTree->Branch("TNBJetsJESUp"  ,        &TNBJetsJESUp,          "TNBJetsJESUp/I"   );
  fMiniTree->Branch("TNBJetsJESDown",        &TNBJetsJESDown,        "TNBJetsJESDown/I" );
  fMiniTree->Branch("TNBJetsJERUp",          &TNBJetsJERUp,          "TNBJetsJERUp/I" );
  fMiniTree->Branch("TNLooseCentral"        ,&NLooseCentral      ,   "TNLooseCentral/I"       );
  fMiniTree->Branch("TNLooseCentralJESUp"   ,&NLooseCentralJESUp  ,  "TNLooseCentralJESUp/I"       );
  fMiniTree->Branch("TNLooseCentralJESDown" ,&NLooseCentralJESDown,  "TNLooseCentralJESDown/I"       );
  fMiniTree->Branch("TNLooseCentralJERUp"   ,&NLooseCentralJERUp  ,  "TNLooseCentralJERUp/I"       );
  fMiniTree->Branch("TNBLooseCentral"       ,&NBLooseCentral      ,  "TNBLooseCentral/I"       );
  fMiniTree->Branch("TNBLooseCentralJESUp"  ,&NBLooseCentralJESUp  , "TNBLooseCentralJESUp/I"       );
  fMiniTree->Branch("TNBLooseCentralJESDown",&NBLooseCentralJESDown, "TNBLooseCentralJESDown/I");
  fMiniTree->Branch("TNBLooseCentralJERUp"  ,&NBLooseCentralJERUp  , "TNBLooseCentralJERUp/I");
  fMiniTree->Branch("TNLooseFwd"            ,&NLooseFwd            , "TNLooseFwd/I"          );
  fMiniTree->Branch("TNLooseFwdJESUp"       ,&NLooseFwdJESUp        ,"TNLooseFwdJESUp/I"     );
  fMiniTree->Branch("TNLooseFwdJESDown"     ,&NLooseFwdJESDown      ,"TNLooseFwdJESDown/I"   );
  fMiniTree->Branch("TNLooseFwdJERUp"       ,&NLooseFwdJERUp        ,"TNLooseFwdJERUp/I"     );
  
  
  fMiniTree->Branch("TLep1_Pt",              &TLep1_Pt,              "TLep1_Pt/F");
  fMiniTree->Branch("TLep1_E",               &TLep1_E,               "TLep1_E/F");
  fMiniTree->Branch("TLep1_Phi",             &TLep1_Phi,             "TLep1_Phi/F");
  fMiniTree->Branch("TLep1_Eta",             &TLep1_Eta,             "TLep1_Eta/F");
  fMiniTree->Branch("TLep2_Pt",              &TLep2_Pt,              "TLep2_Pt/F");
  fMiniTree->Branch("TLep2_E",               &TLep2_E,               "TLep2_E/F");
  fMiniTree->Branch("TLep2_Phi",             &TLep2_Phi,             "TLep2_Phi/F");
  fMiniTree->Branch("TLep2_Eta",             &TLep2_Eta,             "TLep2_Eta/F");
  fMiniTree->Branch("TJet1_Pt",              &TJet1_Pt,              "TJet1_Pt/F");
  fMiniTree->Branch("TJet1_E",               &TJet1_E,               "TJet1_E/F");
  fMiniTree->Branch("TJet1_Phi",             &TJet1_Phi,             "TJet1_Phi/F");
  fMiniTree->Branch("TJet1_Eta",             &TJet1_Eta,             "TJet1_Eta/F");
  fMiniTree->Branch("TJet2_Pt",              &TJet2_Pt,              "TJet2_Pt/F");
  fMiniTree->Branch("TJet2_E",               &TJet2_E,               "TJet2_E/F");
  fMiniTree->Branch("TJet2_Phi",             &TJet2_Phi,             "TJet2_Phi/F");
  fMiniTree->Branch("TJet2_Eta",             &TJet2_Eta,             "TJet2_Eta/F");
  fMiniTree->Branch("TLooseCentral1_Pt",     &TLooseCentral1_Pt,     "TLooseCentral1_Pt/F");
  fMiniTree->Branch("TLooseFwd1_Pt",         &TLooseFwd1_Pt,         "TLooseFwd1_Pt/F");
  fMiniTree->Branch("TMET",                  &TMET,                  "TMET/F");
  fMiniTree->Branch("TMET_Phi",              &TMET_Phi,              "TMET_Phi/F");
  
  fMiniTree->Branch("TLep1Lep2_Pt",          &TLep1Lep2_Pt,          "TLep1Lep2_Pt/F");
  fMiniTree->Branch("TLep1Lep2_M",           &TLep1Lep2_M,           "TLep1Lep2_M/F");
  fMiniTree->Branch("TLep1Lep2_DPhi",        &TLep1Lep2_DPhi,        "TLep1Lep2_DPhi/F");
  fMiniTree->Branch("TLep1Jet1_M" ,          &TLep1Jet1_M,           "TLep1Jet1_M/F");
  fMiniTree->Branch("TLep1Jet1_DPhi",        &TLep1Jet1_DPhi,        "TLep1Jet1_DPhi/F");
  fMiniTree->Branch("TLep2Jet1_M",           &TLep2Jet1_M,           "TLep2Jet1_M/F");
  fMiniTree->Branch("TLep2Jet1_DPhi",        &TLep2Jet1_DPhi,        "TLep2Jet1_DPhi/F");
  fMiniTree->Branch("TLep1Jet2_M",           &TLep1Jet2_M,           "TLep1Jet2_M/F");
  fMiniTree->Branch("TLep1Jet2_DPhi",        &TLep1Jet2_DPhi,        "TLep1Jet2_DPhi/F");
  fMiniTree->Branch("TLep2Jet2_M",           &TLep2Jet2_M,           "TLep2Jet2_M/F");
  fMiniTree->Branch("TLep2Jet2_DPhi",        &TLep2Jet2_DPhi,        "TLep2Jet2_DPhi/F");
  
  fMiniTree->Branch("TLep1Lep2METJet1_Pt",   &Lep1Lep2METJet1_Pt,    "TLep1Lep2METJet1_Pt/F");
  fMiniTree->Branch("TMSys",                 &MSys,                  "TMSys/F");
  fMiniTree->Branch("TLep1Lep2Jet1_Pt",      &Lep1Lep2Jet1_Pt,       "TLep1Lep2Jet1_Pt/F");
  fMiniTree->Branch("TLep1Lep2METJet1_Pz",   &Lep1Lep2METJet1_Pz,    "Lep1Lep2METJet1_Pz/F");
  fMiniTree->Branch("TLep1Lep2METJet1_Eta",  &TLep1Lep2METJet1_Eta,  "TLep1Lep2METJet1_Eta/F");
  fMiniTree->Branch("TMiniMax",              &TMiniMax,              "TMiniMax/F");
  fMiniTree->Branch("TLep1Lep2Jet1_E",       &TLep1Lep2Jet1_E,       "TLep1Lep2Jet1_E/F");
  fMiniTree->Branch("TLep1Lep2METJet1_Mt",   &TLep1Lep2METJet1_Mt,   "TLep1Lep2METJet1_Mt/F");
  fMiniTree->Branch("TLep1Lep2Jet1_M",       &TLep1Lep2Jet1_M,       "TLep1Lep2Jet1_M/F");
  fMiniTree->Branch("THTtot",                &THTtot,                "THTtot/F");
  
  
  // JESUp
  fMiniTree->Branch("TLep1_PtJESUp",              &TLep1_PtJESUp,              "TLep1_PtJESUp/F");
  fMiniTree->Branch("TLep1_EJESUp",               &TLep1_EJESUp,               "TLep1_EJESUp/F");
  fMiniTree->Branch("TLep1_PhiJESUp",             &TLep1_PhiJESUp,             "TLep1_PhiJESUp/F");
  fMiniTree->Branch("TLep1_EtaJESUp",             &TLep1_EtaJESUp,             "TLep1_EtaJESUp/F");
  fMiniTree->Branch("TLep2_PtJESUp",              &TLep2_PtJESUp,              "TLep2_PtJESUp/F");
  fMiniTree->Branch("TLep2_EJESUp",               &TLep2_EJESUp,               "TLep2_EJESUp/F");
  fMiniTree->Branch("TLep2_PhiJESUp",             &TLep2_PhiJESUp,             "TLep2_PhiJESUp/F");
  fMiniTree->Branch("TLep2_EtaJESUp",             &TLep2_EtaJESUp,             "TLep2_EtaJESUp/F");
  fMiniTree->Branch("TJet1_PtJESUp",              &TJet1_PtJESUp,              "TJet1_PtJESUp/F");
  fMiniTree->Branch("TJet1_EJESUp",               &TJet1_EJESUp,               "TJet1_EJESUp/F");
  fMiniTree->Branch("TJet1_PhiJESUp",             &TJet1_PhiJESUp,             "TJet1_PhiJESUp/F");
  fMiniTree->Branch("TJet1_EtaJESUp",             &TJet1_EtaJESUp,             "TJet1_EtaJESUp/F");
  fMiniTree->Branch("TJet2_PtJESUp",              &TJet2_PtJESUp,              "TJet2_PtJESUp/F");
  fMiniTree->Branch("TJet2_EJESUp",               &TJet2_EJESUp,               "TJet2_EJESUp/F");
  fMiniTree->Branch("TJet2_PhiJESUp",             &TJet2_PhiJESUp,             "TJet2_PhiJESUp/F");
  fMiniTree->Branch("TJet2_EtaJESUp",             &TJet2_EtaJESUp,             "TJet2_EtaJESUp/F");
  fMiniTree->Branch("TLooseCentral1_PtJESUp",     &TLooseCentral1_PtJESUp,     "TLooseCentral1_PtJESUp/F");
  fMiniTree->Branch("TLooseFwd1_PtJESUp",         &TLooseFwd1_PtJESUp,         "TLooseFwd1_PtJESUp/F");
  fMiniTree->Branch("TMETJESUp",                  &TMETJESUp,                  "TMETJESUp/F");
  fMiniTree->Branch("TMET_PhiJESUp",              &TMET_PhiJESUp,              "TMET_PhiJESUp/F");
  
  fMiniTree->Branch("TLep1Lep2_PtJESUp",          &TLep1Lep2_PtJESUp,          "TLep1Lep2_PtJESUp/F");
  fMiniTree->Branch("TLep1Lep2_MJESUp",           &TLep1Lep2_MJESUp,           "TLep1Lep2_MJESUp/F");
  fMiniTree->Branch("TLep1Lep2_DPhiJESUp",        &TLep1Lep2_DPhiJESUp,        "TLep1Lep2_DPhiJESUp/F");
  fMiniTree->Branch("TLep1Jet1_MJESUp",           &TLep1Jet1_MJESUp,           "TLep1Jet1_MJESUp/F");
  fMiniTree->Branch("TLep1Jet1_DPhiJESUp",        &TLep1Jet1_DPhiJESUp,        "TLep1Jet1_DPhiJESUp/F");
  fMiniTree->Branch("TLep2Jet1_MJESUp",           &TLep2Jet1_MJESUp,           "TLep2Jet1_MJESUp/F");
  fMiniTree->Branch("TLep2Jet1_DPhiJESUp",        &TLep2Jet1_DPhiJESUp,        "TLep2Jet1_DPhiJESUp/F");
  fMiniTree->Branch("TLep1Jet2_MJESUp",           &TLep1Jet2_MJESUp,           "TLep1Jet2_MJESUp/F");
  fMiniTree->Branch("TLep1Jet2_DPhiJESUp",        &TLep1Jet2_DPhiJESUp,        "TLep1Jet2_DPhiJESUp/F");
  fMiniTree->Branch("TLep2Jet2_MJESUp",           &TLep2Jet2_MJESUp,           "TLep2Jet2_MJESUp/F");
  fMiniTree->Branch("TLep2Jet2_DPhiJESUp",        &TLep2Jet2_DPhiJESUp,        "TLep2Jet2_DPhiJESUp/F");
  
  fMiniTree->Branch("TLep1Lep2METJet1_PtJESUp",   &Lep1Lep2METJet1_PtJESUp,    "TLep1Lep2METJet1_PtJESUp/F");
  fMiniTree->Branch("TMSysJESUp",                 &MSysJESUp,                  "TMSysJESUp/F");
  fMiniTree->Branch("TLep1Lep2Jet1_PtJESUp",      &Lep1Lep2Jet1_PtJESUp,       "TLep1Lep2Jet1_PtJESUp/F");
  fMiniTree->Branch("TLep1Lep2METJet1_PzJESUp",   &Lep1Lep2METJet1_PzJESUp,    "Lep1Lep2METJet1_PzJESUp/F");
  fMiniTree->Branch("TLep1Lep2METJet1_EtaJESUp",  &TLep1Lep2METJet1_EtaJESUp,  "TLep1Lep2METJet1_EtaJESUp/F");
  fMiniTree->Branch("TMiniMaxJESUp",              &TMiniMaxJESUp,              "TMiniMaxJESUp/F");
  fMiniTree->Branch("TLep1Lep2Jet1_EJESUp",       &TLep1Lep2Jet1_EJESUp,       "TLep1Lep2Jet1_EJESUp/F");
  fMiniTree->Branch("TLep1Lep2METJet1_MtJESUp",   &TLep1Lep2METJet1_MtJESUp,   "TLep1Lep2METJet1_MtJESUp/F");
  fMiniTree->Branch("TLep1Lep2Jet1_MJESUp",       &TLep1Lep2Jet1_MJESUp,       "TLep1Lep2Jet1_MJESUp/F");
  fMiniTree->Branch("THTtotJESUp",                &THTtotJESUp,                "THTtotJESUp/F");
  
  
  // JESDown
  fMiniTree->Branch("TLep1_PtJESDown",              &TLep1_PtJESDown,              "TLep1_PtJESDown/F");
  fMiniTree->Branch("TLep1_EJESDown",               &TLep1_EJESDown,               "TLep1_EJESDown/F");
  fMiniTree->Branch("TLep1_PhiJESDown",             &TLep1_PhiJESDown,             "TLep1_PhiJESDown/F");
  fMiniTree->Branch("TLep1_EtaJESDown",             &TLep1_EtaJESDown,             "TLep1_EtaJESDown/F");
  fMiniTree->Branch("TLep2_PtJESDown",              &TLep2_PtJESDown,              "TLep2_PtJESDown/F");
  fMiniTree->Branch("TLep2_EJESDown",               &TLep2_EJESDown,               "TLep2_EJESDown/F");
  fMiniTree->Branch("TLep2_PhiJESDown",             &TLep2_PhiJESDown,             "TLep2_PhiJESDown/F");
  fMiniTree->Branch("TLep2_EtaJESDown",             &TLep2_EtaJESDown,             "TLep2_EtaJESDown/F");
  fMiniTree->Branch("TJet1_PtJESDown",              &TJet1_PtJESDown,              "TJet1_PtJESDown/F");
  fMiniTree->Branch("TJet1_EJESDown",               &TJet1_EJESDown,               "TJet1_EJESDown/F");
  fMiniTree->Branch("TJet1_PhiJESDown",             &TJet1_PhiJESDown,             "TJet1_PhiJESDown/F");
  fMiniTree->Branch("TJet1_EtaJESDown",             &TJet1_EtaJESDown,             "TJet1_EtaJESDown/F");
  fMiniTree->Branch("TJet2_PtJESDown",              &TJet2_PtJESDown,              "TJet2_PtJESDown/F");
  fMiniTree->Branch("TJet2_EJESDown",               &TJet2_EJESDown,               "TJet2_EJESDown/F");
  fMiniTree->Branch("TJet2_PhiJESDown",             &TJet2_PhiJESDown,             "TJet2_PhiJESDown/F");
  fMiniTree->Branch("TJet2_EtaJESDown",             &TJet2_EtaJESDown,             "TJet2_EtaJESDown/F");
  fMiniTree->Branch("TLooseCentral1_PtJESDown",     &TLooseCentral1_PtJESDown,     "TLooseCentral1_PtJESDown/F");
  fMiniTree->Branch("TLooseFwd1_PtJESDown",         &TLooseFwd1_PtJESDown,         "TLooseFwd1_PtJESDown/F");
  fMiniTree->Branch("TMETJESDown",                  &TMETJESDown,                  "TMETJESDown/F");
  fMiniTree->Branch("TMET_PhiJESDown",              &TMET_PhiJESDown,              "TMET_PhiJESDown/F");
  
  fMiniTree->Branch("TLep1Lep2_PtJESDown",          &TLep1Lep2_PtJESDown,          "TLep1Lep2_PtJESDown/F");
  fMiniTree->Branch("TLep1Lep2_MJESDown",           &TLep1Lep2_MJESDown,           "TLep1Lep2_MJESDown/F");
  fMiniTree->Branch("TLep1Lep2_DPhiJESDown",        &TLep1Lep2_DPhiJESDown,        "TLep1Lep2_DPhiJESDown/F");
  fMiniTree->Branch("TLep1Jet1_MJESDown",           &TLep1Jet1_MJESDown,           "TLep1Jet1_MJESDown/F");
  fMiniTree->Branch("TLep1Jet1_DPhiJESDown",        &TLep1Jet1_DPhiJESDown,        "TLep1Jet1_DPhiJESDown/F");
  fMiniTree->Branch("TLep2Jet1_MJESDown",           &TLep2Jet1_MJESDown,           "TLep2Jet1_MJESDown/F");
  fMiniTree->Branch("TLep2Jet1_DPhiJESDown",        &TLep2Jet1_DPhiJESDown,        "TLep2Jet1_DPhiJESDown/F");
  fMiniTree->Branch("TLep1Jet2_MJESDown",           &TLep1Jet2_MJESDown,           "TLep1Jet2_MJESDown/F");
  fMiniTree->Branch("TLep1Jet2_DPhiJESDown",        &TLep1Jet2_DPhiJESDown,        "TLep1Jet2_DPhiJESDown/F");
  fMiniTree->Branch("TLep2Jet2_MJESDown",           &TLep2Jet2_MJESDown,           "TLep2Jet2_MJESDown/F");
  fMiniTree->Branch("TLep2Jet2_DPhiJESDown",        &TLep2Jet2_DPhiJESDown,        "TLep2Jet2_DPhiJESDown/F");
  
  fMiniTree->Branch("TLep1Lep2METJet1_PtJESDown",   &Lep1Lep2METJet1_PtJESDown,    "TLep1Lep2METJet1_PtJESDown/F");
  fMiniTree->Branch("TMSysJESDown",                 &MSysJESDown,                  "TMSysJESDown/F");
  fMiniTree->Branch("TLep1Lep2Jet1_PtJESDown",      &Lep1Lep2Jet1_PtJESDown,       "TLep1Lep2Jet1_PtJESDown/F");
  fMiniTree->Branch("TLep1Lep2METJet1_PzJESDown",   &Lep1Lep2METJet1_PzJESDown,    "Lep1Lep2METJet1_PzJESDown/F");
  fMiniTree->Branch("TLep1Lep2METJet1_EtaJESDown",  &TLep1Lep2METJet1_EtaJESDown,  "TLep1Lep2METJet1_EtaJESDown/F");
  fMiniTree->Branch("TMiniMaxJESDown",              &TMiniMaxJESDown,              "TMiniMaxJESDown/F");
  fMiniTree->Branch("TLep1Lep2Jet1_EJESDown",       &TLep1Lep2Jet1_EJESDown,       "TLep1Lep2Jet1_EJESDown/F");
  fMiniTree->Branch("TLep1Lep2METJet1_MtJESDown",   &TLep1Lep2METJet1_MtJESDown,   "TLep1Lep2METJet1_MtJESDown/F");
  fMiniTree->Branch("TLep1Lep2Jet1_MJESDown",       &TLep1Lep2Jet1_MJESDown,       "TLep1Lep2Jet1_MJESDown/F");
  fMiniTree->Branch("THTtotJESDown",                &THTtotJESDown,                "THTtotJESDown/F");
  
  
  // JERUp
  fMiniTree->Branch("TLep1_PtJERUp",              &TLep1_PtJERUp,              "TLep1_PtJERUp/F");
  fMiniTree->Branch("TLep1_EJERUp",               &TLep1_EJERUp,               "TLep1_EJERUp/F");
  fMiniTree->Branch("TLep1_PhiJERUp",             &TLep1_PhiJERUp,             "TLep1_PhiJERUp/F");
  fMiniTree->Branch("TLep1_EtaJERUp",             &TLep1_EtaJERUp,             "TLep1_EtaJERUp/F");
  fMiniTree->Branch("TLep2_PtJERUp",              &TLep2_PtJERUp,              "TLep2_PtJERUp/F");
  fMiniTree->Branch("TLep2_EJERUp",               &TLep2_EJERUp,               "TLep2_EJERUp/F");
  fMiniTree->Branch("TLep2_PhiJERUp",             &TLep2_PhiJERUp,             "TLep2_PhiJERUp/F");
  fMiniTree->Branch("TLep2_EtaJERUp",             &TLep2_EtaJERUp,             "TLep2_EtaJERUp/F");
  fMiniTree->Branch("TJet1_PtJERUp",              &TJet1_PtJERUp,              "TJet1_PtJERUp/F");
  fMiniTree->Branch("TJet1_EJERUp",               &TJet1_EJERUp,               "TJet1_EJERUp/F");
  fMiniTree->Branch("TJet1_PhiJERUp",             &TJet1_PhiJERUp,             "TJet1_PhiJERUp/F");
  fMiniTree->Branch("TJet1_EtaJERUp",             &TJet1_EtaJERUp,             "TJet1_EtaJERUp/F");
  fMiniTree->Branch("TJet2_PtJERUp",              &TJet2_PtJERUp,              "TJet2_PtJERUp/F");
  fMiniTree->Branch("TJet2_EJERUp",               &TJet2_EJERUp,               "TJet2_EJERUp/F");
  fMiniTree->Branch("TJet2_PhiJERUp",             &TJet2_PhiJERUp,             "TJet2_PhiJERUp/F");
  fMiniTree->Branch("TJet2_EtaJERUp",             &TJet2_EtaJERUp,             "TJet2_EtaJERUp/F");
  fMiniTree->Branch("TLooseCentral1_PtJERUp",     &TLooseCentral1_PtJERUp,     "TLooseCentral1_PtJERUp/F");
  fMiniTree->Branch("TLooseFwd1_PtJERUp",         &TLooseFwd1_PtJERUp,         "TLooseFwd1_PtJERUp/F");
  fMiniTree->Branch("TMETJERUp",                  &TMETJERUp,                  "TMETJERUp/F");
  fMiniTree->Branch("TMET_PhiJERUp",              &TMET_PhiJERUp,              "TMET_PhiJERUp/F");
  
  fMiniTree->Branch("TLep1Lep2_PtJERUp",          &TLep1Lep2_PtJERUp,          "TLep1Lep2_PtJERUp/F");
  fMiniTree->Branch("TLep1Lep2_MJERUp",           &TLep1Lep2_MJERUp,           "TLep1Lep2_MJERUp/F");
  fMiniTree->Branch("TLep1Lep2_DPhiJERUp",        &TLep1Lep2_DPhiJERUp,        "TLep1Lep2_DPhiJERUp/F");
  fMiniTree->Branch("TLep1Jet1_MJERUp",           &TLep1Jet1_MJERUp,           "TLep1Jet1_MJERUp/F");
  fMiniTree->Branch("TLep1Jet1_DPhiJERUp",        &TLep1Jet1_DPhiJERUp,        "TLep1Jet1_DPhiJERUp/F");
  fMiniTree->Branch("TLep2Jet1_MJERUp",           &TLep2Jet1_MJERUp,           "TLep2Jet1_MJERUp/F");
  fMiniTree->Branch("TLep2Jet1_DPhiJERUp",        &TLep2Jet1_DPhiJERUp,        "TLep2Jet1_DPhiJERUp/F");
  fMiniTree->Branch("TLep1Jet2_MJERUp",           &TLep1Jet2_MJERUp,           "TLep1Jet2_MJERUp/F");
  fMiniTree->Branch("TLep1Jet2_DPhiJERUp",        &TLep1Jet2_DPhiJERUp,        "TLep1Jet2_DPhiJERUp/F");
  fMiniTree->Branch("TLep2Jet2_MJERUp",           &TLep2Jet2_MJERUp,           "TLep2Jet2_MJERUp/F");
  fMiniTree->Branch("TLep2Jet2_DPhiJERUp",        &TLep2Jet2_DPhiJERUp,        "TLep2Jet2_DPhiJERUp/F");
  
  fMiniTree->Branch("TLep1Lep2METJet1_PtJERUp",   &Lep1Lep2METJet1_PtJERUp,    "TLep1Lep2METJet1_PtJERUp/F");
  fMiniTree->Branch("TMSysJERUp",                 &MSysJERUp,                  "TMSysJERUp/F");
  fMiniTree->Branch("TLep1Lep2Jet1_PtJERUp",      &Lep1Lep2Jet1_PtJERUp,       "TLep1Lep2Jet1_PtJERUp/F");
  fMiniTree->Branch("TLep1Lep2METJet1_PzJERUp",   &Lep1Lep2METJet1_PzJERUp,    "Lep1Lep2METJet1_PzJERUp/F");
  fMiniTree->Branch("TLep1Lep2METJet1_EtaJERUp",  &TLep1Lep2METJet1_EtaJERUp,  "TLep1Lep2METJet1_EtaJERUp/F");
  fMiniTree->Branch("TMiniMaxJERUp",              &TMiniMaxJERUp,              "TMiniMaxJERUp/F");
  fMiniTree->Branch("TLep1Lep2Jet1_EJERUp",       &TLep1Lep2Jet1_EJERUp,       "TLep1Lep2Jet1_EJERUp/F");
  fMiniTree->Branch("TLep1Lep2METJet1_MtJERUp",   &TLep1Lep2METJet1_MtJERUp,   "TLep1Lep2METJet1_MtJERUp/F");
  fMiniTree->Branch("TLep1Lep2Jet1_MJERUp",       &TLep1Lep2Jet1_MJERUp,       "TLep1Lep2Jet1_MJERUp/F");
  fMiniTree->Branch("THTtotJERUp",                &THTtotJERUp,                "THTtotJERUp/F");
  
  
  // VAINAS DE GENERACION
  // PARTICLE LEVEL: DRESS
  // PARTON LEVEL: PART
  // SHARED: GEN
  fMiniTree->Branch("TPassDress",               &TPassDress,             "TPassDress/B");
  fMiniTree->Branch("TPassPart",                &TPassPart,              "TPassPart/B");
  fMiniTree->Branch("TGenChannel",              &GenChannel,             "GenChannel/I");
  fMiniTree->Branch("TGenIsSS",                 &TGenIsSS,               "TGenIsSS/B");
  
  // Particle level variables
  fMiniTree->Branch("TDressNJets",              &DressNJets,             "TDressNJets/I");
  fMiniTree->Branch("TDressNBJets",             &DressNBJets,            "TDressNBJets/I");
  fMiniTree->Branch("TDressNLooseCentral",      &DressNLooseCentral,     "TDressNLooseCentral/I");
  fMiniTree->Branch("TDressNBLooseCentral",     &DressNBLooseCentral,    "TDressNBLooseCentral/I");
  fMiniTree->Branch("TDressNLooseFwd",          &DressNLooseFwd,         "TDressNLooseFwd/I");
  
  fMiniTree->Branch("TDressLep1_Pt",            &TDressLep1_Pt,          "TDressLep1_Pt/F");
  fMiniTree->Branch("TDressLep1_E",             &TDressLep1_E,           "TDressLep1_E/F");
  fMiniTree->Branch("TDressLep1_Phi",           &TDressLep1_Phi,         "TDressLep1_Phi/F");
  fMiniTree->Branch("TDressLep1_Eta",           &TDressLep1_Eta,         "TDressLep1_Eta/F");
  fMiniTree->Branch("TDressLep2_Pt",            &TDressLep2_Pt,          "TDressLep2_Pt/F");
  fMiniTree->Branch("TDressLep2_E",             &TDressLep2_E,           "TDressLep2_E/F");
  fMiniTree->Branch("TDressLep2_Phi",           &TDressLep2_Phi,         "TDressLep2_Phi/F");
  fMiniTree->Branch("TDressLep2_Eta",           &TDressLep2_Eta,         "TDressLep2_Eta/F");
  fMiniTree->Branch("TDressJet1_Pt",            &TDressJet1_Pt,          "TDressJet1_Pt/F");
  fMiniTree->Branch("TDressJet1_E",             &TDressJet1_E,           "TDressJet1_E/F");
  fMiniTree->Branch("TDressJet1_Phi",           &TDressJet1_Phi,         "TDressJet1_Phi/F");
  fMiniTree->Branch("TDressJet1_Eta",           &TDressJet1_Eta,         "TDressJet1_Eta/F");
  fMiniTree->Branch("TDressJet2_Pt",            &TDressJet2_Pt,          "TDressJet2_Pt/F");
  fMiniTree->Branch("TDressJet2_E",             &TDressJet2_E,           "TDressJet2_E/F");
  fMiniTree->Branch("TDressJet2_Phi",           &TDressJet2_Phi,         "TDressJet2_Phi/F");
  fMiniTree->Branch("TDressJet2_Eta",           &TDressJet2_Eta,         "TDressJet2_Eta/F");
  fMiniTree->Branch("TDressLooseCentral1_Pt",   &TDressLooseCentral1_Pt, "TDressLooseCentral1_Pt/F");
  fMiniTree->Branch("TDressLooseCentral1_E",    &TDressLooseCentral1_E,  "TDressLooseCentral1_E/F");
  fMiniTree->Branch("TDressLooseCentral1_Phi",  &TDressLooseCentral1_Phi,"TDressLooseCentral1_Phi/F");
  fMiniTree->Branch("TDressLooseCentral1_Eta",  &TDressLooseCentral1_Eta,"TDressLooseCentral1_Eta/F");
  fMiniTree->Branch("TDressLooseFwd1_Pt",       &TDressLooseFwd1_Pt,     "TDressLooseFwd1_Pt/F");
  fMiniTree->Branch("TDressLooseFwd1_E",        &TDressLooseFwd1_E,      "TDressLooseFwd1_E/F");
  fMiniTree->Branch("TDressLooseFwd1_Phi",      &TDressLooseFwd1_Phi,    "TDressLooseFwd1_Phi/F");
  fMiniTree->Branch("TDressLooseFwd1_Eta",      &TDressLooseFwd1_Eta,    "TDressLooseFwd1_Eta/F");
  
  fMiniTree->Branch("TDressMET",                &TDressMET,              "TDressMET/F");
  fMiniTree->Branch("TDressMET_Phi",            &TDressMET_Phi,          "TDressMET_Phi/F");
  fMiniTree->Branch("TDressM_Lep1Jet1" ,        &TDressM_Lep1Jet1,       "TDressM_Lep1Jet1/F");
  fMiniTree->Branch("TDressM_Lep2Jet1",         &TDressM_Lep2Jet1,       "TDressM_Lep2Jet1/F");
  fMiniTree->Branch("TDressM_Lep1Jet2" ,        &TDressM_Lep1Jet2,       "TDressM_Lep1Jet2/F");
  fMiniTree->Branch("TDressM_Lep2Jet2",         &TDressM_Lep2Jet2,       "TDressM_Lep2Jet2/F");
  fMiniTree->Branch("TDressE_Lep1Lep2Jet1",     &TDressE_Lep1Lep2Jet1,   "TDressE_Lep1Lep2Jet1/F");
  fMiniTree->Branch("TDressMiniMax",            &TDressMiniMax,          "TDressMiniMax/F");
  fMiniTree->Branch("TDressLep1Lep2METJet1_Mt", &TDressLep1Lep2METJet1_Mt,"TDressLep1Lep2METJet1_Mt/F");
  fMiniTree->Branch("TDressLep1Lep2Jet1_M",     &TDressLep1Lep2Jet1_M,   "TDressLep1Lep2Jet1_M/F");
  fMiniTree->Branch("TDressLep1Lep2_Pt",        &TDressLep1Lep2_Pt,      "TDressLep1Lep2_Pt/F");
  fMiniTree->Branch("TDressLep1Lep2Jet1_Pt",    &TDressLep1Lep2Jet1_Pt,  "TDressLep1Lep2Jet1_Pt/F");
  fMiniTree->Branch("TDressLep1Lep2METJet1_Pt", &TDressLep1Lep2METJet1_Pt,"TDressLep1Lep2METJet1_Pt/F");
  fMiniTree->Branch("TDressHTtot",              &TDressHTtot,            "TDressHTtot/F");
  fMiniTree->Branch("TDressLep1Lep2METJet1_Pz", &TDressLep1Lep2METJet1_Pz"TDressLep1Lep2METJet1_Pz/F");
  fMiniTree->Branch("TDressLep1Lep2METJet1_Eta",&TDressLep1Lep2METJet1_Eta,"TDressLep1Lep2METJet1_Eta/F");
  fMiniTree->Branch("TDressMSys",               &TDressMSys,             "TDressMSys/F");
  fMiniTree->Branch("TDressLep1Lep2_M",         &TDressLep1Lep2_M,       "TDressLep1Lep2_M/F");
  fMiniTree->Branch("TDressLep1Lep2_DPhi",      &TDressLep1Lep2_DPhi,    "TDressLep1Lep2_DPhi/F");
  fMiniTree->Branch("TDressLep1Jet1_DPhi",      &TDressLep1Jet1_DPhi,    "TDressLep1Jet1_DPhi/F");
  fMiniTree->Branch("TDressLep2Jet1_DPhi",      &TDressLep2Jet1_DPhi,    "TDressLep2Jet1_DPhi/F");
  
  // Parton level variables
  fMiniTree->Branch("TPartMET",                &TPartMET,              "TPartMET/F");
  fMiniTree->Branch("TPartMET_Phi",            &TPartMET_Phi,          "TPartMET_Phi/F");
}


void TWTTbarAnalysis::ResetTWTTbarVariables() {
  hasTW = false;
  
  Lep1Lep2METJet1_Pt  = -99;
  Lep1METJetPt   = -99;
  DPtDilep_JetMET= -99;
  DPtDilep_MET   = -99;
  DPtLep1_MET    = -99;
  Lep1Lep2METJet1_Pz = -99;
  Lep1Lep2METJet1_PzJESUp = -99;
  Lep1Lep2METJet1_PzJESDown = -99;
  Lep1Lep2METJet1_PzJERUp = -99;
  NLooseCentral  = -99;
  NLooseCentralJESUp  = -99;
  NLooseCentralJESDown  = -99;
  NLooseCentralJERUp  = -99;
  NLooseFwd        = -99;
  NLooseFwdJESUp   = -99;
  NLooseFwdJESDown = -99;
  NLooseFwdJERUp   = -99;
  NBLooseCentral = -99;
  NBLooseCentralJESUp = -99;
  NBLooseCentralJESDown   = -99;
  NBLooseCentralJERUp = -99;
  TJet2_CSV       = -99;
  MSys           = -99;
  MSysJESUp      = -99;
  MSysJESDown    = -99;
  MSysJERUp      = -99;
  TLep1Lep2_M           = -99;
  TLep1Lep2_MJESUp      = -99;
  TLep1Lep2_MJESDown    = -99;
  TLep1Lep2_MJERUp      = -99;
  TLooseCentral1_Pt = -99;
  Lep1Lep2Jet1_Pt     = -99;
  
  TJet1_E    = -99;
  TJet1_Pt   = -99;
  TJet1_Phi  = -99;
  TJet1_Eta  = -99;
  TJet1_EJESUp    = -99;
  TJet1_PtJESUp   = -99;
  TJet1_PhiJESUp  = -99;
  TJet1_EtaJESUp  = -99;
  TJet1_EJESDown    = -99;
  TJet1_PtJESDown   = -99;
  TJet1_PhiJESDown  = -99;
  TJet1_EtaJESDown  = -99;
  TJet1_EJERUp    = -99;
  TJet1_PtJERUp   = -99;
  TJet1_PhiJERUp  = -99;
  TJet1_EtaJERUp  = -99;
  TLep1_E    = -99;
  TLep1_Pt   = -99;
  TLep1_Phi  = -99;
  TLep1_Eta  = -99;
  TLep1_EJESUp    = -99;
  TLep1_PtJESUp   = -99;
  TLep1_PhiJESUp  = -99;
  TLep1_EtaJESUp  = -99;
  TLep1_EJESDown    = -99;
  TLep1_PtJESDown   = -99;
  TLep1_PhiJESDown  = -99;
  TLep1_EtaJESDown  = -99;
  TLep1_EJERUp    = -99;
  TLep1_PtJERUp   = -99;
  TLep1_PhiJERUp  = -99;
  TLep1_EtaJERUp  = -99;
  TLep2_E    = -99;
  TLep2_Pt   = -99;
  TLep2_Phi  = -99;
  TLep2_Eta  = -99;
  TLep2_EJESUp    = -99;
  TLep2_PtJESUp   = -99;
  TLep2_PhiJESUp  = -99;
  TLep2_EtaJESUp  = -99;
  TLep2_EJESDown    = -99;
  TLep2_PtJESDown   = -99;
  TLep2_PhiJESDown  = -99;
  TLep2_EtaJESDown  = -99;
  TLep2_EJERUp    = -99;
  TLep2_PtJERUp   = -99;
  TLep2_PhiJERUp  = -99;
  TLep2_EtaJERUp  = -99;
  TLep1Lep2_Pt        = -99;
  TLep1Jet1_M     = -99;
  TLep2Jet1_M  = -99;
  TLep1Jet2_M        = -99;
  TLep2Jet2_M     = -99;
  TMiniMax = -99;
  TLep1Lep2Jet1_E          = -99;
  TLep1Lep2METJet1_Mt      = -99;
  TLep1Lep2Jet1_M          = -99;
  TLep1Lep2METJet1_Eta      = -99;
  TLep1Jet1_MJESUp     = -99;
  TLep2Jet1_MJESUp  = -99;
  TLep1Lep2Jet1_EJESUp          = -99;
  TLep1Lep2METJet1_MtJESUp      = -99;
  TLep1Lep2Jet1_MJESUp          = -99;
  TLep1Lep2_PtJESUp        = -99;
  TLep1Lep2METJet1_EtaJESUp      = -99;
  TLep1Jet1_MJESDown     = -99;
  TLep2Jet1_MJESDown  = -99;
  TLep1Lep2Jet1_EJESDown          = -99;
  TLep1Lep2METJet1_MtJESDown      = -99;
  TLep1Lep2Jet1_MJESDown          = -99;
  TLep1Lep2_PtJESDown        = -99;
  TLep1Lep2METJet1_EtaJESDown      = -99;
  TLep1Jet1_MJERUp     = -99;
  TLep2Jet1_MJERUp  = -99;
  TLep1Lep2Jet1_EJERUp          = -99;
  TLep1Lep2METJet1_MtJERUp      = -99;
  TLep1Lep2Jet1_MJERUp          = -99;
  TLep1Lep2_PtJERUp        = -99;
  TLep1Lep2METJet1_EtaJERUp      = -99;
  
  TLep1Lep2_DPhi              = -99.;
  TLep1Lep2_DPhiJESUp         = -99.;
  TLep1Lep2_DPhiJESDown       = -99.;
  TLep1Lep2_DPhiJERUp         = -99.;
  TLep1Jet1_DPhi         = -99.;
  TLep1Jet1_DPhiJESUp    = -99.;
  TLep1Jet1_DPhiJESDown  = -99.;
  TLep1Jet1_DPhiJERUp    = -99.;
  TLep2Jet1_DPhi      = -99.;
  TLep2Jet1_DPhiJESUp = -99.;
  TLep2Jet1_DPhiJESDown= -99.;
  TLep2Jet1_DPhiJERUp = -99.;
  
  TDressJet1_E    = -99;
  TDressJet1_Pt   = -99;
  TDressJet1_Phi  = -99;
  TDressJet1_Eta  = -99;
  TDressLooseFwd1_Pt   = -99.;
  TDressLooseFwd1_E    = -99.;
  TDressLooseFwd1_Phi  = -99.;
  TDressLooseFwd1_Eta  = -99.;
  TDressLooseCentral2_Pt  = -99.;
  TDressLooseCentral2_E   = -99.;
  TDressLooseCentral2_Phi = -99.;
  TDressLooseCentral2_Eta = -99.;
  TDressLep1_E    = -99;
  TDressLep1_Pt   = -99;
  TDressLep1_Phi  = -99;
  TDressLep1_Eta  = -99;
  TDressLep2_E    = -99;
  TDressLep2_Pt   = -99;
  TDressLep2_Phi  = -99;
  TDressLep2_Eta  = -99;
  TDressLep1Jet1_M     = -99;
  TDressLep2Jet1_M  = -99;
  TDressLep1Lep2Jet1_E          = -99;
  TDressLep1Jet2_M   = -99;
  TDressLep2Jet2_M= -99;
  TDressMiniMax=-99;
  TGenLep1Lep2METJet1_Mt      = -99;
  TGenLep1Lep2Jet1_M          = -99;
  TGenLep1Lep2Jet1_Pt     = -99;
  TGenLep1Lep2METJet1_Pt  = -99;
  TGenHTtot          = -99;
  TGenLep1Lep2METJet1_Pz = -99;
  TGenLep1Lep2METJet1_Eta      = -99;
  TGenMSys           = -99;
  TGenDPhiLL         = -99.;
  TGenDPhiLeadJet    = -99.;
  TGenDPhiSubLeadJet = -99.;
  TPassReco          = 0;
  TPassRecoJESUp     = 0;
  TPassRecoJESDown   = 0;
  TPassRecoJERUp     = 0;
  TPassGen           = 0;
  TGenIsSS           = 0;
  GenJets.clear();
  GenLooseCentralJets.clear();
  GeNLooseFwdJets.clear();
  GenLeps.clear();
  GenMET.SetXYZT(0, 0, 0, 0);
  tpL.SetXYZT(0, 0, 0, 0);
  tpJ.SetXYZT(0, 0, 0, 0);
  tMET.SetXYZT(0, 0, 0, 0);
  TIsFid          = false;
  GenChannel      = -1;
  TGenMET         = -99;
  TGenMET_Phi     = -99;
  nGenJets     = -99;
  nGenbJets    = -99;
  nGenLooseCentralJets= -99;
  nGenLooseCentralbJets= -99;
  nGeNLooseFwdJets= -99;
  nGenLeps     = -99;
  nGenJets  = -99;
  nGenLeps  = -99;
  nGenMET   = -99;
  TMET            = -99;
  TMET_Phi        = -99;
  TMETJESUp       = -99;
  TMET_PhiJESUp   = -99;
  TMETJESDown     = -99;
  TMET_PhiJESDown = -99;
  TMETJERUp       = -99;
  TMET_PhiJERUp   = -99;
  TWeight_normal  = -99;
}


void TWTTbarAnalysis::CalculateTWTTbarVariables() {
  if (hasTW) return;
  hasTW = true;

  get20Jets();
  TLorentzVector met;
  if (TNJets == 2 && TNBJets == 2) {
    met.SetPtEtaPhiE(TMET,0,TMET_Phi,TMET);
    
    Lep1Lep2METJet1_Pt     =  getLep1Lep2METJet1_Pt()   ;
    Lep1Lep2Jet1_Pt        =  getLep1Lep2Jet1_Pt()      ;
    Lep1METJetPt      =  getLep1METJetPt()    ;
    DPtDilep_JetMET   =  getDPtDilep_JetMET() ;
    DPtLep1_MET       =  getDPtLep1_MET()     ;
    Lep1Lep2METJet1_Pz    =  getLep1Lep2METJet1_Pz()  ;

    TJet1_pt          = selJets.at(0).Pt();
    MSys              = getSysM();
    THTtot            = getHTtot();
    DilepmetjetOverHT = Lep1Lep2METJet1_Pt/THTtot;
    TLep1Lep2_Pt          = (selLeptons.at(0).p + selLeptons.at(1).p).Pt();
    TLep1Lep2_M              = (selLeptons.at(0).p + selLeptons.at(1).p).M();
    
    TLep1Jet1_M       = (selJets.at(0).p + selLeptons.at(0).p).M();
    TLep2Jet1_M    = (selJets.at(0).p + selLeptons.at(1).p).M();
    TLep1Jet2_M     = (selJets.at(1).p + selLeptons.at(0).p).M();
    TLep2Jet2_M  = (selJets.at(1).p + selLeptons.at(1).p).M();
    if(GreaterThan(TLep1Jet1_M,TLep2Jet2_M)) {
        if(GreaterThan(TLep2Jet1_M,TLep1Jet2_M)) {
            if(!GreaterThan(TLep1Jet1_M,TLep2Jet1_M))       {TMiniMax = TLep1Jet1_M;}
            else                                               {TMiniMax = TLep2Jet1_M;}
        }
        else {
            if(!GreaterThan(TLep1Jet1_M,TLep1Jet2_M))       {TMiniMax = TLep1Jet1_M; }
            else                                               {TMiniMax = TLep1Jet2_M;}
        }
    }
    else {
        if(GreaterThan(TLep2Jet1_M,TLep1Jet2_M)) {
            if(!GreaterThan(TLep2Jet2_M,TLep2Jet1_M)) {TMiniMax = TLep2Jet2_M;}
            else                                               {TMiniMax = TLep2Jet1_M;}
        }
        else {
            if(!GreaterThan(TLep2Jet2_M,TLep1Jet2_M)) {TMiniMax = TLep2Jet2_M;}
            else                                               {TMiniMax = TLep1Jet2_M;}
        }
    }
    TLep1Lep2Jet1_E            = (selJets.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p).E();
    TLep1Lep2METJet1_Mt        = (selJets.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p + met).Mt();
    TLep1Lep2Jet1_M            = (selJets.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p).M();
    TLep1Lep2METJet1_Eta        = (selJets.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p + met).Eta();
    TLep1Lep2_DPhi           = getDPhiLep1_Lep2();
    TLep1Jet1_DPhi      = getDPhiLep1_Jet1();
    TLep2Jet1_DPhi   = getDPhiLep2_Jet1();
  }
  else {
    Lep1Lep2METJet1_Pt    = -99.;
    Lep1Lep2Jet1_Pt       = -99.;
    Lep1METJetPt     = -99.;
    DPtDilep_JetMET  = -99.;
    DPtLep1_MET      = -99.;
    Lep1Lep2METJet1_Pz   = -99.;
    C_jll            = -99.;
    TJet1_pt         = -99.;
    DilepmetjetOverHT= -99.;
    MSys             = -99.;
    THTtot           = -99.;
    TLep1Lep2_Pt         = -99.;
    TLep1Lep2_M             = -99.;
    TLep1Jet1_M      = -99.;
    TLep2Jet1_M   = -99.;
    TLep1Lep2Jet1_E           = -99.;
    TLep1Jet2_M    = -99.;
    TLep2Jet2_M = -99.;
    TMiniMax= -99.;
    TLep1Lep2METJet1_Mt       = -99.;
    TLep1Lep2Jet1_M           = -99.;
    TLep1Lep2METJet1_Eta       = -99.;
    TLep1Lep2_DPhi          = -99.;
    TLep1Jet1_DPhi     = -99.;
    TLep2Jet1_DPhi  = -99.;
  }

  if (!gIsData) {
    if (TNJetsJESUp == 1 && TNBJetsJESUp == 1) {
      // cout << "1j1b Up" << endl;
      met.SetPtEtaPhiE(TMETJESUp,0,TMET_PhiJESUp,TMETJESUp);
      
      Lep1Lep2METJet1_PtJESUp   =  getLep1Lep2METJet1_Pt(1)   ;
      Lep1Lep2Jet1_PtJESUp      =  getLep1Lep2Jet1_Pt("JESUp")      ;
      Lep1METJetPtJESUp    =  getLep1METJetPt("JESUp")    ;
      DPtDilep_JetMETJESUp =  getDPtDilep_JetMET("JESUp") ;
      DPtDilep_METJESUp    =  getDPtDilep_MET("JESUp")    ;
      DPtLep1_METJESUp     =  getDPtLep1_MET("JESUp")     ;
      Lep1Lep2METJet1_PzJESUp  =  getLep1Lep2METJet1_Pz("JESUp")  ;
      
      C_jllJESUp = (selJetsJecUp.at(0).p.Et() + selLeptons.at(0).p.Et() + selLeptons.at(1).p.Et()) / (selJetsJecUp.at(0).p.E() + selLeptons.at(0).p.E() + selLeptons.at(1).p.E());
      TJet1_ptJESUp           = selJetsJecUp.at(0).Pt();
      MSysJESUp               = getSysM("JESUp");
      THTtotJESUp             = getHTtot("JESUp");
      DilepmetjetOverHTJESUp  = Lep1Lep2METJet1_PtJESUp/THTtotJESUp          ;
      TLep1Lep2_PtJESUp           = (selLeptons.at(0).p + selLeptons.at(1).p).Pt();
      
      TLep1Jet1_MJESUp        = (selJetsJecUp.at(0).p + selLeptons.at(0).p).M();
      TLep2Jet1_MJESUp     = (selJetsJecUp.at(0).p + selLeptons.at(1).p).M();
      TLep1Lep2Jet1_EJESUp             = (selJetsJecUp.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p).E();
      TLep1Lep2METJet1_MtJESUp         = (selJetsJecUp.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p + met).Mt();
      TLep1Lep2Jet1_MJESUp             = (selJetsJecUp.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p).M();
      TLep1Lep2METJet1_EtaJESUp         = (selJetsJecUp.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p + met).Eta();
      TLep1Lep2_MJESUp               = (selLeptons.at(0).p + selLeptons.at(1).p).M();
      TLep1Lep2_DPhiJESUp            = getDPhiLep1_Lep2();
      TLep1Jet1_DPhiJESUp       = getDPhiLep1_Jet1("JESUp");
      TLep2Jet1_DPhiJESUp    = getDPhiLep2_Jet1("JESUp");
    }
    else {
      Lep1Lep2METJet1_PtJESUp    = -99.;
      Lep1Lep2Jet1_PtJESUp       = -99.;
      Lep1METJetPtJESUp     = -99.;
      DPtDilep_JetMETJESUp  = -99.;
      DPtLep1_METJESUp      = -99.;
      Lep1Lep2METJet1_PzJESUp   = -99.;
      C_jllJESUp            = -99.;
      TJet1_ptJESUp         = -99.;
      DilepmetjetOverHTJESUp= -99.;
      MSysJESUp             = -99.;
      TLep1Lep2_PtJESUp         = -99.;
      THTtotJESUp           = -99.;
      TLep1Jet1_MJESUp      = -99.;
      TLep2Jet1_MJESUp   = -99.;
      TLep1Lep2Jet1_EJESUp           = -99.;
      TLep1Lep2METJet1_MtJESUp       = -99.;
      TLep1Lep2Jet1_MJESUp           = -99.;
      TLep1Lep2METJet1_EtaJESUp       = -99.;
      TLep1Lep2_MJESUp             = -99.;
      TLep1Lep2_DPhiJESUp          = -99.;
      TLep1Jet1_DPhiJESUp     = -99.;
      TLep2Jet1_DPhiJESUp  = -99.;
    }
    
    if (TNJetsJESDown == 1 && TNBJetsJESDown == 1) {
      // cout << "1j1b Down"<< endl;
      met.SetPtEtaPhiE(TMETJESDown,0,TMET_PhiJESDown,TMETJESDown);
      
      Lep1Lep2METJet1_PtJESDown   =  getLep1Lep2METJet1_Pt(-1)   ;
      Lep1Lep2Jet1_PtJESDown      =  getLep1Lep2Jet1_Pt("JESDown")      ;
      Lep1METJetPtJESDown    =  getLep1METJetPt("JESDown")    ;
      DPtDilep_JetMETJESDown =  getDPtDilep_JetMET("JESDown") ;
      DPtLep1_METJESDown     =  getDPtLep1_MET("JESDown")     ;
      Lep1Lep2METJet1_PzJESDown  =  getLep1Lep2METJet1_Pz("JESDown")  ;
      
      C_jllJESDown = (selJetsJecDown.at(0).p.Et() + selLeptons.at(0).p.Et() + selLeptons.at(1).p.Et()) / (selJetsJecDown.at(0).p.E() + selLeptons.at(0).p.E() + selLeptons.at(1).p.E());
      TJet1_ptJESDown           = selJetsJecDown.at(0).p.Pt();
      THTtotJESDown             = getHTtot("JESDown");
      DilepmetjetOverHTJESDown  = Lep1Lep2METJet1_PtJESDown/THTtotJESDown          ;
      TLep1Lep2_PtJESDown           = (selLeptons.at(0).p + selLeptons.at(1).p).Pt();
      MSysJESDown               = getSysM("JESDown");
      
      TLep1Jet1_MJESDown        = (selJetsJecDown.at(0).p + selLeptons.at(0).p).M();
      TLep2Jet1_MJESDown     = (selJetsJecDown.at(0).p + selLeptons.at(1).p).M();
      TLep1Lep2Jet1_EJESDown             = (selJetsJecDown.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p).E();
      TLep1Lep2METJet1_MtJESDown         = (selJetsJecDown.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p + met).Mt();
      TLep1Lep2Jet1_MJESDown             = (selJetsJecDown.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p).M();
      TLep1Lep2METJet1_EtaJESDown         = (selJetsJecDown.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p + met).Eta();
      TLep1Lep2_MJESDown               = (selLeptons.at(0).p + selLeptons.at(1).p).M();
      TLep1Lep2_DPhiJESDown            = getDPhiLep1_Lep2();
      TLep1Jet1_DPhiJESDown       = getDPhiLep1_Jet1("JESDown");
      TLep2Jet1_DPhiJESDown    = getDPhiLep2_Jet1("JESDown");
    }
    else {
      Lep1Lep2METJet1_PtJESDown    = -99.;
      Lep1Lep2Jet1_PtJESDown       = -99.;
      Lep1METJetPtJESDown     = -99.;
      DPtDilep_JetMETJESDown  = -99.;
      DPtLep1_METJESDown      = -99.;
      Lep1Lep2METJet1_PzJESDown   = -99.;
      C_jllJESDown            = -99.;
      TJet1_ptJESDown         = -99.;
      DilepmetjetOverHTJESDown= -99.;
      MSysJESDown             = -99.;
      TLep1Lep2_PtJESDown         = -99.;
      THTtotJESDown           = -99.;
      TLep1Jet1_MJESDown      = -99.;
      TLep2Jet1_MJESDown   = -99.;
      TLep1Lep2Jet1_EJESDown           = -99.;
      TLep1Lep2METJet1_MtJESDown       = -99.;
      TLep1Lep2Jet1_MJESDown           = -99.;
      TLep1Lep2METJet1_EtaJESDown       = -99.;
      TLep1Lep2_MJESDown             = -99.;
      TLep1Lep2_DPhiJESDown          = -99.;
      TLep1Jet1_DPhiJESDown     = -99.;
      TLep2Jet1_DPhiJESDown  = -99.;
    }
    
    if (TNJetsJERUp == 1 && TNBJetsJERUp == 1) {
      met.SetPtEtaPhiE(TMETJERUp,0,TMET_PhiJERUp,TMETJERUp);
      
      Lep1Lep2METJet1_PtJERUp   =  getLep1Lep2METJet1_Pt(-2)   ; // -2 is JERUp :D
      Lep1Lep2Jet1_PtJERUp      =  getLep1Lep2Jet1_Pt("JER")      ;
      Lep1METJetPtJERUp    =  getLep1METJetPt("JER")    ;
      DPtDilep_JetMETJERUp =  getDPtDilep_JetMET("JER") ;
      DPtDilep_METJERUp    =  getDPtDilep_MET("JER")    ;
      DPtLep1_METJERUp     =  getDPtLep1_MET("JER")     ;
      Lep1Lep2METJet1_PzJERUp  =  getLep1Lep2METJet1_Pz("JER")  ;
      
      C_jllJERUp = (selJetsJER.at(0).p.Et() + selLeptons.at(0).p.Et() + selLeptons.at(1).p.Et()) / (selJetsJER.at(0).p.E() + selLeptons.at(0).p.E() + selLeptons.at(1).p.E());
      TJet1_ptJERUp           = selJetsJER.at(0).p.Pt();
      THTtotJERUp             = getHTtot("JER");
      TLep1Lep2_PtJERUp           = (selLeptons.at(0).p + selLeptons.at(1).p).Pt();
      DilepmetjetOverHTJERUp  = Lep1Lep2METJet1_PtJERUp/THTtotJERUp          ;
      MSysJERUp               = getSysM("JER");

      TLep1Jet1_MJERUp        = (selJetsJER.at(0).p + selLeptons.at(0).p).M();
      TLep2Jet1_MJERUp     = (selJetsJER.at(0).p + selLeptons.at(1).p).M();
      TLep1Lep2Jet1_EJERUp             = (selJetsJER.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p).E();
      TLep1Lep2METJet1_MtJERUp         = (selJetsJER.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p + met).Mt();
      TLep1Lep2Jet1_MJERUp             = (selJetsJER.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p).M();
      TLep1Lep2METJet1_EtaJERUp         = (selJetsJER.at(0).p + selLeptons.at(0).p + selLeptons.at(1).p + met).Eta();
      TLep1Lep2_MJERUp               = (selLeptons.at(0).p + selLeptons.at(1).p).M();
      TLep1Lep2_DPhiJERUp            = getDPhiLep1_Lep2();
      TLep1Jet1_DPhiJERUp       = getDPhiLep1_Jet1("JER");
      TLep2Jet1_DPhiJERUp    = getDPhiLep2_Jet1("JER");
    }
    else {
      Lep1Lep2METJet1_PtJERUp    = -99.;
      Lep1Lep2Jet1_PtJERUp       = -99.;
      Lep1METJetPtJERUp     = -99.;
      DPtDilep_JetMETJERUp  = -99.;
      DPtDilep_METJERUp     = -99.;
      DPtLep1_METJERUp      = -99.;
      Lep1Lep2METJet1_PzJERUp   = -99.;
      C_jllJERUp            = -99.;
      TJet1_ptJERUp         = -99.;
      DilepmetjetOverHTJERUp= -99.;
      MSysJERUp             = -99.;
      THTtotJERUp           = -99.;
      TLep1Lep2_PtJERUp         = -99.;
      TLep1Jet1_MJERUp      = -99.;
      TLep2Jet1_MJERUp   = -99.;
      TLep1Lep2Jet1_EJERUp           = -99.;
      TLep1Lep2METJet1_MtJERUp       = -99.;
      TLep1Lep2Jet1_MJERUp           = -99.;
      TLep1Lep2METJet1_EtaJERUp       = -99.;
      TLep1Lep2_MJERUp             = -99.;
      TLep1Lep2_DPhiJERUp          = -99.;
      TLep1Jet1_DPhiJERUp     = -99.;
      TLep2Jet1_DPhiJERUp  = -99.;
    }
  }

  if (TNJets == 2 && TNBJets == 1){
    TDilepMETPt     = getDilepMETPt();
    THtRejJ2        = selLeptons.at(0).Pt() + selLeptons.at(1).Pt() +  selJets.at(0).Pt() + TMET;
    TDPtL1_L2       = selLeptons.at(0).Pt() - selLeptons.at(1).Pt();
    TDPtJ2_L2       = selJets.at(1).Pt()    - selLeptons.at(1).Pt();
    TDR_L1_J1       = getDeltaRLep1_Jet1();
    TDR_L1L2_J1J2   = getDeltaRDilep_Jets12();
    TDR_L1L2_J1J2MET= getDeltaRDilep_METJets12();
    LeadingLeptPt_    = selLeptons.at(0).Pt();
    LeadingLeptEta_   = selLeptons.at(0).Eta();
    jetPtSubLeading_  = selJets.at(1).Pt();
    jetEtaSubLeading_ = selJets.at(1).Eta();
  }
  else{
    TDilepMETPt      = -99;
    TETSys           = -99.;
    TET_LLJetMET     = -99.;
    TDPtL1_L2        = -99.;
    TDPtJ2_L2        = -99.;
    TDR_L1_J1        = -99.;
    TDR_L1L2_J1J2    = -99.;
    TDR_L1L2_J1J2MET = -99.;
    THtRejJ2         = -99.;
    THTtot2j         = -99.;
  }

  if (!gIsData){
    if (TNJetsJESUp == 2 && TNBJetsJESUp == 1){
      TDilepMETPtJESUp       = getDilepMETPt(+1);
      THtRejJ2JESUp          = selLeptons.at(0).Pt() + selLeptons.at(1).Pt() +  selJetsJecUp.at(0).Pt() + TMETJESUp;
      TDPtL1_L2JESUp         = selLeptons.at(0).Pt() - selLeptons.at(1).Pt();
      TDPtJ2_L2JESUp         = selJetsJecUp.at(1).Pt()    - selLeptons.at(1).Pt();
      TDR_L1_J1JESUp         = getDeltaRLep1_Jet1(1);
      TDR_L1L2_J1J2JESUp     = getDeltaRDilep_Jets12(1);
      TDR_L1L2_J1J2METJESUp  = getDeltaRDilep_METJets12(1);
      LeadingLeptPt_JESUp    = selLeptons.at(0).Pt();
      LeadingLeptEta_JESUp   = selLeptons.at(0).Eta();
      jetPtSubLeading_JESUp  = selJetsJecUp.at(1).Pt();
      jetEtaSubLeading_JESUp = selJetsJecUp.at(1).Eta();
    }
    else{
      TDilepMETPtJESUp       = -99.;
      TETSysJESUp            = -99.;
      TET_LLJetMETJESUp      = -99.;
      THtRejJ2JESUp          = -99.;
      TDPtL1_L2JESUp         = -99.;
      TDPtJ2_L2JESUp         = -99.;
      TDR_L1_J1JESUp         = -99.;
      TDR_L1L2_J1J2JESUp     = -99.;
      TDR_L1L2_J1J2METJESUp  = -99.;
      LeadingLeptPt_JESUp    = -99.;
      LeadingLeptEta_JESUp   = -99.;
      jetPtSubLeading_JESUp  = -99.;
      jetEtaSubLeading_JESUp = -99.;
      THTtot2jJESUp          = -99.;
    }
    
    if (TNJetsJESDown == 2 && TNBJetsJESDown == 1){
      TDilepMETPtJESDown       = getDilepMETPt(-1);
      THtRejJ2JESDown          = selLeptons.at(0).Pt() + selLeptons.at(1).Pt() +  selJetsJecDown.at(0).Pt() + TMETJESDown;
      TDPtL1_L2JESDown         = selLeptons.at(0).Pt() - selLeptons.at(1).Pt();
      TDPtJ2_L2JESDown         = selJetsJecDown.at(1).Pt()    - selLeptons.at(1).Pt();
      TDR_L1_J1JESDown         = getDeltaRLep1_Jet1(-1);
      TDR_L1L2_J1J2JESDown     = getDeltaRDilep_Jets12(-1);
      TDR_L1L2_J1J2METJESDown  = getDeltaRDilep_METJets12(-1);
      LeadingLeptPt_JESDown    = selLeptons.at(0).Pt();
      LeadingLeptEta_JESDown   = selLeptons.at(0).Eta();
      jetPtSubLeading_JESDown  = selJetsJecDown.at(1).Pt();
      jetEtaSubLeading_JESDown = selJetsJecDown.at(1).Eta();
    }
    else{
      TDilepMETPtJESDown       = -99.;
      TETSysJESDown            = -99.;
      TET_LLJetMETJESDown      = -99.;
      THtRejJ2JESDown          = -99.;
      TDPtL1_L2JESDown         = -99.;
      TDPtJ2_L2JESDown         = -99.;
      TDR_L1_J1JESDown         = -99.;
      TDR_L1L2_J1J2JESDown     = -99.;
      TDR_L1L2_J1J2METJESDown  = -99.;
      LeadingLeptPt_JESDown    = -99.;
      LeadingLeptEta_JESDown   = -99.;
      jetPtSubLeading_JESDown  = -99.;
      jetEtaSubLeading_JESDown = -99.;
      THTtot2jJESDown          = -99.;
    }

    if (TNJetsJERUp == 2 && TNBJetsJERUp == 1) {
      TDilepMETPtJERUp       = getDilepMETPt(-2);
      TETSysJERUp            = getETSys(-2);
      THtRejJ2JERUp          = selLeptons.at(0).Pt() + selLeptons.at(1).Pt() +  selJetsJER.at(0).Pt() + TMETJERUp;
      TDPtL1_L2JERUp         = selLeptons.at(0).Pt() - selLeptons.at(1).Pt();
      TDPtJ2_L2JERUp         = selJetsJER.at(1).Pt() - selLeptons.at(1).Pt();
      TDR_L1_J1JERUp         = getDeltaRLep1_Jet1(-2);
      TDR_L1L2_J1J2JERUp     = getDeltaRDilep_Jets12(-2);
      TDR_L1L2_J1J2METJERUp  = getDeltaRDilep_METJets12(-2);
      LeadingLeptPt_JERUp    = selLeptons.at(0).Pt();
      LeadingLeptEta_JERUp   = selLeptons.at(0).Eta();
      jetPtSubLeading_JERUp  = selJetsJER.at(1).Pt();
      jetEtaSubLeading_JERUp = selJetsJER.at(1).Eta();
    }
    else{
      TDilepMETPtJERUp       = -99.;
      TETSysJERUp            = -99.;
      TET_LLJetMETJERUp      = -99.;
      THtRejJ2JERUp          = -99.;
      TDPtL1_L2JERUp         = -99.;
      TDPtJ2_L2JERUp         = -99.;
      TDR_L1_J1JERUp         = -99.;
      TDR_L1L2_J1J2JERUp     = -99.;
      TDR_L1L2_J1J2METJERUp  = -99.;
      LeadingLeptPt_JERUp    = -99.;
      LeadingLeptEta_JERUp   = -99.;
      jetPtSubLeading_JERUp  = -99.;
      jetEtaSubLeading_JERUp = -99.;
      THTtot2jJERUp          = -99.;
    }
  }

  if (TNJets > 1)        TJet2_Pt        = selJets.at(1).Pt();
  else                   TJet2_Pt        = -99.;

  if (!gIsData){
    if (TNJetsJESUp > 1)   TJet2_PtJESUp   = selJetsJecUp.at(1).Pt();
    else                   TJet2_PtJESUp   = -99.;
    if (TNJetsJESDown > 1) TJet2_PtJESDown = selJetsJecDown.at(1).Pt();
    else                   TJet2_PtJESDown = -99.;
    if (TNJetsJERUp > 1)   TJet2_PtJERUp   = selJetsJER.at(1).Pt();
    else                   TJet2_PtJERUp   = -99.;
  }
  return;
}


void TWTTbarAnalysis::get20Jets() {
  vector<Float_t> looseJetPt;
  vector<Float_t> looseJetPtJESUp;
  vector<Float_t> looseJetPtJESDown;
  vector<Float_t> looseJetPtJER;

  vector<Float_t> looseJetCentralPt;
  vector<Float_t> looseJetCentralPtJESUp;
  vector<Float_t> looseJetCentralPtJESDown;
  vector<Float_t> looseJetCentralPtJER;
  
  vector<Float_t> looseJetFwdPt;
  vector<Float_t> looseJetFwdPtJESUp;
  vector<Float_t> looseJetFwdPtJESDown;
  vector<Float_t> looseJetFwdPtJER;
  
  nLooseCentral  = 0.; NLooseCentralJESUp  = 0.; NLooseCentralJESDown  = 0.; NLooseCentralJERUp  = 0.;
  NLooseFwd      = 0.; NLooseFwdJESUp      = 0.; NLooseFwdJESDown      = 0.; NLooseFwdJERUp      = 0.;
  NBLooseCentral = 0.; NBLooseCentralJESUp = 0.; NBLooseCentralJESDown = 0.; NBLooseCentralJERUp = 0.;
  TJet2_CSV       = 0.; TJet2_CSVJESUp       = 0.; TJet2_CSVJESDown       = 0.; TJet2_CSVJERUp       = 0.;
  TJetLoosept    = 0.; TJetLooseptJESUp    = 0.; TJetLooseptJESDown    = 0.; TJetLooseptJERUp    = 0.;
  TLooseCentral1_Pt        = 0.; TLooseCentral1_PtJESUp = 0.;
  TLooseCentral1_PtJESDown = 0.; TLooseCentral1_PtJERUp = 0.;
  TJetLooseFwdpt        = 0.; TJetLooseFwdptJESUp = 0.;
  TJetLooseFwdptJESDown = 0.; TJetLooseFwdptJERUp = 0.;

  for (unsigned int j = 0; j < vetoJets.size(); ++j){
    if (vetoJets.at(j).p.Pt() > 20.){
      looseJetPt.push_back( vetoJets.at(j).p.Pt() );
      if (TMath::Abs(vetoJets.at(j).p.Eta()) < 2.4){
        nLooseCentral++;
        looseJetCentralPt.push_back(vetoJets.at(j).p.Pt());
      }
      else {
        NLooseFwd++;
        looseJetFwdPt.push_back(vetoJets.at(j).p.Pt());
      }
      if (vetoJets.at(j).isBtag){
    if (TMath::Abs(vetoJets.at(j).p.Eta()) < 2.4) NBLooseCentral++;
      }
    }

    if (vetoJets.at(j).pTJESUp > 20.){
      looseJetPtJESUp.push_back( vetoJets.at(j).pTJESUp );
      if (TMath::Abs(vetoJets.at(j).p.Eta()) < 2.4){
        NLooseCentralJESUp++;
        looseJetCentralPtJESUp.push_back(vetoJets.at(j).pTJESUp);
      }
      else {
        NLooseFwdJESUp++;
        looseJetFwdPtJESUp.push_back(vetoJets.at(j).pTJESUp);
      }
      if (vetoJets.at(j).isBtag){
        if (TMath::Abs(vetoJets.at(j).p.Eta()) < 2.4) NBLooseCentralJESUp++;
      }
    }

    if (vetoJets.at(j).pTJESDown > 20.){
      looseJetPtJESDown.push_back( vetoJets.at(j).pTJESDown );
      if (TMath::Abs(vetoJets.at(j).p.Eta()) < 2.4){
        NLooseCentralJESDown++;
        looseJetCentralPtJESDown.push_back(vetoJets.at(j).pTJESDown);
      }
      else {
        NLooseFwdJESDown++;
        looseJetFwdPtJESDown.push_back(vetoJets.at(j).pTJESDown);
      }
      if (vetoJets.at(j).isBtag){
        if (TMath::Abs(vetoJets.at(j).p.Eta()) < 2.4) NBLooseCentralJESDown++;
      }
    }

    if (vetoJets.at(j).pTJERUp > 20.){
      looseJetPtJER.push_back( vetoJets.at(j).pTJERUp );
      if (TMath::Abs(vetoJets.at(j).p.Eta()) < 2.4){
        NLooseCentralJERUp++;
        looseJetCentralPtJER.push_back(vetoJets.at(j).pTJERUp);
      }
      else {
        NLooseFwdJERUp++;
        looseJetFwdPtJER.push_back(vetoJets.at(j).pTJERUp);
      }
      if (vetoJets.at(j).isBtag){
        if (TMath::Abs(vetoJets.at(j).p.Eta()) < 2.4) NBLooseCentralJERUp++;
      }
    }
  }

  std::sort( looseJetPt.begin()       , looseJetPt.end()       , GreaterThan);
  std::sort( looseJetPtJESUp.begin()  , looseJetPtJESUp.end()  , GreaterThan);
  std::sort( looseJetPtJESDown.begin(), looseJetPtJESDown.end(), GreaterThan);
  std::sort( looseJetPtJER.begin()    , looseJetPtJER.end()    , GreaterThan);

  std::sort( looseJetCentralPt.begin()       , looseJetCentralPt.end()       , GreaterThan);
  std::sort( looseJetCentralPtJESUp.begin()  , looseJetCentralPtJESUp.end()  , GreaterThan);
  std::sort( looseJetCentralPtJESDown.begin(), looseJetCentralPtJESDown.end(), GreaterThan);
  std::sort( looseJetCentralPtJER.begin()    , looseJetCentralPtJER.end()    , GreaterThan);

  if (NLooseFwd + nLooseCentral > 1){
    TJet2_CSV = TMath::Max( vetoJets.at(1).csv   , (Float_t) 0.);
    TJetLoosept = looseJetPt.at(1);
  }
  if (NLooseFwdJESUp + NLooseCentralJESUp > 1){
    TJet2_CSVJESUp = TMath::Max( vetoJets.at(1).csv   , (Float_t) 0.);
    TJetLooseptJESUp = looseJetPtJESUp.at(1);
  }
  if (NLooseFwdJESDown + NLooseCentralJESDown > 1){
    TJet2_CSVJESDown = TMath::Max( vetoJets.at(1).csv   , (Float_t) 0.);
    TJetLooseptJESDown = looseJetPtJESDown.at(1);
  }
  if (NLooseFwdJERUp + NLooseCentralJERUp > 1){
    TJet2_CSVJERUp = TMath::Max( vetoJets.at(1).csv   , (Float_t) 0.);
    TJetLooseptJERUp = looseJetPtJER.at(1);
  }

  if (nLooseCentral > 1)
    TLooseCentral1_Pt = looseJetCentralPt.at(1);

  if (NLooseCentralJESUp > 1)
    TLooseCentral1_PtJESUp = looseJetCentralPtJESUp.at(1);

  if (NLooseCentralJESDown > 1)
    TLooseCentral1_PtJESDown = looseJetCentralPtJESDown.at(1);

  if (NLooseCentralJERUp > 1)
    TLooseCentral1_PtJERUp = looseJetCentralPtJER.at(1);

  if (NLooseFwd > 1)
    TJetLooseFwdpt = looseJetFwdPt.at(1);

  if (NLooseFwdJESUp > 1)
    TJetLooseFwdptJESUp = looseJetFwdPtJESUp.at(1);

  if (NLooseFwdJESDown > 1)
    TJetLooseFwdptJESDown = looseJetFwdPtJESDown.at(1);

  if (NLooseFwdJERUp > 1)
    TJetLooseFwdptJERUp = looseJetFwdPtJER.at(1);

  return;
}


Double_t TWTTbarAnalysis::getHTtot(const TString& sys) {
  if (sys == "Norm")
    return selLeptons.at(0).p.Pt() + selLeptons.at(1).p.Pt() + TMET + selJets.at(0).p.Pt();
  else if (sys == "JESUp")
    return selLeptons.at(0).p.Pt() + selLeptons.at(1).p.Pt() + TMETJESUp + selJetsJecUp.at(0).p.Pt();
  else if (sys == "JESDown")
    return selLeptons.at(0).p.Pt() + selLeptons.at(1).p.Pt() + TMETJESDown + selJetsJecDown.at(0).p.Pt();
  else if (sys == "JER")
    return selLeptons.at(0).p.Pt() + selLeptons.at(1).p.Pt() + TMETJERUp     + selJetsJER.at(0).p.Pt();
  else{
    cout << "Unset systematic" << endl;
    return -99.;
  }
}


Double_t TWTTbarAnalysis::getLep1Lep2METJet1_Pt(int sys) {
  TLorentzVector vSystem[4];
  vSystem[0] = selLeptons.at(0).p;                               // lep1
  vSystem[1] = selLeptons.at(1).p;                               // lep2
  if (sys == 0){
    vSystem[2].SetPtEtaPhiE(TMET,0,TMET_Phi,TMET); // met
    vSystem[3] = selJets.at(0).p;                                    // jet 1
  }
  else if (sys > 0){
    vSystem[2].SetPtEtaPhiE(TMETJESUp,0,TMET_PhiJESUp,TMETJESUp); // met
    vSystem[3] = selJetsJecUp.at(0).p;
  }
  else if (sys == -1){
    vSystem[2].SetPtEtaPhiE(TMETJESDown,0,TMET_PhiJESDown,TMETJESDown); // met
    vSystem[3] = selJetsJecDown.at(0).p;
  }
  else if (sys == -2){
    vSystem[2].SetPtEtaPhiE(TMETJERUp,0,TMET_Phi,TMETJERUp); // met
    vSystem[3] = selJetsJER.at(0).p;
  }
  else{
    cout << "Wrong systematic in TWTTbarAnalysis::getLep1Lep2METJet1_Pt" << endl;
    return -9999.;
  }
  return getPtSys(vSystem,4);
}


Double_t TWTTbarAnalysis::getDilepMETPt(int sys) {
  TLorentzVector vSystem[3];
  vSystem[0] = selLeptons.at(0).p;                               // lep1
  vSystem[1] = selLeptons.at(1).p;                               // lep2
  if (sys == 0){
    vSystem[2].SetPtEtaPhiE(TMET,0,TMET_Phi,TMET); // met
  }
  else if (sys > 0){
    vSystem[2].SetPtEtaPhiE(TMETJESUp,0,TMET_PhiJESUp,TMETJESUp); // met
  }
  else if (sys == -1){
    vSystem[2].SetPtEtaPhiE(TMETJESDown,0,TMET_PhiJESDown,TMETJESDown); // met
  }
  else if (sys == -2){
    vSystem[2].SetPtEtaPhiE(TMETJERUp,0,TMET_PhiJERUp,TMETJERUp); // met
  }
  else{
    cout << "Wrong systematic in TWTTbarAnalysis::getLep1Lep2METJet1_Pt" << endl;
    return -9999.;
  }

  return getPtSys(vSystem,3);
}


Double_t TWTTbarAnalysis::getLep1Lep2Jet1_Pt(const TString& sys) {
  TLorentzVector vSystem[3];
  vSystem[0] = selLeptons.at(0).p;                               // lep1
  vSystem[1] = selLeptons.at(1).p;                               // lep2
  if (sys == "Norm")
    vSystem[2] = selJets.at(0).p;                                    // jet 1
  else if (sys == "JESUp")
    vSystem[2] = selJetsJecUp.at(0).p;                                    // jet 1
  else if (sys == "JESDown")
    vSystem[2] = selJetsJecDown.at(0).p;                                    // jet 1
  else if (sys == "JER")
    vSystem[2] = selJetsJER.at(0).p;                                    // jet 1
  else{
    cout << "Wrong systematic in TWTTbarAnalysis::getLep1Lep2Jet1_Pt" << endl;
    return -9999.;
  }

  return getPtSys(vSystem,3);
}


Double_t TWTTbarAnalysis::getLep1METJetPt(const TString& sys) {
  TLorentzVector vSystem[3];
  vSystem[0] = selLeptons.at(0).p;
  if (sys == "Norm"){
    vSystem[1].SetPtEtaPhiE(TMET,0,TMET_Phi,TMET);
    vSystem[2] = selJets.at(0).p;
  }
  else if (sys == "JESUp"){
    vSystem[1].SetPtEtaPhiE(TMETJESUp,0,TMET_PhiJESUp,TMETJESUp); // met
    vSystem[2] = selJetsJecUp.at(0).p;
  }
  else if (sys == "JESDown"){
    vSystem[1].SetPtEtaPhiE(TMETJESDown,0,TMET_PhiJESDown,TMETJESDown); // met
    vSystem[2] = selJetsJecDown.at(0).p;
  }  
  else if (sys == "JER"){
    vSystem[1].SetPtEtaPhiE(TMETJERUp,0,TMET_PhiJERUp,TMETJERUp); // met
    vSystem[2] = selJetsJER.at(0).p;
  }  
  else{
    cout << "Wrong systematic in TWTTbarAnalysis::getLep1METJetPt" << endl;
    return -9999.;
  }
  
  return getPtSys(vSystem,3);
}


Double_t TWTTbarAnalysis::getPtSys(TLorentzVector* vSystem, int nSystem) {
  TLorentzVector vSyst;
  vSyst = vSystem[0];
  for (int i = 1; i < nSystem; ++i){
    vSyst += vSystem[i];
  }
  return vSyst.Pt();
}


Double_t TWTTbarAnalysis::getLep1Lep2METJet1_Pz(const TString& sys) {
  TLorentzVector vSystem[4];
  vSystem[0] = selLeptons.at(0).p;                               // lep1
  vSystem[1] = selLeptons.at(1).p;                               // lep2
  if (sys == "Norm"){
    vSystem[2].SetPtEtaPhiE(TMET,0,TMET_Phi,TMET); // met
    vSystem[3] = selJets.at(0).p;                                    // jet 1
  }
  else if (sys == "JESUp"){
    vSystem[2].SetPtEtaPhiE(TMETJESUp,0,TMET_PhiJESUp,TMETJESUp); // met
    vSystem[3] = selJetsJecUp.at(0).p;                                    // jet 1
  }
  else if (sys == "JESDown"){
    vSystem[2].SetPtEtaPhiE(TMETJESDown,0,TMET_PhiJESDown,TMETJESDown); // met
    vSystem[3] = selJetsJecDown.at(0).p;                                    // jet 1
  }
  else if (sys == "JER"){
    vSystem[2].SetPtEtaPhiE(TMETJERUp,0,TMET_PhiJERUp,TMETJERUp); // met
    vSystem[3] = selJetsJER.at(0).p;                                    // jet 1
  }
  else{
    cout << "Wrong systematic in TWTTbarAnalysis::getLep1Lep2METJet1_Pz" << endl;
    return -9999.;
  }

  return getPzSys(vSystem,4);
}


Double_t TWTTbarAnalysis::getPzSys(TLorentzVector* vSystem, int nSystem) {
 TLorentzVector vSyst;
  vSyst = vSystem[0];
  for (int i = 1; i < nSystem; ++i){
    vSyst += vSystem[i];
  }
  return vSyst.Pz();
}


Double_t TWTTbarAnalysis::getDeltaPt(vector<TLorentzVector> col1, vector<TLorentzVector> col2){
  TLorentzVector v1, v2;
  if (col1.size()*col2.size() == 0){
    cout << "[TAT::getDeltaPt] One of the collections is empty" << endl;
    return -99.;
  }
  v1 = col1[0];
  v2 = col2[0];
  for (unsigned int i = 1; i < col1.size(); ++i){
    v1 += col1[i];
  }
  for (unsigned int i = 1; i < col2.size(); ++i){
    v2 += col2[i];
  }
  return (v1-v2).Pt();
}


Double_t TWTTbarAnalysis::getDPhiLep1_Lep2()
{
  vector<TLorentzVector>  col1;
  vector<TLorentzVector>  col2;
  col1.push_back(selLeptons.at(0).p);
  col2.push_back(selLeptons.at(1).p);
  return getDeltaPhi(col1,col2);
}


Double_t TWTTbarAnalysis::getDPhiLep1_Jet1(const TString& sys)
{
  vector<TLorentzVector>  col1;
  vector<TLorentzVector>  col2;
  col1.push_back(selLeptons.at(0).p);
  if (sys == "Norm")
    col2.push_back(selJets.at(0).p);
  else if (sys == "JESUp")
    col2.push_back(selJetsJecUp.at(0).p);
  else if (sys == "JESDown")
    col2.push_back(selJetsJecDown.at(0).p);
  else if (sys == "JER")
    col2.push_back(selJetsJER.at(0).p);
  else{
    cout << "Wrong systematic in TWTTbarAnalysis::getDPhiLep1_Jet1" << endl;
    return -9999.;
  }
  return getDeltaPhi(col1,col2);
}


Double_t TWTTbarAnalysis::getDPhiLep2_Jet1(const TString& sys) {
  vector<TLorentzVector>  col1;
  vector<TLorentzVector>  col2;
  col1.push_back(selLeptons.at(1).p);
  if (sys == "Norm")
    col2.push_back(selJets.at(0).p);
  else if (sys == "JESUp")
    col2.push_back(selJetsJecUp.at(0).p);
  else if (sys == "JESDown")
    col2.push_back(selJetsJecDown.at(0).p);
  else if (sys == "JER")
    col2.push_back(selJetsJER.at(0).p);
  else{
    cout << "Wrong systematic in TWTTbarAnalysis::getDPhiLep1_Jet1" << endl;
    return -9999.;
  }
  return getDeltaPhi(col1,col2);
}


Double_t TWTTbarAnalysis::getDeltaPhi(vector<TLorentzVector> col1, vector<TLorentzVector> col2) {
  TLorentzVector v1, v2;
  if (col1.size()*col2.size() == 0){
    cout << "[TAT::getDeltaPt] One of the collections is empty" << endl;
    return -99.;
  }
  v1 = col1[0];
  v2 = col2[0];
  for (unsigned int i = 1; i < col1.size(); ++i){
    v1 += col1[i];
  }
  for (unsigned int i = 1; i < col2.size(); ++i){
    v2 += col2[i];
  }
  return v1.DeltaPhi(v2);
}


Double_t TWTTbarAnalysis::getSysM(const TString& sys) {
  vector<TLorentzVector> col;
  TLorentzVector met;
  if (sys == "Norm"){
    met.SetPtEtaPhiE(TMET,0,TMET_Phi,TMET);
    col.push_back(selJets.at(0).p);
  }
  else if (sys == "JESUp"){
    met.SetPtEtaPhiE(TMETJESUp,0,TMET_PhiJESUp,TMETJESUp);
    col.push_back(selJetsJecUp.at(0).p);
  }
  else if (sys == "JESDown"){
    met.SetPtEtaPhiE(TMETJESDown,0,TMET_PhiJESDown,TMETJESDown);
    col.push_back(selJetsJecDown.at(0).p);
  }
  else if (sys == "JER"){
    met.SetPtEtaPhiE(TMETJERUp,0,TMET_PhiJERUp,TMETJERUp);
    col.push_back(selJetsJER.at(0).p);
  }
  else{
    return -9999.;
  }

  col.push_back(met);
  col.push_back(selLeptons.at(0).p);
  col.push_back(selLeptons.at(1).p);
  return getM(col);
}


Double_t TWTTbarAnalysis::getM(vector<TLorentzVector> col) {
  if (col.size() == 0){
    return -99.;
  }
  TLorentzVector v;
  v  = col[0];
  for (unsigned int i = 1; i < col.size(); ++i){
    v += col[i];
  }
  return v.M();
}

Double_t TWTTbarAnalysis::getDilepPt() {
  TLorentzVector vSystem[2];
  vSystem[0] = selLeptons.at(0).p;
  vSystem[1] = selLeptons.at(1).p;
  return getPtSys(vSystem,2);
}

Double_t TWTTbarAnalysis::getDeltaR(vector<TLorentzVector> col1, vector<TLorentzVector> col2) {
  TLorentzVector v1, v2;
  if (col1.size()*col2.size() == 0){
    cout << "[TAT::getDeltaR] One of the collections is empty" << endl;
    return -99.;
  }
  v1 = col1[0];
  v2 = col2[0];
  for (unsigned int i = 1; i < col1.size(); ++i){
    v1 += col1[i];
  }
  for (unsigned int i = 1; i < col2.size(); ++i){
    v2 += col2[i];
  }
  return v1.DeltaR(v2);
}
