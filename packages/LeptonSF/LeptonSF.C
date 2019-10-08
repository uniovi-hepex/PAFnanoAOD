#include "LeptonSF.h"

#include "PAF.h"

#include <iostream>

LeptonSF::LeptonSF(TString path, Int_t year, TString options):
  fMuonTrackerSF(0),  // Muon Reco
  fMuonIdSF(0), // Muon Id
  fMuonIdSF_BCDEF(0), // Muon Id
  fMuonIdSF_GH(0),    // Muon Id
  fMuonIdSFSUSY(0),   // Muon Id
  fMuonIsoSF(0),// Muon Iso
  fMuonIsoSF_BCDEF(0),// Muon Iso
  fMuonIsoSF_GH(0),   // Muon Iso
  fMuonIsoSFSUSY(0),  // Muon Iso
  fMuonIP2DSF(0),     // Muon IP2d
  fMuonSIP3DSF(0),    // Muon SIP
  fMuonFastSim(0),
  fMuonIsoFastSimStop(0),
  fMuonIdFastSimStop(0),
  fMuonlepMVA2lSSttH(0),    // Muon ttH LepMVA 2lSS
  fMuonlepMVA3l4lttH(0),    // Muon ttH LepMVA 3l, 4l
  fMuonLooseTracksttH(0),   // Muon ttH Loose Tracks
  fMuonLooseMiniIsottH(0),  // Muon ttH Loose MiniIso
  fMuonTightIP2DttH(0),     // Muon ttH Tight IP2D


  fElecTrackerSF(0),  // Electron Reco
  fElecTrackerSFF(0),  // Electron Reco 2017
  fElecIdSF(0),       // Electron Id (+Iso)
  fElecIsoSF(0),      // Electron Iso
  fElecIP2DSF(0),     // Electron IP2D
  fElecSIP3DSF(0),    // Electron SIP3D
  fElecFastSim(0),
  fElecIsoFastSimStop(0),
  fElecIdFastSimStop(0),
  fEleclepMVA2lSSttH(0),    // Elec ttH LepMVA 2lSS
  fEleclepMVA3l4lttH(0),    // Elec ttH LepMVA 3l, 4l
  fElecTightIP2DM17ttH(0),  // Elec ttH Tight IP2D
  fElecMini4M17ttH(0),      // Elec ttH Mini4 
  fElecConvVetoM17ttH(0),   // Elec ttH ConvVeto 

  fDoubleMuSF(0),     // Trigger Double Muon
  fDoubleElSF(0),     // Trigger Double Elec
  fMuEGSF(0),// Trigger Elec-Muon
  fSingleMuonTrigSF(0),
  fSingleElecTrigSF(0)
   {         

  era = -1;
  path_to_SF_histos = path;
  gIs2016 = false; gIs2017 = false; gIs2018 = false;
  if     (year == 2017) gIs2017 = true;
  else if(year == 2018) gIs2018 = true;
  else if(year == 2016) gIs2016 = true;

  if     (options.Contains("eraB")) era = runB;
  else if(options.Contains("eraC")) era = runC;
  else if(options.Contains("eraD")) era = runD;
  else if(options.Contains("eraE")) era = runE;
  else if(options.Contains("eraF")) era = runF;
};





void LeptonSF::loadHisto(Int_t iHisto, Int_t wp){
  TString filename = ""; TString histoname = "";
// >>>>>>>>>>>>>>>>>>>>>>>>>>>>> Muons 
  if(iHisto == iMuonReco){
    filename = "Tracking_EfficienciesAndSF_BCDEFGH"; histoname = "ratio_eff_eta3_dr030e030_corr";
    fMuonTrackerSF = LoadTrackerMuonSF(path_to_SF_histos + filename + ".root", histoname);
  }
  else if(iHisto == iMuonId){
    if(gIs2018){
      filename = "Run2018ABCD_SF_ID";
      if     (wp == iLoose)   histoname = "NUM_LooseID_DEN_TrackerMuons_pt_abseta";
      else if(wp == iMedium)  histoname = "NUM_MediumID_DEN_TrackerMuons_pt_abseta";
      else if(wp == iTight)   histoname = "NUM_TightID_DEN_TrackerMuons_pt_abseta";
      fMuonIdSF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonIdSF"); 
    }
    else if (gIs2017){
      filename = "RunBCDEF_SF_ID";
      if     (era==runB || era==runC) filename = "MuonSF_RunBC_ID";
      else if(era==runD || era==runE) filename = "MuonSF_RunDE_ID";
      else if(era==runF             ) filename = "MuonSF_RunF_ID";
      if     (wp == iLoose)   histoname = "NUM_LooseID_DEN_genTracks_pt_abseta";
      else if(wp == iMedium)  histoname = "NUM_MediumID_DEN_genTracks_pt_abseta";
      else if(wp == iTight)   histoname = "NUM_TightID_DEN_genTracks_pt_abseta";
      if(era!=-1) histoname = "NUM_TightID_DEN_genTracks";
      fMuonIdSF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonIdSF"); 
    }
    else{ // 2016
      TString fBCDEF = "2016LegacyReReco_RunBCDEF_SF_ID";
      TString fGH    = "2016LegacyReReco_RunGH_SF_ID";
      if(wp == iLoose){
        histoname = "NUM_LooseID_DEN_genTracks_eta_pt";
        fMuonIdSF_BCDEF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + fBCDEF + ".root", histoname, "fMuonIdSF_BCDEF"); 
        fMuonIdSF_GH    = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + fGH    + ".root", histoname, "fMuonIdSF_GH"); 
      }
      else if(wp == iMedium){
        histoname = "NUM_MediumID_DEN_genTracks_eta_pt";
        fMuonIdSF_BCDEF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + fBCDEF + ".root", histoname, "fMuonIdSF_BCDEF"); 
        fMuonIdSF_GH    = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + fGH    + ".root", histoname, "fMuonIdSF_GH"); 
      }
      else if(wp == iTight){
        histoname = "NUM_TightID_DEN_genTracks_eta_pt";
        fMuonIdSF_BCDEF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + fBCDEF + ".root", histoname, "fMuonIdSF_BCDEF"); 
        fMuonIdSF_GH    = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + fGH    + ".root", histoname, "fMuonIdSF_GH"); 
      }
    }
  }
  else if(iHisto == iMuonIsoTightId){ // Tight ISO over tight ID
    if(gIs2018){
      filename = "Run2018ABCD_SF_ISO";
      histoname = "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta";
      fMuonIsoSF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonIsoSF");
    }
    else if (gIs2017){
      filename = "RunBCDEF_SF_ISO";
      if     (era==runB || era==runC) filename = "MuonSF_RunBC_ISO";
      else if(era==runD || era==runE) filename = "MuonSF_RunDE_ISO";
      else if(era==runF             ) filename = "MuonSF_RunF_ISO";
      if(wp == iTight){
        histoname = "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta";
        if(era!=-1) histoname = "NUM_TightRelIso_DEN_TightIDandIPCut";
      }
      fMuonIsoSF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonIsoSF"); 
    }
    else{
      TString fBCDEF = "2016LegacyReReco_RunBCDEF_SF_ISO";
      TString fGH    = "2016LegacyReReco_RunGH_SF_ISO";
      if(wp == iLoose){
        histoname = "NUM_LooseRelIso_DEN_TightIDandIPCut_eta_pt";
        fMuonIsoSF_BCDEF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + fBCDEF + ".root", histoname, "fMuonIsoSF_BCDEF");
        fMuonIsoSF_GH    = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + fGH    + ".root", histoname, "fMuonIsoSF_GH");
      }
      else if(wp == iTight){
        histoname = "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt";
        fMuonIsoSF_BCDEF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + fBCDEF + ".root", histoname, "fMuonIdSF_BCDEF");
        fMuonIsoSF_GH    = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + fGH    + ".root", histoname, "fMuonIdSF_GH");
      }
    }
  }
  else if(iHisto == iMuonIsoMediumId){
    if(wp == iLoose){
      filename = "MuonSFIso_BCDEF"; histoname = "LooseISO_MediumID_pt_eta/pt_abseta_ratio";
      fMuonIsoSF_BCDEF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonIsoSF_BCDEF"); 
      filename = "MuonSFIso_GH"; histoname = "LooseISO_MediumD_pt_eta/pt_abseta_ratio";
      fMuonIsoSF_GH    = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonIsoSF_GH"); 
    }
    if(wp == iTight){
      filename = "MuonSFIso_BCDEF"; histoname = "TightISO_MediumID_pt_eta/pt_abseta_ratio";
      fMuonIsoSF_BCDEF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonIsoSF_BCDEF"); 
      filename = "MuonSFIso_GH"; histoname = "TightISO_MediumID_pt_eta/pt_abseta_ratio";
      fMuonIsoSF_GH    = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonIsoSF_GH"); 
    }
  }
  else if(iHisto == iMuonIsoSUSY){
    if     (wp == iLoose){  filename = ""; histoname = "";}
    else if(wp == iMedium){ filename = ""; histoname = "";}
    else if(wp == iTight){  filename = ""; histoname = "";}
    else if(wp == iVeryTight){  filename = "SUS_MuonVTMultiIsovMediumM17"; histoname = "SF";}
    else if(wp == iWPforStop){  filename = "SUS_MuonRelIso03b012"; histoname = "pt_abseta_ratio";}
    fMuonIsoSFSUSY = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonIsoSFSUSY"); 
  }
  else if(iHisto == iMuonIP2D){
    if     (wp == iLoose){  filename = ""; histoname = "";}
    else if(wp == iTight){  filename = "SUS_MuonTIP2DvMediumM17"; histoname = "SF";}
    fMuonIP2DSF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonIP2DSF"); 
  }
  else if(iHisto == iMuonSIP3D){
    filename = "SUS_MuonSIPb4vMediumM17"; histoname = "SF";
    fMuonSIP3DSF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonSIP3DSF"); 
  }
  else if(iHisto == iMuonFastSim){
    filename = "FullSimFastSim_Muon"; histoname = "SF";
    fMuonFastSim = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonFastSim"); 
  }
  else if(iHisto == iMuonIsoFastSimStop){
    filename = "Fast_vs_Full_mu_mediumID_iso012"; histoname = "histo2D";
    fMuonIsoFastSimStop = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonIsoFastSimStop"); 
  }
  else if(iHisto == iMuonIdFastSimStop){
    filename = "Fast_vs_Full_mu_mediumID"; histoname = "histo2D";
    fMuonIdFastSimStop = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonIdFastSimSrop"); 
  }
  else if(iHisto == iMuonlepMVA2lSSttH){
    filename = "lepMVAEffSF_m_2lss"; histoname = "sf";
    fMuonlepMVA2lSSttH    = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonlepMVA2lSSttH");
  }
  else if(iHisto == iMuonlepMVA3l4lttH){
    filename = "lepMVAEffSF_m_3l"; histoname = "sf";
    fMuonlepMVA3l4lttH    = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonlepMVA3l4lttH");
  }
  else if(iHisto == iMuonLooseTracksttH){
    filename = "MuonSFTTH_TnP_NUM_LooseID_DENOM_generalTracks_VAR_map_pt_eta"; histoname = "SF";
    fMuonLooseTracksttH   = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonLooseTracksttH");
  }
  else if(iHisto == iMuonLooseMiniIsottH){
    filename = "MuonSFTTH_TnP_NUM_MiniIsoLoose_DENOM_LooseID_VAR_map_pt_eta"; histoname = "SF";
    fMuonLooseMiniIsottH  = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonLooseMiniIsottH");
  }
  else if(iHisto == iMuonTightIP2DttH){
    filename = "MuonSFTTH_TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta"; histoname = "SF";
    fMuonTightIP2DttH     = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fMuonTightIP2DttH");
  }


  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>> Electrons
  else if(iHisto == iElecReco){
    // https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/EGM2D_BtoH_GT20GeV_RecoSF_Legacy2016.root
    // https://twiki.cern.ch/twiki/pub/CMS/Egamma2017DataRecommendations/egammaEffi.txt_EGM2D_runBCDEF_passingRECO.root
    // https://twiki.cern.ch/twiki/pub/CMS/EgammaIDRecipesRun2/egammaEffi.txt_EGM2D_updatedAll.root
    histoname = "EGamma_SF2D";
    if      (gIs2018) filename = "2018_ElecRECO";
    else if (gIs2017) filename = "2017_ElecRECO";
    else if (gIs2016) filename = "2016ReReco_ElecRECO";
    fElecTrackerSFF = (TH2F*) GetHistogramFromFileF(path_to_SF_histos + filename + ".root", histoname, "fElecTrackerSFF");
  }
  else if(iHisto == iElecId){
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgammaRunIIRecommendations
    // cutBasedElectronID-Summer16-80X-V1-tight
    // cutBasedElectronID-Fall17-94X-V1-tight
    // cutBasedElectronID-Fall17-94X-V2-tight
    histoname = "EGamma_SF2D";
    if(wp == iTight){
      if     (gIs2018) filename = "2018_ElectronTight";
      else if(gIs2017) filename = "2017_ElectronTight";
      else if(gIs2016) filename = "2016LegacyReReco_ElectronTight_Fall17V2";
    }
    else{
      cout << "Elec ID wp not implemented!" << endl;
    }
    fElecIdSF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fElecIdSF"); 
  }
  else if(iHisto == iElecIdSUSY){
    filename = "SUS_electonSFs";
    if(wp == iWPforStop) histoname = "GsfElectronToCutBasedSpring15T";
    fElecIdSF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fElecIdSF"); 
  }
  else if(iHisto == iElecIsoSUSY){
    filename = "SUS_electonSFs"; 
    if(wp == iWPforStop) histoname = "CutBasedStopsDileptonToRelIso012";
    fElecIsoSF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fElecIsoSF"); 
  }
  else if(iHisto == iElecIP2D){
    if     (wp == iLoose){  filename = ""; histoname = "";}
    else if(wp == iTight){  filename = ""; histoname = "";}
    fElecIP2DSF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fElecIP2DSF"); 
  }
  else if(iHisto == iElecSIP3D){
    filename = ""; histoname = "";
    fElecSIP3DSF = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fElecSIP3DSF"); 
  }
  else if(iHisto == iElecFastSim){
    filename = "FullSimFastSim_Elec"; histoname = "SF";
    fElecFastSim = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fElecFastSim"); 
  }
  else if(iHisto == iElecIsoFastSimStop){
    filename = "Fast_vs_Full_el_tightCB_iso012"; histoname = "histo2D";
    fElecIsoFastSimStop = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fElecIsoFastSimStop"); 
  }
  else if(iHisto == iElecIdFastSimStop){
    filename = "Fast_vs_Full_el_tightCB"; histoname = "histo2D";
    fElecIdFastSimStop = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fElecIdFastSimStop"); 
  }
  else if(iHisto == iEleclepMVA2lSSttH){
    filename = "lepMVAEffSF_e_2lss"; histoname = "sf";
    fEleclepMVA2lSSttH    = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fEleclepMVA2lSSttH");
  }
  else if(iHisto == iEleclepMVA3l4lttH){
    filename = "lepMVAEffSF_e_3l"; histoname = "sf";
    fEleclepMVA3l4lttH    = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fEleclepMVA3l4lttH");
  }
  else if(iHisto == iElecTightIP2DM17ttH){
    filename = "ElecSFTTH_Moriond17"; histoname = "GsfElectronToMVAVLooseFOIDEmuTightIP2D";
    fElecTightIP2DM17ttH  = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fElecTightIP2DM17ttH");
  }
  else if(iHisto == iElecMini4M17ttH){    
    filename = "ElecSFTTH_Moriond17"; histoname = "MVAVLooseElectronToMini4";
    fElecMini4M17ttH      = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fElecMini4M17ttH");
  }
  else if(iHisto == iElecConvVetoM17ttH){
    filename = "ElecSFTTH_Moriond17"; histoname = "MVAVLooseElectronToConvVetoIHit1";
    fElecConvVetoM17ttH   = (TH2D*) GetHistogramFromFileD(path_to_SF_histos + filename + ".root", histoname, "fElecConvVetoM17ttH");
  }


// >>>>>>>>>>>>>>>>>>>>>>>>>>>>> Triggers
  else if(iHisto == iTrigDoubleMuon){
    filename ="TriggerSFfromReza";
    if(gIs2016){
      histoname = "trigSF2016_Muon"; 
    }
    else if(gIs2017){
      histoname = "trigSF2017_Muon"; 
    }
    else{ // 2018
      histoname = "trigSF2018_Muon"; 
    }
    fDoubleMuSF = (TH2F*) GetHistogramFromFileF(path_to_SF_histos + filename + ".root", histoname, "fDoubleMuSF"); 
  }
  else if(iHisto == iTrigDoubleElec){
    filename ="TriggerSFfromReza";
    if(gIs2016){
      histoname = "trigSF2016_Elec";
    }
    else if(gIs2017){
      histoname = "trigSF2017_Elec";
    }
    else{ // 2018
      histoname = "trigSF2018_Elec";
    }
    fDoubleElSF = (TH2F*) GetHistogramFromFileF(path_to_SF_histos + filename + ".root", histoname, "fDoubleElSF"); 
  }
  else if(iHisto == iTrigElMu){
    filename = "TriggerSFfromReza";
    if(gIs2016){
      histoname = "trigSF2016_ElMu";
    }
    else if(gIs2017){
      histoname = "trigSF2017_ElMu";
    }
    else{ // 2018
      histoname = "trigSF2018_ElMu";
    }
    fMuEGSF = (TH2F*) GetHistogramFromFileF(path_to_SF_histos + filename + ".root", histoname, "fMuEGSF"); 
  }

  PAF_INFO("Lepton SF", Form("Loaded histogram %s from file %s%s.root", histoname.Data(), path_to_SF_histos.Data(), filename.Data()));
  loadedHistos.push_back(iHisto);
}

// Get global Lepton SF
Float_t LeptonSF::GetLeptonSF(Float_t pt, Float_t ieta, Int_t type){
  Float_t eta = TMath::Abs(ieta);
  Int_t nSFs = loadedHistos.size();
  Float_t SF = 1; Int_t id; Float_t pr = 1;
  //if(type == 1) cout << "================================= new Elec! [pt,eta] = [" << pt << "," << eta << "]..." << endl;
  //if(type == 0) cout << "================================= new Muon! [pt,eta] = [" << pt << "," << eta << "]..." << endl;
  if (pt <= 20) pt = 20.01;

  for (Int_t i = 0; i < nSFs; i++){
    pr = 1;
    id = loadedHistos[i];
    if(type == 0) {
      if(pt > 120) pt = 119;
      if (((id == iMuonlepMVA2lSSttH) || (id == iMuonlepMVA3l4lttH)) && (pt >= 100)) pt = 99;
      if     (id == iMuonReco)    pr = GetTrackerMuonSF(eta); 
      else if(id == iMuonIdSUSY)  pr = fMuonIdSFSUSY      ->GetBinContent(fMuonIdSFSUSY     ->FindBin(pt,eta));
      else if(id == iMuonId){
        if(gIs2017 || gIs2018){
          pr = fMuonIdSF->GetBinContent(fMuonIdSF->FindBin(pt,eta));
        }
        else{ // 2016
          pr = (fMuonIdSF_BCDEF->GetBinContent(fMuonIdSF_BCDEF->FindBin(eta,pt))*lumiBCDEF + fMuonIdSF_GH->GetBinContent(fMuonIdSF_GH->FindBin(eta,pt))*lumiGH)/(lumiBCDEF+lumiGH);
        }
      }
      else if(id == iMuonIsoSUSY) pr = fMuonIsoSFSUSY     ->GetBinContent(fMuonIsoSFSUSY    ->FindBin(pt,eta));
      else if(id == iMuonIsoMediumId || id == iMuonIsoTightId){
        if(gIs2017 || gIs2018){
          pr = fMuonIsoSF->GetBinContent(fMuonIsoSF->FindBin(pt,eta));
        }
        else{ // 2016
          pr = (fMuonIsoSF_BCDEF->GetBinContent(fMuonIsoSF_BCDEF->FindBin(eta,pt))*lumiBCDEF + fMuonIsoSF_GH->GetBinContent(fMuonIsoSF_GH->FindBin(eta,pt))*lumiGH)/(lumiBCDEF+lumiGH);
        }
      }
      else if (id == iMuonIP2D)            pr = fMuonIP2DSF          ->GetBinContent(fMuonIP2DSF          ->FindBin(pt,eta));
      else if (id == iMuonSIP3D)           pr = fMuonSIP3DSF         ->GetBinContent(fMuonSIP3DSF         ->FindBin(pt,eta));
      else if (id == iMuonIdFastSimStop)   pr = fMuonIdFastSimStop   ->GetBinContent(fMuonIdFastSimStop   ->FindBin(pt,eta));
      else if (id == iMuonIsoFastSimStop)  pr = fMuonIsoFastSimStop  ->GetBinContent(fMuonIsoFastSimStop  ->FindBin(pt,eta));
      else if (id == iMuonlepMVA2lSSttH)   pr = fMuonlepMVA2lSSttH   ->GetBinContent(fMuonlepMVA2lSSttH   ->FindBin(pt,eta));
      else if (id == iMuonlepMVA3l4lttH)   pr = fMuonlepMVA3l4lttH   ->GetBinContent(fMuonlepMVA3l4lttH   ->FindBin(pt,eta));
      else if (id == iMuonLooseTracksttH)  pr = fMuonLooseTracksttH  ->GetBinContent(fMuonLooseTracksttH  ->FindBin(pt,eta));
      else if (id == iMuonLooseMiniIsottH) pr = fMuonLooseMiniIsottH ->GetBinContent(fMuonLooseMiniIsottH ->FindBin(pt,eta));
      else if (id == iMuonTightIP2DttH)    pr = fMuonTightIP2DttH    ->GetBinContent(fMuonTightIP2DttH    ->FindBin(pt,eta));
    }
    else if (type == 1) { // For electrons
      if (pt > 500) pt = 499;
      if (((id == iEleclepMVA2lSSttH) || (id == iEleclepMVA3l4lttH)) && (pt >= 100)) pt = 99;
      if (id == iElecReco) {
        eta = ieta;
        pr = fElecTrackerSFF->GetBinContent(fElecTrackerSFF->FindBin(eta,pt));
      }
      else if (id == iElecIdSUSY)          pr = fElecIdSF            ->GetBinContent(fElecIdSF            ->FindBin(pt,eta));
      else if (id == iElecIsoSUSY)         pr = fElecIsoSF           ->GetBinContent(fElecIsoSF           ->FindBin(pt,eta));
      else if (id == iElecId)              pr = fElecIdSF            ->GetBinContent(fElecIdSF            ->FindBin(eta,pt));
      else if (id == iElecIP2D)            pr = fElecIP2DSF          ->GetBinContent(fElecIP2DSF          ->FindBin(eta,pt));
      else if (id == iElecSIP3D)           pr = fElecSIP3DSF         ->GetBinContent(fElecSIP3DSF         ->FindBin(eta,pt));
      else if (id == iElecIsoFastSimStop)  pr = fElecIsoFastSimStop  ->GetBinContent(fElecIsoFastSimStop  ->FindBin(pt, eta));
      else if (id == iElecIdFastSimStop)   pr = fElecIdFastSimStop   ->GetBinContent(fElecIdFastSimStop   ->FindBin(pt, eta));
      else if (id == iEleclepMVA2lSSttH)   pr = fEleclepMVA2lSSttH   ->GetBinContent(fEleclepMVA2lSSttH   ->FindBin(pt,eta));
      else if (id == iEleclepMVA3l4lttH)   pr = fEleclepMVA3l4lttH   ->GetBinContent(fEleclepMVA3l4lttH   ->FindBin(pt,eta));
      else if (id == iElecTightIP2DM17ttH) pr = fElecTightIP2DM17ttH ->GetBinContent(fElecTightIP2DM17ttH ->FindBin(pt,eta));
      else if (id == iElecMini4M17ttH)     pr = fElecMini4M17ttH     ->GetBinContent(fElecMini4M17ttH     ->FindBin(pt,eta));
      else if (id == iElecConvVetoM17ttH)  pr = fElecConvVetoM17ttH  ->GetBinContent(fElecConvVetoM17ttH  ->FindBin(pt,eta));
    }
    SF *= pr;
    //cout << "   " << SFTString[id] << ", SF = " << pr << ", total SF = " << SF << endl;
  }
  return SF;
}

// Get global error for Lepton SF
Float_t LeptonSF::GetLeptonSFerror(Float_t pt, Float_t ieta, Int_t type){
  Float_t eta = TMath::Abs(ieta);
  float t = 0;
  Int_t nSFs = loadedHistos.size();
  Float_t err = 0; Int_t id;

  if (pt <= 20) pt = 20.01;
  for (Int_t i = 0; i < nSFs; i++) {
    id = loadedHistos[i];
    if(type == 0) { // For muons
      if(pt > 120) pt = 119;
      if (((id == iMuonlepMVA2lSSttH) || (id == iMuonlepMVA3l4lttH)) && (pt >= 100)) pt = 99;
      if      (id == iMuonReco)   err += 0;
      else if (id == iMuonIdSUSY) err += p2(fMuonIdSFSUSY->GetBinError(fMuonIdSFSUSY->FindBin(pt,eta)));
      else if (id == iMuonId) {
        if (gIs2017 || gIs2018) {
          err += p2( fMuonIdSF->GetBinError(fMuonIdSF->FindBin(pt,eta)));
        }
        else { // 2016
          err += p2( (fMuonIdSF_BCDEF->GetBinError(fMuonIdSF_BCDEF->FindBin(eta,pt))*lumiBCDEF + fMuonIdSF_GH->GetBinError(fMuonIdSF_GH->FindBin(eta,pt))*lumiGH)/(lumiBCDEF+lumiGH) );
        }
      }
      else if (id == iMuonIsoSUSY) err += p2(fMuonIsoSFSUSY ->GetBinError(fMuonIsoSFSUSY->FindBin(pt,eta)));
      else if (id == iMuonIsoMediumId || id == iMuonIsoTightId) {
        if (gIs2017 || gIs2018) {
          err += p2(fMuonIsoSF->GetBinError(fMuonIsoSF->FindBin(pt,eta)));
        }
        else{ // 2016
          err += p2( (fMuonIsoSF_BCDEF->GetBinError(fMuonIsoSF_BCDEF->FindBin(pt,eta))*lumiBCDEF + fMuonIsoSF_GH->GetBinError(fMuonIsoSF_GH->FindBin(pt,eta))*lumiGH)/(lumiBCDEF+lumiGH) );
        }
      }
      else if (id == iMuonIP2D)            err += p2(fMuonIP2DSF          ->GetBinError(fMuonIP2DSF         ->FindBin(pt,eta)));
      else if (id == iMuonSIP3D)           err += p2(fMuonSIP3DSF         ->GetBinError(fMuonSIP3DSF        ->FindBin(pt,eta)));
      else if (id == iMuonIdFastSimStop)   err += p2(fMuonIdFastSimStop   ->GetBinError(fMuonIdFastSimStop  ->FindBin(pt,eta)));
      else if (id == iMuonIsoFastSimStop)  err += p2(fMuonIsoFastSimStop  ->GetBinError(fMuonIsoFastSimStop ->FindBin(pt,eta)));
      else if (id == iMuonlepMVA2lSSttH)   err += p2(fMuonlepMVA2lSSttH   ->GetBinError(fMuonlepMVA2lSSttH  ->FindBin(pt,eta)));
      else if (id == iMuonlepMVA3l4lttH)   err += p2(fMuonlepMVA3l4lttH   ->GetBinError(fMuonlepMVA3l4lttH  ->FindBin(pt,eta)));
      else if (id == iMuonLooseTracksttH)  err += p2(fMuonLooseTracksttH  ->GetBinError(fMuonLooseTracksttH ->FindBin(pt,eta)));
      else if (id == iMuonLooseMiniIsottH) err += p2(fMuonLooseMiniIsottH ->GetBinError(fMuonLooseMiniIsottH->FindBin(pt,eta)));
      else if (id == iMuonTightIP2DttH)    err += p2(fMuonTightIP2DttH    ->GetBinError(fMuonTightIP2DttH   ->FindBin(pt,eta)));
    }
    else if (type == 1) { // For electrons
      if (pt > 500) pt = 499;
      if (((id == iEleclepMVA2lSSttH) || (id == iEleclepMVA3l4lttH)) && (pt >= 100)) pt = 99;
      if (id == iElecReco) {
        eta = ieta;
        err += p2(fElecTrackerSFF->GetBinError(fElecTrackerSFF->FindBin(eta, pt)));
      }
      else if(id == iElecIdSUSY)          err += p2(fElecIdSF            ->GetBinError(fElecIdSF            ->FindBin(pt,eta)));
      else if(id == iElecIsoSUSY)         err += p2(fElecIsoSF           ->GetBinError(fElecIsoSF           ->FindBin(pt,eta)));
      else if(id == iElecId)              err += p2(fElecIdSF            ->GetBinError(fElecIdSF            ->FindBin(eta,pt)));
      else if(id == iElecIP2D)            err += p2(fElecIP2DSF          ->GetBinError(fElecIP2DSF          ->FindBin(eta,pt)));
      else if(id == iElecSIP3D)           err += p2(fElecSIP3DSF         ->GetBinError(fElecSIP3DSF         ->FindBin(eta,pt)));
      else if(id == iElecIsoFastSimStop)  err += p2(fElecIsoFastSimStop  ->GetBinError(fElecIsoFastSimStop  ->FindBin(pt, eta)));
      else if(id == iElecIdFastSimStop)   err += p2(fElecIdFastSimStop   ->GetBinError(fElecIdFastSimStop   ->FindBin(pt, eta)));
      else if(id == iEleclepMVA2lSSttH)   err += p2(fEleclepMVA2lSSttH   ->GetBinError(fEleclepMVA2lSSttH   ->FindBin(pt,eta)));
      else if(id == iEleclepMVA3l4lttH)   err += p2(fEleclepMVA3l4lttH   ->GetBinError(fEleclepMVA3l4lttH   ->FindBin(pt,eta)));
      else if(id == iElecTightIP2DM17ttH) err += p2(fElecTightIP2DM17ttH ->GetBinError(fElecTightIP2DM17ttH ->FindBin(pt,eta)));
      else if(id == iElecMini4M17ttH)     err += p2(fElecMini4M17ttH     ->GetBinError(fElecMini4M17ttH     ->FindBin(pt,eta)));
      else if(id == iElecConvVetoM17ttH)  err += p2(fElecConvVetoM17ttH  ->GetBinError(fElecConvVetoM17ttH  ->FindBin(pt,eta)));
    }
  }
  return TMath::Sqrt(err);
}



// Functions to load the histograms
TH2* LeptonSF::GetHistogramFromFileD(const char* filename, const char* histoname, const char* newhname) {
  TFile* file  = TFile::Open(filename);
  if (!file) {
    PAF_FATAL("LeptonSF", Form("Could not load file\"%s\"", filename));
    return 0;
  }
  TH2* h = (TH2*) file->Get(histoname)->Clone(newhname);
  if (!h) {
    PAF_ERROR("LeptonSF", Form("Could not find histogram \"%s\"", histoname));
    return 0;
  }
  h->SetDirectory(0);
  file->Close();
  return h;
}

TH2* LeptonSF::GetHistogramFromFileF(const char* filename, const char* histoname, const char* newhname) const {
  TFile* file  = TFile::Open(filename);
  if (!file) {
    PAF_FATAL("LeptonSF", Form("Could not load file\"%s\"", filename));
    return 0;
  }
  TH2* h = (TH2*) file->Get(histoname)->Clone(newhname);
  if (!h) {
    PAF_ERROR("LeptonSF", Form("Could not find histogram \"%s\"", histoname));
    return 0;
    }
  h->SetDirectory(0);
  file->Close();
 
  return h;
}

// -----------------------------------------
TGraphAsymmErrors* LeptonSF::LoadTrackerMuonSF(const char* file, const char* histo){
 TFile *f = TFile::Open(file);
 f->GetObject(histo, fMuonTrackerSF);
 return fMuonTrackerSF;
}

Float_t LeptonSF::GetTrackerMuonSF(Float_t eta){
  Float_t val = 0; Float_t xlow = 0; Float_t xhigh = 0;
  for(Int_t i = 0; i < fMuonTrackerSF->GetN(); i++){
    xlow  = fMuonTrackerSF->GetX()[i] - fMuonTrackerSF->GetErrorXlow(i);
    xhigh = fMuonTrackerSF->GetX()[i] + fMuonTrackerSF->GetErrorXhigh(i);
    if(xlow <= eta && xhigh > eta) val = fMuonTrackerSF->GetY()[i];
  }
  return val;
}

///////////////////////////////////////////////////
// Trigger SFs
///////////////////////////////////////////////////

Float_t LeptonSF::GetTrigDoubleMuSF(Float_t eta1, Float_t eta2) const { // binned in eta1, eta2
  eta1 = TMath::Abs(eta1);
  eta2 = TMath::Abs(eta2);
  float pt1 = eta1; if(pt1 > 200) pt1 = 199;
  float pt2 = eta2; if(pt2 > 200) pt2 = 199;
  return fDoubleMuSF->GetBinContent(fDoubleMuSF->FindBin(pt1, pt2));
  return 1;
}

Float_t LeptonSF::GetTrigDoubleElSF(Float_t eta1, Float_t eta2) const { // binned in eta1, eta2
  eta1 = TMath::Abs(eta1);
  eta2 = TMath::Abs(eta2);
  float pt1 = eta1; if(pt1 > 200) pt1 = 199;
  float pt2 = eta2; if(pt2 > 200) pt2 = 199;
  return fDoubleElSF->GetBinContent(fDoubleElSF->FindBin(pt1, pt2));
}

Float_t LeptonSF::GetTrigElMuSF(Float_t eta1, Float_t eta2) const { // binned in eta1, eta2
  eta1 = TMath::Abs(eta1);
  eta2 = TMath::Abs(eta2);
  float pt1 = eta1; if(pt1 > 200) pt1 = 199;
  float pt2 = eta2; if(pt2 > 200) pt2 = 199;
  return fMuEGSF->GetBinContent(fMuEGSF->FindBin(pt1, pt2));
}

Float_t LeptonSF::GetTrigSingleMuonSF(Float_t pt, Float_t eta) const{
  if(pt > 500) pt = 499;
  float aeta = TMath::Abs(eta);
  return fSingleMuonTrigSF->GetBinContent(fSingleMuonTrigSF->FindBin(aeta,pt));
}

// Trigger SF errors
Float_t LeptonSF::GetTrigDoubleMuSF_err(Float_t eta1, Float_t eta2) const { // binned in eta1, eta2
  if(eta1 > 200) eta1 = 199;
  if(eta2 > 200) eta2 = 199;
  return fDoubleMuSF->GetBinError(fDoubleMuSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));
}
Float_t LeptonSF::GetTrigDoubleElSF_err(Float_t eta1, Float_t eta2) const { // binned in eta1, eta2
  if(eta1 > 200) eta1 = 199;
  if(eta2 > 200) eta2 = 199;
  return fDoubleElSF->GetBinError(fDoubleElSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));
}
Float_t LeptonSF::GetTrigElMuSF_err(Float_t eta1, Float_t eta2) const { // binned in eta1, eta2
  if(eta1 > 200) eta1 = 199;
  if(eta2 > 200) eta2 = 199;
  return fMuEGSF->GetBinError(fMuEGSF->FindBin(TMath::Abs(eta1),TMath::Abs(eta2)));
}

///////////////////////////////////////////////////
// FullSim/FastSim SF
///////////////////////////////////////////////////
Float_t LeptonSF::GetFSSF(Float_t pt, Float_t eta, Int_t id){
  Float_t SF = 1;
  if(pt >= 200) pt = 199; eta = TMath::Abs(eta); id  = TMath::Abs(id);
  if(id == 13) SF = fMuonFastSim   ->GetBinContent(fMuonFastSim  ->FindBin(pt,eta));
  else         SF = fElecFastSim ->GetBinContent(fElecFastSim  ->FindBin(pt, eta));
  return SF;
}

Float_t LeptonSF::GetFSSFerr(Float_t pt, Float_t eta, Int_t id){
  Float_t SFerr = 0;
  if(pt >= 200) pt = 199; eta = TMath::Abs(eta); id  = TMath::Abs(id);
  if(id == 13) SFerr = fMuonFastSim   ->GetBinError(fMuonFastSim  ->FindBin(pt,eta));
  else         SFerr = fElecFastSim ->GetBinError(fElecFastSim  ->FindBin(pt, eta));
  return SFerr;
}
