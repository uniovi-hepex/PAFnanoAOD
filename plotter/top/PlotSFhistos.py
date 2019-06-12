from ROOT import *
path = '/nfs/fanae/user/juanr/nanoAOD/InputFiles/'
out = '/nfs/fanae/user/juanr/www/plots2017/SF/withErrors/'

gROOT.SetBatch(1)
gROOT.LoadMacro('/nfs/fanae/user/juanr/nanoAOD/InputFiles/PrintHistos.C')


#PrintHistos(fileName, hName,outname = "histo2D", min = 0.95, max = 1.05, tit = "", xtit = "", ytit = "", doErrors = false){
PrintHistos(path + "triggerSummary_mumu", "scalefactor_eta2d_with_syst", out + "trigDoubleMuon", 0.95, 1.03, "Double #mu trigger","","",1)
PrintHistos(path + "triggerSummary_ee",   "scalefactor_eta2d_with_syst", out + "trigDoubleElec", 0.95, 1.03, "Double e trigger","","",1)
PrintHistos(path + "triggerSummary_emu",  "scalefactor_eta2d_with_syst", out + "trigElMu",       0.95, 1.03, "e#mu trigger","","",1)
PrintHistos(path + "Elec2017_runBCDEF_passingRECO", "EGamma_SF2D", out + "ElecReco2017", 0.92, 1.03, "Electron reconstruction","","",1)
PrintHistos(path + "RunBCDEF_SF_ID", "NUM_TightID_DEN_genTracks_pt_abseta", out + "MuonTightID2017", 0.96, 1.02, "Muon tight ID", "p_{T} (GeV)", "|#eta|",1)
PrintHistos(path + "RunBCDEF_SF_ISO", "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta", out + "MuonTightISO2017", 0.96, 1.02, "Muon tight RelIso", "p_{T} (GeV)", "|#eta|",1)
PrintHistos(path + "ElecTightCBid94X", "EGamma_SF2D", out + "ElecTightCBIdIso2017", 0.88, 1.03, "Electron tight cut-based ID","","",1)
