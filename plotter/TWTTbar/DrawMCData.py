import ROOT     as r
import varList  as vl
import sys, os, copy
from multiprocessing import Pool
from array import array
print "===== Minitrees MC/DATA graphs plotting\n"
vl.SetUpWarnings()

NameOfTree  = vl.treename;
StandardCut = "TPassReco == 1";
ControlCut  = "TIsSS == 0 && TNJets == 1  && TNBtags == 1 && TNLooseCentral > 1";
StandardDressCut = "TPassDress == 1";
#systlist    = vl.GiveMeTheExpNamesWOJER(vl.varList["Names"]["ExpSysts"])
systlist    = ""
systdresslist = ""
#labelsignal = "+2j2b+0j_{loose}"
#labelcontrol= "+2j2b+>0j_{loose}"
labelsignal = "+2j2b"
labelcontrol= "+2j2b"
#labelsignal = ""
#labelcontrol= ""
chandir        = {}
chandir["All"] = "\\ell_{1}^{\\pm}\\ell_{2}^{\\mp}"
chandir["ElMu"]= "e^{#pm}#mu^{#mp}"
chandir["Elec"]= "e^{#pm}e^{#mp}"
chandir["Muon"]= "#mu^{#pm}#mu^{#mp}"
#legtxtsize  = 0.028
legtxtsize  = 0.055
labelpos    = (0.275, 0.89)
doPrefChecks= False
doNorm      = False
pathToTree  = vl.minipath
nCores      = 1

if (len(sys.argv) > 1):
    nCores      = int(sys.argv[1])
    if nCores != 1: print ('> Parallelization will be done with ' + str(nCores) + ' cores')
    else:           print ('> Sequential execution mode chosen')
else:
    print '> Sequential execution mode chosen'

if (len(sys.argv) > 2):
    if sys.argv[2] == 'last':
        pathToTree    = vl.GetLastFolder(vl.storagepath)
    else:
        pathToTree    = vl.storagepath + sys.argv[2] + "/"
print "> Minitrees will be read from:", pathToTree, "\n"

r.gROOT.SetBatch(True)
r.gROOT.LoadMacro('../Histo.C+')
r.gROOT.LoadMacro('../Looper.C+')
r.gROOT.LoadMacro('../TResultsTable.C+')
r.gROOT.LoadMacro('../Plot.C+')
r.gROOT.LoadMacro('../PlotToPy.C+')
r.gROOT.LoadMacro('../PlotToPyC.C+')
r.gROOT.LoadMacro('../Datacard.C+')
r.gROOT.LoadMacro('../PDFunc.C+')


def plotvariable(tsk):
    var, cut, chan = tsk
    nbins    = int(20) if "ndescbins" not in vl.varList[var] else int(vl.varList[var]['ndescbins'])
    lowedge  = float(vl.varList[var]['recobinning'][0])  if "descbinning" not in vl.varList[var] else float(vl.varList[var]['descbinning'][0])
    highedge = float(vl.varList[var]['recobinning'][-1]) if "descbinning" not in vl.varList[var] else float(vl.varList[var]['descbinning'][1])
    width    = (highedge - lowedge)/nbins
    
    p = r.PlotToPy(r.TString(vl.varList[var]['var']), r.TString(StandardCut) if cut == "signal" else r.TString(ControlCut), r.TString(chan), nbins, lowedge, highedge, r.TString(var), r.TString(vl.varList[var]['xaxis']))
    p.SetPath(pathToTree);
    p.SetTreeName(NameOfTree);
    p.SetPathSignal(pathToTree);
    p.SetTitleY("Events / " + str(int(round(width, 0))) + " GeV" if ("Eta" not in var and "DPhi" not in var) else "Events / bin" )
    p.SetLumi(vl.Lumi)
    p.verbose  = False;
    p.verbose  = True;
    p.SetChLabel(chandir[chan] + (labelsignal if cut == "signal" else labelcontrol))
    p.SetChLabelPos(labelpos[0], labelpos[1], -1)
    
    p.SetCanvasHeight(600)
    p.SetCanvasWidth(600)
    
    p.SetWeight('TWeight_normal') # TEMPORAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAL *********"*"*"*"**"*"*"*"*"*"*

    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    
    #p.AddSample("TTbar_PowhegSemi",             "Non-W/Z",      r.itBkg, 413, systlist)
    #p.AddSample("WJetsToLNu_MLM",               "Non-W/Z",      r.itBkg, 413, systlist)
    
    #p.AddSample("WZ",                           "VV+t#bar{t}V", r.itBkg, 390, systlist);
    #p.AddSample("WW",                           "VV+t#bar{t}V", r.itBkg, 390, systlist);
    #p.AddSample("ZZ",                           "VV+t#bar{t}V", r.itBkg, 390, systlist);
    #p.AddSample("TTWToLNu",                     "VV+t#bar{t}V", r.itBkg, 390, systlist);
    #p.AddSample("TTWToQQ" ,                     "VV+t#bar{t}V", r.itBkg, 390, systlist);
    #p.AddSample("TTZToQQ" ,                     "VV+t#bar{t}V", r.itBkg, 390, systlist);
    #p.AddSample("TTZToLLNuNu",                  "VV+t#bar{t}V", r.itBkg, 390, systlist);

    #p.AddSample("DYJetsToLL_M10to50_aMCatNLO",  "DY",          r.itBkg, 852, systlist);
    #p.AddSample("DYJetsToLL_M50_aMCatNLO",      "DY",          r.itBkg, 852, systlist);
    #p.AddSample("DYJetsToLL_M5to50_MLM",        "DY",           r.itBkg, 852, systlist);
    #p.AddSample("DYJetsToLL_M50_MLM",           "DY",           r.itBkg, 852, systlist);
    
    p.AddSample("TTTo2L2Nu",                 "t#bar{t}",     r.itBkg, 633, systlist)
    
    #specialweight = vl.n_ttbar/vl.sigma_ttbar/(vl.n_ttbar/vl.sigma_ttbar + vl.n_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample('TTbar_Powheg',          't#bar{t}',    r.itBkg, 633, systlist)
    #specialweight = vl.n_dilep/vl.sigma_dilep/(vl.n_ttbar/vl.sigma_ttbar + vl.n_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample('TTbar2L_powheg',        't#bar{t}',    r.itBkg, 633, systlist)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight')

    p.AddSample("tbarW_noFullHad",              "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    p.AddSample("tW_noFullHad",                 "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systlist);
    
    #specialweight = vl.n_tw/vl.sigma_tw/(vl.n_tw/vl.sigma_tw + vl.n_twnohad/vl.sigma_twnohad)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample('TW',                     'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    #specialweight = vl.n_twnohad/vl.sigma_twnohad/(vl.n_tw/vl.sigma_tw + vl.n_twnohad/vl.sigma_twnohad)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample('TW_noFullyHadr',         'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    #specialweight = vl.n_tbarw/vl.sigma_tw/(vl.n_tbarw/vl.sigma_tw + vl.n_tbarwnohad/vl.sigma_twnohad)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample('TbarW',                  'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    #specialweight = vl.n_tbarwnohad/vl.sigma_twnohad/(vl.n_tbarw/vl.sigma_tw + vl.n_tbarwnohad/vl.sigma_twnohad)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample('TbarW_noFullyHadr',      'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight')
    
    p.AddSample("b_bbar_4l",            "WWbb",      r.itSignal, r.kBlue)
    
    
    #p.AddSample("MuonEG",                       "Data",         r.itData);
    #p.AddSample("SingleMuon",                   "Data",         r.itData);
    #p.AddSample("SingleElec",                   "Data",         r.itData);
    
    #p.AddSample("TW"                         ,  "tW",           r.itSys, 1, "JERUp");
    #p.AddSample("TW_noFullyHadr_isrUp"       ,  "tW",           r.itSys, 1, "isrUp");
    #p.AddSample("TW_noFullyHadr_isrDown"     ,  "tW",           r.itSys, 1, "isrDown");
    #p.AddSample("TW_noFullyHadr_fsrUp"       ,  "tW",           r.itSys, 1, "fsrUp");
    #p.AddSample("TW_noFullyHadr_fsrDown"     ,  "tW",           r.itSys, 1, "fsrDown");
    #p.AddSample("TW_noFullyHadr_MEscaleUp"   ,  "tW",           r.itSys, 1, "tW_MEUp");
    #p.AddSample("TW_noFullyHadr_MEscaleDown" ,  "tW",           r.itSys, 1, "tW_MEDown");

    #p.AddSample("TbarW"                        ,"tW",           r.itSys, 1, "JERUp");
    #p.AddSample("TbarW_noFullyHadr_isrUp"      ,"tW",           r.itSys, 1, "isrUp");
    #p.AddSample("TbarW_noFullyHadr_isrDown"    ,"tW",           r.itSys, 1, "isrDown");
    #p.AddSample("TbarW_noFullyHadr_fsrUp"      ,"tW",           r.itSys, 1, "fsrUp");
    #p.AddSample("TbarW_noFullyHadr_fsrDown"    ,"tW",           r.itSys, 1, "fsrDown");
    #p.AddSample("TbarW_noFullyHadr_MEscaleUp"  ,"tW",           r.itSys, 1, "tW_MEUp");
    #p.AddSample("TbarW_noFullyHadr_MEscaleDown","tW",           r.itSys, 1, "tW_MEDown");

    #p.AddSample("TW_noFullyHadr_DS",            "tW",           r.itSys, 1, "DSUp");
    #p.AddSample("TbarW_noFullyHadr_DS",         "tW",           r.itSys, 1, "DSUp");
    #p.AddSymmetricHisto("tW",  "DSUp");
    #p.AddSymmetricHisto("tW",  "JERUp");
    
    #p.AddSample("TTbar_Powheg",                 "t#bar{t}",     r.itSys, 1, "JERUp");
    
    #specialweight = vl.nUEUp_ttbar/vl.sigma_ttbar/(vl.nUEUp_ttbar/vl.sigma_ttbar + vl.nUEUp_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar_Powheg_ueUp",            "t#bar{t}",     r.itSys, 1, "ueUp");
    #specialweight = vl.nUEUp_dilep/vl.sigma_dilep/(vl.nUEUp_ttbar/vl.sigma_ttbar + vl.nUEUp_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar2L_Powheg_ueUp",          "t#bar{t}",     r.itSys, 1, "ueUp");
    #specialweight = vl.nUEDown_ttbar/vl.sigma_ttbar/(vl.nUEDown_ttbar/vl.sigma_ttbar + vl.nUEDown_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar_Powheg_ueDown",          "t#bar{t}",     r.itSys, 1, "ueDown");
    #specialweight = vl.nUEDown_dilep/vl.sigma_dilep/(vl.nUEDown_ttbar/vl.sigma_ttbar + vl.nUEDown_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar2L_Powheg_ueDown",        "t#bar{t}",     r.itSys, 1, "ueDown");
    #specialweight = vl.nhDampUp_ttbar/vl.sigma_ttbar/(vl.nhDampUp_ttbar/vl.sigma_ttbar + vl.nhDampUp_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar_Powheg_hdampUp",         "t#bar{t}",     r.itSys, 1, "hdampUp");
    #specialweight = vl.nhDampUp_dilep/vl.sigma_dilep/(vl.nhDampUp_ttbar/vl.sigma_ttbar + vl.nhDampUp_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar2L_Powheg_hdampUp",       "t#bar{t}",     r.itSys, 1, "hdampUp");
    #specialweight = vl.nhDampDown_ttbar/vl.sigma_ttbar/(vl.nhDampDown_ttbar/vl.sigma_ttbar + vl.nhDampDown_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar_Powheg_hdampDown",       "t#bar{t}",     r.itSys, 1, "hdampDown");
    #specialweight = vl.nhDampDown_dilep/vl.sigma_dilep/(vl.nhDampDown_ttbar/vl.sigma_ttbar + vl.nhDampDown_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar2L_Powheg_hdampDown",     "t#bar{t}",     r.itSys, 1, "hdampDown");
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight')
    
    #p.AddSample("TTbar_Powheg_isrUp"          , "t#bar{t}",     r.itSys, 1, "isrUp");
    #p.AddSample("TTbar_Powheg_isrDown"        , "t#bar{t}",     r.itSys, 1, "isrDown");
    #p.AddSample("TTbar_Powheg_fsrUp"          , "t#bar{t}",     r.itSys, 1, "fsrUp");
    #p.AddSample("TTbar_Powheg_fsrDown"        , "t#bar{t}",     r.itSys, 1, "fsrDown");
    
    #specialweight = vl.nGluonMoveCRTune_ttbar/vl.sigma_ttbar/(vl.nGluonMoveCRTune_ttbar/vl.sigma_ttbar + vl.nGluonMoveCRTune_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight *' + str(specialweight))
    #p.AddSample('TTbar_GluonMoveCRTune',        't#bar{t}',     r.itSys, 1, "GluonMoveCRTune")
    #specialweight = vl.nGluonMoveCRTune_dilep/vl.sigma_dilep/(vl.nGluonMoveCRTune_ttbar/vl.sigma_ttbar + vl.nGluonMoveCRTune_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight *' + str(specialweight))
    #p.AddSample('TTTo2L2Nu_GluonMoveCRTune',    't#bar{t}',     r.itSys, 1, "GluonMoveCRTune")
    #specialweight = vl.nPowhegerdON_ttbar/vl.sigma_ttbar/(vl.nPowhegerdON_ttbar/vl.sigma_ttbar + vl.nPowhegerdON_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight *' + str(specialweight))
    #p.AddSample('TTbar_Powheg_erdON',           't#bar{t}',     r.itSys, 1, "Powheg_erdON")
    #specialweight = vl.nPowhegerdON_dilep/vl.sigma_dilep/(vl.nPowhegerdON_ttbar/vl.sigma_ttbar + vl.nPowhegerdON_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight *' + str(specialweight))
    #p.AddSample('TTTo2L2Nu_Powheg_erdON',       't#bar{t}',     r.itSys, 1, "Powheg_erdON")
    #specialweight = vl.nQCDbasedCRTuneerdON_ttbar/vl.sigma_ttbar/(vl.nQCDbasedCRTuneerdON_ttbar/vl.sigma_ttbar + vl.nQCDbasedCRTuneerdON_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight *' + str(specialweight))
    #p.AddSample('TTbar_GluonMoveCRTune_erdON',  't#bar{t}',     r.itSys, 1, "QCDbasedCRTune_erdON")
    #specialweight = vl.nQCDbasedCRTuneerdON_dilep/vl.sigma_dilep/(vl.nQCDbasedCRTuneerdON_ttbar/vl.sigma_ttbar + vl.nQCDbasedCRTuneerdON_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight *' + str(specialweight))
    #p.AddSample('TTTo2L2Nu_QCDbasedCRTune_erdON','t#bar{t}',    r.itSys, 1, "QCDbasedCRTune_erdON")
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight')
    
    #p.AddSample("TTbar_GluonMoveCRTune_erdON" , "t#bar{t}",     r.itSys, 1, "GluonMoveCRTune_erdON");
    #p.UseEnvelope("t#bar{t}", "GluonMoveCRTune,GluonMoveCRTune_erdON,Powheg_erdON,QCDbasedCRTune_erdON", "ColorReconnection");
    #p.AddSymmetricHisto("t#bar{t}",  "JERUp");
    
    #pdf     = r.PDFToPy(r.TString(pathToTree), r.TString("TTbar_Powheg"), r.TString(NameOfTree), r.TString(StandardCut) if cut == "signal" else r.TString(ControlCut), r.TString("All"), r.TString(vl.varList[var]['var']), nbins, lowedge, highedge);
    #pdf.verbose = False
    #pdf.SetLumi(vl.Lumi * 1000)
    #if doPrefChecks: pdf.SetWeight("TWeight * (1 - prefWeight)")
    #hPDFUp  = pdf.GetSystHisto("up","pdf").CloneHisto();
    #hPDFDown= pdf.GetSystHisto("Down","pdf").CloneHisto();
    #hMEUp   = pdf.GetSystHisto("up","ME").CloneHisto();
    #hMEDown = pdf.GetSystHisto("Down","ME").CloneHisto();
    #p.PrepareHisto(hPDFUp,   "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarPDFUp");
    #p.PrepareHisto(hPDFDown, "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarPDFDown");
    #p.PrepareHisto(hMEUp,    "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarMEUp");
    #p.PrepareHisto(hMEDown,  "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarMEDown");
    #p.AddToSystematicLabels("ttbarPDF");
    #p.AddToSystematicLabels("ttbarME");
    #del pdf
    
    # Other settings
    p.SetDataStyle("psameE1")
    #p.SetSignalStyle("SM")
    p.SetSignalStyle("scan")
    p.doUncInLegend = True;
    p.SetRatioMin( 0.6 );
    p.SetRatioMax( 1.4 );
    p.SetPadPlotMargins(vl.margins)
    p.SetPadRatioMargins(vl.marginsratio)
    p.SetTexChanSize(0.06)
    p.SetTextLumiPosX(0.69)

    #p.SetYratioOffset(0.35)
    p.SetYratioOffset(0.45)
    p.SetRatioYtitle("WWbb/tW+t#bar{t}")
    p.SetRatioStyle("S/B")
    if doNorm:
        p.SetRatioMin(0.7)
        p.SetRatioMax(1.3)
    else:
        p.SetRatioMin(0)
        p.SetRatioMax(0.2)

    p.SetCenterYAxis(False)
    p.SetXaxisOffset(1.1)
    p.ObliterateXErrorBars()
    if 'ncols' in vl.varList[var]: p.SetNColumns(vl.varList[var]['ncols'])
    
    p.SetCMSlabel("CMS");
    if vl.doPre: p.SetCMSmodeLabel("Preliminary");
    else:        p.SetCMSmodeLabel("");
    
    if 'legposdesc' in vl.varList[var]: thepos = vl.varList[var]['legposdesc']
    else:                               thepos = vl.legpos
    p.SetLegendPosition(thepos[0], thepos[1], thepos[2], thepos[3])
    p.SetLegendTextSize(legtxtsize)
    p.SetPlotFolder("./results/MCData/" if cut == 'signal' else "./results/MCData/control/");
    p.doYieldsInLeg = False;
    if "doLogY" in vl.varList[var]: p.doSetLogy = vl.varList[var]['doLogY']
    else:                           p.doSetLogy = False
    if doNorm:                      p.doSetLogy = False
    #p.doData        = False;
    #p.doSignal      = False;
    p.doSignal      = True;
    
    p.NoShowVarName = True;
    p.SetOutputName(vl.varList[var]['var_response'] * ("ATLAS" not in var) + var * ("ATLAS" in var) + "_" + chan + "_norm" * doNorm);
    #p.DrawStack();
    p.DrawStack("", doNorm);
    p.PrintSystematics()
    p.PrintYields("", "", "", "")
    p.PrintSystYields()
    del p
    #del pdf


def plotdressvariable(tsk):
    var, cut, chan = tsk
    nbins    = int(20) if "ndescbins" not in vl.varList[var] else int(vl.varList[var]['ndescbins'])
    lowedge  = float(vl.varList[var]['genbinning'][0])  if "descbinning" not in vl.varList[var] else float(vl.varList[var]['descbinning'][0])
    highedge = float(vl.varList[var]['genbinning'][-1]) if "descbinning" not in vl.varList[var] else float(vl.varList[var]['descbinning'][1])
    width    = (highedge - lowedge)/nbins

    p = r.PlotToPy(r.TString(vl.varList[var]['var_gen']), r.TString(StandardDressCut) if cut == "signal" else r.TString(ControlDressCut), r.TString(chan), nbins, lowedge, highedge, r.TString(var), r.TString(vl.varList[var]['xaxis'] if "xaxis_dress" not in vl.varList[var] else vl.varList[var]['xaxis_dress']))
    p.SetPath(pathToTree);
    p.SetTreeName(NameOfTree);
    p.SetPathSignal(pathToTree);
    p.SetTitleY("Events / " + str(int(round(width, 0))) + " GeV" if ("Eta" not in var and "DPhi" not in var) else "Events / bin" )
    p.SetLumi(vl.Lumi)
    p.verbose  = False;
    p.verbose  = True;
    p.SetChLabel(chandir[chan] + (labelsignal if cut == "signal" else labelcontrol))
    p.SetChLabelPos(labelpos[0], labelpos[1], -1)

    p.SetCanvasHeight(600)
    p.SetCanvasWidth(600)

    p.SetWeight('TWeight_normal')

    #p.AddSample("TTbar_PowhegSemi",             "Non-W/Z",      r.itBkg, 413, systdresslist)
    #p.AddSample("WJetsToLNu_MLM",               "Non-W/Z",      r.itBkg, 413, systdresslist)

    #p.AddSample("WZ",                           "VV+t#bar{t}V", r.itBkg, 390, systdresslist);
    #p.AddSample("WW",                           "VV+t#bar{t}V", r.itBkg, 390, systdresslist);
    #p.AddSample("ZZ",                           "VV+t#bar{t}V", r.itBkg, 390, systdresslist);
    #p.AddSample("TTWToLNu",                     "VV+t#bar{t}V", r.itBkg, 390, systdresslist);
    #p.AddSample("TTWToQQ" ,                     "VV+t#bar{t}V", r.itBkg, 390, systdresslist);
    #p.AddSample("TTZToQQ" ,                     "VV+t#bar{t}V", r.itBkg, 390, systdresslist);
    #p.AddSample("TTZToLLNuNu",                  "VV+t#bar{t}V", r.itBkg, 390, systdresslist);

    #p.AddSample("DYJetsToLL_M10to50_aMCatNLO",  "DY",           r.itBkg, 852, systdresslist);
    #p.AddSample("DYJetsToLL_M50_aMCatNLO",      "DY",           r.itBkg, 852, systdresslist);
    #p.AddSample("DYJetsToLL_M5to50_MLM",        "DY",           r.itBkg, 852, systdresslist);
    #p.AddSample("DYJetsToLL_M50_MLM",           "DY",           r.itBkg, 852, systdresslist);

    p.AddSample("TTTo2L2Nu",                 "t#bar{t}",     r.itBkg, 633, systdresslist)

    #specialweight = vl.n_ttbar/vl.sigma_ttbar/(vl.n_ttbar/vl.sigma_ttbar + vl.n_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample('TTbar_Powheg',          't#bar{t}',    r.itBkg, 633, systlist)
    #specialweight = vl.n_dilep/vl.sigma_dilep/(vl.n_ttbar/vl.sigma_ttbar + vl.n_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample('TTbar2L_powheg',        't#bar{t}',    r.itBkg, 633, systlist)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight')

    p.AddSample("tbarW_noFullHad",              "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systdresslist)
    p.AddSample("tW_noFullHad",                 "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systdresslist);

    #specialweight = vl.n_tw/vl.sigma_tw/(vl.n_tw/vl.sigma_tw + vl.n_twnohad/vl.sigma_twnohad)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample('TW',                     'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    #specialweight = vl.n_twnohad/vl.sigma_twnohad/(vl.n_tw/vl.sigma_tw + vl.n_twnohad/vl.sigma_twnohad)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample('TW_noFullyHadr',         'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    #specialweight = vl.n_tbarw/vl.sigma_tw/(vl.n_tbarw/vl.sigma_tw + vl.n_tbarwnohad/vl.sigma_twnohad)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample('TbarW',                  'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    #specialweight = vl.n_tbarwnohad/vl.sigma_twnohad/(vl.n_tbarw/vl.sigma_tw + vl.n_tbarwnohad/vl.sigma_twnohad)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample('TbarW_noFullyHadr',      'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight')

    p.AddSample("b_bbar_4l",            "WWbb",      r.itSignal, r.kBlue)

    #p.AddSample("MuonEG",                       "Data",         r.itData);
    #p.AddSample("SingleMuon",                   "Data",         r.itData);
    #p.AddSample("SingleElec",                   "Data",         r.itData);

    #p.AddSample("TW_noFullyHadr_isrUp"       ,  "tW",           r.itSys, 1, "isrUp");
    #p.AddSample("TW_noFullyHadr_isrDown"     ,  "tW",           r.itSys, 1, "isrDown");
    #p.AddSample("TW_noFullyHadr_fsrUp"       ,  "tW",           r.itSys, 1, "fsrUp");
    #p.AddSample("TW_noFullyHadr_fsrDown"     ,  "tW",           r.itSys, 1, "fsrDown");
    #p.AddSample("TW_noFullyHadr_MEscaleUp"   ,  "tW",           r.itSys, 1, "tW_MEUp");
    #p.AddSample("TW_noFullyHadr_MEscaleDown" ,  "tW",           r.itSys, 1, "tW_MEDown");

    #p.AddSample("TbarW_noFullyHadr_isrUp"      ,"tW",           r.itSys, 1, "isrUp");
    #p.AddSample("TbarW_noFullyHadr_isrDown"    ,"tW",           r.itSys, 1, "isrDown");
    #p.AddSample("TbarW_noFullyHadr_fsrUp"      ,"tW",           r.itSys, 1, "fsrUp");
    #p.AddSample("TbarW_noFullyHadr_fsrDown"    ,"tW",           r.itSys, 1, "fsrDown");
    #p.AddSample("TbarW_noFullyHadr_MEscaleUp"  ,"tW",           r.itSys, 1, "tW_MEUp");
    #p.AddSample("TbarW_noFullyHadr_MEscaleDown","tW",           r.itSys, 1, "tW_MEDown");

    #p.AddSample("TW_noFullyHadr_DS",            "tW",           r.itSys, 1, "DSUp");
    #p.AddSample("TbarW_noFullyHadr_DS",         "tW",           r.itSys, 1, "DSUp");
    #p.AddSymmetricHisto("tW",  "DSUp");
    #p.AddSymmetricHisto("tW",  "JERUp");


    #specialweight = vl.nUEUp_ttbar/vl.sigma_ttbar/(vl.nUEUp_ttbar/vl.sigma_ttbar + vl.nUEUp_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar_Powheg_ueUp",            "t#bar{t}",     r.itSys, 1, "ueUp");
    #specialweight = vl.nUEUp_dilep/vl.sigma_dilep/(vl.nUEUp_ttbar/vl.sigma_ttbar + vl.nUEUp_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar2L_Powheg_ueUp",          "t#bar{t}",     r.itSys, 1, "ueUp");
    #specialweight = vl.nUEDown_ttbar/vl.sigma_ttbar/(vl.nUEDown_ttbar/vl.sigma_ttbar + vl.nUEDown_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar_Powheg_ueDown",          "t#bar{t}",     r.itSys, 1, "ueDown");
    #specialweight = vl.nUEDown_dilep/vl.sigma_dilep/(vl.nUEDown_ttbar/vl.sigma_ttbar + vl.nUEDown_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar2L_Powheg_ueDown",        "t#bar{t}",     r.itSys, 1, "ueDown");
    #specialweight = vl.nhDampUp_ttbar/vl.sigma_ttbar/(vl.nhDampUp_ttbar/vl.sigma_ttbar + vl.nhDampUp_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar_Powheg_hdampUp",         "t#bar{t}",     r.itSys, 1, "hdampUp");
    #specialweight = vl.nhDampUp_dilep/vl.sigma_dilep/(vl.nhDampUp_ttbar/vl.sigma_ttbar + vl.nhDampUp_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar2L_Powheg_hdampUp",       "t#bar{t}",     r.itSys, 1, "hdampUp");
    #specialweight = vl.nhDampDown_ttbar/vl.sigma_ttbar/(vl.nhDampDown_ttbar/vl.sigma_ttbar + vl.nhDampDown_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar_Powheg_hdampDown",       "t#bar{t}",     r.itSys, 1, "hdampDown");
    #specialweight = vl.nhDampDown_dilep/vl.sigma_dilep/(vl.nhDampDown_ttbar/vl.sigma_ttbar + vl.nhDampDown_dilep/vl.sigma_dilep);
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight*' + str(specialweight))
    #p.AddSample("TTbar2L_Powheg_hdampDown",     "t#bar{t}",     r.itSys, 1, "hdampDown");
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight')

    #p.AddSample("TTbar_Powheg_isrUp"          , "t#bar{t}",     r.itSys, 1, "isrUp");
    #p.AddSample("TTbar_Powheg_isrDown"        , "t#bar{t}",     r.itSys, 1, "isrDown");
    #p.AddSample("TTbar_Powheg_fsrUp"          , "t#bar{t}",     r.itSys, 1, "fsrUp");
    #p.AddSample("TTbar_Powheg_fsrDown"        , "t#bar{t}",     r.itSys, 1, "fsrDown");

    #specialweight = vl.nGluonMoveCRTune_ttbar/vl.sigma_ttbar/(vl.nGluonMoveCRTune_ttbar/vl.sigma_ttbar + vl.nGluonMoveCRTune_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight *' + str(specialweight))
    #p.AddSample('TTbar_GluonMoveCRTune',        't#bar{t}',     r.itSys, 1, "GluonMoveCRTune")
    #specialweight = vl.nGluonMoveCRTune_dilep/vl.sigma_dilep/(vl.nGluonMoveCRTune_ttbar/vl.sigma_ttbar + vl.nGluonMoveCRTune_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight *' + str(specialweight))
    #p.AddSample('TTTo2L2Nu_GluonMoveCRTune',    't#bar{t}',     r.itSys, 1, "GluonMoveCRTune")
    #specialweight = vl.nPowhegerdON_ttbar/vl.sigma_ttbar/(vl.nPowhegerdON_ttbar/vl.sigma_ttbar + vl.nPowhegerdON_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight *' + str(specialweight))
    #p.AddSample('TTbar_Powheg_erdON',           't#bar{t}',     r.itSys, 1, "Powheg_erdON")
    #specialweight = vl.nPowhegerdON_dilep/vl.sigma_dilep/(vl.nPowhegerdON_ttbar/vl.sigma_ttbar + vl.nPowhegerdON_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight *' + str(specialweight))
    #p.AddSample('TTTo2L2Nu_Powheg_erdON',       't#bar{t}',     r.itSys, 1, "Powheg_erdON")
    #specialweight = vl.nQCDbasedCRTuneerdON_ttbar/vl.sigma_ttbar/(vl.nQCDbasedCRTuneerdON_ttbar/vl.sigma_ttbar + vl.nQCDbasedCRTuneerdON_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight *' + str(specialweight))
    #p.AddSample('TTbar_GluonMoveCRTune_erdON',  't#bar{t}',     r.itSys, 1, "QCDbasedCRTune_erdON")
    #specialweight = vl.nQCDbasedCRTuneerdON_dilep/vl.sigma_dilep/(vl.nQCDbasedCRTuneerdON_ttbar/vl.sigma_ttbar + vl.nQCDbasedCRTuneerdON_dilep/vl.sigma_dilep)
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight *' + str(specialweight))
    #p.AddSample('TTTo2L2Nu_QCDbasedCRTune_erdON','t#bar{t}',    r.itSys, 1, "QCDbasedCRTune_erdON")
    #if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    #else:            p.SetWeight('TWeight')

    #p.AddSample("TTbar_GluonMoveCRTune_erdON" , "t#bar{t}",     r.itSys, 1, "GluonMoveCRTune_erdON");
    #p.UseEnvelope("t#bar{t}", "GluonMoveCRTune,GluonMoveCRTune_erdON,Powheg_erdON,QCDbasedCRTune_erdON", "ColorReconnection");
    #p.AddSymmetricHisto("t#bar{t}",  "JERUp");

    #pdf     = r.PDFToPy(r.TString(pathToTree), r.TString("TTbar_Powheg"), r.TString(NameOfTree), r.TString(StandardCut) if cut == "signal" else r.TString(ControlCut), r.TString("All"), r.TString(vl.varList[var]['var']), nbins, lowedge, highedge);
    #pdf.verbose = False
    #pdf.SetLumi(vl.Lumi * 1000)
    #if doPrefChecks: pdf.SetWeight("TWeight * (1 - prefWeight)")
    #hPDFUp  = pdf.GetSystHisto("up","pdf").CloneHisto();
    #hPDFDown= pdf.GetSystHisto("Down","pdf").CloneHisto();
    #hMEUp   = pdf.GetSystHisto("up","ME").CloneHisto();
    #hMEDown = pdf.GetSystHisto("Down","ME").CloneHisto();
    #p.PrepareHisto(hPDFUp,   "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarPDFUp");
    #p.PrepareHisto(hPDFDown, "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarPDFDown");
    #p.PrepareHisto(hMEUp,    "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarMEUp");
    #p.PrepareHisto(hMEDown,  "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarMEDown");
    #p.AddToSystematicLabels("ttbarPDF");
    #p.AddToSystematicLabels("ttbarME");
    #del pdf

    # Other settings
    p.SetDataStyle("psameE1")
    #p.SetSignalStyle("SM")
    p.SetSignalStyle("scan")
    p.doUncInLegend = True;
    p.SetRatioMin( 0.6 );
    p.SetRatioMax( 1.4 );
    p.SetPadPlotMargins(vl.margins)
    p.SetPadRatioMargins(vl.marginsratio)
    p.SetTexChanSize(0.06)
    p.SetTextLumiPosX(0.69)

    #p.SetYratioOffset(0.35)
    p.SetYratioOffset(0.45)
    p.SetRatioYtitle("WWbb/tW+t#bar{t}")
    p.SetRatioStyle("S/B")
    if doNorm:
        p.SetRatioMin(0.7)
        p.SetRatioMax(1.3)
    else:
        p.SetRatioMin(0)
        p.SetRatioMax(0.2)

    p.SetCenterYAxis(False)
    p.SetXaxisOffset(1.1)
    p.ObliterateXErrorBars()
    if 'ncols' in vl.varList[var]: p.SetNColumns(vl.varList[var]['ncols'])

    p.SetCMSlabel("CMS");
    if vl.doPre: p.SetCMSmodeLabel("Preliminary");
    else:        p.SetCMSmodeLabel("");

    if 'legposdesc' in vl.varList[var]: thepos = vl.varList[var]['legposdesc']
    else:                               thepos = vl.legpos
    p.SetLegendPosition(thepos[0], thepos[1], thepos[2], thepos[3])
    p.SetLegendTextSize(legtxtsize)
    p.SetPlotFolder("./results/MCData/" if cut == 'signal' else "./results/MCData/control/");
    p.doYieldsInLeg = False;
    if "doLogY" in vl.varList[var]: p.doSetLogy = vl.varList[var]['doLogY']
    else:                           p.doSetLogy = False
    if doNorm:                      p.doSetLogy = False
    #p.doData        = False;
    #p.doSignal      = False;
    p.doSignal      = True;

    p.NoShowVarName = True;
    p.SetOutputName("Dress" + vl.varList[var]['var_response'] * ("ATLAS" not in var) + var * ("ATLAS" in var) + "_" + chan + "_norm" * doNorm);
    #p.DrawStack();
    p.DrawStack("", doNorm);
    p.PrintSystematics()
    p.PrintYields("", "", "", "")
    p.PrintSystYields()
    del p
    #del pdf


def plotcustomvariable(tsk):
    var, cut, chan = tsk
    binning = array('f', vl.varList[var]['recobinning']) # For some reason, ROOT requires that you create FIRST this object, then put it inside the PlotToPyC.
    p = r.PlotToPyC(r.TString(vl.varList[var]['var']), r.TString(StandardCut) if cut == "signal" else r.TString(ControlCut), r.TString('All'), int(len(vl.varList[var]['recobinning']) - 1), binning, r.TString(var), r.TString(vl.varList[var]['xaxis']))
    p.SetPath(pathToTree);
    p.SetTreeName(NameOfTree);
    p.SetPathSignal(pathToTree);
    p.SetTitleY("Events")
    p.SetLumi(vl.Lumi)
    p.verbose  = False;
    p.verbose  = True;
    p.SetChLabel(labelsignal if cut == "signal" else labelcontrol)
    p.SetChLabelPos(labelpos[0], labelpos[1], -1)
    
    p.SetCanvasHeight(600)
    p.SetCanvasWidth(600)
    
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    
    p.AddSample("TTbar_PowhegSemi",             "Non-W/Z",      r.itBkg, 413, systlist)
    p.AddSample("WJetsToLNu_MLM",               "Non-W/Z",      r.itBkg, 413, systlist)
    
    p.AddSample("WZ",                           "VV+t#bar{t}V", r.itBkg, 390, systlist);
    p.AddSample("WW",                           "VV+t#bar{t}V", r.itBkg, 390, systlist);
    p.AddSample("ZZ",                           "VV+t#bar{t}V", r.itBkg, 390, systlist);
    p.AddSample("TTWToLNu",                     "VV+t#bar{t}V", r.itBkg, 390, systlist);
    p.AddSample("TTWToQQ" ,                     "VV+t#bar{t}V", r.itBkg, 390, systlist);
    p.AddSample("TTZToQQ" ,                     "VV+t#bar{t}V", r.itBkg, 390, systlist);
    p.AddSample("TTZToLLNuNu",                  "VV+t#bar{t}V", r.itBkg, 390, systlist);

    p.AddSample("DYJetsToLL_M10to50_aMCatNLO",  "DY",          r.itBkg, 852, systlist);
    p.AddSample("DYJetsToLL_M50_aMCatNLO",      "DY",          r.itBkg, 852, systlist);
    #p.AddSample("DYJetsToLL_M5to50_MLM",        "DY",           r.itBkg, 852, systlist);
    #p.AddSample("DYJetsToLL_M50_MLM",           "DY",           r.itBkg, 852, systlist);
    
    #p.AddSample("TTbar_Powheg",                 "t#bar{t}",     r.itBkg, 633, systlist)
    
    specialweight = vl.n_ttbar/vl.sigma_ttbar/(vl.n_ttbar/vl.sigma_ttbar + vl.n_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight *' + str(specialweight))
    p.AddSample('TTbar_Powheg',          't#bar{t}',    r.itBkg, 633, systlist)
    specialweight = vl.n_dilep/vl.sigma_dilep/(vl.n_ttbar/vl.sigma_ttbar + vl.n_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight *' + str(specialweight))
    p.AddSample('TTbar2L_powheg',        't#bar{t}',    r.itBkg, 633, systlist)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight')
    
    #p.AddSample("TW",                           "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    #p.AddSample("TbarW",                        "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systlist);
    
    specialweight = vl.n_tw/vl.sigma_tw/(vl.n_tw/vl.sigma_tw + vl.n_twnohad/vl.sigma_twnohad)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight *' + str(specialweight))
    p.AddSample('TW',                     'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    specialweight = vl.n_twnohad/vl.sigma_twnohad/(vl.n_tw/vl.sigma_tw + vl.n_twnohad/vl.sigma_twnohad)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight *' + str(specialweight))
    p.AddSample('TW_noFullyHadr',         'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    specialweight = vl.n_tbarw/vl.sigma_tw/(vl.n_tbarw/vl.sigma_tw + vl.n_tbarwnohad/vl.sigma_twnohad)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight *' + str(specialweight))
    p.AddSample('TbarW',                  'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    specialweight = vl.n_tbarwnohad/vl.sigma_twnohad/(vl.n_tbarw/vl.sigma_tw + vl.n_tbarwnohad/vl.sigma_twnohad)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight *' + str(specialweight))
    p.AddSample('TbarW_noFullyHadr',      'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight')
    
    p.AddSample("MuonEG",                       "Data",         r.itData);
    p.AddSample("SingleMuon",                   "Data",         r.itData);
    p.AddSample("SingleElec",                   "Data",         r.itData);
    
    p.AddSample("TW"                         ,  "tW",           r.itSys, 1, "JERUp");
    p.AddSample("TW_noFullyHadr_isrUp"       ,  "tW",           r.itSys, 1, "isrUp");
    p.AddSample("TW_noFullyHadr_isrDown"     ,  "tW",           r.itSys, 1, "isrDown");
    p.AddSample("TW_noFullyHadr_fsrUp"       ,  "tW",           r.itSys, 1, "fsrUp");
    p.AddSample("TW_noFullyHadr_fsrDown"     ,  "tW",           r.itSys, 1, "fsrDown");
    p.AddSample("TW_noFullyHadr_MEscaleUp"   ,  "tW",           r.itSys, 1, "tW_MEUp");
    p.AddSample("TW_noFullyHadr_MEscaleDown" ,  "tW",           r.itSys, 1, "tW_MEDown");

    p.AddSample("TbarW"                        ,"tW",           r.itSys, 1, "JERUp");
    p.AddSample("TbarW_noFullyHadr_isrUp"      ,"tW",           r.itSys, 1, "isrUp");
    p.AddSample("TbarW_noFullyHadr_isrDown"    ,"tW",           r.itSys, 1, "isrDown");
    p.AddSample("TbarW_noFullyHadr_fsrUp"      ,"tW",           r.itSys, 1, "fsrUp");
    p.AddSample("TbarW_noFullyHadr_fsrDown"    ,"tW",           r.itSys, 1, "fsrDown");
    p.AddSample("TbarW_noFullyHadr_MEscaleUp"  ,"tW",           r.itSys, 1, "tW_MEUp");
    p.AddSample("TbarW_noFullyHadr_MEscaleDown","tW",           r.itSys, 1, "tW_MEDown");

    p.AddSample("TW_noFullyHadr_DS",            "tW",           r.itSys, 1, "DSUp");
    p.AddSample("TbarW_noFullyHadr_DS",         "tW",           r.itSys, 1, "DSUp");
    p.AddSymmetricHisto("tW",  "DSUp");
    p.AddSymmetricHisto("tW",  "JERUp");
    
    p.AddSample("TTbar_Powheg",                 "t#bar{t}",     r.itSys, 1, "JERUp");
    
    specialweight = vl.nUEUp_ttbar/vl.sigma_ttbar/(vl.nUEUp_ttbar/vl.sigma_ttbar + vl.nUEUp_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar_Powheg_ueUp",            "t#bar{t}",     r.itSys, 1, "ueUp");
    specialweight = vl.nUEUp_dilep/vl.sigma_dilep/(vl.nUEUp_ttbar/vl.sigma_ttbar + vl.nUEUp_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar2L_Powheg_ueUp",          "t#bar{t}",     r.itSys, 1, "ueUp");
    specialweight = vl.nUEDown_ttbar/vl.sigma_ttbar/(vl.nUEDown_ttbar/vl.sigma_ttbar + vl.nUEDown_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar_Powheg_ueDown",          "t#bar{t}",     r.itSys, 1, "ueDown");
    specialweight = vl.nUEDown_dilep/vl.sigma_dilep/(vl.nUEDown_ttbar/vl.sigma_ttbar + vl.nUEDown_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar2L_Powheg_ueDown",        "t#bar{t}",     r.itSys, 1, "ueDown");
    specialweight = vl.nhDampUp_ttbar/vl.sigma_ttbar/(vl.nhDampUp_ttbar/vl.sigma_ttbar + vl.nhDampUp_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar_Powheg_hdampUp",         "t#bar{t}",     r.itSys, 1, "hdampUp");
    specialweight = vl.nhDampUp_dilep/vl.sigma_dilep/(vl.nhDampUp_ttbar/vl.sigma_ttbar + vl.nhDampUp_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar2L_Powheg_hdampUp",       "t#bar{t}",     r.itSys, 1, "hdampUp");
    specialweight = vl.nhDampDown_ttbar/vl.sigma_ttbar/(vl.nhDampDown_ttbar/vl.sigma_ttbar + vl.nhDampDown_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar_Powheg_hdampDown",       "t#bar{t}",     r.itSys, 1, "hdampDown");
    specialweight = vl.nhDampDown_dilep/vl.sigma_dilep/(vl.nhDampDown_ttbar/vl.sigma_ttbar + vl.nhDampDown_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar2L_Powheg_hdampDown",     "t#bar{t}",     r.itSys, 1, "hdampDown");
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight')
    
    p.AddSample("TTbar_Powheg_isrUp"          , "t#bar{t}",     r.itSys, 1, "isrUp");
    p.AddSample("TTbar_Powheg_isrDown"        , "t#bar{t}",     r.itSys, 1, "isrDown");
    p.AddSample("TTbar_Powheg_fsrUp"          , "t#bar{t}",     r.itSys, 1, "fsrUp");
    p.AddSample("TTbar_Powheg_fsrDown"        , "t#bar{t}",     r.itSys, 1, "fsrDown");
    
    specialweight = vl.nGluonMoveCRTune_ttbar/vl.sigma_ttbar/(vl.nGluonMoveCRTune_ttbar/vl.sigma_ttbar + vl.nGluonMoveCRTune_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample('TTbar_GluonMoveCRTune',        't#bar{t}',     r.itSys, 1, "GluonMoveCRTune")
    specialweight = vl.nGluonMoveCRTune_dilep/vl.sigma_dilep/(vl.nGluonMoveCRTune_ttbar/vl.sigma_ttbar + vl.nGluonMoveCRTune_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample('TTTo2L2Nu_GluonMoveCRTune',    't#bar{t}',     r.itSys, 1, "GluonMoveCRTune")
    specialweight = vl.nPowhegerdON_ttbar/vl.sigma_ttbar/(vl.nPowhegerdON_ttbar/vl.sigma_ttbar + vl.nPowhegerdON_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample('TTbar_Powheg_erdON',           't#bar{t}',     r.itSys, 1, "Powheg_erdON")
    specialweight = vl.nPowhegerdON_dilep/vl.sigma_dilep/(vl.nPowhegerdON_ttbar/vl.sigma_ttbar + vl.nPowhegerdON_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample('TTTo2L2Nu_Powheg_erdON',       't#bar{t}',     r.itSys, 1, "Powheg_erdON")
    specialweight = vl.nQCDbasedCRTuneerdON_ttbar/vl.sigma_ttbar/(vl.nQCDbasedCRTuneerdON_ttbar/vl.sigma_ttbar + vl.nQCDbasedCRTuneerdON_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample('TTbar_GluonMoveCRTune_erdON',  't#bar{t}',     r.itSys, 1, "QCDbasedCRTune_erdON")
    specialweight = vl.nQCDbasedCRTuneerdON_dilep/vl.sigma_dilep/(vl.nQCDbasedCRTuneerdON_ttbar/vl.sigma_ttbar + vl.nQCDbasedCRTuneerdON_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample('TTTo2L2Nu_QCDbasedCRTune_erdON','t#bar{t}',    r.itSys, 1, "QCDbasedCRTune_erdON")
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight')
    
    p.AddSample("TTbar_GluonMoveCRTune_erdON" , "t#bar{t}",     r.itSys, 1, "GluonMoveCRTune_erdON");
    p.UseEnvelope("t#bar{t}", "GluonMoveCRTune,GluonMoveCRTune_erdON,Powheg_erdON,QCDbasedCRTune_erdON", "ColorReconnection");
    p.AddSymmetricHisto("t#bar{t}",  "JERUp");
    
    pdf     = r.PDFToPyC(r.TString(pathToTree), r.TString("TTbar_Powheg"), r.TString(NameOfTree), r.TString(StandardCut) if cut == "signal" else r.TString(ControlCut), r.TString("All"), r.TString(vl.varList[var]['var']), len(vl.varList[var]['recobinning']) - 1, binning, r.TString(''));
    pdf.verbose = False
    #pdf.verbose = True
    pdf.SetLumi(vl.Lumi * 1000)
    if doPrefChecks: pdf.SetWeight("TWeight * (1 - prefWeight)")
    
    hPDFUp  = pdf.GetSystHisto("up","pdf").CloneHisto();
    hPDFDown= pdf.GetSystHisto("Down","pdf").CloneHisto();
    hMEUp   = pdf.GetSystHisto("up","ME").CloneHisto();
    hMEDown = pdf.GetSystHisto("Down","ME").CloneHisto();
    p.PrepareHisto(hPDFUp,   "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarPDFUp");
    p.PrepareHisto(hPDFDown, "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarPDFDown");
    p.PrepareHisto(hMEUp,    "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarMEUp");
    p.PrepareHisto(hMEDown,  "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarMEDown");
    p.AddToSystematicLabels("ttbarPDF");
    p.AddToSystematicLabels("ttbarME");
    del pdf
    
    # Other settings
    p.SetDataStyle("psameE1")
    #p.SetStackErrorStyle(3144)
    p.doUncInLegend = True;
    p.SetRatioMin( 0.6 );
    p.SetRatioMax( 1.4 );
    p.SetPadPlotMargins(vl.margins)
    p.SetPadRatioMargins(vl.marginsratio)
    p.SetTexChanSize(0.06)
    p.SetYratioOffset(0.35)
    p.SetCenterYAxis(False)
    if "equalbinsfol" in vl.varList[var]: p.ObliterateXErrorBars()
    if 'ncols' in vl.varList[var]: p.SetNColumns(vl.varList[var]['ncols'])
    
    p.SetCMSlabel("CMS");
    if vl.doPre: p.SetCMSmodeLabel("Preliminary");
    else:        p.SetCMSmodeLabel("");
    if 'legpos' in vl.varList[var]: thepos = vl.varList[var]['legpos']
    else:                           thepos = vl.legpos
    p.SetLegendPosition(thepos[0], thepos[1], thepos[2], thepos[3])
    p.SetLegendTextSize(legtxtsize)
    p.SetPlotFolder("/nfs/fanae/user/vrbouza/www/TFM/1j1t/" if cut == 'signal' else "/nfs/fanae/user/vrbouza/www/TFM/1j1t/control/");
    p.doYieldsInLeg = False;
    if "doLogY" in vl.varList[var]: p.doSetLogy = vl.varList[var]['doLogY']
    else:                           p.doSetLogy = False
    #p.doData        = False;
    p.doSignal      = False;
    #p.SetTitleY(r.TString(vl.varList[var]['yaxis']))
    p.SetOutputName("Custom");
    if "abs" in vl.varList[var]['var'] or "min" in vl.varList[var]['var']:
        p.NoShowVarName = True;
        p.SetOutputName('Custom_' + vl.varList[var]['var_response'] * ("ATLAS" not in var) + var * ("ATLAS" in var) + "_" + chan);
    print "JOJOJO"
    p.DrawStack();
    print "JEJEJE"
    p.PrintSystematics()
    p.PrintYields("", "", "", "")
    p.PrintSystYields()
    del p


def plotthenumberofjets(tsk):
    var, cut, chan = tsk
    binning = array('f', vl.varList[var]['recobinning']) # For some reason, ROOT requires that you create FIRST this object, then put it inside the PlotToPyC.
    p = r.PlotToPyC(r.TString(vl.varList[var]['var'] + "- 1"), r.TString("TIsSS == 0 && TNJets == 1  && TNBtags == 1"), r.TString('All'), int(len(vl.varList[var]['recobinning']) - 1), binning, r.TString(var), r.TString(vl.varList[var]['xaxis']))
    p.SetPath(pathToTree);
    p.SetTreeName(NameOfTree);
    p.SetPathSignal(pathToTree);
    p.SetTitleY("Events / bin")
    p.SetLumi(vl.Lumi)
    p.verbose  = False;
    p.verbose  = True;
    p.SetChLabel("e^{#pm}#mu^{#mp}+1j1b")
    p.SetChLabelPos(labelpos[0], labelpos[1], -1)
    
    p.SetCanvasHeight(600)
    p.SetCanvasWidth(600)
    
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    
    p.AddSample("TTbar_PowhegSemi",             "Non-W/Z",      r.itBkg, 413, systlist)
    p.AddSample("WJetsToLNu_MLM",               "Non-W/Z",      r.itBkg, 413, systlist)
    
    p.AddSample("WZ",                           "VV+t#bar{t}V", r.itBkg, 390, systlist);
    p.AddSample("WW",                           "VV+t#bar{t}V", r.itBkg, 390, systlist);
    p.AddSample("ZZ",                           "VV+t#bar{t}V", r.itBkg, 390, systlist);
    p.AddSample("TTWToLNu",                     "VV+t#bar{t}V", r.itBkg, 390, systlist);
    p.AddSample("TTWToQQ" ,                     "VV+t#bar{t}V", r.itBkg, 390, systlist);
    p.AddSample("TTZToQQ" ,                     "VV+t#bar{t}V", r.itBkg, 390, systlist);
    p.AddSample("TTZToLLNuNu",                  "VV+t#bar{t}V", r.itBkg, 390, systlist);

    p.AddSample("DYJetsToLL_M10to50_aMCatNLO",  "DY",          r.itBkg, 852, systlist);
    p.AddSample("DYJetsToLL_M50_aMCatNLO",      "DY",          r.itBkg, 852, systlist);
    #p.AddSample("DYJetsToLL_M5to50_MLM",        "DY",           r.itBkg, 852, systlist);
    #p.AddSample("DYJetsToLL_M50_MLM",           "DY",           r.itBkg, 852, systlist);
    
    #p.AddSample("TTbar_Powheg",                 "t#bar{t}",     r.itBkg, 633, systlist)
    
    specialweight = vl.n_ttbar/vl.sigma_ttbar/(vl.n_ttbar/vl.sigma_ttbar + vl.n_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample('TTbar_Powheg',          't#bar{t}',    r.itBkg, 633, systlist)
    specialweight = vl.n_dilep/vl.sigma_dilep/(vl.n_ttbar/vl.sigma_ttbar + vl.n_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample('TTbar2L_powheg',        't#bar{t}',    r.itBkg, 633, systlist)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight')
    
    #p.AddSample("TW",                           "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    #p.AddSample("TbarW",                        "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systlist);
    
    specialweight = vl.n_tw/vl.sigma_tw/(vl.n_tw/vl.sigma_tw + vl.n_twnohad/vl.sigma_twnohad)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample('TW',                     'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    specialweight = vl.n_twnohad/vl.sigma_twnohad/(vl.n_tw/vl.sigma_tw + vl.n_twnohad/vl.sigma_twnohad)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample('TW_noFullyHadr',         'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    specialweight = vl.n_tbarw/vl.sigma_tw/(vl.n_tbarw/vl.sigma_tw + vl.n_tbarwnohad/vl.sigma_twnohad)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample('TbarW',                  'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    specialweight = vl.n_tbarwnohad/vl.sigma_twnohad/(vl.n_tbarw/vl.sigma_tw + vl.n_tbarwnohad/vl.sigma_twnohad)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample('TbarW_noFullyHadr',      'tW',      r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight')
    
    p.AddSample("MuonEG",                       "Data",         r.itData);
    p.AddSample("SingleMuon",                   "Data",         r.itData);
    p.AddSample("SingleElec",                   "Data",         r.itData);
    
    p.AddSample("TW"                         ,  "tW",           r.itSys, 1, "JERUp");
    p.AddSample("TW_noFullyHadr_isrUp"       ,  "tW",           r.itSys, 1, "isrUp");
    p.AddSample("TW_noFullyHadr_isrDown"     ,  "tW",           r.itSys, 1, "isrDown");
    p.AddSample("TW_noFullyHadr_fsrUp"       ,  "tW",           r.itSys, 1, "fsrUp");
    p.AddSample("TW_noFullyHadr_fsrDown"     ,  "tW",           r.itSys, 1, "fsrDown");
    p.AddSample("TW_noFullyHadr_MEscaleUp"   ,  "tW",           r.itSys, 1, "tW_MEUp");
    p.AddSample("TW_noFullyHadr_MEscaleDown" ,  "tW",           r.itSys, 1, "tW_MEDown");

    p.AddSample("TbarW"                        ,"tW",           r.itSys, 1, "JERUp");
    p.AddSample("TbarW_noFullyHadr_isrUp"      ,"tW",           r.itSys, 1, "isrUp");
    p.AddSample("TbarW_noFullyHadr_isrDown"    ,"tW",           r.itSys, 1, "isrDown");
    p.AddSample("TbarW_noFullyHadr_fsrUp"      ,"tW",           r.itSys, 1, "fsrUp");
    p.AddSample("TbarW_noFullyHadr_fsrDown"    ,"tW",           r.itSys, 1, "fsrDown");
    p.AddSample("TbarW_noFullyHadr_MEscaleUp"  ,"tW",           r.itSys, 1, "tW_MEUp");
    p.AddSample("TbarW_noFullyHadr_MEscaleDown","tW",           r.itSys, 1, "tW_MEDown");

    p.AddSample("TW_noFullyHadr_DS",            "tW",           r.itSys, 1, "DSUp");
    p.AddSample("TbarW_noFullyHadr_DS",         "tW",           r.itSys, 1, "DSUp");
    p.AddSymmetricHisto("tW",  "DSUp");
    p.AddSymmetricHisto("tW",  "JERUp");
    
    p.AddSample("TTbar_Powheg",                 "t#bar{t}",     r.itSys, 1, "JERUp");
    
    specialweight = vl.nUEUp_ttbar/vl.sigma_ttbar/(vl.nUEUp_ttbar/vl.sigma_ttbar + vl.nUEUp_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar_Powheg_ueUp",            "t#bar{t}",     r.itSys, 1, "ueUp");
    specialweight = vl.nUEUp_dilep/vl.sigma_dilep/(vl.nUEUp_ttbar/vl.sigma_ttbar + vl.nUEUp_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar2L_Powheg_ueUp",          "t#bar{t}",     r.itSys, 1, "ueUp");
    specialweight = vl.nUEDown_ttbar/vl.sigma_ttbar/(vl.nUEDown_ttbar/vl.sigma_ttbar + vl.nUEDown_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar_Powheg_ueDown",          "t#bar{t}",     r.itSys, 1, "ueDown");
    specialweight = vl.nUEDown_dilep/vl.sigma_dilep/(vl.nUEDown_ttbar/vl.sigma_ttbar + vl.nUEDown_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar2L_Powheg_ueDown",        "t#bar{t}",     r.itSys, 1, "ueDown");
    specialweight = vl.nhDampUp_ttbar/vl.sigma_ttbar/(vl.nhDampUp_ttbar/vl.sigma_ttbar + vl.nhDampUp_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar_Powheg_hdampUp",         "t#bar{t}",     r.itSys, 1, "hdampUp");
    specialweight = vl.nhDampUp_dilep/vl.sigma_dilep/(vl.nhDampUp_ttbar/vl.sigma_ttbar + vl.nhDampUp_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar2L_Powheg_hdampUp",       "t#bar{t}",     r.itSys, 1, "hdampUp");
    specialweight = vl.nhDampDown_ttbar/vl.sigma_ttbar/(vl.nhDampDown_ttbar/vl.sigma_ttbar + vl.nhDampDown_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar_Powheg_hdampDown",       "t#bar{t}",     r.itSys, 1, "hdampDown");
    specialweight = vl.nhDampDown_dilep/vl.sigma_dilep/(vl.nhDampDown_ttbar/vl.sigma_ttbar + vl.nhDampDown_dilep/vl.sigma_dilep);
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) * ' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight*' + str(specialweight))
    p.AddSample("TTbar2L_Powheg_hdampDown",     "t#bar{t}",     r.itSys, 1, "hdampDown");
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight')
    
    p.AddSample("TTbar_Powheg_isrUp"          , "t#bar{t}",     r.itSys, 1, "isrUp");
    p.AddSample("TTbar_Powheg_isrDown"        , "t#bar{t}",     r.itSys, 1, "isrDown");
    p.AddSample("TTbar_Powheg_fsrUp"          , "t#bar{t}",     r.itSys, 1, "fsrUp");
    p.AddSample("TTbar_Powheg_fsrDown"        , "t#bar{t}",     r.itSys, 1, "fsrDown");
    
    specialweight = vl.nGluonMoveCRTune_ttbar/vl.sigma_ttbar/(vl.nGluonMoveCRTune_ttbar/vl.sigma_ttbar + vl.nGluonMoveCRTune_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight *' + str(specialweight))
    p.AddSample('TTbar_GluonMoveCRTune',        't#bar{t}',     r.itSys, 1, "GluonMoveCRTune")
    specialweight = vl.nGluonMoveCRTune_dilep/vl.sigma_dilep/(vl.nGluonMoveCRTune_ttbar/vl.sigma_ttbar + vl.nGluonMoveCRTune_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight *' + str(specialweight))
    p.AddSample('TTTo2L2Nu_GluonMoveCRTune',    't#bar{t}',     r.itSys, 1, "GluonMoveCRTune")
    specialweight = vl.nPowhegerdON_ttbar/vl.sigma_ttbar/(vl.nPowhegerdON_ttbar/vl.sigma_ttbar + vl.nPowhegerdON_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight *' + str(specialweight))
    p.AddSample('TTbar_Powheg_erdON',           't#bar{t}',     r.itSys, 1, "Powheg_erdON")
    specialweight = vl.nPowhegerdON_dilep/vl.sigma_dilep/(vl.nPowhegerdON_ttbar/vl.sigma_ttbar + vl.nPowhegerdON_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight *' + str(specialweight))
    p.AddSample('TTTo2L2Nu_Powheg_erdON',       't#bar{t}',     r.itSys, 1, "Powheg_erdON")
    specialweight = vl.nQCDbasedCRTuneerdON_ttbar/vl.sigma_ttbar/(vl.nQCDbasedCRTuneerdON_ttbar/vl.sigma_ttbar + vl.nQCDbasedCRTuneerdON_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight *' + str(specialweight))
    p.AddSample('TTbar_GluonMoveCRTune_erdON',  't#bar{t}',     r.itSys, 1, "QCDbasedCRTune_erdON")
    specialweight = vl.nQCDbasedCRTuneerdON_dilep/vl.sigma_dilep/(vl.nQCDbasedCRTuneerdON_ttbar/vl.sigma_ttbar + vl.nQCDbasedCRTuneerdON_dilep/vl.sigma_dilep)
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight) *' + str(specialweight))   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight *' + str(specialweight))
    p.AddSample('TTTo2L2Nu_QCDbasedCRTune_erdON','t#bar{t}',    r.itSys, 1, "QCDbasedCRTune_erdON")
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    else:            p.SetWeight('TWeight')
    
    p.AddSample("TTbar_GluonMoveCRTune_erdON" , "t#bar{t}",     r.itSys, 1, "GluonMoveCRTune_erdON");
    p.UseEnvelope("t#bar{t}", "GluonMoveCRTune,GluonMoveCRTune_erdON,Powheg_erdON,QCDbasedCRTune_erdON", "ColorReconnection");
    p.AddSymmetricHisto("t#bar{t}",  "JERUp");
    
    pdf     = r.PDFToPyC(r.TString(pathToTree), r.TString("TTbar_Powheg"), r.TString(NameOfTree), r.TString("TIsSS == 0 && TNJets == 1  && TNBtags == 1"), r.TString("All"), r.TString(vl.varList[var]['var'] + "- 1"), len(vl.varList[var]['recobinning']) - 1, binning, r.TString(''));
    pdf.verbose = False
    #pdf.verbose = True
    pdf.SetLumi(vl.Lumi * 1000)
    if doPrefChecks: pdf.SetWeight("TWeight * (1 - prefWeight)")
    
    hPDFUp  = pdf.GetSystHisto("up","pdf").CloneHisto();
    hPDFDown= pdf.GetSystHisto("Down","pdf").CloneHisto();
    hMEUp   = pdf.GetSystHisto("up","ME").CloneHisto();
    hMEDown = pdf.GetSystHisto("Down","ME").CloneHisto();
    p.PrepareHisto(hPDFUp,   "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarPDFUp");
    p.PrepareHisto(hPDFDown, "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarPDFDown");
    p.PrepareHisto(hMEUp,    "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarMEUp");
    p.PrepareHisto(hMEDown,  "TTbar_Powheg", "t#bar{t}", r.itSys, 0, "ttbarMEDown");
    p.AddToSystematicLabels("ttbarPDF");
    p.AddToSystematicLabels("ttbarME");
    del pdf
    
    # Other settings
    p.SetDataStyle("psameE1")
    p.doUncInLegend = True
    p.SetRatioMin( 0.6 );
    p.SetRatioMax( 1.4 );
    p.SetPadPlotMargins(vl.margins)
    p.SetPadRatioMargins(vl.marginsratio)
    p.SetTexChanSize(0.06)
    p.SetTextLumiPosX(0.69)
    p.SetYratioOffset(0.35)
    p.SetCenterYAxis(False)
    p.SetXaxisDivisions(010)
    p.ObliterateXErrorBars()
    p.SetBinLabels("0,1,2,3,#geq4")
    p.SetXaxisLabelSize(0.16 / 0.66666) # This is because ROOT automatically changes the size of labels when transforming them to alphanumeric because of reasons
    #r.gPad.Update()
    if 'ncols' in vl.varList[var]: p.SetNColumns(vl.varList[var]['ncols'])
    
    p.SetCMSlabel("CMS");
    if vl.doPre: p.SetCMSmodeLabel("Preliminary");
    else:        p.SetCMSmodeLabel("");
    if 'legpos' in vl.varList[var]: thepos = vl.varList[var]['legpos']
    else:                           thepos = vl.legpos
    p.SetLegendPosition(thepos[0], thepos[1], thepos[2], thepos[3])
    p.SetLegendTextSize(legtxtsize)
    p.SetPlotFolder("/nfs/fanae/user/vrbouza/www/TFM/1j1t/" if cut == 'signal' else "/nfs/fanae/user/vrbouza/www/TFM/1j1t/control/");
    p.doYieldsInLeg = False;
    if "doLogY" in vl.varList[var]: p.doSetLogy = vl.varList[var]['doLogY']
    else:                           p.doSetLogy = False
    p.doSignal      = False;
    #p.SetTitleY(r.TString(vl.varList[var]['yaxis']))
    p.SetOutputName("Custom");
    if "abs" in vl.varList[var]['var'] or "min" in vl.varList[var]['var']:
        p.NoShowVarName = True;
        p.SetOutputName('Custom_' + vl.varList[var]['var_response'] + "_" + chan);
    p.DrawStack();
    p.PrintSystematics()
    p.PrintYields("", "", "", "")
    p.PrintSystYields()
    del p


def lazyoptimisation(tsk):
    var, reg, chn, bnng, lvl = tsk
    if   lvl  == "particle":      return plotdressvariable(tsk[:-2])
    elif var  == "nLooseCentral": return plotthenumberofjets(tsk[:-2])
    elif bnng == "custom":        return plotcustomvariable(tsk[:-2])
    else:                         return plotvariable(tsk[:-2])
    return


if __name__ == '__main__':
    tasks = []
    #tasks.append( ("nLooseCentral", "signal", "wololo") )
    #for v in vl.varList["Names"]["Variables"]:
        #if "Fiducial" in v: continue
        #for ct in ['signal', 'control']:
            #for bnng in ['custom', 'descriptive']:
                #for chn in ["ElMu", "Muon", "Elec", "All"]:
                    #tasks.append( (v, ct, chn, bnng) )
    
    for chn in ["ElMu", "Muon", "Elec", "All"]:
      for lvl in ["detector", "particle"]:
        tasks.append(("Lep1_Pt",       "signal", chn, "descriptive", lvl))
        tasks.append(("Lep1_Eta",      "signal", chn, "descriptive", lvl))
        tasks.append(("Lep2_Pt",       "signal", chn, "descriptive", lvl))
        tasks.append(("Lep2_Eta",      "signal", chn, "descriptive", lvl))
        tasks.append(("LepMuon_Pt",    "signal", chn, "descriptive", lvl))
        tasks.append(("LepMuon_Eta",   "signal", chn, "descriptive", lvl))
        tasks.append(("LepElec_Pt",    "signal", chn, "descriptive", lvl))
        tasks.append(("LepElec_Eta",   "signal", chn, "descriptive", lvl))
        tasks.append(("Jet1_Pt",       "signal", chn, "descriptive", lvl))
        tasks.append(("Jet1_Eta",      "signal", chn, "descriptive", lvl))
        tasks.append(("Jet2_Pt",       "signal", chn, "descriptive", lvl))
        tasks.append(("Jet2_Eta",      "signal", chn, "descriptive", lvl))
        tasks.append(("MiniMax",       "signal", chn, "descriptive", lvl))
        tasks.append(("Lep1Lep2_DPhi", "signal", chn, "descriptive", lvl))
    
    if (nCores > 1):
      print "\n> Launching plotting processes..."
      pool = Pool(nCores)
      pool.map(lazyoptimisation, tasks)
      pool.close()
      pool.join()
      del pool
    else:
      print "\n> Initiating sequential execution..."
      for tsk in tasks: lazyoptimisation(tsk)

    print "> Done!", "\n"
