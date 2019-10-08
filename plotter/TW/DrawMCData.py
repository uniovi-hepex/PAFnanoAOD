import ROOT     as r
import varList  as vl
import sys, os, copy
from multiprocessing import Pool
from array import array
print "===== Minitrees MC/DATA graphs plotting\n"
vl.SetUpWarnings()

cutdir            = {}
cutdir["SR"]      = "TPassReco == 1"
cutdir["CR"]      = "TNJets == 1  && TNBJets == 1 && TNLooseCentral > 1"
cutdir["2j1b"]    = "TNJets == 2  && TNBJets == 1"
cutdir["Dress2j1b"]="TDressNJets == 2 && TDressNBJets == 1 "
cutdir["DressSR"] = "TPassDress == 1"
cutdir["1j1b"]    = "TNJets == 1  && TNBJets == 1"
#StandardCut       = "TPassReco == 1";
#ControlCut        = "TNJets == 1  && TNBtags == 1 && TNLooseCentral > 1";
#StandardDressCut  = "TPassDress == 1";

#systlist          = vl.GiveMeTheExpNamesWOJER(vl.varList["Names"]["ExpSysts"])
systlist          = ""
systdresslist     = ""

labeldir            = {}
#labeldir["SR"]      ="+1j1b+0j_{loose}"
labeldir["SR"]      = "+1j1b"
labeldir["CR"]      = "+1j1b+>0j_{loose}"
labeldir["1j1b"]    = "+1j1b"
labeldir["2j1b"]    = "+2j1b"
labeldir["DressSR"] = "+1j1b"
labeldir["DressCR"] = "+1j1b+>0j_{loose}"
labeldir["Dress1j1b"] = "+1j1b"
labeldir["Dress2j1b"] = "+2j1b"

chandir           = {}
chandir["All"]    = "\\ell_{1}^{\\pm}\\ell_{2}^{\\mp}"
chandir["ElMu"]   = "e^{#pm}#mu^{#mp}"
chandir["Elec"]   = "e^{#pm}e^{#mp}"
chandir["Muon"]   = "#mu^{#pm}#mu^{#mp}"

folderdir         = {}
folderdir["SR"]   = "./results/MCData/"
folderdir["CR"]   = "./results/MCData/control/"
folderdir["1j1b"] = "./results/MCData/"
folderdir["2j1b"] = "./results/MCData/2j1b/"
folderdir["DressSR"]   = "./results/MCData/"
folderdir["DressCR"]   = "./results/MCData/control/"
folderdir["Dress2j1b"] = "./results/MCData/2j1b/"

#legtxtsize  = 0.028
legtxtsize  = 0.055
labelpos    = (0.275, 0.89)
doPrefChecks= False
doNorm      = False
verbosity   = False
NameOfTree  = vl.treename
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


def addnominalsamples(plotobj, year, adddata = True):
    if year == "2016":
        plotobj.AddSample("TTToSemiLeptonic_*",           "Non-W/Z",      r.itBkg, 413, systlist)
        plotobj.AddSample("WJetsToLNu_MLM_*",             "Non-W/Z",      r.itBkg, 413, systlist)

        plotobj.AddSample("WZ_*",                         "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("WW_*",                         "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("ZZ_*",                         "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTWJetsToLNu_*",               "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTWJetsToQQ_*",                "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTZToQQ_*",                    "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTZToLLNuNu_M_10_*",           "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTZToLL_M_1to10_*",            "VV+t#bar{t}V", r.itBkg, 390, systlist)

        #plotobj.AddSample("DYJetsToLL_M10to50_aMCatNLO",  "DY",          r.itBkg, 852, systlist);
        #plotobj.AddSample("DYJetsToLL_M50_aMCatNLO",      "DY",          r.itBkg, 852, systlist);
        plotobj.AddSample("DYJetsToLL_M_10to50_MLM_*",    "DY",           r.itBkg, 852, systlist);
        plotobj.AddSample("DYJetsToLL_M_50_MLM_*",        "DY",           r.itBkg, 852, systlist);

        plotobj.AddSample("TTTo2L2Nu_*",                  "t#bar{t}",     r.itBkg, 633, systlist)

        plotobj.AddSample("tbarW_noFullHad_*",            "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
        plotobj.AddSample("tW_noFullHad_*",               "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systlist);

        if adddata:
            plotobj.AddSample("MuonEG_*",                 "Data",         r.itData);
            plotobj.AddSample("SingleMuon_*",             "Data",         r.itData);
            plotobj.AddSample("SingleElec_*",             "Data",         r.itData);

    elif year == "2017":
        plotobj.AddSample("TTToSemiLeptonic_*",           "Non-W/Z",      r.itBkg, 413, systlist)
        plotobj.AddSample("WJetsToLNu_MLM_*",             "Non-W/Z",      r.itBkg, 413, systlist)

        plotobj.AddSample("WZ_*",                         "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("WWTo2L2Nu_*",                  "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("ZZTo2L2Nu_*",                  "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTWJetsToLNu_*",               "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTWJetsToQQ_*",                "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTZToQQ_*",                    "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTZToLLNuNu_M_10_*",           "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTZToLL_M_1to10_*",            "VV+t#bar{t}V", r.itBkg, 390, systlist)

        #plotobj.AddSample("DYJetsToLL_M_10to50_aMCatNLO_*","DY",          r.itBkg, 852, systlist);
        plotobj.AddSample("DYJetsToLL_M_50_aMCatNLO_*",   "DY",           r.itBkg, 852, systlist);
        plotobj.AddSample("DYJetsToLL_M_10to50_MLM_*",    "DY",           r.itBkg, 852, systlist);
        #plotobj.AddSample("DYJetsToLL_M_50_MLM_*",        "DY",           r.itBkg, 852, systlist);

        plotobj.AddSample("TTTo2L2Nu_*",                  "t#bar{t}",     r.itBkg, 633, systlist)

        plotobj.AddSample("tbarW_noFullHad_*",            "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
        plotobj.AddSample("tW_noFullHad_*",               "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systlist);

        if adddata:
            plotobj.AddSample("MuonEG_*",                 "Data",         r.itData);
            plotobj.AddSample("SingleMuon_*",             "Data",         r.itData);
            plotobj.AddSample("SingleElec_*",             "Data",         r.itData);

    else:
        plotobj.AddSample("TTToSemiLeptonic_*",           "Non-W/Z",      r.itBkg, 413, systlist)
        plotobj.AddSample("WJetsToLNu_MLM_*",             "Non-W/Z",      r.itBkg, 413, systlist)

        plotobj.AddSample("WZ_*",                         "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("WW_*",                         "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("ZZ_*",                         "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTWJetsToLNu_*",               "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTWJetsToQQ_*",                "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTZToQQ_*",                    "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTZToLLNuNu_M_10_*",           "VV+t#bar{t}V", r.itBkg, 390, systlist)
        plotobj.AddSample("TTZToLL_M_1to10_*",            "VV+t#bar{t}V", r.itBkg, 390, systlist)

        #plotobj.AddSample("DYJetsToLL_M10to50_aMCatNLO",  "DY",          r.itBkg, 852, systlist);
        #plotobj.AddSample("DYJetsToLL_M50_aMCatNLO",      "DY",          r.itBkg, 852, systlist);
        plotobj.AddSample("DYJetsToLL_M_10to50_MLM_*",    "DY",           r.itBkg, 852, systlist);
        plotobj.AddSample("DYJetsToLL_M_50_MLM_*",        "DY",           r.itBkg, 852, systlist);

        plotobj.AddSample("TTTo2L2Nu_*",                  "t#bar{t}",     r.itBkg, 633, systlist)

        plotobj.AddSample("tbarW_noFullHad_*",            "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systlist)
        plotobj.AddSample("tW_noFullHad_*",               "tW",           r.itBkg, r.TColor.GetColor("#ffcc33"), systlist);

        if adddata:
            plotobj.AddSample("MuonEG_*",                 "Data",         r.itData);
            plotobj.AddSample("SingleMuon_*",             "Data",         r.itData);
            plotobj.AddSample("EGamma_*",                 "Data",         r.itData);
    return


def setstylesettings(plotobj, var):
    plotobj.SetCanvasHeight(600)
    plotobj.SetCanvasWidth(600)

    plotobj.doSignal = False;

    plotobj.SetDataStyle("psameE1")
    #plotobj.SetSignalStyle("SM")
    #plotobj.SetSignalStyle("scan")

    plotobj.doUncInLegend = True;

    plotobj.SetRatioMin( 0.6 );
    plotobj.SetRatioMax( 1.4 );
    plotobj.SetPadPlotMargins(vl.margins)
    plotobj.SetPadRatioMargins(vl.marginsratio)
    plotobj.SetTexChanSize(0.06)
    plotobj.SetTextLumiPosX(0.69)

    #plotobj.SetYratioOffset(0.35)
    plotobj.SetYratioOffset(0.45)
    #plotobj.SetRatioStyle("S/B")
    if doNorm:
        plotobj.SetRatioMin(0.7)
        plotobj.SetRatioMax(1.3)
    #else:
        #plotobj.SetRatioMin(0)
        #plotobj.SetRatioMax(0.2)

    plotobj.SetCenterYAxis(False)
    plotobj.SetXaxisOffset(1.1)
    plotobj.ObliterateXErrorBars()
    if 'ncols' in vl.varList[var]: plotobj.SetNColumns(vl.varList[var]['ncols'])

    plotobj.SetCMSlabel("CMS");
    if vl.doPre: plotobj.SetCMSmodeLabel("Preliminary");
    else:        plotobj.SetCMSmodeLabel("");

    if 'legposdesc' in vl.varList[var]: thepos = vl.varList[var]['legposdesc']
    else:                               thepos = vl.legpos
    plotobj.SetLegendPosition(thepos[0], thepos[1], thepos[2], thepos[3])
    plotobj.SetLegendTextSize(legtxtsize)
    plotobj.doYieldsInLeg = False;

    if "doLogY" in vl.varList[var]: plotobj.doSetLogy = vl.varList[var]['doLogY']
    else:                           plotobj.doSetLogy = False
    if doNorm:                      plotobj.doSetLogy = False
    return


def plotvariable(tsk):
    var, cut, chan, year = tsk
    nbins    = int(20) if "ndescbins" not in vl.varList[var] else int(vl.varList[var]['ndescbins'])
    lowedge  = float(vl.varList[var]['recobinning'][0])  if "descbinning" not in vl.varList[var] else float(vl.varList[var]['descbinning'][0])
    highedge = float(vl.varList[var]['recobinning'][-1]) if "descbinning" not in vl.varList[var] else float(vl.varList[var]['descbinning'][1])
    width    = (highedge - lowedge)/nbins
    lumi_    = vl.Lumi16 if (year == "2016") else vl.Lumi17 if (year == "2017") else vl.Lumi18 if (year == "2018") else vl.Lumi
    
    p = r.PlotToPy(r.TString(vl.varList[var]['var']), r.TString(cutdir[cut]), r.TString(chan), nbins, lowedge, highedge, r.TString(var), r.TString(vl.varList[var]['xaxis']))
    p.SetPath(pathToTree + year + "/");
    p.SetPathSignal(pathToTree + year + "/");
    p.SetTreeName(NameOfTree);
    p.SetPlotFolder(folderdir[cut] + year + "/");
    p.SetTitleY("Events / " + str(int(round(width, 0))) + " GeV" if ("Eta" not in var and "DPhi" not in var and "DT" not in var) else "Events / bin" )
    p.SetLumi(lumi_)
    p.verbose  = verbosity
    
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    
    addnominalsamples(p, year)

    # Labels
    p.SetChLabel(chandir[chan] + labeldir[cut])
    p.SetChLabelPos(labelpos[0], labelpos[1], -1)

    # General style settings
    setstylesettings(p, var)

    # Output name
    p.NoShowVarName = True;
    p.SetOutputName(vl.varList[var]['var_response'] + "_" + chan + "_norm" * doNorm);

    # Printing
    p.DrawStack("", doNorm);
    p.PrintSystematics()
    p.PrintYields("", "", "", "")
    p.PrintSystYields()
    del p
    #del pdf


def plotdressvariable(tsk):
    var, cut, chan, year = tsk
    nbins    = int(20) if "ndescbins" not in vl.varList[var] else int(vl.varList[var]['ndescbins'])
    lowedge  = float(vl.varList[var]['genbinning'][0])  if "descbinning" not in vl.varList[var] else float(vl.varList[var]['descbinning'][0])
    highedge = float(vl.varList[var]['genbinning'][-1]) if "descbinning" not in vl.varList[var] else float(vl.varList[var]['descbinning'][1])
    width    = (highedge - lowedge)/nbins

    p = r.PlotToPy(r.TString(vl.varList[var]['var_gen']), r.TString(cutdir["Dress" + cut]), r.TString(chan), nbins, lowedge, highedge, r.TString(var), r.TString(vl.varList[var]['xaxis'] if "xaxis_dress" not in vl.varList[var] else vl.varList[var]['xaxis_dress']))
    p.SetPath(pathToTree + year + "/");
    p.SetPathSignal(pathToTree + year + "/");
    p.SetTreeName(NameOfTree);
    p.SetPlotFolder(folderdir[cut] + year + "/");
    p.SetTitleY("Events / " + str(int(round(width, 0))) + " GeV" if ("Eta" not in var and "DPhi" not in var and "DT" not in var) else "Events / bin" )
    p.SetLumi(vl.Lumi)
    p.verbose  = verbosity

    p.SetWeight('TWeight_normal')

    # Samples
    addnominalsamples(p, year, False)

    # Labels
    p.SetChLabel(chandir[chan] + labeldir["Dress" + cut])
    p.SetChLabelPos(labelpos[0], labelpos[1], -1)

    # General style settings
    setstylesettings(p, var)

    # Output name
    p.NoShowVarName = True;
    p.SetOutputName("Dress" + vl.varList[var]['var_response'] + "_" + chan + "_norm" * doNorm);

    # Printing
    p.DrawStack("", doNorm);
    p.PrintSystematics()
    p.PrintYields("", "", "", "")
    p.PrintSystYields()
    del p


def plotthenumberofjets(tsk):
    var, cut, chan, year = tsk
    binning = array('f', vl.varList[var]['recobinning']) # For some reason, ROOT requires that you create FIRST this object, then put it inside the PlotToPyC.
    p = r.PlotToPyC(r.TString(vl.varList[var]['var'] + "- 1"), r.TString(cutdir[cut]), r.TString(chan), int(len(vl.varList[var]['recobinning']) - 1), binning, r.TString(var), r.TString(vl.varList[var]['xaxis']))
    p.SetPath(pathToTree + year + "/");
    p.SetPathSignal(pathToTree + year + "/");
    p.SetTreeName(NameOfTree);
    p.SetPlotFolder(folderdir[cut] + year + "/");
    p.SetTitleY("Events / bin")
    p.SetLumi(vl.Lumi)
    p.verbose  = verbosity
    
    if doPrefChecks: p.SetWeight('TWeight * (1 - prefWeight)')   # FOR PREFIRING CHECKS
    
    # Samples
    addnominalsamples(p, year)

    # Labels
    p.SetChLabel(chandir[chan] + labeldir[cut])
    p.SetChLabelPos(labelpos[0], labelpos[1], -1)

    # General style settings
    setstylesettings(p, var)

    # Other settings
    #p.SetDataStyle("psameE1")
    #p.doUncInLegend = True
    #p.SetRatioMin(0.6);
    #p.SetRatioMax(1.4);
    #p.SetPadPlotMargins(vl.margins)
    #p.SetPadRatioMargins(vl.marginsratio)
    #p.SetTexChanSize(0.06)
    #p.SetTextLumiPosX(0.69)
    ##p.SetYratioOffset(0.35)
    #p.SetYratioOffset(0.45)
    ##p.SetRatioYtitle("WWbb/tW+t#bar{t}")
    ##p.SetRatioStyle("S/B")
    #if doNorm:
        #p.SetRatioMin(0.7)
        #p.SetRatioMax(1.3)
    ##else:
        ##p.SetRatioMin(0)
        ##p.SetRatioMax(0.2)
    #p.SetCenterYAxis(False)
    #p.ObliterateXErrorBars()
    #if 'ncols' in vl.varList[var]: p.SetNColumns(vl.varList[var]['ncols'])
    
    #p.SetCMSlabel("CMS");
    #if vl.doPre: p.SetCMSmodeLabel("Preliminary");
    #else:        p.SetCMSmodeLabel("");
    #if 'legpos' in vl.varList[var]: thepos = vl.varList[var]['legpos']
    #else:                           thepos = vl.legpos
    #p.SetLegendPosition(thepos[0], thepos[1], thepos[2], thepos[3])
    #p.SetLegendTextSize(legtxtsize)
    #p.SetPlotFolder(folderdir[cut]);
    #p.doYieldsInLeg = False;
    #if "doLogY" in vl.varList[var]: p.doSetLogy = vl.varList[var]['doLogY']
    #else:                           p.doSetLogy = False
    #p.doSignal      = False;

    # X axis labels
    p.SetXaxisDivisions(010)
    p.SetBinLabels("0,1,2,3,#geq4")
    p.SetXaxisLabelSize(0.16 / 0.66666) # This is because ROOT automatically changes the size of labels when transforming them to alphanumeric because of reasons

    # Output name
    p.NoShowVarName = True;
    p.SetOutputName(vl.varList[var]['var'] + "_" + chan);

    # Printing
    p.DrawStack();
    p.PrintSystematics()
    p.PrintYields("", "", "", "")
    p.PrintSystYields()
    del p


def lazyoptimisation(tsk):
    var, reg, chn, year, bnng, lvl = tsk
    if   lvl  == "particle":      return plotdressvariable(tsk[:-2])
    elif var  == "NLooseCentral" or var == "NBLooseCentral": return plotthenumberofjets(tsk[:-2])
    elif bnng == "custom":        return plotcustomvariable(tsk[:-2])
    else:                         return plotvariable(tsk[:-2])
    return



if __name__ == '__main__':
    vl.SetUpWarnings()
    tasks = []

    #for chn in ["ElMu", "Muon", "Elec", "All"]:
    for chn in ["ElMu"]:
        for yr in ["2016", "2017", "2018"]:
            for lvl in ["detector", "particle"]:
                tasks.append(("Lep1_Pt",       "SR",   chn, yr, "descriptive", lvl))
                tasks.append(("Lep1_Eta",      "SR",   chn, yr, "descriptive", lvl))
                tasks.append(("Lep2_Pt",       "SR",   chn, yr, "descriptive", lvl))
                tasks.append(("Lep2_Eta",      "SR",   chn, yr, "descriptive", lvl))
                tasks.append(("LepMuon_Pt",    "SR",   chn, yr, "descriptive", lvl))
                tasks.append(("LepMuon_Eta",   "SR",   chn, yr, "descriptive", lvl))
                tasks.append(("LepElec_Pt",    "SR",   chn, yr, "descriptive", lvl))
                tasks.append(("LepElec_Eta",   "SR",   chn, yr, "descriptive", lvl))
                tasks.append(("Jet1_Pt",       "SR",   chn, yr, "descriptive", lvl))
                tasks.append(("Jet1_Eta",      "SR",   chn, yr, "descriptive", lvl))
                tasks.append(("Lep1Lep2_DPhi", "SR",   chn, yr, "descriptive", lvl))
                tasks.append(("Sys_Pt",        "SR",   chn, yr, "descriptive", lvl))
                tasks.append(("Sys_M",         "SR",   chn, yr, "descriptive", lvl))
                tasks.append(("Jet2_Pt",       "2j1b", chn, yr, "descriptive", lvl))
            tasks.append(("HTtot",         "SR",   chn, yr, "descriptive", "detector"))
            tasks.append(("C_jll",         "SR",   chn, yr, "descriptive", "detector"))
            tasks.append(("Lep1Lep2_HTOverHTtot","SR",chn,yr,"descriptive","detector"))
            tasks.append(("LooseCentral1_Pt","SR", chn, yr, "descriptive", "detector"))
            tasks.append(("Sys_PtOverHTtot","SR",  chn, yr, "descriptive", "detector"))
            tasks.append(("Lep1Lep2Jet1_Pt","SR",  chn, yr, "descriptive", "detector"))
            tasks.append(("BDT",           "SR",   chn, yr, "descriptive", "detector"))
            tasks.append(("NLooseCentral", "1j1b", chn, yr, "descriptive", "detector"))
            tasks.append(("NBLooseCentral","1j1b", chn, yr, "descriptive", "detector"))
            tasks.append(("Lep1Jet1_DR",   "2j1b", chn, yr, "descriptive", "detector"))
            tasks.append(("Lep12Jet12_DR", "2j1b", chn, yr, "descriptive", "detector"))
            tasks.append(("Lep12Jet12MET_DR","2j1b",chn,yr, "descriptive", "detector"))
            tasks.append(("BDT2j1b",       "2j1b", chn, yr, "descriptive", "detector"))
            tasks.append(("Jet1_Eta",      "SR",   chn, yr, "descriptive", "detector"))

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
