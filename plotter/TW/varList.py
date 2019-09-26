'''
Library of all the variables so everything is centralised and there are not confusions
'''

import ROOT as r
import warnings as wr
import os

# === ESSENTIAL PARAMETERS OF THE ANALYSIS. CHANGING THIS APPLIES TO EVERYTHING. ===
nuncs       = 4         # Number of uncs. shown in the relative uncertainty plots
asimov      = False     # Use of Asimov dataset or data
doxsec      = True      # Show events or diff. cross section in final results
doPre       = True      # Show "Preliminary" (or not) after CMS in the plots

# === OTHER IMPORTANT DEFINITIONS ===
Lumi16      = 35.92    # In femtobarns
Lumi17      = 41.53    # In femtobarns
Lumi18      = 59.74    # In femtobarns
uncLumi     = 0.025    # In %

Lumi        = Lumi18   # TEMPORAAAAAAAAAAL

sigma_ttbar   = 831.76
sigma_dilep   = 88.28769753
sigma_tw      = 35.85
sigma_twnohad = 19.4674104

n_ttbar                    = 77229341
n_dilep                    = 79092400

nUEUp_ttbar                = 58953660
nUEUp_dilep                = 9907537
nUEDown_ttbar              = 58338240
nUEDown_dilep              = 9691700
nhDampUp_ttbar             = 58858606
nhDampUp_dilep             = 9672473
nhDampDown_ttbar           = 58163976
nhDampDown_dilep           = 9634312

nGluonMoveCRTune_ttbar     = 59037234
nGluonMoveCRTune_dilep     = 9862990
nPowhegerdON_ttbar         = 59882210
nPowhegerdON_dilep         = 9537400
nQCDbasedCRTuneerdON_ttbar = 59620206
nQCDbasedCRTuneerdON_dilep = 9816448

n_tw                       = 6952830
n_twnohad                  = 11345619
n_tbarw                    = 6933094
n_tbarwnohad               = 11408144

plotlimits   = tuple([float(i) for i in "0.00, 0.25, 1.00, 1.00".split(',')]) # xlow, ylow, xup, yup
ratiolimits  = tuple([float(i) for i in "0.00, 0.00, 1.00, 0.25".split(',')]) # xlow, ylow, xup, yup
margins      = "0.06, 0.1, 0.04, 0.1" # top, bottom, right, left
marginsratio = "0.03, 0.4, 0.04, 0.1" # top, bottom, right, left
legpos       = (0.82, 0.65, 0.93, 0.93)

if asimov: labellegend = 'Pseudodata'
else:      labellegend = 'Data'

storagepath = "/pool/ciencias/userstorage/vrbouza/proyectos/TW_inclusivo_run2/MiniTrees/"
minipath    = "../../temp_TW/"
treename    = "fMiniTree"

def GetLastFolder(stpth):
    savefolders   = next(os.walk(stpth))[1]
    saveyears     = map(int, [i[:4]  for i in savefolders])
    savefolders   = [i for i in savefolders if int(i[:4]) == max(saveyears)]
    savemonths    = map(int, [i[5:7] for i in savefolders])
    savefolders   = [i for i in savefolders if int(i[5:7]) == max(savemonths)]
    savedays      = map(int, [i[8:]  for i in savefolders])
    savefolders   = [i for i in savefolders if int(i[8:]) == max(savedays)]
    return (stpth + savefolders[0] + "/")


def SetUpWarnings():
    wr.simplefilter("always", UserWarning)
    wr.filterwarnings(action = 'ignore', category = RuntimeWarning, message = 'TClass::Init:0: RuntimeWarning: no dictionary for class')
    return

def mean(numbers):
    return float(sum(numbers)) / max(len(numbers), 1)

def GiveMeTheExpNamesWOJER(inl):
    l = []
    for el in inl:
        if "JER" in el or "Down" in el: continue
        l.append(el.replace("Up", ""))
    return ",".join(l)



# === DICTIONARIES ===
# var           := name of the variable in the tree to make cards
# var_response  := name of the variable in the response matrix without the M
varList = {}
varList['LCurve'] = {
    'xaxis'       : 'log L_{1}',
    'yaxis'       : 'log \\frac{L_{2}}{\\tau^{2}}',
}

varList['Names'] = {
    'Variables'   : ["Lep1_Pt", "Jet1_Pt"],
    'ExpSysts'    : ["JESUp", "JESDown", "JERUp", "ElecEffUp", "ElecEffDown", "MuonEffUp", #   DO NOT MOVE THE FIRST THREE TO OTHER
                     "MuonEffDown", "TrigUp", "TrigDown", "PUUp", "PUDown", "BtagUp",      # POSITION: it will affect the calculus
                     #"BtagDown", "MistagUp", "MistagDown"],                                # of the response matrices.
                     "BtagDown", "MistagUp", "MistagDown", "TopPtUp", "TopPtDown"],        # of the response matrices.  WITH TOP PT REW
    #'ExpSysts'    : ["JESUp", "JESDown", "JERUp", "ElecEffUp", "ElecEffDown", "MuonEffUp", # WITH TOP PT REW.
                     #"MuonEffDown", "TrigUp", "TrigDown", "PUUp", "PUDown", "BtagUp",
                     #"BtagDown", "MistagUp", "MistagDown", "TopPtUp", "TopPtDown"],
    'ttbarSysts'  : ["ttbarMEUp", "ttbarMEDown", "pdfUp", "pdfDown", "hDampUp", "hDampDown", "UEUp", "UEDown"],
    'specialSysts': ["JERDown", "DSDown"],
    'colorSysts'  : ["GluonMoveCRTuneerdON", "PowhegerdON", "QCDbasedCRTuneerdON", "GluonMoveCRTune"],
    'NormSysts'   : ["ttbarUp", "ttbarDown", "Non-WorZUp", "Non-WorZDown", "DYUp", "DYDown", "VVttbarVUp", "VVttbarVDown"],
}

varList['E_LLB'] = {
    'xaxis'       : 'E(\\ell_{1}, \\ell_{2}, j) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0, 190, 330, 550, 700],
    'recobinning' : [0, 160, 220, 280, 340, 400, 450, 550, 700],
    'var'         : 'TE_LLB',
    'var_response': 'ELLB',
    'var_gen'     : 'TGenE_LLB',
}
varList['E_LLBuncertainties'] = {
    'xaxis'       : varList['E_LLB']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}

varList['LeadingJetE'] = {
    'xaxis'       : 'E(j) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0, 75, 275, 400],
    'recobinning' : [0., 40., 70., 120., 175., 275., 400.],
    'var'         : 'TLeadingJetE',
    'var_response': 'LeadingJetE',
    'var_gen'     : 'TGenLeadingJetE',
}
varList['LeadingJetEuncertainties'] = {
    'xaxis'       : varList['LeadingJetE']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}

varList['MT_LLMETB'] = {
    #'xaxis'       : 'm_{T}(\\ell_{1}, \\ell_{2},\\slash{E}_{T}, j) (GeV)',
    'xaxis'       : '#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j}) (GeV)',
    'yaxis'       : 'd#sigma/d(#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j})) (pb)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j})) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j})) (1/GeV)',
    #'genbinning'  : [0., 225., 325., 425., 500.],                         # usado en drafts, en nota, etc. (4 bins)
    #'recobinning' : [0., 225., 250., 275., 300., 325., 350., 425., 500.], # usado en drafts, en nota, etc. (4 bins)

    'genbinning'  : [0., 220., 300., 325., 400., 500.], # propuesta (5 bins)
    'recobinning' : [0., 220., 250., 260., 270., 280., 290., 300., 325., 400., 500.], # propuesta (5 bins) SELECCIONADA

    #'genbinning'  : [0., 220., 270., 300., 325., 400., 500.], # propuesta (6 bins)
    #'recobinning' : [0., 220., 240., 250., 260., 270., 280., 290., 300., 310., 325., 400., 500.], # propuesta (6 bins)

    #'genbinning'  : [0., 220., 260., 280., 300., 325., 400., 500.], # propuesta (7 bins)
    #'recobinning' : [0., 220., 240., 250., 260., 270., 275., 280., 285., 290., 300., 310., 325., 400., 500.], # propuesta (7 bins)

    #'genbinning'  : [0., 220., 240., 270., 300., 320., 350., 430., 500.], # propuesta (8 bins)
    #'recobinning' : [0., 190., 220., 230., 240., 250., 260., 270., 280., 290., 300., 310., 320., 330., 350., 430., 500.], # propuesta (8 bins)

    #'genbinning'  : [0., 220., 240., 250., 260., 270., 280., 320., 400., 500.], # propuesta (9 bins)
    #'recobinning' : [0., 220., 230., 240., 250., 255., 260., 265., 270., 275., 280., 285., 290., 300., 310., 320., 330., 400., 500.], # propuesta (9 bins)

    #'genbinning'  : [0., 220., 240., 250., 260., 270., 280., 300., 320., 400., 500.], # propuesta (10 bins)
    #'recobinning' : [0., 220., 230., 240., 250., 255., 260., 265., 270., 275., 280., 285., 290., 295., 300., 305., 310., 320., 330., 400., 500.], # propuesta (10 bins)

    #'descbinning' : [0., 800.], ## TEMPORAL
    #'ndescbins'   : 20, ## TEMPORAL
    'var'         : 'min(TMT_LLMETB, 499.)',
    'var_response': 'MTLLMETB',
    'var_gen'     : 'min(TGenMT_LLMETB, 499.)',
    'legpos'      : (0.51, 0.55, 0.71, 0.93),
    'legposdesc'  : (0.15, 0.425, 0.35, 0.81),
    'legpos_foldas':"BL",
    'legpos_fold' : "BL",
    'legpos_fid'  : "TL",
    "legpos_unf"   : (.18, .65, .38, .40),
    'legpos_unfas': "TL",
    'uncleg_fold' : "TL",
    'uncleg_unf'  : "TC",
    'uncleg_fid'  : "TL",
    'legpos_fidbin':"ML",
    'uncleg_fidbin':"TL",
    'resptxtsize' : 0.9,
    'covtxtsizefol': 0.75,
    'covtxtsizeunf': 1.5,
    "covtxtsizefidnorm": 1.4,
    'covtxtangleunffidnorm': "45",
}
varList['MT_LLMETBuncertainties'] = {
    'xaxis'       : varList['MT_LLMETB']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}

varList['MT_LLMETBATLAS'] = {
    #'xaxis'       : 'm_{T}(\\ell_{1}, \\ell_{2},\\slash{E}_{T}, j) (GeV)',
    'xaxis'       : '#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j}) (GeV)',
    'yaxis'       : 'd#sigma/d(#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j})) (pb)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j})) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(#it{m}_{T}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{p}_{T}^{miss}, #it{j})) (1/GeV)',
    'genbinning'  : [0., 275., 375., 500., 1000.],
    'recobinning' : [0., 200., 235., 275., 300., 325., 400., 500., 1000.],
    'var'         : 'TMT_LLMETB',
    'var_response': 'MTLLMETBATLAS',
    'var_gen'     : 'TGenMT_LLMETB',
    'legpos'      : (0.51, 0.55, 0.71, 0.93),
    'legposdesc'  : (0.15, 0.425, 0.35, 0.81),
    'legpos_foldas':"BL",
    'legpos_fold' : "BL",
    'legpos_fid'  : "TL",
    "legpos_unf"  : (.18, .65, .38, .40),
    'legpos_unfas': "TL",
    'uncleg_fold' : "TL",
    'uncleg_unf'  : "TC",
    'uncleg_fid'  : "TL",
    'legpos_fidbin':"ML",
    'uncleg_fidbin':"TL",
    'resptxtsize' : 0.9,
    'covtxtsizefol': 0.75,
    'covtxtsizeunf': 1.5,
    "covtxtsizefidnorm": 1.4,
    'covtxtangleunffidnorm': "45",
}
varList['MT_LLMETBATLASuncertainties'] = {
    'xaxis'       : varList['MT_LLMETB']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)',
}

varList['M_LLB'] = {
    #'xaxis'       : 'm(#ell_{1}, #ell_{2}, j) (GeV)',
    'xaxis'       : '#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j}) (GeV)',
    'yaxis'       : 'd#sigma/d(#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (pb)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(#it{m}(#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (1/GeV)',
    #'genbinning'  : [0., 125., 250., 350., 400.],                         # usado en drafts, en nota, etc. (4 bins)
    #'recobinning' : [0., 125., 170., 180., 200., 240., 290., 350., 400.], # usado en drafts, en nota, etc. (4 bins)

    #'genbinning'  : [0., 130., 170., 200., 325., 400.], # propuesta (5 bins)
    #'recobinning' : [0., 130., 140., 150., 160., 170., 180., 190., 200., 325., 400.], # propuesta (5 bins)

    'genbinning'  : [0., 130., 170., 200., 250., 325., 400.], # propuesta (6 bins)
    'recobinning' : [0., 130., 140., 150., 160., 170., 180., 190., 200., 210., 220., 250., 325., 400.], # propuesta (6 bins) SELECCIONADA

    #'genbinning'  : [0., 130., 150., 170., 200., 250., 325., 400.], # propuesta (7 bins)
    #'recobinning' : [0., 130., 140., 150., 155., 160., 165., 170., 180., 190., 200., 210., 220., 250., 325., 400.], # propuesta (7 bins)

    #'genbinning'  : [0., 125., 170., 180., 200., 240., 290., 350., 400.], # propuesta (8 bins)
    #'recobinning' : [0., 125., 150., 160., 170., 180., 190., 200., 210., 220., 230., 240., 250., 275., 290., 350., 400.], # propuesta (8 bins)

    #'genbinning'  : [0., 130., 150., 160., 170., 180., 200., 250., 325., 400.], # propuesta (9 bins)
    #'recobinning' : [0., 130., 140., 145., 150., 155., 160., 165., 170., 175., 180., 185., 190., 195., 200., 225., 250., 325., 400.], # propuesta (9 bins)

    #'genbinning'  : [0., 130., 140., 150., 160., 170., 180., 200., 250., 325., 400.], # propuesta (10 bins)
    #'recobinning' : [0., 130., 140., 145., 150., 155., 160., 165., 170., 175., 180., 185., 190., 195., 200., 205., 210., 220., 250., 325., 400.], # propuesta (10 bins)

    'descbinning' : [0., 400.],
    'ndescbins'   : 16,
    'legpos'      : (0.75, 0.55, 0.95, 0.93),
    #'legposdesc'  : (0.11, 0.435, 0.31, 0.82),
    'legposdesc'  : (0.65, 0.55, 0.85, 0.93),
    'var'         : 'min(TM_LLB, 399.)',
    'var_response': 'MLLB',
    'var_gen'     : 'min(TGenM_LLB, 399.)',
    'uncleg_fold' : "TC",
    'uncleg_fid'  : "TL",
    'legpos_fidbin':"ML",
    "legpos_unf"   : (.18, .4, .38, .15),
    'covtxtsizefol': 0.5,
    'covtxtsizeunf': 1.3,
    'covtxtsizefidnorm': 1.25,
    'covtxtangleunffidnorm': "45",
}
varList['M_LLBuncertainties'] = {
    'xaxis'       : varList['M_LLB']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['BDT'] = {
    'xaxis'       : 'BDT disc. (adim.)',
    'yaxis'       : 'd#sigma/d(BDT disc.) (pb)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(BDT disc.) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(BDT disc.) (1/GeV)',
    'genbinning'  : [-1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.],
    'recobinning' : [-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.],
    'descbinning' : [-1., 1.],
    'ndescbins'   : 20,
    'legpos'      : (0.75, 0.55, 0.95, 0.93),
    #'legposdesc'  : (0.11, 0.435, 0.31, 0.82),
    'legposdesc'  : (0.65, 0.55, 0.85, 0.93),
    'var'         : 'TBDT',
    'var_response': 'BDT',
    'var_gen'     : 'TGenBDT',
    'uncleg_fold' : "TC",
    'uncleg_fid'  : "TL",
    'legpos_fidbin':"ML",
    'covtxtsizefol': 0.5,
    'covtxtsizeunf': 1.3,
}
varList['BDTuncertainties'] = {
    'xaxis'       : varList['BDT']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['BDT2j1b'] = {
    'xaxis'       : 'BDT disc. (adim.)',
    'yaxis'       : 'd#sigma/d(BDT disc.) (pb)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(BDT disc.) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(BDT disc.) (1/GeV)',
    'genbinning'  : [-1., -0.8, -0.6, -0.4, -0.2, 0., 0.2, 0.4, 0.6, 0.8, 1.],
    'recobinning' : [-1., -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.],
    'descbinning' : [-1., 1.],
    'ndescbins'   : 20,
    'legpos'      : (0.75, 0.55, 0.95, 0.93),
    #'legposdesc'  : (0.11, 0.435, 0.31, 0.82),
    'legposdesc'  : (0.65, 0.55, 0.85, 0.93),
    'var'         : 'TBDT2j1b',
    'var_response': 'BDT2j1b',
    'var_gen'     : 'TDressBDT2j1b',
    'uncleg_fold' : "TC",
    'uncleg_fid'  : "TL",
    'legpos_fidbin':"ML",
    'covtxtsizefol': 0.5,
    'covtxtsizeunf': 1.3,
}
varList['BDT2j1buncertainties'] = {
    'xaxis'       : varList['BDT']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}


varList['M_LeadingB'] = {
    'xaxis'       : 'm(\\ell_{1}, j) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0., 75., 175., 275., 400.],
    'recobinning' : [0., 75., 95., 115., 135., 175., 225., 275., 400.],
    'var'         : 'TM_LeadingB',
    'var_response': 'MLeadingB',
    'var_gen'     : 'TGenM_LeadingB',
}
varList['M_LeadingBuncertainties'] = {
    'xaxis'       : varList['M_LeadingB']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}

varList['M_SubLeadingB'] = {
    'xaxis'       : 'm(\\ell_{2}, j) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0., 60., 100., 150., 300.],
    'recobinning' : [0., 60., 70., 80., 90., 100., 125., 150., 300.],
    'var'         : 'TM_SubLeadingB',
    'var_response': 'MSubLeadingB',
    'var_gen'     : 'TGenM_SubLeadingB',
}
varList['M_SubLeadingBuncertainties'] = {
    'xaxis'       : varList['M_SubLeadingB']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}

varList['MET'] = {
    'xaxis'       : '\\slash{E}_{T} (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0, 50, 140, 200],
    'recobinning' : [0., 20., 35., 50., 70., 140., 200.],
    'var'         : 'TMET',
    'var_response': 'MET',
    'var_gen'     : 'TGenMET'
}
varList['METuncertainties'] = {
    'xaxis'       : varList['MET']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}

varList['MET_Phi'] = {
    'xaxis'       : '\\varphi(\\slash{E}_{T}) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [-r.TMath.Pi(), -1.5, 0, 1.5, r.TMath.Pi()],
    'recobinning' : [-r.TMath.Pi(), -2.25, -1.5, -.75, 0, .75, 1.5, 2.25, r.TMath.Pi()],
    'var'         : 'TMET_Phi',
    'var_response': 'METPhi',
    'var_gen'     : 'TGenMET_Phi',
}
varList['MET_Phiuncertainties'] = {
    'xaxis'       : varList['MET_Phi']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}


varList['Jet1_Pt'] = {
    'xaxis'       : 'Leading jet #it{p}_{T} (GeV)',
    'yaxis'       : 'd#sigma/d(leading jet #it{p}_{T}) (GeV)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(leading jet #it{p}_{T}) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(leading jet #it{p}_{T}) (1/GeV)',
    #'genbinning'  : [0, 75, 200, 300],
    #'recobinning' : [0., 50., 75., 110., 150., 200., 300.],
   #'genbinning'  : [0., 75., 140., 200., 300.],                            # binning presentado en singletop
   #'recobinning' : [0., 60., 80., 105., 120., 140., 170., 210., 300.],
#    'genbinning'  : [0., 60., 110., 150., 300.],
#    'recobinning' : [0., 60., 75., 90., 110., 125., 150., 175., 300.],
    #'genbinning'  : [0., 60., 90., 120., 150.],
    #'recobinning' : [0., 60., 70., 80, 90., 100., 110., 120., 150.],
    #'genbinning'  : [30., 60., 90., 120., 150.],                        # usado en drafts, en nota, etc. (4 bins)
    #'recobinning' : [30., 60., 70., 80, 90., 100., 110., 120., 150.],   # usado en drafts, en nota, etc. (4 bins)

    'genbinning'  : [30, 50, 70, 90, 125, 150], # propuesta (5 bins)
    'recobinning' : [30, 40, 45, 50, 55, 60, 65, 70, 90, 125, 150], # propuesta (5 bins) SELECCTIONADA

    #'genbinning'  : [30., 40., 50., 70., 90., 125., 150.], # propuesta (6 bins)
    #'recobinning' : [30., 35., 40., 45., 50., 55., 60., 65., 70., 90., 105., 125., 150.], # propuesta (6 bins)

    #'genbinning'  : [30., 40., 50., 60., 70., 90., 125., 150.], # propuesta (7 bins)
    #'recobinning' : [30., 35., 37.5,40., 42.5, 45., 50., 55.,  60., 65., 70., 80., 105., 125., 150.], # propuesta (7 bins)

    #'genbinning'  : [30., 50., 60., 70., 80., 90., 105., 125., 150.], # propuesta (8 bins)
    #'recobinning' : [30., 40., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 115., 125., 150.], # propuesta (8 bins)

    #'genbinning'  : [30., 40., 45., 50., 55., 60., 70., 80., 125., 150.], # propuesta (9 bins)
    #'recobinning' : [30., 37.5, 40.,42.5, 45.,47.5, 50.,52.5,  55., 57.5, 60., 65., 70., 75., 80., 90., 105., 125., 150.], # propuesta (9 bins)

    #'genbinning'  : [30., 40., 45., 50., 55., 60., 70., 80., 90., 125., 150.], # propuesta (10 bins)
    #'recobinning' : [30., 35.,37.5, 40.,42.5, 45.,47.5, 50.,52.5,  55., 57.5, 60., 65., 70., 75., 80., 85., 90., 105., 125., 150.], # propuesta (10 bins)

    #'descbinning' : [10., 150.],
    'descbinning' : [10., 200.],
    'ndescbins'   : 14,
    'legpos'      : (0.7, 0.55, 0.90, 0.93),
    'legposdesc'  : (0.7, 0.55, 0.90, 0.93),
    #"legpos_unf"  : "TC",
    #'var'         : 'min(TJet1_Pt, 149.)',
    'var'         : 'min(TJet1_Pt, 199.)',
    'var_response': 'Jet1_Pt',
    #'var_gen'     : 'min(TDressJet1_Pt, 149.)',
    'var_gen'     : 'min(TDressJet1_Pt, 199.)',
    'uncleg_fold' : "TL",
    'uncleg_unf'  : "TL",
    'uncleg_fid'  : "TL",
    'uncleg_fidbin':"TL",
    'legpos_fidbin':"TC",
    'resptxtsize' : 1.5,
    'covtxtsizeunf': 2,
    'covtxtsizefol': 1.25,
    "covtxtsizefidnorm": 1.8,
}
varList['Jet1_Ptuncertainties'] = {
    'xaxis'       : varList['Jet1_Pt']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseJet1_Pt'] = {
    'xaxis'       : 'Gen leading jet #it{p}_{T} (GeV)',
    'yaxis'       : 'Leading jet #it{p}_{T} (GeV)'
}


varList['Jet2_Pt'] = {
    'xaxis'       : 'Subleading jet #it{p}_{T} (GeV)',
    'yaxis'       : 'd#sigma/d(subleading jet #it{p}_{T}) (GeV)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(subleading jet #it{p}_{T}) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(subleading jet #it{p}_{T}) (1/GeV)',
    #'genbinning'  : [0, 75, 200, 300],
    #'recobinning' : [0., 50., 75., 110., 150., 200., 300.],
   #'genbinning'  : [0., 75., 140., 200., 300.],                            # binning presentado en singletop
   #'recobinning' : [0., 60., 80., 105., 120., 140., 170., 210., 300.],
#    'genbinning'  : [0., 60., 110., 150., 300.],
#    'recobinning' : [0., 60., 75., 90., 110., 125., 150., 175., 300.],
    #'genbinning'  : [0., 60., 90., 120., 150.],
    #'recobinning' : [0., 60., 70., 80, 90., 100., 110., 120., 150.],
    #'genbinning'  : [30., 60., 90., 120., 150.],                        # usado en drafts, en nota, etc. (4 bins)
    #'recobinning' : [30., 60., 70., 80, 90., 100., 110., 120., 150.],   # usado en drafts, en nota, etc. (4 bins)

    'genbinning'  : [30, 50, 70, 90, 125, 150], # propuesta (5 bins)
    'recobinning' : [30, 40, 45, 50, 55, 60, 65, 70, 90, 125, 150], # propuesta (5 bins) SELECCTIONADA

    #'genbinning'  : [30., 40., 50., 70., 90., 125., 150.], # propuesta (6 bins)
    #'recobinning' : [30., 35., 40., 45., 50., 55., 60., 65., 70., 90., 105., 125., 150.], # propuesta (6 bins)

    #'genbinning'  : [30., 40., 50., 60., 70., 90., 125., 150.], # propuesta (7 bins)
    #'recobinning' : [30., 35., 37.5,40., 42.5, 45., 50., 55.,  60., 65., 70., 80., 105., 125., 150.], # propuesta (7 bins)

    #'genbinning'  : [30., 50., 60., 70., 80., 90., 105., 125., 150.], # propuesta (8 bins)
    #'recobinning' : [30., 40., 50., 55., 60., 65., 70., 75., 80., 85., 90., 95., 100., 105., 115., 125., 150.], # propuesta (8 bins)

    #'genbinning'  : [30., 40., 45., 50., 55., 60., 70., 80., 125., 150.], # propuesta (9 bins)
    #'recobinning' : [30., 37.5, 40.,42.5, 45.,47.5, 50.,52.5,  55., 57.5, 60., 65., 70., 75., 80., 90., 105., 125., 150.], # propuesta (9 bins)

    #'genbinning'  : [30., 40., 45., 50., 55., 60., 70., 80., 90., 125., 150.], # propuesta (10 bins)
    #'recobinning' : [30., 35.,37.5, 40.,42.5, 45.,47.5, 50.,52.5,  55., 57.5, 60., 65., 70., 75., 80., 85., 90., 105., 125., 150.], # propuesta (10 bins)

    #'descbinning' : [10., 150.],
    'descbinning' : [10., 200.],
    'ndescbins'   : 14,
    'legpos'      : (0.7, 0.55, 0.90, 0.93),
    'legposdesc'  : (0.7, 0.55, 0.90, 0.93),
    #"legpos_unf"  : "TC",
    #'var'         : 'min(TJet2_Pt, 149.)',
    'var'         : 'min(TJet2_Pt, 199.)',
    'var_response': 'Jet2_Pt',
    #'var_gen'     : 'min(TDressJet2_Pt, 149.)',
    'var_gen'     : 'min(TDressJet2_Pt, 199.)',
    'uncleg_fold' : "TL",
    'uncleg_unf'  : "TL",
    'uncleg_fid'  : "TL",
    'uncleg_fidbin':"TL",
    'legpos_fidbin':"TC",
    'resptxtsize' : 1.5,
    'covtxtsizeunf': 2,
    'covtxtsizefol': 1.25,
    "covtxtsizefidnorm": 1.8,
}
varList['Jet2_Ptuncertainties'] = {
    'xaxis'       : varList['Jet2_Pt']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseJet2_Pt'] = {
    'xaxis'       : 'Gen subleading jet #it{p}_{T} (GeV)',
    'yaxis'       : 'Subleading jet #it{p}_{T} (GeV)'
}


varList['Jet1_Eta'] = {
    'xaxis'       : 'Leading jet #eta (adim.)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [-2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4],
    'recobinning' : [-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4],
    'var'         : 'TJet1_Eta',
    'var_response': 'Jet1_Eta',
    'var_gen'     : 'TDressJet1_Eta',
}
varList['Jet1_Etauncertainties'] = {
    'xaxis'       : varList['Jet1_Eta']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseJet1_Eta'] = {
    'xaxis'       : 'Dressed leading jet #eta (adim.)',
    'yaxis'       : 'Leading jet #eta (adim.)'
}


varList['Jet2_Eta'] = {
    'xaxis'       : 'Subleading jet #eta (adim.)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [-2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4],
    'recobinning' : [-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4],
    'var'         : 'TJet2_Eta',
    'var_response': 'Jet2_Eta',
    'var_gen'     : 'TDressJet2_Eta',
}
varList['Jet2_Etauncertainties'] = {
    'xaxis'       : varList['Jet2_Eta']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseJet2_Eta'] = {
    'xaxis'       : 'Dressed subleading jet #eta (adim.)',
    'yaxis'       : 'Subleading jet #eta (adim.)'
}


varList['LeadingJetPhi'] = {
    'xaxis'       : '\\varphi(j) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [-r.TMath.Pi(), -1.5, 0, 1.5, r.TMath.Pi()],
    'recobinning' : [-r.TMath.Pi(), -2.25, -1.5, -.75, 0, .75, 1.5, 2.25, r.TMath.Pi()],
    'var'         : 'TLeadingJetPhi',
    'var_response': 'LeadingJetPhi',
    'var_gen'     : 'TGenLeadingJetPhi',
}
varList['LeadingJetPhiuncertainties'] = {
    'xaxis'       : varList['LeadingJetPhi']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}


varList['LooseCentral1_Pt'] = {
    'xaxis'       : 'Leading loose jet #it{p}_{T} (GeV)',
    'yaxis'       : 'd#sigma/d(leading jet #it{p}_{T}) (GeV)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(leading loose jet #it{p}_{T}) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(leading loose jet #it{p}_{T}) (1/GeV)',
    'genbinning'  : [30, 50, 70, 90, 125, 150],
    'recobinning' : [30, 40, 45, 50, 55, 60, 65, 70, 90, 125, 150],

    'descbinning' : [10., 200.],
    'ndescbins'   : 14,
    'legpos'      : (0.7, 0.55, 0.90, 0.93),
    'legposdesc'  : (0.7, 0.55, 0.90, 0.93),
    #"legpos_unf"  : "TC",
    #'var'         : 'min(TLooseCentral1_Pt, 149.)',
    'var'         : 'min(TLooseCentral1_Pt, 199.)',
    'var_response': 'LooseCentral1_Pt',
    #'var_gen'     : 'min(TDressLooseCentral1_Pt, 149.)',
    'var_gen'     : 'min(TDressLooseCentral1_Pt, 199.)',
    'uncleg_fold' : "TL",
    'uncleg_unf'  : "TL",
    'uncleg_fid'  : "TL",
    'uncleg_fidbin':"TL",
    'legpos_fidbin':"TC",
    'resptxtsize' : 1.5,
    'covtxtsizeunf': 2,
    'covtxtsizefol': 1.25,
    "covtxtsizefidnorm": 1.8,
}
varList['LooseCentral1_Ptuncertainties'] = {
    'xaxis'       : varList['LooseCentral1_Pt']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseLooseCentral1_Pt'] = {
    'xaxis'       : 'Gen leading loose jet #it{p}_{T} (GeV)',
    'yaxis'       : 'Leading loose jet #it{p}_{T} (GeV)'
}



varList['LeadingLepE'] = {
    'xaxis'       : 'E(\\ell_{1}) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0, 70, 120, 250, 350],
    'recobinning' : [0., 40., 60., 80., 100., 120., 150., 250., 350.],
    'var'         : 'TLeadingLepE',
    'var_response': 'LeadingLepE',
    'var_gen'     : 'TGenLeadingLepE',
}
varList['LeadingLepEuncertainties'] = {
    'xaxis'       : varList['LeadingLepE']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['C_jll'] = {
    'xaxis'       : 'C(\\ell_{1},\\ell_{2}, j) (adim.)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0., 0.2, 0.4, 0.6, 0.8, 1.],
    'recobinning' : [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    'var'         : 'TC_jll',
    'var_response': 'C_jll',
    'var_gen'     : 'TDressC_jll',
}
varList['C_jlluncertainties'] = {
    'xaxis'       : varList['C_jll']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['Lep1Lep2_HTOverHTtot'] = {
    'xaxis'       : '(p_{T}(\\ell_{1}) + p_{T}(\\ell_{2}))/H_{T} (adim.)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0., 0.2, 0.4, 0.6, 0.8, 1.],
    'recobinning' : [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    'var'         : 'TLep1Lep2_HTOverHTtot',
    'var_response': 'Lep1Lep2_HTOverHTtot',
    'var_gen'     : 'TDressLep1Lep2_HTOverHTtot',
}
varList['Lep1Lep2_HTOverHTtotuncertainties'] = {
    'xaxis'       : varList['Lep1Lep2_HTOverHTtot']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['Lep1_Pt'] = {
    #'xaxis'       : 'p_{T}(\\ell_{1}) (GeV)',
    'xaxis'       : 'Leading lepton #it{p}_{T} (GeV)',
    'yaxis'       : 'd#sigma/d(leading lepton #it{p}_{T}) (GeV)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(leading lepton #it{p}_{T}) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(leading lepton #it{p}_{T}) (1/GeV)',
    #'genbinning'  : [0, 50, 120, 190, 250],
    #'recobinning' : [0, 50, 65, 85, 97, 110, 145, 180, 250],
   #'genbinning'  : [0, 50, 120, 160, 250],                            # binning presentado en singletop
   #'recobinning' : [0, 50, 65, 85, 97, 110, 145, 180, 250],
#    'genbinning'  : [0., 50., 90., 135., 250.],
#    'recobinning' : [0., 50., 65., 80., 95., 110., 130., 150., 250.],
    #'genbinning'  : [0., 50., 90., 135., 250.],
    #'recobinning' : [0., 50., 60., 75., 90., 105., 115., 135., 250.],
    #'genbinning'  : [0., 50., 90., 125., 150.],
    #'recobinning' : [0., 50., 60., 70., 80., 90., 105., 125., 150.], # antes de 15-10-2018
    #'genbinning'  : [25., 50., 90., 125., 150.],                      # usado en drafts, en nota, etc. (4 bins)
    #'recobinning' : [25., 50., 60., 70., 80., 90., 105., 125., 150.], # usado en drafts, en nota, etc. (4 bins)

    #'genbinning'  : [25., 45., 65., 90., 125., 150.], # propuesta (5 bins)
    #'recobinning' : [25., 35., 45., 50.,  55., 60., 65., 75., 90., 125., 150.], # propuesta (5 bins)

    #'genbinning'  : [25., 45., 65., 85., 100., 125., 150.], # propuesta (6 bins)
    #'recobinning' : [25., 35., 45., 50., 55., 60., 65., 75., 85., 100., 115., 125., 150.], # propuesta (6 bins)

    #'genbinning'  : [25., 45., 55., 65., 85., 100., 125., 150.], # propuesta (7 bins)
    #'recobinning' : [25., 35., 45., 50., 55.,  60.,  65.,  70., 75., 80., 85., 100., 115., 125., 150.], # propuesta (7 bins)

    'genbinning'  : [25, 50, 60, 70, 80, 90, 105, 125, 150], # propuesta (8 bins)
    'recobinning' : [25, 35, 50, 55, 60, 65,  70,  75,  80, 85, 90, 95, 100, 105, 115, 125, 150], # propuesta (8 bins) SELECCIONADA

    #'genbinning'  : [25., 45., 50., 55., 60., 65., 85., 100., 125., 150.], # propuesta (9 bins)
    #'recobinning' : [25., 40., 42.5, 45.,47.5,50.,52.5,  55., 57.5,  60., 62.5, 65., 70., 75., 85., 100., 115., 125., 150.], # propuesta (9 bins)

    #'genbinning'  : [25., 35., 45., 50., 55., 60., 65., 85., 100., 125., 150.], # propuesta (10 bins)
    #'recobinning' : [25., 35., 40., 42.5, 45., 47.5, 50., 52.5, 55., 57.5, 60., 62.5, 65., 67.5, 70., 75., 85., 100., 115., 125., 150.], # propuesta (10 bins)

    #'descbinning' : [10., 150.],
    'descbinning' : [0., 200.],
    #'ndescbins'   : 14,
    'legpos'      : (0.7, 0.55, 0.90, 0.93),
    'legposdesc'  : (0.58, 0.55, 0.78, 0.93),
    #'var'         : 'min(TLep1_Pt, 149.)',
    'var'         : 'min(TLep1_Pt, 199.)',
    'var_response': 'Lep1_Pt',
    #'var_gen'     : 'min(TDressLep1_Pt, 149.)',
    'var_gen'     : 'min(TDressLep1_Pt, 199.)',
    'uncleg_fold' : "TL",
    'uncleg_fid'  : "TL",
    "legpos_unf"   : (.18, .3, .32, .05),
    'legpos_fidbin':"TC",
    'uncleg_fidbin':"TL",
    'resptxtsize' : 1.5,
    'covtxtsizefol': 1.2,
    'covtxtsizeunf': 2,
    "covtxtsizefidnorm": 1.45,
}
varList['Lep1_Ptuncertainties'] = {
    'xaxis'       : varList['Lep1_Pt']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseLeadingLepPt'] = {
    'xaxis'       : 'Gen Leading lep #it{p}_{T} (GeV)',
    'yaxis'       : 'Leading lep #it{p}_{T} (GeV)'
}


varList['Lep2_Pt'] = {
    #'xaxis'       : 'p_{T}(\\ell_{2}) (GeV)',
    'xaxis'       : 'Subleading lepton #it{p}_{T} (GeV)',
    'yaxis'       : 'd#sigma/d(subleading lepton #it{p}_{T}) (GeV)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(subleading lepton #it{p}_{T}) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(subleading lepton #it{p}_{T}) (1/GeV)',
    #'genbinning'  : [0, 50, 120, 190, 250],
    #'recobinning' : [0, 50, 65, 85, 97, 110, 145, 180, 250],
   #'genbinning'  : [0, 50, 120, 160, 250],                            # binning presentado en singletop
   #'recobinning' : [0, 50, 65, 85, 97, 110, 145, 180, 250],
#    'genbinning'  : [0., 50., 90., 135., 250.],
#    'recobinning' : [0., 50., 65., 80., 95., 110., 130., 150., 250.],
    #'genbinning'  : [0., 50., 90., 135., 250.],
    #'recobinning' : [0., 50., 60., 75., 90., 105., 115., 135., 250.],
    #'genbinning'  : [0., 50., 90., 125., 150.],
    #'recobinning' : [0., 50., 60., 70., 80., 90., 105., 125., 150.], # antes de 15-10-2018
    #'genbinning'  : [25., 50., 90., 125., 150.],                      # usado en drafts, en nota, etc. (4 bins)
    #'recobinning' : [25., 50., 60., 70., 80., 90., 105., 125., 150.], # usado en drafts, en nota, etc. (4 bins)

    #'genbinning'  : [25., 45., 65., 90., 125., 150.], # propuesta (5 bins)
    #'recobinning' : [25., 35., 45., 50.,  55., 60., 65., 75., 90., 125., 150.], # propuesta (5 bins)

    #'genbinning'  : [25., 45., 65., 85., 100., 125., 150.], # propuesta (6 bins)
    #'recobinning' : [25., 35., 45., 50., 55., 60., 65., 75., 85., 100., 115., 125., 150.], # propuesta (6 bins)

    #'genbinning'  : [25., 45., 55., 65., 85., 100., 125., 150.], # propuesta (7 bins)
    #'recobinning' : [25., 35., 45., 50., 55.,  60.,  65.,  70., 75., 80., 85., 100., 115., 125., 150.], # propuesta (7 bins)

    'genbinning'  : [25, 50, 60, 70, 80, 90, 105, 125, 150], # propuesta (8 bins)
    'recobinning' : [25, 35, 50, 55, 60, 65,  70,  75,  80, 85, 90, 95, 100, 105, 115, 125, 150], # propuesta (8 bins) SELECCIONADA

    #'genbinning'  : [25., 45., 50., 55., 60., 65., 85., 100., 125., 150.], # propuesta (9 bins)
    #'recobinning' : [25., 40., 42.5, 45.,47.5,50.,52.5,  55., 57.5,  60., 62.5, 65., 70., 75., 85., 100., 115., 125., 150.], # propuesta (9 bins)

    #'genbinning'  : [25., 35., 45., 50., 55., 60., 65., 85., 100., 125., 150.], # propuesta (10 bins)
    #'recobinning' : [25., 35., 40., 42.5, 45., 47.5, 50., 52.5, 55., 57.5, 60., 62.5, 65., 67.5, 70., 75., 85., 100., 115., 125., 150.], # propuesta (10 bins)

    #'descbinning' : [10., 150.],
    'descbinning' : [10., 200.],
    'ndescbins'   : 14,
    'legpos'      : (0.7, 0.55, 0.90, 0.93),
    'legposdesc'  : (0.58, 0.55, 0.78, 0.93),
    #'var'         : 'min(TLep2_Pt, 149.)',
    'var'         : 'min(TLep2_Pt, 199.)',
    'var_response': 'Lep2_Pt',
    #'var_gen'     : 'min(TDressLep2_Pt, 149.)',
    'var_gen'     : 'min(TDressLep2_Pt, 199.)',
    'uncleg_fold' : "TL",
    'uncleg_fid'  : "TL",
    "legpos_unf"   : (.18, .3, .32, .05),
    'legpos_fidbin':"TC",
    'uncleg_fidbin':"TL",
    'resptxtsize' : 1.5,
    'covtxtsizefol': 1.2,
    'covtxtsizeunf': 2,
    "covtxtsizefidnorm": 1.45,
}
varList['Lep1_Ptuncertainties'] = {
    'xaxis'       : varList['Lep2_Pt']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseLeadingLepPt'] = {
    'xaxis'       : 'Gen Subleading lep #it{p}_{T} (GeV)',
    'yaxis'       : 'Subleading lep #it{p}_{T} (GeV)'
}


varList['LepMuon_Pt'] = {
    #'xaxis'       : 'p_{T}(\\ell_{1}) (GeV)',
    'xaxis'       : 'Muon #it{p}_{T} (GeV)',
    'yaxis'       : 'd#sigma/d(muon #it{p}_{T}) (GeV)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(muon #it{p}_{T}) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(muon #it{p}_{T}) (1/GeV)',
    #'genbinning'  : [0, 50, 120, 190, 250],
    #'recobinning' : [0, 50, 65, 85, 97, 110, 145, 180, 250],
   #'genbinning'  : [0, 50, 120, 160, 250],                            # binning presentado en singletop
   #'recobinning' : [0, 50, 65, 85, 97, 110, 145, 180, 250],
#    'genbinning'  : [0., 50., 90., 135., 250.],
#    'recobinning' : [0., 50., 65., 80., 95., 110., 130., 150., 250.],
    #'genbinning'  : [0., 50., 90., 135., 250.],
    #'recobinning' : [0., 50., 60., 75., 90., 105., 115., 135., 250.],
    #'genbinning'  : [0., 50., 90., 125., 150.],
    #'recobinning' : [0., 50., 60., 70., 80., 90., 105., 125., 150.], # antes de 15-10-2018
    #'genbinning'  : [25., 50., 90., 125., 150.],                      # usado en drafts, en nota, etc. (4 bins)
    #'recobinning' : [25., 50., 60., 70., 80., 90., 105., 125., 150.], # usado en drafts, en nota, etc. (4 bins)

    #'genbinning'  : [25., 45., 65., 90., 125., 150.], # propuesta (5 bins)
    #'recobinning' : [25., 35., 45., 50.,  55., 60., 65., 75., 90., 125., 150.], # propuesta (5 bins)

    #'genbinning'  : [25., 45., 65., 85., 100., 125., 150.], # propuesta (6 bins)
    #'recobinning' : [25., 35., 45., 50., 55., 60., 65., 75., 85., 100., 115., 125., 150.], # propuesta (6 bins)

    #'genbinning'  : [25., 45., 55., 65., 85., 100., 125., 150.], # propuesta (7 bins)
    #'recobinning' : [25., 35., 45., 50., 55.,  60.,  65.,  70., 75., 80., 85., 100., 115., 125., 150.], # propuesta (7 bins)

    'genbinning'  : [25, 50, 60, 70, 80, 90, 105, 125, 150], # propuesta (8 bins)
    'recobinning' : [25, 35, 50, 55, 60, 65,  70,  75,  80, 85, 90, 95, 100, 105, 115, 125, 150], # propuesta (8 bins) SELECCIONADA

    #'genbinning'  : [25., 45., 50., 55., 60., 65., 85., 100., 125., 150.], # propuesta (9 bins)
    #'recobinning' : [25., 40., 42.5, 45.,47.5,50.,52.5,  55., 57.5,  60., 62.5, 65., 70., 75., 85., 100., 115., 125., 150.], # propuesta (9 bins)

    #'genbinning'  : [25., 35., 45., 50., 55., 60., 65., 85., 100., 125., 150.], # propuesta (10 bins)
    #'recobinning' : [25., 35., 40., 42.5, 45., 47.5, 50., 52.5, 55., 57.5, 60., 62.5, 65., 67.5, 70., 75., 85., 100., 115., 125., 150.], # propuesta (10 bins)

    #'descbinning' : [10., 150.],
    'descbinning' : [0., 200.],
    #'ndescbins'   : 14,
    'legpos'      : (0.7, 0.55, 0.90, 0.93),
    'legposdesc'  : (0.58, 0.55, 0.78, 0.93),
    #'var'         : 'min(TLepMuon_Pt, 149.)',
    'var'         : 'min(TLepMuon_Pt, 199.)',
    'var_response': 'LepMuon_Pt',
    #'var_gen'     : 'min(TDressLepMuon_Pt, 149.)',
    'var_gen'     : 'min(TDressLepMuon_Pt, 199.)',
    'uncleg_fold' : "TL",
    'uncleg_fid'  : "TL",
    "legpos_unf"   : (.18, .3, .32, .05),
    'legpos_fidbin':"TC",
    'uncleg_fidbin':"TL",
    'resptxtsize' : 1.5,
    'covtxtsizefol': 1.2,
    'covtxtsizeunf': 2,
    "covtxtsizefidnorm": 1.45,
}
varList['LepMuon_Ptuncertainties'] = {
    'xaxis'       : varList['LepMuon_Pt']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseLeadingLepPt'] = {
    'xaxis'       : 'Dressed muon #it{p}_{T} (GeV)',
    'yaxis'       : 'Muon #it{p}_{T} (GeV)'
}


varList['LepElec_Pt'] = {
    #'xaxis'       : 'p_{T}(\\ell_{2}) (GeV)',
    'xaxis'       : 'Electron #it{p}_{T} (GeV)',
    'yaxis'       : 'd#sigma/d(electron #it{p}_{T}) (GeV)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(electron #it{p}_{T}) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(electron #it{p}_{T}) (1/GeV)',
    #'genbinning'  : [0, 50, 120, 190, 250],
    #'recobinning' : [0, 50, 65, 85, 97, 110, 145, 180, 250],
   #'genbinning'  : [0, 50, 120, 160, 250],                            # binning presentado en singletop
   #'recobinning' : [0, 50, 65, 85, 97, 110, 145, 180, 250],
#    'genbinning'  : [0., 50., 90., 135., 250.],
#    'recobinning' : [0., 50., 65., 80., 95., 110., 130., 150., 250.],
    #'genbinning'  : [0., 50., 90., 135., 250.],
    #'recobinning' : [0., 50., 60., 75., 90., 105., 115., 135., 250.],
    #'genbinning'  : [0., 50., 90., 125., 150.],
    #'recobinning' : [0., 50., 60., 70., 80., 90., 105., 125., 150.], # antes de 15-10-2018
    #'genbinning'  : [25., 50., 90., 125., 150.],                      # usado en drafts, en nota, etc. (4 bins)
    #'recobinning' : [25., 50., 60., 70., 80., 90., 105., 125., 150.], # usado en drafts, en nota, etc. (4 bins)

    #'genbinning'  : [25., 45., 65., 90., 125., 150.], # propuesta (5 bins)
    #'recobinning' : [25., 35., 45., 50.,  55., 60., 65., 75., 90., 125., 150.], # propuesta (5 bins)

    #'genbinning'  : [25., 45., 65., 85., 100., 125., 150.], # propuesta (6 bins)
    #'recobinning' : [25., 35., 45., 50., 55., 60., 65., 75., 85., 100., 115., 125., 150.], # propuesta (6 bins)

    #'genbinning'  : [25., 45., 55., 65., 85., 100., 125., 150.], # propuesta (7 bins)
    #'recobinning' : [25., 35., 45., 50., 55.,  60.,  65.,  70., 75., 80., 85., 100., 115., 125., 150.], # propuesta (7 bins)

    'genbinning'  : [25, 50, 60, 70, 80, 90, 105, 125, 150], # propuesta (8 bins)
    'recobinning' : [25, 35, 50, 55, 60, 65,  70,  75,  80, 85, 90, 95, 100, 105, 115, 125, 150], # propuesta (8 bins) SELECCIONADA

    #'genbinning'  : [25., 45., 50., 55., 60., 65., 85., 100., 125., 150.], # propuesta (9 bins)
    #'recobinning' : [25., 40., 42.5, 45.,47.5,50.,52.5,  55., 57.5,  60., 62.5, 65., 70., 75., 85., 100., 115., 125., 150.], # propuesta (9 bins)

    #'genbinning'  : [25., 35., 45., 50., 55., 60., 65., 85., 100., 125., 150.], # propuesta (10 bins)
    #'recobinning' : [25., 35., 40., 42.5, 45., 47.5, 50., 52.5, 55., 57.5, 60., 62.5, 65., 67.5, 70., 75., 85., 100., 115., 125., 150.], # propuesta (10 bins)

    'descbinning' : [10., 200.],
    'ndescbins'   : 14,
    'legpos'      : (0.7, 0.55, 0.90, 0.93),
    'legposdesc'  : (0.58, 0.55, 0.78, 0.93),
    'var'         : 'min(TLepElec_Pt, 199.)',
    'var_response': 'LepElec_Pt',
    'var_gen'     : 'min(TDressLepElec_Pt, 199.)',
    'uncleg_fold' : "TL",
    'uncleg_fid'  : "TL",
    "legpos_unf"   : (.18, .3, .32, .05),
    'legpos_fidbin':"TC",
    'uncleg_fidbin':"TL",
    'resptxtsize' : 1.5,
    'covtxtsizefol': 1.2,
    'covtxtsizeunf': 2,
    "covtxtsizefidnorm": 1.45,
}
varList['LepElec_Ptuncertainties'] = {
    'xaxis'       : varList['LepElec_Pt']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseLeadingLepPt'] = {
    'xaxis'       : 'Dressed electron #it{p}_{T} (GeV)',
    'yaxis'       : 'Electron #it{p}_{T} (GeV)'
}


varList['Fiducial'] = {
    'xaxis'       : 'a.u.',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [25., 150.],
    'recobinning' : [25., 150.],
    #'recobinning' : [25., 80., 150.],
    'var'         : 'min(TLep1_Pt, 149.)',
    'var_response': 'Fiducial',
    'var_gen'     : 'min(TDressLep1_Pt, 149.)',
}
varList['Fiducialuncertainties'] = {
    'xaxis'       : 'a.u.',
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseFiducial'] = {
    'xaxis'       : 'a.u.',
    'yaxis'       : 'a.u.'
}

varList['FiducialtWttbar'] = {
    'xaxis'       : 'a.u.',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [25., 150.],
    'recobinning' : [25., 150.],
    'var'         : 'min(TLep1_Pt, 149.)',
    'var_response': 'FiducialtWttbar',
    'var_gen'     : 'min(TDressLep1_Pt, 149.)',
}
varList['FiducialtWttbaruncertainties'] = {
    'xaxis'       : 'a.u.',
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseFiducialtWttbar'] = {
    'xaxis'       : 'a.u.',
    'yaxis'       : 'a.u.'
}

varList['LeadingLepPhi'] = {
    'xaxis'       : '\\varphi(\\ell_{1}) (rad)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [-r.TMath.Pi(), -1.5, 0, 1.5, r.TMath.Pi()],
    'recobinning' : [-r.TMath.Pi(), -2.25, -1.5, -.75, 0, .75, 1.5, 2.25, r.TMath.Pi()],
    'var'         : 'TLeadingLepPhi',
    'var_response': 'LeadingLepPhi',
    'var_gen'     : 'TGenLeadingLepPhi',
}
varList['LeadingLepPhiuncertainties'] = {
    'xaxis'       : varList['LeadingLepPhi']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}


varList['Lep1_Eta'] = {
    'xaxis'       : 'Leading lepton #eta (adim.)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [-2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4],
    'recobinning' : [-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4],
    #'genbinning'  : [0., 0.5, 1., 1.6, 2.4],
    #'recobinning' : [0., 0.25, 0.5, 0.75, 1., 1.3, 1.6, 2., 2.4],
    'var'         : 'TLep1_Eta',
    'var_response': 'Lep1_Eta',
    'var_gen'     : 'TDressLep1_Eta',
}
varList['Lep1_Etauncertainties'] = {
    'xaxis'       : varList['Lep1_Eta']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseLep1_Eta'] = {
    'xaxis'       : 'Dressed leading lepton |#eta| (adim.)',
    'yaxis'       : 'Leading lepton |#eta| (adim.)'
}


varList['LepMuon_Eta'] = {
    'xaxis'       : 'Muon #eta (adim.)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [-2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4],
    'recobinning' : [-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4],
    #'genbinning'  : [0., 0.5, 1., 1.6, 2.4],
    #'recobinning' : [0., 0.25, 0.5, 0.75, 1., 1.3, 1.6, 2., 2.4],
    'var'         : 'TLepMuon_Eta',
    'var_response': 'LepMuon_Eta',
    'var_gen'     : 'TDressLepMuon_Eta',
}
varList['LepMuon_Etauncertainties'] = {
    'xaxis'       : varList['LepMuon_Eta']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseLepMuon_Eta'] = {
    'xaxis'       : 'Dressed muon |#eta| (adim.)',
    'yaxis'       : 'Muon |#eta| (adim.)'
}


varList['LepElec_Eta'] = {
    'xaxis'       : 'Elec #eta (adim.)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [-2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4],
    'recobinning' : [-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4],
    #'genbinning'  : [0., 0.5, 1., 1.6, 2.4],
    #'recobinning' : [0., 0.25, 0.5, 0.75, 1., 1.3, 1.6, 2., 2.4],
    'var'         : 'TLepElec_Eta',
    'var_response': 'LepElec_Eta',
    'var_gen'     : 'TDressLepElec_Eta',
}
varList['LepElec_Etauncertainties'] = {
    'xaxis'       : varList['LepElec_Eta']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}
varList['ResponseLepElec_Eta'] = {
    'xaxis'       : 'Dressed muon |#eta| (adim.)',
    'yaxis'       : 'Elec |#eta| (adim.)'
}


varList['SubLeadingLepE'] = {
    'xaxis'       : 'E(\\ell_{2}) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0, 50, 100, 175, 250],
    'recobinning' : [0., 30., 50., 70., 90., 115., 140., 175., 250.],
    'var'         : 'TSubLeadingLepE',
    'var_response': 'SubLeadingLepE',
    'var_gen'     : 'TGenSubLeadingLepE',
}
varList['SubLeadingLepEuncertainties'] = {
    'xaxis'       : varList['SubLeadingLepE']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}


varList['SubLeadingLepPt'] = {
    'xaxis'       : 'p_{T}(\\ell_{2}) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0, 30, 60, 100, 150],
    'recobinning' : [0., 30., 40., 50., 58., 68., 78., 100., 150.],
    'var'         : 'TSubLeadingLepPt',
    'var_response': 'SubLeadingLepPt',
    'var_gen'     : 'TGenSubLeadingLepPt',
}
varList['SubLeadingLepPtuncertainties'] = {
    'xaxis'       : varList['SubLeadingLepPt']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}

varList['SubLeadingLepPhi'] = {
    'xaxis'       : '\\varphi(\\ell_{2}) (rad)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [-r.TMath.Pi(), -1.5, 0, 1.5, r.TMath.Pi()],
    'recobinning' : [-r.TMath.Pi(), -2.25, -1.5, -.75, 0, .75, 1.5, 2.25, r.TMath.Pi()],
    'var'         : 'TSubLeadingLepPhi',
    'var_response': 'SubLeadingLepPhi',
    'var_gen'     : 'TGenSubLeadingLepPhi',
}
varList['SubLeadingLepPhiuncertainties'] = {
    'xaxis'       : varList['SubLeadingLepPhi']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}


varList['Lep2_Eta'] = {
    'xaxis'       : 'Subleading lepton #eta',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [-2.4, -1.8, -1.2, -0.6, 0., 0.6, 1.2, 1.8, 2.4],
    'recobinning' : [-2.4, -2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0., 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1, 2.4],
    'var'         : 'TLep2_Eta',
    'var_response': 'Lep2_Eta',
    'var_gen'     : 'TDressLep2_Eta',
}
varList['Lep2_Etauncertainties'] = {
    'xaxis'       : varList['Lep2_Eta']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}


varList['Lep1Lep2_Pt'] = {
    'xaxis'       : 'p_{T}(\\ell_{1}, \\ell_{2}) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0, 50, 100, 150, 200],
    'recobinning' : [0., 40., 60., 80., 90., 110., 130., 150., 200.],
    'var'         : 'TLep1Lep2_Pt',
    'var_response': 'Lep1Lep2_Pt',
    'var_gen'     : 'TDressLep1Lep2_Pt',
}
varList['Lep1Lep2_Ptuncertainties'] = {
    'xaxis'       : varList['Lep1Lep2_Pt']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['Lep1Lep2Jet1_Pt'] = {
    'xaxis'       : 'p_{T}(\\ell_{1}, \\ell_{2}, j) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0, 50, 100, 150, 200],
    'recobinning' : [0., 40., 60., 80., 90., 110., 130., 150., 200.],
    'var'         : 'TLep1Lep2Jet1_Pt',
    'var_response': 'Lep1Lep2Jet1_Pt',
    'var_gen'     : 'TDressLep1Lep2Jet1_Pt',
}
varList['Lep1Lep2_Ptuncertainties'] = {
    'xaxis'       : varList['Lep1Lep2Jet1_Pt']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['DilepJetPt'] = {
    'xaxis'       : 'p_{T}(\\ell_{1}, \\ell_{2}, j) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0., 40., 80., 120., 200.],
    'recobinning' : [0., 20., 30., 40., 50., 60., 80., 120., 200.],
    'var'         : 'TDilepJetPt',
    'var_response': 'DilepJetPt',
    'var_gen'     : 'TGenDilepJetPt',
}
varList['DilepJetPtuncertainties'] = {
    'xaxis'       : varList['DilepJetPt']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}

varList['DilepMETJetPt'] = {
    'xaxis'       : 'p_{T}(\\ell_{1} ,\\ell_{2}, j,\\slash{E}_{T}) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0., 20., 40., 70., 150.],
    'recobinning' : [0., 10., 20., 30., 40., 50., 60., 70., 150.],
    'var'         : 'TDilepMETJetPt',
    'var_response': 'DilepMETJetPt',
    'var_gen'     : 'TGenDilepMETJetPt',
}
varList['DilepMETJetPtuncertainties'] = {
    'xaxis'       : varList['DilepMETJetPt']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}

varList['HTtot'] = {
    'xaxis'       : 'H_{T}(\\ell_{1}, \\ell_{2}, j,\\slash{E}_{T}) (GeV)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0, 200, 300, 450, 600],
    'recobinning' : [0., 150., 200., 250., 300., 350., 400., 450., 600.],
    'var'         : 'THTtot',
    'var_response': 'HTtot',
    'var_gen'     : 'TDressHTtot',
}
varList['HTtotuncertainties'] = {
    'xaxis'       : varList['HTtot']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}


varList['DilepMETJet1Pz'] = {
    'xaxis'       : '#it{p}_{Z} (#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j}) (GeV)',
    'yaxis'       : 'd#sigma/d(#it{p}_{Z} (#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (pb)',
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#it{p}_{Z} (#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(#it{p}_{Z} (#it{e}^{#pm}, #it{#mu}^{#mp}, #it{j})) (1/GeV)',
    #'genbinning'  : [0., 100., 200., 350., 450.],                        # usado en drafts, en nota, etc. (4 bins)
    #'recobinning' : [0., 50., 100., 150., 200., 250., 300., 350., 450.], # usado en drafts, en nota, etc. (4 bins)

    #'genbinning'  : [0., 40., 100., 175., 300., 450.],  # propuesta (5 bins)
    #'recobinning' : [0., 20., 40., 60., 80., 100., 125., 150., 175., 350., 450.], # propuesta (5 bins)

    #'genbinning'  : [0., 40., 100., 175., 250., 350., 450.],  # propuesta (6 bins)
    #'recobinning' : [0., 20., 40., 60., 80., 100., 125., 150., 175., 200., 250., 350., 450.], # propuesta (6 bins)

    #'genbinning'  : [0., 40., 60., 100., 175., 250., 350., 450.],  # propuesta (7 bins)
    #'recobinning' : [0., 20., 30.,  40.,  50.,  60.,  80., 100., 125., 150., 175., 200., 250., 350., 450.], # propuesta (7 bins)

    'genbinning'  : [0., 50., 100., 150., 200., 250., 300., 350., 450.],                        # propuesta (8 bins)
    'recobinning' : [0., 20., 40., 60., 80., 100., 125., 150., 175., 200., 225., 250., 275., 300., 350., 400., 450.], # propuesta (8 bins) SELECCIONADA

    #'genbinning'  : [0., 20., 40., 60., 80., 100., 175., 250., 350., 450.],  # propuesta (9 bins)
    #'recobinning' : [0., 5.,  10., 15., 20.,  25.,  30.,  40.,  50.,  60.,  70.,  80.,  90., 100., 150., 175., 250., 350., 450.], # propuesta (9 bins)

    #'genbinning'  : [0., 20., 40., 60., 80., 100., 125., 175., 250., 350., 450.],  # propuesta (10 bins)
    #'recobinning' : [0., 5.,  10., 15., 20.,  25.,  30.,  40.,  50.,  60.,  70.,  80.,  90., 100., 125., 150., 175., 200., 250., 350., 450.], # propuesta (10 bins)

    'descbinning' : [0., 450.],
    'ndescbins'   : 18,
    'var'         : 'min(abs(TDilepMETJet1Pz), 449.)',
    'var_response': 'DilepMETJet1Pz',
    'var_gen'     : 'min(abs(TGenDilepMETJet1Pz), 449.)',
    'legpos'      : (0.51, 0.55, 0.71, 0.93),
    'legposdesc'  : (0.57, 0.55, 0.78, 0.93),
    #'legpos_fold':  "TC",
    'legpos_fold' : (.5, .685, .725, .93),
    'legpos_foldas':"BL",
    'uncleg_fold' : "TL",
    'uncleg_fid'  : "TL",
    'legpos_fidbin':"BL",
    "legpos_unf"  : "TC",
    'resptxtsize' : 1.7,
    'covtxtsizeunf': 2,
    'covtxtsizefol': 1.4,
    "covtxtsizefidnorm": 1.6,
}
varList['DilepMETJet1Pzuncertainties'] = {
    'xaxis'       : varList['DilepMETJet1Pz']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}

varList['LLMETBEta'] = {
    'xaxis'       : '\\eta(\\ell_{1}, \\ell_{2}, j)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0., 1.25, 2.5, 3.75, 5.],
    'recobinning' : [0., 0.75, 1.25, 1.75, 2.25, 2.75, 3.25, 3.75, 5.],
    'var'         : 'abs(TLLMETBEta)',
    'var_response': 'LLMETBEta',
    'var_gen'     : 'abs(TGenLLMETBEta)',
}
varList['LLMETBEtauncertainties'] = {
    'xaxis'       : varList['LLMETBEta']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['Sys_M'] = {
    'xaxis'       : 'm(\\ell_{1}, \\ell_{2}, j,\\slash{E}_{T})',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0., 225., 325., 425., 700.],
    'recobinning' : [0., 225., 250., 275., 300., 325., 350., 425., 700.],
    'var'         : 'TSys_M',
    'var_response': 'Sys_M',
    'var_gen'     : 'TDressSys_M',
}
varList['Sys_Muncertainties'] = {
    'xaxis'       : varList['Sys_M']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['Sys_Pt'] = {
    'xaxis'       : 'p_{T}(\\ell_{1}, \\ell_{2}, j,\\slash{E}_{T})',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0., 20., 40., 70., 150.],
    'recobinning' : [0., 10., 20., 30., 40., 50., 60., 70., 150.],
    'var'         : 'TSys_Pt',
    'var_response': 'Sys_Pt',
    'var_gen'     : 'TDressSys_Pt',
}
varList['Sys_Ptuncertainties'] = {
    'xaxis'       : varList['Sys_Pt']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['Sys_PtOverHTtot'] = {
    'xaxis'       : 'p_{T}(\\ell_{1}, \\ell_{2}, j,\\slash{E}_{T})/H_{T}',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0., 0.2, 0.4, 0.6, 0.8, 1.],
    'recobinning' : [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    'var'         : 'TSys_PtOverHTtot',
    'var_response': 'Sys_PtOverHTtot',
    'var_gen'     : 'TDressSys_PtOverHTtot',
}
varList['Sys_PtOverHTtotuncertainties'] = {
    'xaxis'       : varList['Sys_PtOverHTtot']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['Mll'] = {
    'xaxis'       : 'm(\\ell_{1}, \\ell_{2})',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0., 50., 100., 150., 300.],
    'recobinning' : [0., 25., 45., 60., 75., 100., 125., 150., 300.],
    'var'         : 'TMll',
    'var_response': 'Mll',
    'var_gen'     : 'TGenMll',
}
varList['Mlluncertainties'] = {
    'xaxis'       : varList['Mll']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}


varList['Lep1Lep2_DPhi'] = {
    #'xaxis'       : '\\Delta \\varphi(\\ell_{1}, \\ell_{2}) (rad)',
    'xaxis'       : "#Delta#it{#varphi}(#it{e}^{#pm}, #it{#mu}^{#mp})/#it{#pi} (adim.)",
    'yaxis'       : "d#sigma/d(#Delta#it{#varphi}(#it{e}^{#pm}, #it{#mu}^{#mp})/#it{#pi}) (pb)",
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#Delta#it{#varphi}(#it{e}^{#pm}, #it{#mu}^{#mp})/#it{#pi}) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(#Delta#it{#varphi}(#it{e}^{#pm}, #it{#mu}^{#mp})/#it{#pi}) (adim.)',
    #'genbinning'  : [0., .25, .50, 0.75, 1.0],                        # usado en drafts, en nota, etc. (4 bins)
    #'recobinning' : [0., .125, .25, .375, .50, .625, .75, .875, 1.0], # usado en drafts, en nota, etc. (4 bins)
    #'genbinning'  : [0., .75, 1.5, 2.25, r.TMath.Pi()],
    #'recobinning' : [0., .35, .85, 1.25, 1.65, 2.05, 2.45, 2.85, r.TMath.Pi()],

    #'genbinning'  : [0., .2, .4, .6, .8, 1.0],                        # propuesta (5 bins)
    #'recobinning' : [0., 0.1, .2, 0.3, 0.4, 0.5, .6, 0.7, 0.8, 0.9, 1.0], # propuesta (5 bins)

    'genbinning'  : [0., .17, .33, .5, .67, .83, 1.0],                        # propuesta (6 bins)
    'recobinning' : [0., 0.083, .17, 0.25, 0.33, 0.417, .5, 0.583, 0.67, 0.75, 0.83, 0.917, 1.0], # propuesta (6 bins)

    #'genbinning'  : [0., .14, .29, .42, .57, .71, .86, 1.0],                        # propuesta (7 bins)
    #'recobinning' : [0., 0.07, .14, 0.21, 0.29, 0.36, .43, 0.5, 0.57, .64, 0.71, .79, 0.86, 0.93, 1.0], # propuesta (7 bins)

    #'genbinning'  : [0., .125, .25, .375, .50, .625, .75, .875, 1.0],                        # propuesta (8 bins)
    #'recobinning' : [0., 0.0625, .125, 0.1875, 0.25, 0.3125, .375, 0.4375, 0.50, 0.5625, 0.625, 0.6875, 0.75, 0.8125, 0.875, 0.9375, 1.0], # propuesta (8 bins)

    #'genbinning'  : [0., .11, .22, .33, .44, .56, .67, .78, .89, 1.0],                        # propuesta (9 bins)
    #'recobinning' : [0., 0.056, .11, 0.17, 0.22, 0.28, .33, 0.39, 0.44, 0.5, 0.56, .61, .67, .72, .78, .83, .89, .94, 1.0], # propuesta (9 bins)

    #'genbinning'  : [0., .1, .2, .3, .4, .5, .6, .7, .8, .9, 1.0],                        # propuesta (10 bins)
    #'recobinning' : [0., 0.05, .1, 0.15, 0.2, 0.25, .3, 0.35, 0.4, 0.45, 0.5, 0.55, .6, .65, .7, .75, .8, .85, .9, .95, 1.0], # propuesta (10 bins)

    'var'         : 'abs(TLep1Lep2_DPhi)/' + str(r.TMath.Pi()),
    #'var'         : 'abs(TDPhiLL)',
    'var_response': 'Lep1Lep2_DPhi',
    'var_gen'     : 'abs(TDressLep1Lep2_DPhi)/' + str(r.TMath.Pi()),
    #'var_gen'     : 'abs(TGenDPhiLL)',
    #'legpos'      : (0.82, 0.14, 0.93, 0.47),
    'legpos'      : (0.15, 0.425, 0.35, 0.81),
    'legposdesc'  : (0.15, 0.425, 0.35, 0.81),
    'legpos_fold' : "TL",
    'legpos_fid'  : "BR",
    #'legpos_unf'  : "BC",
    "legpos_unf"   : (.36, .255, .56, .04),
    'legpos_foldas':"TL",
    'legpos_unfas': "TL",
    'legpos_fidbin':"BC",
    'resptxtsize' : 2,
    'covtxtsizeunf': 2,
    'covtxtsizefol': 1.75,
    "covtxtsizefidnorm": 1.6,
    "equalbinsunf" : True,
    "equalbinsfol" : True,
}
varList['Lep1Lep2_DPhiuncertainties'] = {
    'xaxis'       : varList['Lep1Lep2_DPhi']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['Lep1Jet1_DR'] = {
    'xaxis'       : "\\Delta R(\\ell_{1}, j) (adim.)",
    'yaxis'       : "d\\sigma/d(\\Delta R(\\ell_{1}, j)) (pb)",
    'yaxisfid'    : '(1/\\sigma_{fid.})d\\sigma/d(\\Delta R(\\ell_{1}, j)) (adim.)',
    'yaxisfidbin' : '(1/\\sigma_{fid.})d\\sigma/d(\\Delta R(\\ell_{1}, j)) (adim.)',

    'genbinning'  : [0., 1.0, 2.0, 3.0, 4.0, 5.0],
    'recobinning' : [0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0],

    'var'         : 'TLep1Jet1_DR',
    'var_response': 'Lep1Jet1_DR',
    'var_gen'     : 'TDressLep1Jet1_DR',
    'legpos'      : (0.15, 0.425, 0.35, 0.81),
    'legposdesc'  : (0.15, 0.425, 0.35, 0.81),
    'legpos_fold' : "TL",
    'legpos_fid'  : "BR",
    "legpos_unf"   : (.36, .255, .56, .04),
    'legpos_foldas':"TL",
    'legpos_unfas': "TL",
    'legpos_fidbin':"BC",
    'resptxtsize' : 2,
    'covtxtsizeunf': 2,
    'covtxtsizefol': 1.75,
    "covtxtsizefidnorm": 1.6,
    "equalbinsunf" : True,
    "equalbinsfol" : True,
}
varList['Lep1Jet1_DRuncertainties'] = {
    'xaxis'       : varList['Lep1Jet1_DR']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}




varList['Lep12Jet12_DR'] = {
    'xaxis'       : "#Delta R(e^{#pm}#mu^{#mp}, j_{1}j_{2}) (adim.)",
    'yaxis'       : "d#sigma/d(#Delta R(e^{#pm}#mu^{#mp}, j_{1}j_{2})) (pb)",
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#Delta R(e^{#pm}#mu^{#mp}, j_{1}j_{2})) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(#Delta R(e^{#pm}#mu^{#mp}, j_{1}j_{2})) (adim.)',

    'genbinning'  : [0., 1.0, 2.0, 3.0, 4.0, 5.0],
    'recobinning' : [0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0],

    'var'         : 'TLep12Jet12_DR',
    'var_response': 'Lep12Jet12_DR',
    'var_gen'     : 'TDressLep12Jet12_DR',
    'legpos'      : (0.15, 0.425, 0.35, 0.81),
    'legposdesc'  : (0.15, 0.425, 0.35, 0.81),
    'legpos_fold' : "TL",
    'legpos_fid'  : "BR",
    "legpos_unf"   : (.36, .255, .56, .04),
    'legpos_foldas':"TL",
    'legpos_unfas': "TL",
    'legpos_fidbin':"BC",
    'resptxtsize' : 2,
    'covtxtsizeunf': 2,
    'covtxtsizefol': 1.75,
    "covtxtsizefidnorm": 1.6,
    "equalbinsunf" : True,
    "equalbinsfol" : True,
}
varList['Lep12Jet12_DRuncertainties'] = {
    'xaxis'       : varList['Lep12Jet12_DR']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['Lep12Jet12MET_DR'] = {
    'xaxis'       : "#Delta R(e^{#pm}#mu^{#mp}, j_{1}j_{2}p_{T}^{miss}) (adim.)",
    'yaxis'       : "d#sigma/d(#Delta R(e^{#pm}#mu^{#mp}, j_{1}j_{2})) (pb)",
    'yaxisfid'    : '(1/#sigma_{fid.})d#sigma/d(#Delta R(e^{#pm}#mu^{#mp}, j_{1}j_{2}p_{T}^{miss})) (adim.)',
    'yaxisfidbin' : '(1/#sigma_{fid.})d#sigma/d(#Delta R(e^{#pm}#mu^{#mp}, j_{1}j_{2}p_{T}^{miss})) (adim.)',

    'genbinning'  : [0., 1.0, 2.0, 3.0, 4.0, 5.0],
    'recobinning' : [0., 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0],

    'var'         : 'TLep12Jet12MET_DR',
    'var_response': 'Lep12Jet12MET_DR',
    'var_gen'     : 'TDressLep12Jet12MET_DR',
    'legpos'      : (0.15, 0.425, 0.35, 0.81),
    'legposdesc'  : (0.15, 0.425, 0.35, 0.81),
    'legpos_fold' : "TL",
    'legpos_fid'  : "BR",
    "legpos_unf"   : (.36, .255, .56, .04),
    'legpos_foldas':"TL",
    'legpos_unfas': "TL",
    'legpos_fidbin':"BC",
    'resptxtsize' : 2,
    'covtxtsizeunf': 2,
    'covtxtsizefol': 1.75,
    "covtxtsizefidnorm": 1.6,
    "equalbinsunf" : True,
    "equalbinsfol" : True,
}
varList['Lep12Jet12MET_DRuncertainties'] = {
    'xaxis'       : varList['Lep12Jet12MET_DR']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['DPhiLeadJet'] = {
    'xaxis'       : '\\Delta \\varphi(\\ell_{1}, j) (rad)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0, 1., 1.75, 2.5, r.TMath.Pi()],
    'recobinning' : [0, .5, 1., 1.5, 1.75, 2., 2.5, 2.85, r.TMath.Pi()],
    'var'         : 'abs(TDPhiLeadJet)',
    'var_response': 'DPhiLeadJet',
    'var_gen'     : 'abs(TGenDPhiLeadJet)',
}
varList['DPhiLeadJetuncertainties'] = {
    'xaxis'       : varList['DPhiLeadJet']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}

varList['DPhiSubLeadJet'] = {
    'xaxis'       : '\\Delta \\varphi(\\ell_{2}, j) (rad)',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [0, 1., 1.75, 2.5, r.TMath.Pi()],
    'recobinning' : [0, .5, 1., 1.5, 1.75, 2., 2.5, 2.85, r.TMath.Pi()],
    'var'         : 'abs(TDPhiSubLeadJet)',
    'var_response': 'DPhiSubLeadJet',
    'var_gen'     : 'abs(TGenDPhiSubLeadJet)',
}
varList['DPhiSubLeadJetuncertainties'] = {
    'xaxis'       : varList['DPhiSubLeadJet']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['NLooseCentral'] = {
    'xaxis'       : 'Number of loose jets',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [-0.5, 1.5, 3.5, 4.5],
    'recobinning' : [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5],
    'legpos'      : (0.70, 0.55, 0.90, 0.93),
    'var'         : 'TNLooseCentral',
    'var_response': 'NLooseCentral',
    'var_gen'     : 'TDressNLooseCentral',
}
varList['NLooseCentraluncertainties'] = {
    'xaxis'       : varList['NLooseCentral']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



varList['NBLooseCentral'] = {
    'xaxis'       : 'Number of b-tagged loose jets',
    'yaxis'       : 'd#sigma (pb)',
    'genbinning'  : [-0.5, 1.5, 3.5, 4.5],
    'recobinning' : [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5],
    'legpos'      : (0.70, 0.55, 0.90, 0.93),
    'var'         : 'TNBLooseCentral',
    'var_response': 'NBLooseCentral',
    'var_gen'     : 'TDressNBLooseCentral',
}
varList['NBLooseCentraluncertainties'] = {
    'xaxis'       : varList['NBLooseCentral']['xaxis'],
    'yaxis'       : 'Relative uncertainty (adim.)'
}



systMap = {
    'fsrUp' : { 'TW'             : 'TW_noFullyHadr_fsrUp',
                'TbarW'          : 'TbarW_noFullyHadr_fsrUp',
                'TTbar_Powheg'   : 'TTbar_Powheg_fsrUp'
                },
    'fsrDown' : { 'TW'           : 'TW_noFullyHadr_fsrDown',
                  'TbarW'        : 'TbarW_noFullyHadr_fsrDown',
                  'TTbar_Powheg' : 'TTbar_Powheg_fsrDown'
                  },
    'isrUp' : { 'TW'             : 'TW_noFullyHadr_isrUp',
                'TbarW'          : 'TbarW_noFullyHadr_isrUp',
                'TTbar_Powheg'   : 'TTbar_Powheg_isrUp'
                },
    'isrDown' : { 'TW'           : 'TW_noFullyHadr_isrDown',
                  'TbarW'        : 'TbarW_noFullyHadr_isrDown',
                  'TTbar_Powheg' : 'TTbar_Powheg_isrDown'
                  },
    'tWMEUp': { 'TW'           : 'TW_noFullyHadr_MEscaleUp',
                  'TbarW'        : 'TbarW_noFullyHadr_MEscaleUp',
                  'TTbar_Powheg' : 'TTbar_Powheg'
                  },
    'tWMEDown': { 'TW'         : 'TW_noFullyHadr_MEscaleDown',
                  'TbarW'        : 'TbarW_noFullyHadr_MEscaleDown',
                  'TTbar_Powheg' : 'TTbar_Powheg'
                  },
    'tWPSUp': { 'TW'           : 'TW_noFullyHadr_PSscaleUp',
                  'TbarW'        : 'TbarW_noFullyHadr_PSscaleUp',
                  'TTbar_Powheg' : 'TTbar_Powheg'
                  },
    'tWPSDown': { 'TW'           : 'TW_noFullyHadr_PSscaleDown',
                    'TbarW'        : 'TbarW_noFullyHadr_PSscaleDown',
                    'TTbar_Powheg' : 'TTbar_Powheg'
                  },
    'DSUp'        : { 'TW'           : 'TW_noFullyHadr_DS',
                    'TbarW'        : 'TbarW_noFullyHadr_DS',
                    'TTbar_Powheg' : 'TTbar_Powheg'
                  },
    'hDampUp'   : { 'TW'           : 'TW',
                    'TbarW'        : 'TbarW',
                    'TTbar_Powheg' : 'TTbar_Powheg_hdampUp'
                  },
    'hDampDown'   : { 'TW'           : 'TW',
                      'TbarW'        : 'TbarW',
                      'TTbar_Powheg' : 'TTbar_Powheg_hdampDown'
                      },
    'UEUp'        : { 'TW'           : 'TW',
                      'TbarW'        : 'TbarW',
                      'TTbar_Powheg' : 'TTbar_Powheg_ueUp'
                      },
    'UEDown'      : { 'TW'           : 'TW',
                      'TbarW'        : 'TbarW',
                      'TTbar_Powheg' : 'TTbar_Powheg_ueDown'
                      },
    'GluonMoveCRTune' : { 'TW'           : 'TW',
                          'TbarW'        : 'TbarW',
                          'TTbar_Powheg' : 'TTbar_GluonMoveCRTune'
                          },
    'GluonMoveCRTuneerdON' : { 'TW'           : 'TW',
                                'TbarW'        : 'TbarW',
                                'TTbar_Powheg' : 'TTbar_GluonMoveCRTune_erdON'
                                },
    'PowhegerdON' : { 'TW'           : 'TW',
                       'TbarW'        : 'TbarW',
                       'TTbar_Powheg' : 'TTbar_Powheg_erdON'
                       },
    'QCDbasedCRTuneerdON' : { 'TW'           : 'TW',
                               'TbarW'        : 'TbarW',
                               'TTbar_Powheg' : 'TTbar_QCDbasedCRTune_erdON'
                               },
}


ColorMap = {
    'JES'                   : r.TColor.GetColor("#a6cee3"),
    'fsr'                   : r.TColor.GetColor("#1f77b4"),
    'isr'                   : r.TColor.GetColor("#b2df8a"),
    'tWME'                  : r.TColor.GetColor("#33a02c"),
    'DS'                    : r.TColor.GetColor("#e31a1c"),
    'hDamp'                 : r.TColor.GetColor("#fdbf6f"),
    'UE'                    : r.TColor.GetColor("#ff7f00"),
    'ttbarME'               : r.TColor.GetColor("#6a3d9a"),
    'pdf'                   : r.kOrange-6,
    'ColorR'                : r.TColor.GetColor("#cab2d6"),
    'JER'                   : r.kAzure-6,
    'ElecEff'               : r.kMagenta,
    'MuonEff'               : r.kGray+2,
    'Trig'                  : r.kBlue,
    'PU'                    : r.kRed,
    'Btag'                  : r.kGreen+4,
    'Mistag'                : r.kYellow-6,
    'ttbar'                 : r.kTeal-7,
    'Non-WorZ'              : r.kViolet-2,
    'DY'                    : r.kPink+1,
    'VVttbarV'              : r.kSpring-9,
    'Stat.'                 : r.TColor.GetColor("#fb9a99"),
    'TopPt'                 : r.kGreen+1,
    'mtop'                  : r.kBlue+3,
    'Lumi'                  : r.kPink-5,
    'asimov'                : r.kTeal,
}


SysNameTranslator = {
    'JES'                   : "JES",
    'fsr'                   : "FSR",
    'isr'                   : "ISR",
    'tWME'                  : "tW #mu_{R}/#mu_{F}",
    'DS'                    : "tW DS-DR",
    'hDamp'                 : "ME/PS",
    'UE'                    : "UE",
    'ttbarME'               : "t#bar{t} #mu_{R}/#mu_{F}",
    'pdf'                   : "PDF",
    'ColorR'                : "Colour rec.",
    'JER'                   : "JER",
    'ElecEff'               : "Electron eff.",
    'MuonEff'               : "Muon eff.",
    'Trig'                  : "Trigger eff.",
    'TopPt'                 : "Top p_{T} rew.",
    'mtop'                  : "m_{top} unc.",
    'PU'                    : "PU",
    'Btag'                  : "B-tagging",
    'Mistag'                : "Mistagging",
    'ttbar'                 : "t#bar{t} norm.",
    'Non-WorZ'              : "Non-W/Z norm.",
    'DY'                    : "DY norm.",
    'VVttbarV'              : "VV+t#bar{t}V norm.",
    'Stat.'                 : "Statistical",
    'Lumi'                  : "Luminosity",
    'asimov'                : "asimov",
}
