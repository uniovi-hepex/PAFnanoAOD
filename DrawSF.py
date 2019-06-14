import os
from ROOT import TCanvas, TH2F, gStyle, TFile, gSystem, gROOT
gROOT.SetBatch(1)
path = 'Inputfiles'

SF = [
['Run2018ABCD_SF_ID', 'NUM_TightID_DEN_TrackerMuons_pt_abseta', 'MuonID2018'],
['RunBCDEF_SF_ID', 'NUM_TightID_DEN_genTracks_pt_abseta', 'MuonID2017'],
['Run2018ABCD_SF_ISO','NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta','MuonIso2018'],
['RunBCDEF_SF_ISO','NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta','MuonIso2017'],
['2018_ElecRECO', 'EGamma_SF2D','ElecReco2018'],
['2017_ElecRECO', 'EGamma_SF2D','ElecReco2017'],
['2016ReReco_ElecRECO', 'EGamma_SF2D','ElecReco2016LegacyReReco'],
['2018_ElectronTight','EGamma_SF2D','ElecIDIso2018'],
['2017_ElectronTight','EGamma_SF2D','ElecIDIso2017'],
['2016LegacyReReco_ElectronTight','EGamma_SF2D','ElecIDIso2016LegacyReReco']
]

def DrawPlot(fname, hname, outname, outpath = "./SFplots/", inpath = 'InputFiles/', doErrors = True):
  if not outpath.endswith('/'    ): outpath += '/'
  if not  inpath.endswith('/'    ):  inpath += '/'
  if not   fname.endswith('.root'):   fname += '.root'
  f = TFile.Open(inpath+fname)
  h = f.Get(hname)
  c = TCanvas("c", "c", 10, 10, 1600, 1200);
  gStyle.SetOptStat(0);
  gStyle.SetPalette(1);
  gStyle.SetPaintTextFormat("1.2f"); 

  if doErrors: h.Draw("colz, text, errors")
  else       : h.Draw("colz, text")
  if not os.path.isdir(outpath): gSystem.mkdir(outpath, 1)
  c.Print(outpath + outname + '.pdf', 'pdf')
  c.Print(outpath + outname + '.png', 'png')

for fname, hname, outname in SF:
  DrawPlot(fname, hname, outname)
