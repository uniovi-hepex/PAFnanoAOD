'''
  Example: 
  python run.py tt -p /pool/ciencias/userstorage/juanr/nanoAODv4/2017/ -s TTTo2L2Nu
'''
import os, sys
from fileReader import getDicFiles, GetAllInfoFromFile

def loadxsecdic(fname):
  xsecdir = {}
  if not os.path.isfile('./'+fname):
    l = filter(lambda x: x[0] == fname, [x.split('.') for x in os.listdir('.')])
    if len(l) == 0:
      print 'ERROR: not found file %s with cross sections...'
      return
    else: l = l[0]
    fname = l[0] + '.' + l[1]
  print ' >> Reading cross section from %s...'%fname
  f = open(fname)
  lines = f.readlines()
  for l in lines:
    l = l.replace(' ', '')
    l = l.replace('\n', '')
    if l.startswith('#'): continue
    if not ':' in l: continue
    if '#' in l: l = l.split('#')[0]
    if l == '': continue
    lst = l.split(':')
    key = lst[0]
    val = lst[1]
    if val == '': val = 1
    xsecdir[key] = float(val)
  return xsecdir 
    
def RunSamplePAF(selection, path, sample, year = 2018, xsec = 1, nSlots = 1, outname = '', outpath = '', nEvents = 0, FirstEvent = 0, prefix = 'Tree', verbose = True, pretend = False, dotest = False):

  if ',' in sample:
    sample.replace(' ', '')
    sample = sample.split(',')
  if not isinstance(sample, list): sample = [sample]
  
  if not path.endswith('/'): path += '/'
  
  if isinstance(xsec, int): xsec = float(xsec)
  if isinstance(xsec, str): 
    xsecdic = loadxsecdic(xsec)
    s = sample[0]
    if not s in xsecdic.keys():
      print 'ERROR: not found xsec value for sample %s'%s
      xsec = 1
    else: xsec = xsecdic[s]
  
  # Get dictionary with all files in the path directory
  dic = getDicFiles(path)
  nfileInPath = len(dic)
  if verbose: print 'Found %i files in path %s'%(nfileInPath, path)
  
  samples = []
  for s in sample: 
    dk = dic.keys()
    if not s in dk: s = prefix+'_'+s
    if not s in dk:
      print 'WARNING: file %s not in path %s'%(s, path)
    else:
      samples += dic[s]
  
  nEventsInTree, nGenEvents, nSumOfWeights, isData = GetAllInfoFromFile([path + x for x in samples])
  isamcatnlo = True if nGenEvents != nSumOfWeights else False
  
  stipe = 'MC' if not isData else 'DATA'
  if verbose: 
    print '## Procesing a %i %s sample...'%(year, stipe)
    print '## Files found in %s'%path
    for f in samples: print ' >> %s'%f
    print '## Total events:      %i'%nEventsInTree
    print '## Generated events:  %i'%nGenEvents
    if isamcatnlo:
      print '## Sum of weights:    %1.2f'%nSumOfWeights
      print '## This sample has gen weights!!'
  
  if pretend: exit()
  
  from ROOT import PAFProject, PAFIExecutionEnvironment ,PAFSequentialEnvironment, PAFPROOFLiteEnvironment, PAFPoDEnvironment
  from ROOT import vector, TString, gSystem
  
  # PAF mode selection (based on number of slots)
  pafmode = PAFSequentialEnvironment();
  if   nSlots <=  1: pafmode = PAFSequentialEnvironment();
  elif nSlots <= 64: pafmode = PAFPROOFLiteEnvironment(nSlots);
  else             : pafmode = PAFPoDEnvironment(nSlots);
  
  myProject = PAFProject(pafmode); 
  
  # Add the input data files
  v = vector(TString)()
  for s in samples: 
    t = TString(s)
    v.push_back(t)
  myProject.AddDataFiles(v); 
  myProject.SetDefaultTreeName("Events");
  
  # Deal with first and last event
  if nEvents > 0    : myProject.SetNEvents(nEvents);
  if FirstEvent > 0 : myProject.SetFirstEvent(FirstEvent);
  
  # Set output file
  if outpath == '': outpath = "./"+selection+"_temp"
  gSystem.mkdir(outpath, 1);
  myProject.SetOutputFile(outpath + "/" + outname + ".root");
  print '## Output to: %s'%(outpath + "/" + outname + ".root")
  
  # Parameters for the analysis
  myProject.SetInputParam("sampleName",        outname);
  myProject.SetInputParam("IsData",            isData    );
  myProject.SetInputParam("weight",            xsec/nSumOfWeights);
  myProject.SetInputParam("IsMCatNLO",         isamcatnlo);
  myProject.SetInputParam("selection",         selection);
  myProject.SetInputParam("WorkingDir",        os.getcwd());
  myProject.SetInputParam("nEntries",          nEventsInTree);
  myProject.SetInputParam("xsec",              xsec);
  myProject.SetInputParam("_options",          options);
  myProject.SetInputParam("Year",              year);
  
  # Name of analysis class
  myProject.AddSelectorPackage("LeptonSelector");
  #myProject.AddSelectorPackage("JetSelector");
  #myProject.AddSelectorPackage("EventBuilder");
  
  # Analysis selector
  if(selection == 'ttbar' or selection == 'TT' or selection == 'ttxsec'): selection = 'tt'
  if(selection == 'TWTT' or selection == 'WWbb' or selection == 'twtt'):  selection = 'tWtt'
  #if    selection == "tt"  :  myProject.AddSelectorPackage("TopAnalysis");
  #elif  selection == "tWtt":  myProject.AddSelectorPackage("TWTTbarAnalysis");
  #else: print "UNKNOWN SELECTOR."
  
  # Additional packages
  myProject.AddPackage("Lepton");
  myProject.AddPackage("Jet");
  myProject.AddPackage("mt2");
  myProject.AddPackage("Functions");
  myProject.AddPackage("LeptonSF");
  myProject.AddPackage("BTagSFUtil");
  #myProject.AddPackage("PUWeight");
  
  myProject.Run();

  
################################################################################
### Execute
################################################################################
import argparse
parser = argparse.ArgumentParser(description='Run with PAF')
parser.add_argument('--verbose','-v'    , action='store_true'  , help = 'Activate the verbosing')
parser.add_argument('--pretend'         , action='store_true'  , help = 'Create the files but not send the jobs')
parser.add_argument('--test','-t'       , action='store_true'  , help = 'Sends only one or two jobs, as a test')
parser.add_argument('selection'                                , help = 'Name of the selector')
parser.add_argument('--path','-p'       , default=''           , help = 'Path to look for nanoAOD')
parser.add_argument('--sample','-s'     , default=''           , help = 'Sample(s) to process')
parser.add_argument('--xsec','-x'       , default='xsec'       , help = 'Cross section')
parser.add_argument('--year','-y'       , default=2018         , help = 'Year')
parser.add_argument('--conf','-c'       , default=''           , help = 'Config file (not yet implemented')
parser.add_argument('--options','-o'    , default=''           , help = 'Options to pass to your analysis')
parser.add_argument('--analysis','-a'   , default='Top'        , help = 'Analysis name')
parser.add_argument('--prefix'          , default='Tree'       , help = 'Prefix of the name...')
parser.add_argument('--outname'         , default=''           , help = 'Name of the output file')
parser.add_argument('--outpath'         , default=''           , help = 'Output path')
parser.add_argument('--firstEvent'      , default=0            , help = 'First event')
parser.add_argument('--nEvents'         , default=0            , help = 'Number of events')
parser.add_argument('--nSlots','-n'     , default=1            , help = 'Number of slots')

args = parser.parse_args()

verbose     = args.verbose
pretend     = args.pretend
dotest      = args.test
sample      = args.sample
path        = args.path
options     = args.options
xsec        = args.xsec
year        = args.year
selection   = args.selection
prefix      = args.prefix
outname     = args.outname
outpath     = args.outpath
nEvents     = args.nEvents
nSlots      = args.nSlots
FirstEvent  = args.firstEvent

# Check if ROOT and PAF is loaded...
import imp
try:
    imp.find_module('ROOT')
    found = True
except ImportError:
  print 'Please, load ROOT... (typically by executing \'root6\')'
  exit()

try:
  from ROOT import PAFProject
except ImportError:
  print 'Please, load PAF... (typically by executing \'source /opt/PAF/PAF_setup.sh\')'
  exit()


# Check if a cfg file is given as first argument
fname = selection
if not os.path.isfile('./'+fname):
  l = filter(lambda x: x[0] == fname, [x.split('.') for x in os.listdir('.')])
  if len(l) != 0: 
    l = l[0]
    fname = l[0] + '.' + l[1]
if os.path.isfile(fname):
  print ' >> Using config file \'%s\'...'%fname
  spl = []
  samplefiles = {}
  f = open(fname)
  lines = f.readlines()
  for l in lines:
    l = l.replace(' ', '')
    l = l.replace('\n', '')
    if l.startswith('#'): continue
    if '#' in l: l = l.split('#')[0]
    if l == '': continue
    if not ':' in l: 
      if   l == 'verbose': verbose = 1
      elif l == 'pretend': pretend = 1
      elif l == 'test'   : dotest  = 1
      else:
        spl.append(l)
        samplefiles[l]=l
    else:
      lst = l.split(':')
      key = lst[0]
      val = lst[1]
      if   key == 'pretend'   : pretend   = 1
      elif key == 'verbose'   : verbose   = 1
      elif key == 'test'      : dotest    = 1
      elif key == 'path'      : path      = val
      elif key == 'sample'    : sample    = val
      elif key == 'options'   : options   = val
      elif key == 'selection' : selection = val
      elif key == 'xsec'      : xsec      = val
      elif key == 'prefix'    : prefix    = val
      elif key == 'outpath'   : outpath   = val
      elif key == 'year'      : year      = int(val)
      elif key == 'nSlots'    : nSlots    = int(val)
      elif key == 'nEvents'   : nEvents   = int(val)
      elif key == 'firstEvent': FirstEvent= int(val)
      else:
        spl.append(key)
        samplefiles[key] = val

  for sname in spl:
    outname = sname
    sample  = samplefiles[sname]
    RunSamplePAF(selection, path, sample, year, xsec, nSlots, outname, outpath, nEvents, FirstEvent, prefix, verbose, pretend, dotest)

else: # no config file...
  RunSamplePAF(selection, path, sample, year, xsec, nSlots, outname, outpath, nEvents, FirstEvent, prefix, verbose, pretend, dotest)
 
