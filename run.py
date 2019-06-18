'''
  Example: 
  python run.py tt -p /pool/ciencias/userstorage/juanr/nanoAODv4/2017/ -s TTTo2L2Nu
  python run.py 2017.cfg
'''

runC = True

# Check if ROOT and PAF is loaded...
import imp, os, sys, time
try:
    imp.find_module('ROOT')
    found = True
except ImportError:
  print 'Please, load ROOT... (typically by executing \'root6\')'
  exit()

from ROOT import gROOT
gROOT.SetBatch(1)

try:
  from ROOT import PAFProject
except ImportError:
  print 'Please, load PAF... (typically by executing \'source /opt/PAF/PAF_setup.sh\')'
  exit()

from fileReader import getDicFiles, GetAllInfoFromFile, IsVarInTree

def ex(command, verbose = False, pretend = False):
  if verbose:
    print 'Executing command: ', command
  if pretend: 
    print 'Pretending...'
    return
  os.system(command)

def loadxsecdic(fname, verbose):
  xsecdir = {}
  if not os.path.isfile('./'+fname):
    l = filter(lambda x: x[0] == fname, [x.split('.') for x in os.listdir('.')])
    if len(l) == 0:
      print 'ERROR: not found file %s with cross sections...'
      return
    else: l = l[0]
    fname = l[0] + '.' + l[1]
  if verbose: print ' >> Reading cross section from %s...'%fname
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

def GetXsec(xsec, s, verbose):
  if isinstance(xsec, int): xsec = float(xsec)
  if isinstance(xsec, str):
    xsecdic = loadxsecdic(xsec, verbose)
    if not s in xsecdic.keys():
      print 'ERROR: not found xsec value for sample %s'%s
      xsec = 1
    else: xsec = xsecdic[s]
  return xsec

def GetSampleList(path, sample):
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
  return samples

def GetOptions(path, sample, options = ""):
  if not path.endswith('/'): path += '/'
  if not sample.endswith(".root"): sample += '.root'
  doPUweight  = 'PUweight,' if IsVarInTree(path+sample, 'puWeight') else ''
  doJECunc    = 'JECunc,'   if IsVarInTree(path+sample, 'Jet_pt_jesTotalUp') else ''
  useJetPtNom = 'JetPtNom,' if IsVarInTree(path+sample, 'Jet_pt_nom') else ''
  useLepGood  = 'LepGood,'  if IsVarInTree(path+sample, 'nLepGood') else ''
  options += doPUweight + doJECunc + useJetPtNom + useLepGood + options
  if options.endswith(','): options = options[:-1]
  return options

def GetTStringVectorSamples(path, samples):
  from ROOT import vector, TString, gSystem
  # Add the input data files
  v = vector(TString)()
  for s in samples: 
    t = TString(path+s)
    v.push_back(t)
  return v
  v = GetTStringVector(samples)


def RunSamplePAF(selection, path, sample, year = 2018, xsec = 1, nSlots = 1, outname = '', outpath = '', options = '', nEvents = 0, FirstEvent = 0, prefix = 'Tree', verbose = False, pretend = False, dotest = False, sendJobs = False, queu = 'short'):
  if ',' in sample:
    sample.replace(' ', '')
    sample = sample.split(',')
  if not isinstance(sample, list): sample = [sample]
  
  if not path.endswith('/'): path += '/'
  
  xsec = GetXsec(xsec, outname, verbose)
  
  # Get dictionary with all files in the path directory
  samples = GetSampleList(path, sample)
  
  nEventsInTree, nGenEvents, nSumOfWeights, isData = GetAllInfoFromFile([path + x for x in samples])
  isamcatnlo = True if nGenEvents != nSumOfWeights else False
  stipe = 'MC' if not isData else 'DATA'
  if verbose: 
    print '## Processing a %i %s sample...'%(year, stipe)
    print '## Files found in %s'%path
    for f in samples: print ' >> %s'%f
    print '## Total events:      %i'%nEventsInTree
    print '## Generated events:  %i'%nGenEvents
    if isamcatnlo:
      print '## Sum of weights:    %1.2f'%nSumOfWeights
      print '## This sample has gen weights!!'

  # Check what is in the tree...
  options = GetOptions(path, samples[0], options)
  workingdir = os.getcwd()
  outpath = workingdir + '/' + outpath
  if verbose:
    print '## Runing with options: %s'%(options)
    print '## Output to: %s'%(outpath + "/" + outname)

  SamplesVector = GetTStringVectorSamples(path,samples)
  SampString    = ''
  for s in samples: SampString += '%s/%s,'%(path,s)
  if SampString.endswith(','): SampString = SampString[:-1]


  command = '\'run.C(\"%s\", \"%s\", %f, %f, %i, \"%s\", %i, \"%s\", \"%s\", %i, %i, %i, %i, \"%s\")\''%(SampString, selection, xsec, nSumOfWeights, year, outname, nSlots, outpath, options, isamcatnlo, isData, nEvents, FirstEvent, workingdir)
  command = 'root -l -b -q ' + command

  if sendJobs:
    tag = "%i"%(int(time.time()*1e6)%1e12)
    if verbose: print 'Creating job with tag: %s'%tag
    jobname = 'job_%s'%tag
    path = os.getcwd()+'/jobs/'
    if not os.path.isdir(path): os.mkdir(path)
    pathjob = path + jobname + '/'
    jobfile = "%s%s.sh"%(pathjob,jobname)
    os.mkdir(pathjob)
    ex('cp -r packages %s'%pathjob, verbose)
    ex('cp xsec.cfg %s'   %pathjob, verbose)
    ex('cp run.C %s'      %pathjob, verbose)
    ex('cp run.py %s'     %pathjob, verbose)
    f = open(jobfile,'w')
    f.write('#!/bin/sh\n\n')
    f.write('cd %s\n\n'%pathjob)
    f.write(command)
    f.close()
    jname = 'PAF%s'%(tag)
    errname = '%sERR%s.out'%(pathjob,tag)
    outname = '%sOUT%s.out'%(pathjob,tag)
    runCommand = "sbatch -p %s -c %i -J %s -e %s -o %s %s"%(queue, nSlots, jname, errname, outname, jobfile)
    ex(runCommand, verbose, pretend)

  elif runC:
    ex(command, verbose, pretend)
  else:
    if pretend: 
      print 'Pretending...'
      return
    RunPAF(SamplesVector, selection, xsec, nSumOfWeights, year, outname, nSlots, outpath, options, isamcatnlo, isData, nEvents, FirstEvent, workingdir)


def RunPAF(samples, selection, xsec, nSumOfWeights, year, outname, nSlots = 1, outpath = '', options = '', isamcatnlo = False, isData = False, nEvents = 0, FirstEvent = 0, workingdir = ''):

  from ROOT import PAFProject, PAFIExecutionEnvironment ,PAFSequentialEnvironment, PAFPROOFLiteEnvironment, PAFPoDEnvironment
  from ROOT import vector, TString, gSystem
  
  # PAF mode selection (based on number of slots)
  pafmode = PAFSequentialEnvironment();
  if   nSlots <=  1: pafmode = PAFSequentialEnvironment();
  else             : pafmode = PAFPROOFLiteEnvironment(nSlots);
  
  myProject = PAFProject(pafmode); 
  
  myProject.AddDataFiles(samples); 
  myProject.SetDefaultTreeName("Events");
  
  # Deal with first and last event
  if nEvents > 0    : myProject.SetNEvents(nEvents);
  if FirstEvent > 0 : myProject.SetFirstEvent(FirstEvent);
  if workingdir == "": workingdir = os.getcwd()
  
  # Set output file
  if outpath == '': outpath = workingdir + "/" + selection + "_temp"
  gSystem.mkdir(outpath, 1);
  myProject.SetOutputFile(outpath + "/" + outname + ".root");
  if verbose: print '## Output to: %s'%(outpath + "/" + outname + ".root")
  
  # Parameters for the analysis
  myProject.SetInputParam("sampleName",        outname);
  myProject.SetInputParam("IsData",            isData    );
  myProject.SetInputParam("weight",            xsec/nSumOfWeights);
  myProject.SetInputParam("IsMCatNLO",         isamcatnlo);
  myProject.SetInputParam("selection",         selection);
  myProject.SetInputParam("WorkingDir",        workingdir);
  myProject.SetInputParam("xsec",              xsec);
  myProject.SetInputParam("_options",          options);
  myProject.SetInputParam("year",              str(year));
  
  # Name of analysis class
  myProject.AddSelectorPackage("LeptonSelector");
  myProject.AddSelectorPackage("JetSelector");
  myProject.AddSelectorPackage("EventBuilder");
  
  # Analysis selector
  if selection != "": myProject.AddSelectorPackage(selection)
  
  # Additional packages
  myProject.AddPackage("Lepton");
  myProject.AddPackage("Jet");
  myProject.AddPackage("mt2");
  myProject.AddPackage("Functions");
  myProject.AddPackage("LeptonSF");
  myProject.AddPackage("BTagSFUtil");
  
  myProject.Run();

  
################################################################################
### Execute
################################################################################
import argparse
parser = argparse.ArgumentParser(description='Run with PAF')
parser.add_argument('--verbose','-v'    , action='store_true'  , help = 'Activate the verbosing')
parser.add_argument('--pretend','-p'    , action='store_true'  , help = 'Create the files but not send the jobs')
parser.add_argument('--test','-t'       , action='store_true'  , help = 'Sends only one or two jobs, as a test')
parser.add_argument('--sendJobs','-j'   , action='store_true'  , help = 'Send jobs!')
parser.add_argument('selection'                                , help = 'Name of the selector')
parser.add_argument('--path'            , default=''           , help = 'Path to look for nanoAOD')
parser.add_argument('--sample','-s'     , default=''           , help = 'Sample(s) to process')
parser.add_argument('--xsec','-x'       , default='xsec'       , help = 'Cross section')
parser.add_argument('--year','-y'       , default=-1           , help = 'Year')
parser.add_argument('--conf','-c'       , default=''           , help = 'Config file (not yet implemented')
parser.add_argument('--options','-o'    , default=''           , help = 'Options to pass to your analysis')
parser.add_argument('--prefix'          , default='Tree'       , help = 'Prefix of the name...')
parser.add_argument('--outname'         , default=''           , help = 'Name of the output file')
parser.add_argument('--outpath'         , default=''           , help = 'Output path')
parser.add_argument('--queue'           , default='short'      , help = 'Queue to send jobs')
parser.add_argument('--firstEvent'      , default=0            , help = 'First event')
parser.add_argument('--nEvents'         , default=0            , help = 'Number of events')
parser.add_argument('--nSlots','-n'     , default=1            , help = 'Number of slots')

args = parser.parse_args()
aarg = sys.argv
selection   = args.selection
verbose     = args.verbose
pretend     = args.pretend
dotest      = args.test
sample      = args.sample
path        = args.path
options     = args.options
xsec        = args.xsec
year        = args.year
prefix      = args.prefix
outname     = args.outname
outpath     = args.outpath
nEvents     = args.nEvents
nSlots      = args.nSlots
FirstEvent  = args.firstEvent
sendJobs    = args.sendJobs
queue       = args.queue
ncores = nSlots

# Check if a cfg file is given as first argument
fname = selection
if not os.path.isfile('./'+fname):
  l = filter(lambda x: x[0] == fname, [x.split('.') for x in os.listdir('.')])
  if len(l) != 0: 
    l = l[0]
    fname = l[0] + '.' + l[1]
if os.path.isfile(fname):
  if verbose: print ' >> Using config file \'%s\'...'%fname
  selection = ''
  spl = []
  samplefiles = {}
  nslots = {}
  f = open(fname)
  lines = f.readlines()
  for l in lines:
    l = l.replace(' ', '')
    l = l.replace('\n', '')
    if l.startswith('#'): continue
    if '#' in l: l = l.split('#')[0]
    if l == '': continue
    if l.endswith(':'): l = l[:-1]
    if not ':' in l: 
      if   l == 'verbose': verbose = 1
      elif l == 'pretend': pretend = 1
      elif l == 'test'   : dotest  = 1
      elif l in ['path', 'sample', 'options', 'selection', 'xsec', 'prefix', 'outpath', 'year', 'nSlots', 'nEvents', 'firstEvent', 'queue']: continue
      else:
        spl.append(l)
        samplefiles[l]=l
        nslots[l]=nSlots
    else:
      lst = l.split(':')
      key = lst[0]
      val = lst[1] if lst[1] != '' else lst[0]
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
      elif key == 'queue'     : queue     = val
      elif key == 'year'      : year      = int(val)
      elif key == 'nSlots'    : nSlots    = int(val)
      elif key == 'nEvents'   : nEvents   = int(val)
      elif key == 'firstEvent': FirstEvent= int(val)
      else:
        ncor = nSlots if len(lst) < 3 else int(lst[2])
        spl.append(key)
        samplefiles[key] = val
        nslots[key] = ncor

  # Re-assign arguments...
  if '--verbose' in aarg or '-v' in aarg : verbose     = args.verbose
  if '--pretend' in aarg or '-p' in aarg : pretend     = args.pretend
  if '--test'    in aarg or '-t' in aarg : dotest      = args.test
  if '--sendJobs'in aarg or '-j' in aarg : sendJobs    = args.sendJobs
  if args.sample     != ''     : sample      = args.sample
  if args.path       != ''     : path        = args.path
  if args.options    != ''     : options     = args.options
  if args.xsec       != 'xsec' : xsec        = args.xsec
  if args.year       != -1     : year        = args.year
  if args.prefix     != 'Tree' : prefix      = args.prefix
  if args.outname    != ''     : outname     = args.outname
  if args.outpath    != ''     : outpath     = args.outpath
  if args.nEvents    != 0      : nEvents     = args.nEvents
  if args.nSlots     != 1      : nSlots      = args.nSlots
  if args.firstEvent != 0      : FirstEvent  = args.firstEvent
  if args.queue      != 0      : queue       = args.queue

  for sname in spl:
    outname = sname
    sample  = samplefiles[sname]
    ncores  = nslots[sname]
    RunSamplePAF(selection, path, sample, year, xsec, ncores, outname, outpath, options, nEvents, FirstEvent, prefix, verbose, pretend, dotest, sendJobs, queue)

else: # no config file...
  RunSamplePAF(selection, path, sample, year, xsec, nSlots, outname, outpath, options, nEvents, FirstEvent, prefix, verbose, pretend, dotest, sendJobs, queue)
