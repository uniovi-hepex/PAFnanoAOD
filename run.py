'''
  Example:
  python run.py tt -p /pool/ciencias/userstorage/juanr/nanoAODv4/2017/ -s TTTo2L2Nu
  python run.py 2017.cfg
'''

runC = True
doSendPath = True

# Check if ROOT and PAF is loaded...
import imp, os, sys, time, re, argparse
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


def ExecOrder(command, verbose = False, pretend = False):
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


def GetXsec(xsec, s, verbose, isdata):
  if isdata: return 1
  if isinstance(xsec, int): xsec = float(xsec)
  if isinstance(xsec, str):
    xsecdic = loadxsecdic(xsec, verbose)
    if not re.sub("_(([0-9])|([1-9][0-9])|([1-9][0-9][0-9]))$", "", s) in xsecdic.keys():
      print 'ERROR: not found xsec value for sample %s'%s
      xsec = 1
    else: xsec = xsecdic[re.sub("_(([0-9])|([1-9][0-9])|([1-9][0-9][0-9]))$", "", s)]
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
  doPrefire   = 'prefire,'  if IsVarInTree(path+sample, 'PrefireWeight') else ''
  doPUweight  = 'PUweight,' if IsVarInTree(path+sample, 'puWeight') else ''
  doPS        = 'PS,'       if IsVarInTree(path+sample, 'nPSWeight') else ''
  doScale     = 'Scale,'    if IsVarInTree(path+sample, 'nLHEScaleWeight') else ''
  doPDF       = 'PDF,'      if IsVarInTree(path+sample, 'nLHEPdfWeight') else ''
  doJECunc    = 'JECunc,'   if IsVarInTree(path+sample, 'Jet_pt_jesTotalUp') else ''
  useJetPtNom = 'JetPtNom,' if IsVarInTree(path+sample, 'Jet_pt_nom') else ''
  useLepGood  = 'LepGood,'  if IsVarInTree(path+sample, 'nLepGood') else ''
  options += doPUweight + doPrefire + doPS + doScale + doPDF + doJECunc + useJetPtNom + useLepGood + options
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


def RunSamplePAF(selection, path, sample, year = 2018, xsec = 1, nSlots = 1, outname = '', outpath = '', options = '', nEvents = 0, FirstEvent = 0, prefix = 'Tree', verbose = False, pretend = False, dotest = False, sendJobs = False, queue = 'short'):
  if ',' in sample:
    sample.replace(' ', '')
    sample = sample.split(',')
  if not isinstance(sample, list): sample = [sample]

  if not path.endswith('/'): path += '/'


  # Get dictionary with all files in the path directory
  samples = GetSampleList(path, sample)

  nEventsInTree, nGenEvents, nSumOfWeights, isData = GetAllInfoFromFile([path + x for x in samples])
  xsec = GetXsec(xsec, outname, verbose, isData) if not dotest else 1
  isamcatnlo = True if nGenEvents != nSumOfWeights else False
  if isData: isamcatnlo = False

  if verbose:
    stipe = 'MC' if not isData else 'DATA'
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
  SampString  = ''
  pathSamples = ''
  if doSendPath:
    pathSamples = path
    for s in samples: SampString += '%s,'%(s)
  else:
    for s in samples: SampString += '%s/%s,'%(path,s)
  if SampString.endswith(','): SampString = SampString[:-1]

  command = '\'run.C(\"%s\", \"%s\", %f, %f, %i, \"%s\", %i, \"%s\", \"%s\", %i, %i, %i, %i, \"%s\", \"%s\", %i)\''%(SampString, selection, xsec, nSumOfWeights, year, outname, nSlots, outpath, options, isamcatnlo, isData, nEvents, FirstEvent, workingdir, pathSamples, verbose)
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
    ExecOrder('cp -r packages %s'%pathjob, verbose)
    ExecOrder('cp xsec.cfg %s'   %pathjob, verbose)
    ExecOrder('cp run.C %s'      %pathjob, verbose)
    ExecOrder('cp run.py %s'     %pathjob, verbose)
    f = open(jobfile,'w')
    f.write('#!/bin/sh\n\n')
    f.write("sleep 2s\n\n")
    f.write('cd %s\n\n'%pathjob)
    f.write(command + "\n\n")
    f.write("sleep 2s")
    f.close()
    jname = 'PAF%s'%(tag)
    errname = '%sERR%s.out'%(pathjob,tag)
    outname = '%sOUT%s.out'%(pathjob,tag)
    runCommand = "sbatch -p %s -c %i -J %s -e %s -o %s %s"%(queue, nSlots, jname, errname, outname, jobfile)
    ExecOrder(runCommand, verbose, pretend)

  elif runC:
    ExecOrder(command, verbose, pretend)
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
  if nEvents > 0:      myProject.SetNEvents(nEvents);
  if FirstEvent > 0:   myProject.SetFirstEvent(FirstEvent);
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
  if selection != "": 
    myProject.AddSelectorPackage(selection)
  #if selection == 'TopAnalysis':
  #  myProject.AddSelectorPackage("LepEffTop")
 

  # Additional packages
  myProject.AddPackage("Lepton");
  myProject.AddPackage("Jet");
  myProject.AddPackage("mt2");
  myProject.AddPackage("Functions");
  myProject.AddPackage("LeptonSF");
  myProject.AddPackage("BTagSFUtil");

  myProject.Run();


def CheckFileEvents(snames, sfiles, path, outpath, chkdir = None):
  ''' Check events in mother file and output file, using the fhDummy histogram '''
  from ROOT import TFile
  outlist = []
  getlist = lambda x : [x] if not ',' in x else x.replace(' ', '').split(',')
  for f in snames:
    if f not in chkdir:
      trypath = '%s/%s.root'%(outpath, f)
      if not os.path.isfile(trypath):
        outlist.append(f)
        print '\033[0;31mSample \033[0;36m%s\033[0;31m not found\033[0m'%f
      else:
        tf = TFile.Open(trypath)
        if not hasattr(tf, 'fhDummy'):
          print '\033[0;31mSample \033[0;36m%s\033[0;31m does not contain fhDummy\033[0m'%f
          outlist.append(f)
          continue
        dummy = tf.fhDummy.GetEntries()
        samples = GetSampleList(path, getlist(sfiles[f]) )
        nEventsInTree, nGenEvents, nSumOfWeights, isData = GetAllInfoFromFile([path + x for x in samples])
        if nEventsInTree == dummy:
          print '\033[0;32mOK  \033[0;34m%s\033[0m'%f
        else:
          print '\033[0;31mBAD \033[0;34m%s \033[0;33m(%1.2f %s)\033[0m'%(f, (nEventsInTree-dummy)/nEventsInTree*100, '%')
          outlist.append(f)
    else:
      samples      = GetSampleList(path, getlist(sfiles[f]))
      nTrueEntries = GetAllInfoFromFile([path + x for x in samples])[0]
      tmpnchs      = int(chunkdir[f])
      tmpnEvents   = nTrueEntries / tmpnchs
      for ich in range(tmpnchs):
        trypath = '{op}/{nm}_{ch}.root'.format(op = outpath, nm = f, ch = ich)
        if not os.path.isfile(trypath):
          outlist.append((f, ich))
          print '\033[0;31mSample \033[0;36m%s_%s\033[0;31m not found\033[0m'%(f, ich)
        else:
          tf = TFile.Open(trypath)
          if not hasattr(tf, 'fhDummy'):
            print '\033[0;31mSample \033[0;36m%s\033[0;31m does not contain fhDummy\033[0m'%f
            outlist.append((f, ich))
            continue
          dummy   = tf.fhDummy.GetEntries()
          tmpFirstEvent = tmpnEvents * ich
          if (ich == tmpnchs - 1): tmpnEvents = nTrueEntries - tmpFirstEvent
          if tmpnEvents == dummy:
            print '\033[0;32mOK  \033[0;34m%s_%s\033[0m'%(f, ich)
          else:
            print '\033[0;31mBAD \033[0;34m%s \033[0;33m(%1.2f %s)\033[0m'%(f, (tmpnEvents - dummy) / tmpnEvents * 100, '%')
            outlist.append((f, ich))

  return outlist



################################################################################
### Execute
################################################################################
if __name__ == "__main__":
  parser = argparse.ArgumentParser(description='Run with PAF')
  parser.add_argument('--verbose',    '-v', action='store_true', help = 'Activate the verbosing')
  parser.add_argument('--pretend',    '-p', action='store_true', help = 'Create the files but not send the jobs')
  parser.add_argument('--test',       '-t', action='store_true', help = 'Sends only one or two jobs, as a test')
  parser.add_argument('--sendJobs',   '-j', action='store_true', help = 'Send jobs!')
  parser.add_argument('selection'                              , help = 'Name of the selector')
  parser.add_argument('--path',             default=''         , help = 'Path to look for nanoAOD')
  parser.add_argument('--sample',     '-s', default=''         , help = 'Sample(s) to process')
  parser.add_argument('--xsec',       '-x', default='xsec'     , help = 'Cross section')
  parser.add_argument('--year',       '-y', default=-1         , help = 'Year')
  parser.add_argument('--conf',       '-c', default=''         , help = 'Config file (not yet implemented')
  parser.add_argument('--options',    '-o', default=''         , help = 'Options to pass to your analysis')
  parser.add_argument('--prefix',           default='Tree'     , help = 'Prefix of the name...')
  parser.add_argument('--outname',          default=''         , help = 'Name of the output file')
  parser.add_argument('--outpath',          default=''         , help = 'Output path')
  parser.add_argument('--queue',      '-q', default='short'    , help = 'Queue to send jobs')
  parser.add_argument('--firstEvent',       default=0          , help = 'First event')
  parser.add_argument('--nEvents',          default=0          , help = 'Number of events')
  parser.add_argument('--nSlots',     '-n', default=-1         , help = 'Number of slots')
  parser.add_argument('--nChunks',    '-N', default=-1         , help = 'Number of chunks')
  parser.add_argument('--fixedchunk', '-f', default=-1         , help = 'Chunk to be produced alone. It requires to specify the number of chunks in the configuration file. IMPORTANT: is the CHUNK number (from 0 to N-1).')
  parser.add_argument('--check'           , action='store_true', help = 'Check the output trees')
  parser.add_argument('--resubmit'        , action='store_true', help = 'Resubmit jobs')

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
  nSlots      = args.nSlots
  FirstEvent  = args.firstEvent
  sendJobs    = args.sendJobs
  queue       = args.queue
  nEvents     = int(args.nEvents)
  fixedchunk  = int(args.fixedchunk)
  ncores      = nSlots
  doCheck     = args.check
  doReSubmit  = args.resubmit
  if doReSubmit: doCheck = True

  # Check if a cfg file is given as first argument
  fname = selection
  if not os.path.isfile('./'+fname):
    l = filter(lambda x: x[0] == fname, [x.split('.') for x in os.listdir('.')])
    if len(l) != 0:
      l = l[0]
      fname = l[0] + '.' + l[1]
  if os.path.isfile(fname):
    if verbose: print ' >> Using config file \'%s\'...'%fname
    selection   = ''
    spl         = []
    samplefiles = {}
    nslots      = {}
    chunkdir    = {}
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
          nslots[key]      = ncor
          if len(lst) > 3: chunkdir[key] = lst[3]

    # Re-assign arguments...
    if '--verbose' in aarg or '-v' in aarg : verbose     = args.verbose
    if '--pretend' in aarg or '-p' in aarg : pretend     = args.pretend
    if '--test'    in aarg or '-t' in aarg : dotest      = args.test
    if '--sendJobs'in aarg or '-j' in aarg : sendJobs    = args.sendJobs
    if args.path       != ''     : path        = args.path
    if args.options    != ''     : options     = args.options
    if args.xsec       != 'xsec' : xsec        = args.xsec
    if args.year       != -1     : year        = args.year
    if args.prefix     != 'Tree' : prefix      = args.prefix
    if args.outname    != ''     : outname     = args.outname
    if args.outpath    != ''     : outpath     = args.outpath
    if args.nEvents    != 0      : nEvents     = int(args.nEvents)
    if args.firstEvent != 0      : FirstEvent  = int(args.firstEvent)
    if args.queue      != "short": queue       = args.queue

    if args.nSlots     != -1:
      nSlots      = int(args.nSlots)
      for k in nslots.keys(): nslots[k] = nSlots
    elif nSlots == -1:
      nSlots = 1

    if args.nChunks != -1:
      nChunks = int(args.nChunks)
      for k in nslots.keys(): chunkdir[k] = nChunks

    if args.sample     != '':
      sample = args.sample
      if not sample in samplefiles.keys():
        print 'WARNING: Sample \'%s\' not found in cfg file!!'%sample
        samplefiles[sample] = sample
        nslots[sample] = nSlots
      spl = [sample]

    if doCheck:
      outlist = CheckFileEvents(spl, samplefiles, path, outpath, chunkdir)
      if doReSubmit:
        spl = outlist
        if len(spl) == 0: print 'Everything went fine!! :)'
        else:
          print 'Files to resubmit: '
          for f in spl:
            if not isinstance(f, tuple): print ' >> %s'%f
            else:                        print ' >> %s_%s'%(f[0], f[1])
      else:
        exit()

    if dotest:
      nEvents = 1000
      queue = 'cpupower'
      spl = spl[0:1]
      nslots[spl[0]] = 1
      #samplefiles[spl[0]] = [samplefiles[spl[0]][0]]
      outname = 'test'

    for sname in spl:
      if not isinstance(sname, tuple): truesname = sname
      else:
        truesname = sname[0]
        fixedchunk = sname[1]

      outname = truesname
      sample  = samplefiles[truesname]
      ncores  = nslots[truesname]

      if truesname in chunkdir:
        tmpsample = sample
        tmppath   = path
        tmpnchs   = int(chunkdir[truesname])
        if verbose: print " >> The sample {smp} is going to be separated into {chs} chunks.".format(smp = truesname, chs = str(chunkdir[truesname]))
        if ',' in tmpsample:
          tmpsample.replace(' ', '')
          tmpsample = tmpsample.split(',')
        if not isinstance(tmpsample, list): tmpsample = [tmpsample]
        if not tmppath.endswith('/'): tmppath += '/'
        tmpsamples   = GetSampleList(tmppath, tmpsample)
        nTrueEntries = GetAllInfoFromFile([tmppath + x for x in tmpsamples])[0]

        if not sendJobs:
          ExecOrder("sleep 2s")
          ExecOrder("resetpaf")
        if fixedchunk < 0:
          for ich in range(tmpnchs):
            tmpoutname    = outname + "_{ch}".format(ch = ich)
            tmpnEvents    = nTrueEntries / tmpnchs
            tmpFirstEvent = tmpnEvents * ich
            if (ich == tmpnchs - 1): tmpnEvents = nTrueEntries - tmpFirstEvent
            RunSamplePAF(selection, path, sample, year, xsec, ncores, tmpoutname, outpath, options, tmpnEvents, tmpFirstEvent, prefix, verbose, pretend, dotest, sendJobs, queue)
            if not sendJobs:
              ExecOrder("resetpaf")
              ExecOrder("sleep 3s")
        else:
          if verbose: print " >> We only are going to produce the chunk with index {chk}.".format(chk = str(fixedchunk))
          tmpoutname    = outname + "_{ch}".format(ch = fixedchunk)
          tmpnEvents    = nTrueEntries / tmpnchs
          tmpFirstEvent = tmpnEvents * fixedchunk
          if (fixedchunk == tmpnchs - 1): tmpnEvents = nTrueEntries - tmpFirstEvent
          RunSamplePAF(selection, path, sample, year, xsec, ncores, tmpoutname, outpath, options, tmpnEvents, tmpFirstEvent, prefix, verbose, pretend, dotest, sendJobs, queue)

      else:
        RunSamplePAF(selection, path, sample, year, xsec, ncores, outname, outpath, options, nEvents, FirstEvent, prefix, verbose, pretend, dotest, sendJobs, queue)
  else: # no config file...
    RunSamplePAF(selection, path, sample, year, xsec, nSlots, outname, outpath, options, nEvents, FirstEvent, prefix, verbose, pretend, dotest, sendJobs, queue)

