import os, sys, argparse
from ROOT import TFile, TTree, TH1F

defaultPath = '/afs/cern.ch/work/j/jrgonzal/public/Run2017G/skim2l'

def isdigit(a):
  ''' Redefinition of str.isdigit() that takes into account negative numbers '''
  if a.isdigit(): return True
  m = a[0]; n = a[1:]
  if m == '-' and n.isdigit(): return True
  return False

def findValidRootfiles(path, sampleName = '', getOnlyNumberedFiles = False, verbose = False, FullPaths = False):
  ''' Find rootfiles in path with a given name '''
  files = []
  if ',' in sampleName:
    sl = sampleName.replace(' ', '').split(',')
    return findValidRootfiles(path, sl, getOnlyNumberedFiles, verbose, FullPaths)
  elif isinstance(sampleName, list):
    for s in sampleName:
      files += findValidRootfiles(path, s, getOnlyNumberedFiles, verbose, FullPaths)
      #if len(files) == 0:
      #  files += findValidRootfiles(path, 'Tree_'+s, getOnlyNumberedFiles, verbose, FullPaths)
    return files
  if not path[-1] == '/': path += '/'
  if verbose: print ' >> Looking for files in path: ' + path
  for f in os.listdir(path):
    if not f[-5:] == '.root': continue
    if not '_' in f: continue
    n = f[:-5].split('_')[-1]
    s = f[:-5].split('_')[:-1]
    if not isdigit(n): s.append(n)
    fname = ''
    for e in s: fname+=e+'_'
    if fname[-1] == '_': fname = fname[:-1]
    if getOnlyNumberedFiles and not n.isdigit(): continue
    if sampleName != '' and fname != sampleName and (fname+'_'+n) != sampleName: continue
    if verbose: print ' >> Adding file: ', f
    files.append(f)
  if FullPaths: files = [path + x for x in files]
  if len(files) == 0: 
    files = findValidRootfiles(path, 'Tree_' + sampleName, getOnlyNumberedFiles, verbose, FullPaths)
  if len(files) == 0: 
    print '[ERROR]: Not files "' + sampleName + '" found in: ' + path
  return files

def GetFiles(path, name, verbose = False):
  ''' Get all rootfiles in path for a given process name'''
  return findValidRootfiles(path, name, False, verbose, FullPaths = True)

def GetNGenEvents(fname):
  ''' Returns number of events from the 'Count' histograms '''
  if isinstance(fname, list):
    c = 0
    for f in fname: c+=GetNGenEvents(f)
    return c
  elif isinstance(fname, str):
    f = TFile.Open(fname)
    h = f.Get('Count')
    return h.GetBinContent(1)
  else: print '[ERROR] [GetNGenEvents]: wrong input' 

def GetHisto(fname, hname):
  ''' Returns a histogram from files fname '''
  if isinstance(fname, list):
    h0 = GetHisto(fname[0], hname)
    for fi in fname[1:]: h0.Add(GetHisto(fi, hname))
    return h0
  else:
    f = TFile.Open(fname)
    h = f.Get(hname)
    h.SetDirectory(0)
    f.Close()
    return h

def GetSumWeights(fname):
  ''' Returns number of events from the 'SumWeights' histograms '''
  if isinstance(fname, list):
    c = 0
    for f in fname: c+=GetSumWeights(f)
    return c
  elif isinstance(fname, str):
    f = TFile.Open(fname)
    h = f.Get('SumWeights')
    return h.GetBinContent(1)
  else: print '[ERROR] [GetSumWeights]: wrong input' 

def GetHistoFromSetOfFiles(fname, histoname):
  ''' Returns the sum of a histo with a name in a list of files '''
  if isinstance(fname, list):
    h = GetHistoFromSetOfFiles(fname[0], histoname)
    for f in fname[1:]: h.Add(GetHistoFromSetOfFiles(f, histoname))
    return h
  elif isinstance(fname, str):
    f = TFile.Open(fname)
    if not hasattr(f, histoname):
      print '[ERROR] [GetHistoFromSetOfFiles] Histogram \'%s\' does not exist in file %s !!'%(hitoname, fnmae)
    h = f.Get(histoname)
    h.SetDirectory(0)
    return h
  else: print '[ERROR] [GetHistoFromSetOfFiles]: wrong input' 

def GetEntries(fname, treeName = 'Events'):
  ''' Returns number of events from the tree 'Events' in a file '''
  if isinstance(fname, list):
    c = 0
    for f in fname: c+=GetEntries(f, treeName)
    return c
  elif isinstance(fname, str):
    f = TFile.Open(fname)
    t = f.Get(treeName)
    return t.GetEntries()
  else: print '[ERROR] [GetEntries]: wrong input' 

def GuessIsData(fname):
  ''' Guess if a tree is data or simulation '''
  if isinstance(fname, list): fname = fname[0] # Assume all files are the same process/dataset
  f = TFile.Open(fname)
  t = f.Get('Events')
  if hasattr(t,'genWeight'): return False
  return True

def guessPathAndName(p):
  ''' Guess path and sample name for a given rootfile '''
  path = ''; n = -1
  while '/' in p:
    path += p[:p.index('/')+1]
    p     = p[p.index('/')+1:]
  if p[-5:] == '.root': p = p[:-5]
  elif os.path.isdir(path + p):
    path = path + p
    p = ''
  if len(path) > 0 and not path[-1] == '/': path += '/'
  if '_' in p: 
    n = p.split('_')[-1]
    s = p.split('_')[:-1]
    if not isdigit(n): 
      s.append(n)
      n = '-1'
    p = ''
    for e in s: p+=e+'_'
    if p[-1] == '_': p = p[:-1]
  return path, p, n

def guessProcessName(fileName):
  ''' Guess the name of the process for a given file '''
  if isinstance(fileName, list): 
    path, name, n = guessPathAndName(fileName[0])
    fileName = name
    if fileName[-5:] == '.root': fileName = fileName[:-5]
    while '/' in fileName: fileName = fileName[fileName.index('/')+1:]
  return fileName

def groupFilesInDic(listOfFiles, name, i=-1, verbose = False):
  ''' Manages a dictionary with sample names and lists of samples '''
  if isinstance(name, list):
    for e in name:
      path, nam, n = guessPathAndName(e)  
      groupFilesInDic(listOfFiles, nam, n)
    return
  fname = name + '_' + str(i) + '.root' if str(i).isdigit() else name + '.root'
  if name in listOfFiles: listOfFiles[name].append(fname)
  else: 
    newList = [fname]
    listOfFiles[name] = newList
    if verbose: print ' >> Sample found: ' + name

def getDicFiles(inFolder):
  ''' Get a dictionary with sample names and lists of files  '''
  listOfFiles = {}
  files = findValidRootfiles(inFolder)
  groupFilesInDic(listOfFiles,files)
  return listOfFiles
  
def GetAllInfoFromFile(fname, treeName = 'Events'):
  ''' Returns a list with all the info of a file ''' 
  if isinstance(fname, list):
    nEvents = 0
    nGenEvents = 0
    nSumOfWeights = 0
    isData = False
    for f in fname: 
      iE, iG, iS, isData = GetAllInfoFromFile(f, treeName)
      nEvents += iE
      nGenEvents += iG
      nSumOfWeights += iS
    return [nEvents, nGenEvents, nSumOfWeights, isData]
  elif isinstance(fname, str):
    f = TFile.Open(fname)
    t = f.Get(treeName)
    hs = f.Get('SumWeights')
    hc = f.Get('Count')
    nEvents = t.GetEntries()
    nGenEvents = hc.GetBinContent(1) if isinstance(hc,TH1F) else 1
    nSumOfWeights = hs.GetBinContent(1) if isinstance(hs,TH1F) else 1
    isData = not hasattr(t,'genWeight')
    return [nEvents, nGenEvents, nSumOfWeights, isData]
  else: print '[ERROR] [GetAllInfoFromFile]: wrong input' 

def GetProcessInfo(path, process='', treeName = 'Events'):
  ''' Prints all info from a process in path '''
  if isinstance(path, list): 
    files = path
    path, process, k = guessPathAndName(files[0])
  else: files = GetFiles(path, process)
  nEvents, nGenEvents, nSumOfWeights, isData = GetAllInfoFromFile(files, treeName)
  fileType = '(Data)' if isData else ('(MC)')
  print '\n##################################################################'
  print ' path: ' + path
  print ' Process:            ' + process + ' ' + fileType
  print ' Number of files:    ' + str(len(files))
  print ' Total entries:      ' + str(nEvents)
  if isData:
    print ' Triggered events:   ' + str(nGenEvents)
  else: 
    print ' Generated events:   ' + str(nGenEvents)
    print ' Sum of gen weights: ' + str(nSumOfWeights)
  print '##################################################################\n'

def IsVarInTree(fname, var, treeName = 'Events'):
  ''' Check if a given file and tree contains a branch '''
  if not os.path.isfile(fname):
    print 'ERROR: %s does not exists!'%fname
    return False
  f = TFile.Open(fname)
  t = f.Get(treeName)
  return hasattr(t, var)

def GetValOfVarInTree(fname, var, treeName = 'Events'):
  ''' Check the value of a var in a tree '''
  if not os.path.isfile(fname):
    print 'ERROR: %s does not exists!'%fname
    return False
  f = TFile.Open(fname)
  t = f.Get(treeName)
  t.GetEntry(1)
  return getattr(t,var)


##################################
# Extra functions to work check .root files from terminal

def addToListOfFiles(listOfFiles, name, i):
  ''' Manages a dictionary with sample names and lists of samples '''
  fname = name + '_' + str(i) + '.root'
  if name in listOfFiles: listOfFiles[name].append(fname)
  else: 
    newList = [fname]
    listOfFiles[name] = newList
    if verbose: print ' >> Sample found: ' + name

def main():
 # Parsing options
 path = defaultPath
 sample = ''
 pr = argparse.ArgumentParser()
 pr.add_argument('path', help='Input folder', type = str, default = defaultPath)
 pr.add_argument('--sample', type = str, default = '')
 pr.add_argument('-p','--inspect', action='store_true', help='Print branches')
 pr.add_argument('-t','--treeName', default='Events', help='Name of the tree')
 args = pr.parse_args()
 if args.sample:  sample = args.sample
 treeName = args.treeName
 printb = args.inspect
 path = args.path
 if os.path.isdir(path) and not path[-1] == '/': path += '/'

 if sample == '':
   origpath = path
   path, sample, n = guessPathAndName(path)

   if sample == '': 
     d = getDicFiles(path)
     for c in d:
       print ' >> ' + c + ': ', d[c]

   else:
     totfile = path + sample + '_' + n + '.root' if int(n) >= 0 else path + sample + '.root'
     if os.path.isfile(totfile): 
       GetProcessInfo([totfile], treeName = treeName)
       exit()
     else:
       GetProcessInfo(path, sample, treeName)
 else:
   GetProcessInfo(path, sample, treeName)
   exit()

if __name__ == '__main__':
  main()
