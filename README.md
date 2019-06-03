# PAFnanoAOD

Download the code
====

Download all the code from github.

    git clone https://github.com/Oviedo-PAF/PAFnanoAOD
    cd PAFnanoAOD





Set the enviroment
====

    root6
    source /opt/PAF/PAF_setup.sh


Change the paths
====

Edit the file RunAnalyserPAF.C and modify:
   - The tab in the spreadsheet with samples (SelectedTab) 
   - The default year (2017 by default... it can be also be introduced as an input)
   - The name of the datasets (SelectedDataset), only to run on data


Run the analysis
====

Execute RunAnalyserPAF(TString samplename, TString selector, int nSlots). Example: 

    root -l -b -q 'RunAnalyserPAF.C("ZZ_ext", "Top", 1)'

Or introduce a path to run on a sample not using the spreadsheet (and with a customized cross section for the sample). Examples:

    root -b -q 'RunAnalyserPAF.C("/pool/ciencias/userstorage/juanr/nanoAODv4/feb22/Tree_DoubleMuon_Run2017*", "Top", 40, 0, 0, 1, "Data2017")'
    root -l -q 'RunAnalyserPAF.C("LocalFile:/pool/ciencias/HeppyTreesSummer16/v2/eq1lepton_nanoAOD/jun20/tbarW_inclusiveDecays*", "TTbarSemilep", 40, 0, 0, 1, "xsec:35.85,2016")'



