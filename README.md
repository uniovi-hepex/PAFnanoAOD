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


Run the analysis
====

Execute RunAnalyserPAF(TString samplename, TString selector, int nSlots). Example: 

    root -l -b -q 'RunAnalyserPAF.C("ZZ_ext", "Top", 1)'

