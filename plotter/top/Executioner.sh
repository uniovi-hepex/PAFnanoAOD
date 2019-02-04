# Preliminary definitions of samples and so on

samples=( #"TTTo2L2Nu" "TTToSemiLeptonic"
#          "tbarW_noFullHad" "tW_noFullHad"
#          "DYJetsToLL_M_10to50_MLM" "DYJetsToLL_M_50" #"DYJetsToLL_M5to50_madgraphMLM"
#          "WJetsToLNu_MLM" "WWTo2L2Nu"
#          "ZZTo4L" "ZZTo2L2Nu"
#          "ZZZ"
#          "WZ" "WZTo2L2Q" "WZTo3LNu"
#          "TTWJetsToQQ" "TTWJetsToLNu"
#          "TTZToQQ" "TTZToLLNuNu_M_10" "TTZToLL_M_1to10"
         "MuonEG" "SingleElectron" "SingleMuon" "DoubleMuon" "DoubleEG" "MET")

samples_syst=( #"TTTo2L2Nu_TuneCP5up" "TTTo2L2Nu_TuneCP5down"
              "TTTo2L2Nu_hdampUP" "TTTo2L2Nu_hdampDOWN")

runsamples=( #"TTTo2L2Nu" "TTToSemiLeptonic"
#             "tbarW_noFullHad_ext1 & tbarW_noFullHad" "tW_noFullHad"
#             "DYJetsToLL_M_10to50_MLM" "DYJetsToLL_M_50 & DYJetsToLL_M_50_ext1" #"DYJetsToLL_M5to50_madgraphMLM"
#             "WJetsToLNu_MLM" "WWTo2L2Nu"
#             "ZZTo4L & ZZTo4L_ext1" "ZZTo2L2Nu"
#             "ZZZ"
#             "WZ" "WZTo2L2Q" "WZTo3LNu"
#             "TTWJetsToQQ" "TTWJetsToLNu"
#             "TTZToQQ" "TTZToLLNuNu_M_10" "TTZToLL_M_1to10"
            "MuonEG" "SingleElectron" "SingleMuon" "DoubleMuon" "DoubleEG" "MET")

runsamples_syst=( #"TTTo2L2Nu_TuneCP5up" "TTTo2L2Nu_TuneCP5down"
                 "TTTo2L2Nu_hdampUP" "TTTo2L2Nu_hdampDOWN")


uplimit=$((${#runsamples[@]}-1))
uplimit_syst=$((${#runsamples_syst[@]}-1))

init="Tree_"
final=".root"
sel="top"

if [ "$3" != "" ]; then
  sel=$3
fi
slash="/"
allok=0
path=""
checker=0
actualsize=0

# Minimum size (in bytes!):
minimumsize=50000

workingpath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
prepath="../.."
plotspath=$workingpath$prepath
if [ "$4" != "" ]; then
  plotspath=$4
fi


if [ "$1" == "an" ]; then
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TOP ANALYSIS EXECUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Setting up the environment"
    source pre_start.sh

    echo "%%%%%> DONE"
    echo ""
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Starting analysis"
    cd ../..
  
    echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Running general-purpose samples..."
    for ((i=0; i<=$uplimit; i++)); do
        root -l -b -q "RunAnalyserPAF.C(\"${runsamples[i]}\", \"$sel\", $2)"
        resetpaf -a
    done

#     echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Running samples for systematic uncertanties..."
#     for ((i=0; i<=$uplimit_syst; i++)); do
#         root -l -b -q "RunAnalyserPAF.C(\"${runsamples_syst[i]}\", \"$sel\", $2)"
#         resetpaf -a
#     done
  
  
elif [ "$1" == "ch" ]; then
  echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TOP CHECKER EXECUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
  echo ""
  echo ""
  echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Running prolegomena..."
  source pre_start.sh
  echo "%%%%%> DONE"
  echo " "
  
  cd ../..
  
  echo "Checking that the respective..."
  echo $uplimit
  echo "...root files of samples..."
  echo $uplimit_syst
  echo "...root files of systematic samples exist in..."
  echo $plotspath
  echo "...with selection..."
  echo $sel
  echo "...with at least..."
  echo $minimumsize
  echo "...bytes of size. If they do not fulfil these requirements, they will be reanalysed with..."
  echo $2
  echo "...cores."
  echo " "

  echo "%%%%%%%%> Checking general-purpose samples..."
  while [ $allok != ${#samples[@]} ]; do
    checker=$(($checker+1))
    allok=0
    for ((i=0; i<=$uplimit; i++)); do
      unset path
      unset actualsize
      
      path=$plotspath$slash$init${samples[i]}$final
      
      if [ ! -e $path ]; then
        echo " "
        echo "%%%% => ROOT file not found. The sample that is missing is:"
        echo ${samples[i]}
        echo "Reanalysing..."
        echo " "
        
        root -l -b -q "RunAnalyserPAF.C(\"${runsamples[i]}\", \"$sel\", $2)"
        resetpaf -a
        
        allok=$(($allok-8))
      fi
      
      if [ -e $path ]; then
        actualsize=$(wc -c <"$path")
        if [ $actualsize -le $minimumsize ]; then
          echo " "
          echo "%%%% => ROOT file with..."
          echo $actualsize
          echo "...bytes of size, which are lower than the minimum. This sample is:"
          echo ${samples[i]}
          echo "Reanalysing..."
          echo " "
          root -l -b -q "RunAnalyserPAF.C(\"${runsamples[i]}\", \"$sel\", $2)"
          resetpaf -a
          
          allok=$(($allok-8))
        fi
      fi
      
      allok=$(($allok+1))
    done
    
    if [ $checker == 10 ]; then
      echo " "
      echo "%%%% => ERROR: limit of iterations (10) reached. There has been a problem with the execution or the general-purpose sample files."
      echo "%%%% => The bash script will now end."
      echo " "
      cd plotter/top
      return
    fi
    sleep 5
  done

  path=""
  checker=0
  actualsize=0
  
#   echo "%%%%%%%%> Checking samples for systematic uncertanties..."
#   while [ $allok != ${#samples_syst[@]} ]; do
#     checker=$(($checker+1))
#     allok=0
#     for ((i=0; i<=$uplimit_syst; i++)); do
#       unset path
#       unset actualsize
#       
#       path=$plotspath$slash$init${samples_syst[i]}$final
#       
#       if [ ! -e $path ]; then
#         echo " "
#         echo "%%%% => ROOT file not found. The sample that is missing is:"
#         echo ${samples_syst[i]}
#         echo "Reanalysing..."
#         echo " "
#         root -l -b -q "RunAnalyserPAF.C(\"${runsamples_syst[i]}\", \"$sel\", $2)"
#         resetpaf -a
#         
#         allok=$(($allok-8))
#       fi
#       
#       if [ -e $path ]; then
#         actualsize=$(wc -c <"$path")
#         if [ $actualsize -le $minimumsize ]; then
#           echo " "
#           echo "%%%% => ROOT file with..."
#           echo $actualsize
#           echo "...bytes of size, which are lower than the minimum. This sample is:"
#           echo ${samples_syst[i]}
#           echo "Reanalysing..."
#           echo " "
#           root -l -b -q "RunAnalyserPAF.C(\"${runsamples_syst[i]}\", \"$sel\", $2)"
#           resetpaf -a
#           
#           allok=$(($allok-8))
#         fi
#       fi
#       
#       allok=$(($allok+1))
#     done
#     if [ $checker == 10 ]; then
#       echo " "
#       echo "%%%% => ERROR: limit of iterations (10) reached. There has been a problem with the execution or the sample files for systematic uncertanties."
#       echo "%%%% => The bash script will now end."
#       echo " "
#       cd plotter/top
#       return
#     fi
#     sleep 5
#   done
  
  echo "%%%% => All expected ROOT files are in the folder"
  
else
  echo "ERROR - No valid arguments given"
  echo "Please, execute this script with a valid argument"
fi

cd plotter/top
