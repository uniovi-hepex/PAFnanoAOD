# This script makes a copy of the full content of ../../PROCESS_temp into ~/Storage/PROCESS/MiniTrees/YYYY_MM_DD
#
#

lowerbar="_"
slash="/"
plots="Plots"
minitrees="MiniTrees"

indstr="index.php"
indexpath="/nfs/fanae/user/vrbouza/www/index.php"

sourcepath=""
storagepath="/pool/cienciasrw/userstorage/vrbouza/proyectos/TW_inclusivo_run2/"
webpath="/nfs/fanae/user/vrbouza/www/Proyectos/tw_inclusivo_run2/results/"
savepath="../../temp_TW/"

if [ "$1" == "w" ] || [ "$1" == "web" ]; then
  echo "===> Copying results to the web!"
  echo " "
  echo "Plots will be copied into..."
  echo $webpath
  
  echo " "
  echo "Copying results..."
  rsync -arvzP results/* $webpath
  echo "Copying index.php to all subfolders..."
  find $webpath -type d -exec cp $indexpath {} \;
  echo " "
  echo "Done!"
  return
elif [ "$1" == "r" ] || [ "$1" == "results" ]; then
  echo "===> Copying results to the userstorage!"
  echo " "

  if [ "$2" != "" ]; then
    savefolder=$2
  else
    d=$(date +%d)
    m=$(date +%m)
    y=$(date +%Y)
    savefolder=$y$lowerbar$m$lowerbar$d
  fi

  savepath=$storagepath$Plots$slash$savefolder
  echo "Results will be copied into..."
  echo $savepath

  echo " "
  echo "Creating folder (if it does not exist)..."
  mkdir -p $savepath

  echo " "
  echo "Copying results..."
  savepath=$savepath$slash
  rsync -arvzP ../../temp_TW/ $savepath

  echo " "
  echo "Done!"
  return
else
  d=$(date +%d)
  m=$(date +%m)
  y=$(date +%Y)
  savefolder=$y$lowerbar$m$lowerbar$d
  savepath=$storagepath$savefolder
fi


echo "===> Copying minitrees to the userstorage!"
echo " "
echo "Files will be copied into..."
echo $savepath

echo " "
echo "Creating folder (if it does not exist)..."
mkdir -p $savepath

echo " "
echo "Copying files (note that this process can last even minutes depending on the size of the minitrees)..."
savepath=$savepath$slash
rsync -arvzP $sourcepath $savepath

echo " "
echo "Done!"

