# This script makes a copy of the full content of ../../PROCESS_temp into ~/Storage/PROCESS/MiniTrees/YYYY_MM_DD
#
#

# TODO: copiar tolos plots
plotspath="/nfs/fanae/user/vrbouza/www/Proyectos/twttbar_run2/results/"
mcdatastr="MCData/"
indstr="index.php"
cnstr="CondN/"
covstr="CovMat/"
indexpath="/nfs/fanae/user/vrbouza/www/index.php"

if [ "$1" == "p" ]; then
  echo "===> Copying unfolding results!"
  echo " "
  echo "Plots will be copied into..."
  echo $plotspath
  
  echo " "
  echo "Copying plots..."
  cp -R results/. $plotspath
  #rsync -avzP results/. $plotspath
  cp $indexpath $plotspath$indstr
  cp $indexpath $plotspath$mcdatastr$indstr
  #cp $indexpath $plotspath$covstr$indstr
  #cp $indexpath $plotspath$cnstr$indstr
  echo " "
  echo "Done!"
  return
elif [ "$1" != "" ]; then
  d=$1
else
  d=$(date +%d)
fi

if [ "$2" != "" ]; then
  m=$2
else
  m=$(date +%m)
fi

if [ "$3" != "" ]; then
  y=$3
else
  y=$(date +%Y)
fi

lowerbar="_"
slash="/"
savefolder=$y$lowerbar$m$lowerbar$d
storagepath="/pool/cienciasrw/userstorage/vrbouza/proyectos/TWTTbar/MiniTrees/"
savepath=$storagepath$savefolder

echo "===> Copying minitrees!"
echo " "
echo "Files will be copied into..."
echo $savepath
mkdir -p $savepath

echo " "
echo "Copying files (note that this process can last even minutes depending on the size of the minitrees)..."
savepath=$savepath$slash
#cp -R ../../temp_TWTTbar/. $savepath
rsync -avzP ../../temp_TWTTbar/ $savepath

echo " "
echo "Done!"

