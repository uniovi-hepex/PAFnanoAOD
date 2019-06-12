#===============================================================================
#
#                         Analysis of the TOP process
#
#===============================================================================

ext="logs/"
mkdir -p $ext
workingpath="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
logpath=$workingpath$ext

echo ""
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TOP ANALYSIS EXECUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo ""
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Creating jobs..."
An=$(qsub -q batch -l nodes=1:ppn=$1 -l walltime=20:00:00 -N "PAF_TOP_EXEC" -o $logpath -e $logpath -d $workingpath -F "an $1 $2 $3" Executioner.sh)
echo $An
qsub -q batch -l nodes=1:ppn=$1 -l walltime=20:00:00 -N "PAF_TOP_CHECK" -o $logpath -e $logpath -d $workingpath -W depend=afterany:$An -F "ch $1 $2 $3" Executioner.sh
echo "Done!!"