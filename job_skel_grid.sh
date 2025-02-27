#!/bin/bash
#SBATCH --job-name=skl2
#SBATCH --partition=pscomp
#SBATCH --mem=100gb
#SBATCH --time=2:00:00
#SBATCH --output=skl2_%j.out
#SBATCH --mail-user=fabien.castillo@etu.unistra.fr
#SBATCH --mail-type=ALL

module () {
  eval $(/usr/bin/modulecmd bash $*)
}

source /etc/profile.d/modules.sh
module purge
module load disperse/0.9.24

path_dat="/data100/fcastillo/RESULT/benchM/"
path="/data100/fcastillo/RESULT/benchM/"

file="4_densite"

persistence=0.1
smoothing=1


/softs/disperse/0.9.24/bin/mse $path_dat$file".fits" -cut $persistence -upSkl -manifolds -outName $path$file -periodicity 0 -forceLoops


/softs/disperse/0.9.24/bin/skelconv $path$file"_c"$persistence".up.NDskl" -outName $path$file"_c"$persistence".up.NDskl" -smooth $smoothing -to NDskl_ascii
#/softs/disperse/0.9.24/bin/skelconv $path$file"_c"$persistence".up.NDskl" -outName $path$file"_c"$persistence".up.NDskl" -smooth $smoothing -to NDskl_ascii

##/softs/disperse/0.9.24/bin/skelconv $path$file"_s"$persistence".up.NDskl" -outName $path$file"_s"$persistence".up.NDskl" -breakdown -smooth $smoothing


##/softs/disperse/0.9.24-CentOS6/bin/mse $path_dat$file".ND" -nthreads 40 -cut $persistence -upSkl -manifolds -outName $path$file -periodicity 0 -forceLoops
##/softs/disperse/0.9.24-CentOS6/bin/skelconv $path$file"_c"$persistence".up.NDskl" -outName $path$file"_c"$persistence".up.NDskl" -breakdown -smooth $smoothing
##/softs/disperse/0.9.24-CentOS6/bin/skelconv $path$file"_c"$persistence".up.NDskl" -outName $path$file"_c"$persistence".up.NDskl" -breakdown
##/softs/disperse/0.9.24-CentOS6/bin/skelconv $path$file"_c"$persistence".up.NDskl" -outName $path$file"_c"$persistence".up.NDskl" -smooth $smoothing
##/softs/disperse/0.9.24-CentOS6/bin/skelconv $path$file"_c"$persistence".up.NDskl" -outName $path$file"_c"$persistence".up.NDskl" -breakdown -smooth $smoothing
##/softs/disperse/0.9.24-CentOS6/bin/skelconv $path$file"_c"$persistence".up.NDskl" -outName $path$file"_c"$persistence".up.NDskl" -breakdown -smooth $smoothing -to vtp
##/softs/disperse/0.9.24-CentOS6/bin/mse $path_dat$file".ND" -loadMSC $path$file".MSC" -cut $persistence -outName $path$file -dumpManifolds J1a -forceLoops



exit 0
