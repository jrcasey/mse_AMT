#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -p glenn
#SBATCH -A snic001-12-249
#SBATCH -o my.stdout
#SBATCH -t 167:00:00
#SBATCH --mail-user gatto@chalmers.se


module load matlab/7.14

MOSEKPLATFORM=linux64x86
export PATH=$PATH:/c3se/users/gatto/Glenn/mosek/6/tools/platform/linux64x86/bin
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/c3se/users/gatto/Glenn/mosek/6/tools/platform/linux64x86/bin
export MOSEKLM_LICENSE_FILE=/c3se/users/gatto/Glenn/mosek/6/licenses/mosek.lic

matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(1,2)" > INITOutputLog12 2> INITErrorLog12 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(2,2)" > INITOutputLog22 2> INITErrorLog22 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(3,2)" > INITOutputLog32 2> INITErrorLog32 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(4,2)" > INITOutputLog42 2> INITErrorLog42 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(5,2)" > INITOutputLog52 2> INITErrorLog52 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(6,2)" > INITOutputLog62 2> INITErrorLog62 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(7,2)" > INITOutputLog72 2> INITErrorLog72 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(8,2)" > INITOutputLog82 2> INITErrorLog82 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(9,2)" > INITOutputLog92 2> INITErrorLog92 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(10,2)" > INITOutputLog102 2> INITErrorLog102 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(11,2)" > INITOutputLog112 2> INITErrorLog112 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(12,2)" > INITOutputLog122 2> INITErrorLog122 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(13,2)" > INITOutputLog132 2> INITErrorLog132 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(14,2)" > INITOutputLog142 2> INITErrorLog142 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(15,2)" > INITOutputLog152 2> INITErrorLog152 &
matlab -nojvm -nodesktop -r "getINITmodelFromArrayDataInServer(16,2)" > INITOutputLog162 2> INITErrorLog162 &
wait

MOSEKPLATFORM=linux64x86
export PATH=$/c3se/users/gatto/Glenn//mosek/7/tools/platform/linux64x86/bin/mosek:$PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/c3se/users/gatto/Glenn/mosek/7/tools/platform/linux64x86/bin
export MOSEKLM_LICENSE_FILE=/c3se/users/gatto/Glenn/mosek/7/licenses/mosek.lic

matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(1,2)" > TINITOutputLog12 2> TINITErrorLog12 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(2,2)" > TINITOutputLog22 2> TINITErrorLog22 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(3,2)" > TINITOutputLog32 2> TINITErrorLog32 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(4,2)" > TINITOutputLog42 2> TINITErrorLog42 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(5,2)" > TINITOutputLog52 2> TINITErrorLog52 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(6,2)" > TINITOutputLog62 2> TINITErrorLog62 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(7,2)" > TINITOutputLog72 2> TINITErrorLog72 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(8,2)" > TINITOutputLog82 2> TINITErrorLog82 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(9,2)" > TINITOutputLog92 2> TINITErrorLog92 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(10,2)" > TINITOutputLog102 2> TINITErrorLog102 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(11,2)" > TINITOutputLog112 2> TINITErrorLog112 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(12,2)" > TINITOutputLog122 2> TINITErrorLog122 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(13,2)" > TINITOutputLog132 2> TINITErrorLog132 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(14,2)" > TINITOutputLog142 2> TINITErrorLog142 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(15,2)" > TINITOutputLog152 2> TINITErrorLog152 &
matlab -nojvm -nodesktop -r "getTINITmodelFromArrayDataInServer(16,2)" > TINITOutputLog162 2> TINITErrorLog162 &
wait