#!/bin/bash
################################################################################
# COMPILE RESULTS
file_list="cmdview_results.sh"
temp_list=""

syslist=""
syslist="${syslist} watersymbo_verlong"
syslist="${syslist} watersymbo_vermult"
syslist="${syslist} watersymbo_magnus"
syslist="${syslist} watersymbo_maglet"
syslist="${syslist} watersymbo_magmult"
syslist="${syslist} watersymbo_verlet"

for system in ${syslist}; do

src_bo=${system}/watersymbo_bo.out
src_eh=${system}/watersymbo_eh.out
result="${system}.tkv"
temp_list="${temp_list} ${result}"


grep 'Etot' ${src_bo} > tempbo0
grep 'Etot' ${src_eh} > tempeh0
head -n -2 tempbo0 > tempbo1
head -n -2 tempeh0 > tempeh1
paste tempbo1 tempeh1 > ${result}
rm -f tempbo0 tempbo1 tempeh0 tempeh1

done

file_list="${file_list} ${temp_list}"


################################################################################
# COPY RESULTS AND CLEAN
targetbox="francisco@logos"
targetbox=${1}
targetdir="/home/francisco/Ehrenfest/testamber/."

echo; echo Exporting to  :  ${targetbox}:${targetdir}
echo; echo List of files :
echo; echo ${file_list} | tr " " "\n"
echo; rsync -rP ${file_list} ${targetbox}:${targetdir}

rm -f ${temp_list}

################################################################################
