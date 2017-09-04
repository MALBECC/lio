#/bin/bash
################################################################################
echo; echo "This test should take around 10 min ..."

echo; echo "Running BO0 ..."
./test_runbo0.sh

echo; echo "Running BO1 ..."
./test_runbo1.sh

echo; echo "Running BO2 ..."
./test_runbo2.sh

echo; echo "Running EH0 ..."
./test_runeh0.sh

echo; echo "Running EH1 ..."
./test_runeh1.sh

echo; echo "Running EH2 ..."
./test_runeh2.sh


filebo0=watersymbo_bo0.mdfrz
filebo1=watersymbo_bo1.mdfrz
filebo2=watersymbo_bo2.mdfrz
filebo3=watersymbo_bo3.mdfrz

fileeh0=watersymbo_eh0.mdfrz
fileeh1=watersymbo_eh1.mdfrz
fileeh2=watersymbo_eh2.mdfrz
fileeh3=watersymbo_eh3.mdfrz


tail -n +2 ${filebo2} > temp
cat ${filebo1} temp >  ${filebo3}
rm temp

tail -n +2 ${fileeh2} > temp
cat ${fileeh1} temp >  ${fileeh3}
rm temp


echo; echo "Comparing forces between eh rst and bo rst (should be similar):"
echo "vimdiff ${fileeh3} ${filebo3}"
sleep 1
vimdiff ${fileeh3} ${filebo3}

echo; echo "Comparing forces between eh rst and eh full (should be similar):"
echo "vimdiff ${fileeh3} ${fileeh0}"
sleep 1
vimdiff ${fileeh3} ${fileeh0}

################################################################################
