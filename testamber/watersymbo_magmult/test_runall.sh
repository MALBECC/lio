#/bin/bash
################################################################################
echo; echo "This test should take around 3 min ..."

echo; echo "Running BO ..."
./test_runbo.sh

echo; echo "Running EH ..."
./test_runeh.sh


file1=watersymbo_bo.mdfrz
file2=watersymbo_eh.mdfrz
echo; echo "Comparing forces (should be similar):"
echo "vimdiff ${file1} ${file2}"
sleep 1
#vimdiff ${file1} ${file2}

################################################################################
