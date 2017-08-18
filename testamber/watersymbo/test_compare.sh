#/bin/bash
################################################################################
file1=watersymbo_bo.mdfrz
file2=watersymbo_eh.mdfrz

echo; echo "Comparing EH and BO forces (should be similar):"
echo "vimdiff ${file1} ${file2}"
sleep 1
vimdiff ${file1} ${file2}


file1=watersymbo_bo.mdfrz
file2=watersymbo_bo.gdfrz

echo; echo "Comparing BO forces and reference (should be similar):"
echo "vimdiff ${file1} ${file2}"
sleep 1
vimdiff ${file1} ${file2}


file1=watersymbo_eh.mdfrz
file2=watersymbo_eh.gdfrz

echo; echo "Comparing EH forces and reference (should be similar):"
echo "vimdiff ${file1} ${file2}"
sleep 1
vimdiff ${file1} ${file2}

################################################################################
