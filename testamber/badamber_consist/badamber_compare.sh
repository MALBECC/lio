#/bin/bash
################################################################################
echo
echo " Comparing coordinates, velocities and forces."
echo " All should be equal but they are not..."
sleep 3

echo; echo " Comparing coordinates - these tend to be equal"
echo " vimdiff waterlong.mdcrd watershrt.mdcrd"
sleep 2
vimdiff waterlong.mdcrd watershrt.mdcrd

echo; echo " Comparing velocities - something wrong happens at the last step"
echo " vimdiff waterlong.mdvel watershrt.mdvel"
sleep 2
vimdiff waterlong.mdvel watershrt.mdvel

echo; echo " Comparing forces - these are out of phase"
echo " vimdiff waterlong.mdfrz watershrt.mdfrz"
sleep 2
vimdiff waterlong.mdfrz watershrt.mdfrz

################################################################################
