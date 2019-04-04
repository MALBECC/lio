#/bin/bash
################################################################################
echo
echo " Intro general message..."
sleep 3

echo; echo " Comparing velocities of scratch VS restart with 0 vels"
echo " vimdiff water_scratch.mdvel water_restart.mdvel"
sleep 2
vimdiff water_scratch.mdvel water_restart.mdvel

echo; echo " Comparing forces of scratch VS restart with 0 vels"
echo " vimdiff water_scratch.mdfrz water_restart.mdfrz"
sleep 2
vimdiff water_scratch.mdfrz water_restart.mdfrz

################################################################################
