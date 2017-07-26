#/bin/bash
################################################################################
echo "vimdiff waterlong.mdcrd watershrt.mdcrd"
sleep 1
vimdiff waterlong.mdcrd watershrt.mdcrd

echo "vimdiff waterlong.mdvel watershrt.mdvel"
sleep 1
vimdiff waterlong.mdvel watershrt.mdvel

echo "vimdiff waterlong.mdfrz watershrt.mdfrz"
sleep 1
vimdiff waterlong.mdfrz watershrt.mdfrz

################################################################################
