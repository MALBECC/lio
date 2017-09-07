#!/bin/bash
################################################################################
rungrace='xmgrace -geometry 1200x900+200'

srcfile1="watersymbo_verlong.tkv"
srcfile2="watersymbo_vermult.tkv"
srcfile3="watersymbo_maglet.tkv"
srcfile4="watersymbo_magnus.tkv"
srcfile5="watersymbo_magmult.tkv"

echo
echo "Standard colors:"
echo " * Black:   born opp"
echo " * Red:     verlet"
echo " * Green:   vermult"
echo " * Blue:    maglet"
echo " * Yellow:  magnus"
echo " * Brown:   magmult"

echo
echo -n "Comparing total energies"
echo -n "."; sleep .5; echo -n "."; sleep .5; echo "."; sleep .5; echo
${rungrace} -block ${srcfile1} -bxy 0:3  -block ${srcfile1} -bxy 0:12 \
            -block ${srcfile2} -bxy 0:12 -block ${srcfile3} -bxy 0:12 \
            -block ${srcfile4} -bxy 0:12 -block ${srcfile5} -bxy 0:12

echo
echo -n "Comparing potential energies"
echo -n "."; sleep .5; echo -n "."; sleep .5; echo "."; sleep .5; echo
${rungrace} -block ${srcfile1} -bxy 0:9  -block ${srcfile1} -bxy 0:18 \
            -block ${srcfile2} -bxy 0:18 -block ${srcfile3} -bxy 0:18 \
            -block ${srcfile4} -bxy 0:18 -block ${srcfile5} -bxy 0:18


################################################################################
