#/bin/bash
################################################################################
for numb in 0 1 2
do

   rm -f watersymbo_bo${numb}.stdo  watersymbo_bo${numb}.out
   rm -f watersymbo_bo${numb}.rst7  watersymbo_bo${numb}.mdcrd
   rm -f watersymbo_bo${numb}.mdfrz watersymbo_bo${numb}.mdvel

   rm -f watersymbo_eh${numb}.stdo  watersymbo_eh${numb}.out
   rm -f watersymbo_eh${numb}.rst7  watersymbo_eh${numb}.mdcrd
   rm -f watersymbo_eh${numb}.mdfrz watersymbo_eh${numb}.mdvel

done

rm -f watersymbo_bo3.mdfrz watersymbo_eh3.mdfrz

rm -f qm.xyz mdinfo Check_posvel.log dipole_moment dipole_moment.dat
################################################################################
