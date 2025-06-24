#!/bin/bash

grep "Current TB electrode           1" currentTB.ok | awk '{print$5}'  > I1ok.dat
grep "Current TB electrode           2" currentTB.ok | awk '{print$5}'  > I2ok.dat
grep "Current DFT part M" currentTB | awk '{print$5}'  > I0ok.dat

grep "Mulliken TB electrode           1" mullikenTB.ok | awk '{print$5}'  > MG1ok.dat
grep "Mulliken TB electrode           2" mullikenTB.ok | awk '{print$5}'  > MG2ok.dat
grep "Mulliken DFT  part M" mullikenTB | awk '{print$5}'  > MG0ok.dat

grep "Current TB electrode           1" currentTB | awk '{print$5}'  > I1.dat
grep "Current TB electrode           2" currentTB | awk '{print$5}'  > I2.dat
grep "Current DFT part M" currentTB | awk '{print$5}'  > I0ok.dat

grep "Mulliken TB electrode           1" mullikenTB | awk '{print$5}'  > MG1.dat
grep "Mulliken TB electrode           2" mullikenTB | awk '{print$5}'  > MG2.dat
grep "Mulliken DFT  part M" mullikenTB | awk '{print$5}'  > MG0.dat

xmgrace I*
xmgrace MG*

rm I*
rm MG*

