#!/bin/bash
# Este script extrae la matriz densidad excitada de GAUSSIAN
# como un vector.

ARCHIVO=$1

formchk $ARCHIVO.chk

awk '$2=="CI",$1=="QEq"' $ARCHIVO.fchk | awk '$3=="Density",$1=="QEq"' > scratch1
sed '1d' scratch1 | sed '$d' > scratch2
awk '{for (i=1;i<=NF;i++) print $i}' scratch2 > gau_dens
rm scratch* $ARCHIVO.fchk
