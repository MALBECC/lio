#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../../liosolo/liosolo
fi

## Chequeamos si puso commo input la cantidad de iteraciones
iteraciones=0
if [ -n "$1" ]
    then
	iteraciones=$1
fi

ENERGIA=energia

# Timestamp for the test folder name
now=$(date +"%Y%m%d%T")

echo "========================="
echo "= STRESS TEST: caroteno_GPU ="
echo "========================="

### El usuario ingresa la cantidad de test a ejecutar
if [ $iteraciones -eq 0 ]
    then
	echo -n "Cantidad de iteraciones:"
	read iteraciones
fi

## Renombramos el achivo de salida
SALIDA=${iteraciones}_salida_${now}
ENERGIA=${iteraciones}_energia_${now}

echo "STRESS TEST: caroteno_GPU" >> $SALIDA

### Ciclo para ejecutar varias veces el mismo test
counter=1
while [ $counter -le $iteraciones ]
do
#    echo "================="
#    echo "= Test $counter/$iteraciones        ="
#    echo $LIOBIN -i caroteno.in -b basis -c caroteno.xyz -v
    echo "=================" >> $SALIDA
    echo "= Test $counter        =" >> $SALIDA
    echo "=================" >> $SALIDA
    $LIOBIN -i carotenox.in -b DZVP -c caroteno.xyz -v >> $SALIDA
    ((counter++))
done

## Ahora filtramos los resultados

echo "========================"
echo "= Filtrando resultados ="
echo "= Energias             ="
echo "grep -E 'COULOMB|ONE ELECTRON|NUCLEAR|EXC.|TOTAL' $SALIDA"
echo "Energias " >> $ENERGIA
grep -E "COULOMB|ONE ELECTRON|NUCLEAR|EXC.|TOTAL" $SALIDA >> $ENERGIA
echo "Formateando la salida"
#echo sed -i 's/\║/,/g' $ENERGIA
sed -i 's/\║/,/g' $ENERGIA

## Separamos los archivos de energias
#COULOMB
echo "Separando COULOMB"
echo "COULOMB" >> ${iteraciones}_coulomb_${now}
grep -E "COULOMB" $ENERGIA | sed 's/    ,   COULOMB        ,/ /g' | sed 's/       ,/ /g' >> ${iteraciones}_coulomb_${now} | sed 's/\./\,/g'
sed -i 's/\./\,/g' ${iteraciones}_coulomb_${now}

#ONE ELECTRON
echo "Separando ONE ELECTRON"
echo "ONE ELECTRON" >> ${iteraciones}_one_electron_${now}
grep -E "ONE ELECTRON" $ENERGIA | sed 's/    ,   ONE ELECTRON   ,/ /g' | sed 's/       ,/ /g' >> ${iteraciones}_one_electron_${now} | sed 's/\./\,/g'
sed -i 's/\./\,/g' ${iteraciones}_one_electron_${now}

#NUCLEAR
echo "Separando NUCLEAR"
echo "NUCLEAR" >> ${iteraciones}_nuclear_${now}
grep -E "NUCLEAR" $ENERGIA | sed 's/    ,   NUCLEAR        ,/ /g' | sed 's/       ,/ /g' >> ${iteraciones}_nuclear_${now} | sed 's/\./\,/g'
sed -i 's/\./\,/g' ${iteraciones}_nuclear_${now}

#EXC. - CORR.
echo "Separando EXC. - CORR."
echo "EXC. - CORR." >> ${iteraciones}_exc_corr_${now}
grep -E "EXC. - CORR." $ENERGIA | sed 's/    ,   EXC. - CORR.   ,/ /g' | sed 's/       ,/ /g' >> ${iteraciones}_exc_corr_${now} | sed 's/\./\,/g'
sed -i 's/\./\,/g' ${iteraciones}_exc_corr_${now}

#TOTAL
echo "Separando TOTAL"
echo "TOTAL" >> ${iteraciones}_total_${now}
grep -E "TOTAL" $ENERGIA | sed 's/    ,   TOTAL          ,/ /g' | sed 's/       ,/ /g' >> ${iteraciones}_total_${now}
sed  -i 's/\./\,/g' ${iteraciones}_total_${now}

#Tiempos del timer (solo para ver la parte en que se usa libxc
echo "Separanto Exchange-correlation Fock de la parte de timers"
echo "Exchange-correlation Fock" >> ${iteraciones}_exc_corr_time_${now}
grep -E "Exchange-correlation Fock" $SALIDA | sed 's/Exchange-correlation Fock/ /g' >> ${iteraciones}_exc_corr_time_${now}
## Le quitamos los porcentajes
sed -i 's/(.*)/ /g' ${iteraciones}_exc_corr_time_${now}
## Por ultimo el character s al final lo quitamos
sed -i 's/s/ /g' ${iteraciones}_exc_corr_time_${now}

#Lio deja los datos con . y los necesito con , para las hojas de calculo
sed -i 's/\./\,/g' ${iteraciones}_exc_corr_time_${now}

echo "Fin"
