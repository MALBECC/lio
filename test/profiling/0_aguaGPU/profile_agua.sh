#! /bin/bash
if [ -z "$LIOBIN" ] ; then
  LIOBIN=../../../liosolo/liosolo
fi
SALIDA=salida
if [ -n "$1" ]
  then
    SALIDA=$1
fi

ENERGIA=energia

# Timestamp for the test folder name
now=$(date +"%Y%m%d%T")

echo "========================="
echo "= PROFILING: AGUA_GPU   ="
echo "========================="

### Antes de ejecutar los test, borramos los archivos de salida
#rm salida*
#rm energias*

### El usuario ingresa la cantidad de test a ejecutar
##echo -n "Cantidad de iteraciones:"
##read iteraciones

## Chequeamos si puso commo input la cantidad de iteraciones
iteraciones=0
if [ -n "$1" ]
    then
	iteraciones=$1
fi

### El usuario ingresa la cantidad de test a ejecutar
if [ $iteraciones -eq 0 ]
    then
	echo -n "Cantidad de iteraciones:"
	read iteraciones
fi

## Renombramos el achivo de salida
SALIDA=${iteraciones}_salida_${now}
ENERGIA=${iteraciones}_energia_${now}

echo "Profiling: AGUA_GPU" >> $SALIDA

### Ciclo para ejecutar varias veces el mismo test
counter=1
while [ $counter -le $iteraciones ]
do
    #echo "================="
    #echo "= Test $counter/$iteraciones        ="
    #echo $LIOBIN -i agua.in -b basis -c agua.xyz -v
    echo "=================" >> $SALIDA
    echo "= Profile $counter        =" >> $SALIDA
    echo "=================" >> $SALIDA
    nvprof --analysis-metrics -o agua.nvprof $LIOBIN -i agua.in -b basis -c agua.xyz -v >> $SALIDA
    ((counter++))
done

echo "Fin"
