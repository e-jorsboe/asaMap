#!/bin/bash

PRG=""


if [ $# -eq 0 ] 
then
    exit 1;
fi

if [ $# -eq 1 ]
then
    PRG=$1
fi


echo "--------------------"
echo "Using PRG: '${PRG}' "
echo "--------------------"

WDIR=`dirname $PRG`

RVAL=0

$WDIR/asaMap -p plink/plink -o plink1 -y plink/plink.phe -Q plink/plink.2.Q -f plink/plink.2.P -c plink/plink.cov -m 0 -r 778
res1new=`cat plink1.res | tail -n +2 | cut -f6`
$WDIR/asaMap -p plink/plink -o plink2 -y plink/plink.phe -Q plink/plink.2.Q -f plink/plink.2.P -c plink/plink.cov -m 1 -r 778
res2new=`cat plink2.res | tail -n +2 | cut -f6`
$WDIR/asaMap -p plink/plink -o plink3 -y plink/plink.phe -Q plink/plink.2.Q -f plink/plink.2.P -c plink/plink.cov -w 0 -r 778
res3new=`cat plink3.res | tail -n +2 | cut -f6`

res1=234.591730

if [ ! "$res1" = "$res1new"  ] ;then
    echo "--------------"
    echo "Problem with test1"
    echo "--------------"
    RVAL=1
fi

res2=235.326567

if [ ! "$res2" = "$res2new"  ] ;then
    echo "--------------"
    echo "Problem with test2"
    echo "--------------"
    RVAL=2
fi

res3=234.607741

if [ ! "$res3" = "$res3new"  ] ;then
    echo "--------------"
    echo "Problem with test3"
    echo "--------------"
    RVAL=3
fi

exit $RVAL

