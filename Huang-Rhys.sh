#!/bin/bash
#in EC file
EPMR=`grep -n ' POSITION' OUTCAR |tail -n 1 |awk -F ':' '{print $1}'`
EPMR=`expr $EPMR + 2`
NIONS=`grep -n ' NIONS' OUTCAR |awk -F '=' '{print $3}'`
EPNR=`expr $EPMR + $NIONS`

GPMR=`grep -n ' position of ions in cartesian coordinates  (Angst):' OUTCAR|awk -F ':' '{print $1}'`
GPMR=`expr $GPMR + 1`
GPNR=`expr $GPMR + $NIONS`

awk "{if(NR>=$EPMR && NR<$EPNR) print \$1,"\t",\$2,"\t",\$3}" OUTCAR > PC.dat
awk "{if(NR>=$GPMR && NR<$GPNR) print \$1,"\t",\$2,"\t",\$3}" OUTCAR > PA.dat
paste PA.dat PC.dat > PAC.dat; rm -rf PC.dat PA.dat
for i in $(seq 1 $N1)
do 
echo "28.084" >> M.dat
done
for i in $(seq 1 $N2)
do 
echo "12.0096" >> M.dat
done
paste M.dat PAC.dat > MPAC.dat

awk '{printf("%12.8f\t%12.8f\t%12.8f\n",$4-$1,$5-$2,$6-$3)}' PAC.dat >  R.dat

#---区域范围的确定----#
PHMR=`grep -n 'dynamical' OUTCAR-phonon |awk -F ':' '{print $1}'`
PHMR=`expr $PHMR + 3`
PHNR=`grep -n 'Finite differences POTIM='  OUTCAR-phonon |awk -F ':' '{print $1}'`
PHNR=`expr $PHNR - 1`


PERIOD=`expr $NIONS + 3` #每个模式有多少行
DEG=`grep -n 'Degree of freedom'  OUTCAR-phonon |awk -F '/' '{print $2}' |head -n 1` #多少种模式
PPHMR=`expr $PHMR + 1` #决定F=所在行数
MOD=$PPHMR%$PERIOD  #决定频率行所在位置
awk "{if(NR>=$PHMR&&NR<=$PHNR&&NR%$PERIOD==$MOD) print \$10}" OUTCAR-phonon > fmode.dat



for i in $(seq 1 $DEG)
do
    for j in $(seq 1 $NIONS)
do
printf "%8i\t%8i\n" $i $j >> 1-NB.dat
awk "NR==$i" fmode.dat >> 2-ffreq.dat
L=`expr $PPHMR + 1 + $j + $PERIOD \* \( $i - 1 \)`
awk -v "L=$L" '{if(NR==L) printf("%12.8f\t%12.8f\t%12.8f",$4,$5,$6)}' OUTCAR-phonon >> 5-umode.dat
done
cat MPAC.dat >> 4-MPACC.dat
done
