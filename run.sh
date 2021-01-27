#!/bin/bash
filelist=`ls /public/home/weiyd/postdoc/amor`
for file in $filelist
do
 mkdir ./POSCARS/${file}v
 mv ./POSCARS/${file} ./POSCARS/${file}v/POSCAR
 cp {INCAR,POTCAR,KPOINTS,run} ./POSCARS/${file}v
done
#for nn in 1..206
#filelist2=`find ./POSCARS -type d -name "9*"`
#filelist2=`cat log`
#for file in $filelist2
#do
# cd $file
# yhbatch -N 1 -n 24 vasp-weiyd
# cp vasprun.xml ../${file}_vasprun.xml
# cd ../
#done
