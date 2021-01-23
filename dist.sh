#!/bin/bash
	for m in {1..11}
		do
		dist=`echo "scale=2;${m}*0.1+2.4" |bc`
		var1=`echo "scale=8;2*${dist}+2.68" |bc`
		var2=`echo "scale=8;1.34/(2.68+${dist})" |bc`
		var3=`echo "scale=8;${var2}+0.5"|bc`
		mkdir dist-${dist};cd dist-${dist}
                cp ../{POTCAR,POSCAR,INCAR,KPOINTS,vasp.pbs,vdw*} .
		sed -i 's/var1/'$var1'/g' POSCAR
		sed -i 's/var2/'$var2'/g' POSCAR
		sed -i 's/var3/'$var3'/g' POSCAR
		echo 'dist='${dist} 'c='${var1} 'As1z'=${var2} 'As2z'=${var3}
		sed -i 's/abcd/'${dist}'/g' vasp.pbs
		qsub vasp.pbs
		cd ..
		done 
