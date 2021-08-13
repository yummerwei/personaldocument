grep 'number of k points' band.out |awk '{print $5}'  #得到总k点数

echo `grep -n 'number of k points' band.out | awk -F ':' '{print $1}'` + 2 |bc   #得到cart起始行数

echo `grep -n ' cryst. coord.' band.out | awk -F ':' '{print $1}'` + 1 |bc #得到分数坐标起始行数
