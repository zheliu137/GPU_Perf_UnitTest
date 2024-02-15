rm -rf time.dat
for i in `seq 20 40`
do
  ./main.x $i > out
  time=`grep "avg time per matmul" out `
  echo $i $time >> time.dat
done
