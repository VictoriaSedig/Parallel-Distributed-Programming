f="result_10_cores"
touch $f
rm $f
touch $f
for i in {1..10}
do
        echo "computing "
        echo $i
        mpirun --bind-to none -np $i stencil 3628800 100 >> $f
done

