f="result_32_cores"
touch $f
rm $f
touch $f
for i in 2 4 8 16 32 64
do
        echo "computing "
        echo $i
        mpirun -np $i stencil 1048576 100 >> $f
done

