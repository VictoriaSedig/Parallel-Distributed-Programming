f="result_16_cores"
fref="result_16_cores_ref"
touch $f
rm $f
touch $f
touch $fref
rm $fref
touch $fref

i=$(( 128 ))
while [ $i -le 10000000 ]
do
        echo "computing "
        echo $i
        mpirun --bind-to none -np 16 stencil $i 100 >> $f
        ./stencil_ref $i 100 >> $fref
        i=$(( $i *  2 ))
done

