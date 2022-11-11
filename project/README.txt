The code is built by typing make.
The code has one argument inputfile. The first number in the inputfile is an int specifying the size n of a n*n matrix. n has to be evenly didvided by the number of processors. After that there are the values to the matrix row-wise seperated by a space.  

The program is run as following:

```
mpirun -np p shear [inputfile.txt]
```

where p is the number of processors.

The time measurement will be printed. 
The results will be written to output.txt will be created if it does not already exist.