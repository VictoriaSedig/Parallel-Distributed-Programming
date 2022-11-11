#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ROOT 0


void matrixMultiply(double *A, double *B, int chunk, double *C) {

	///// Matrix multiplication for the blocks
		for (int k = 0; k < chunk; k++){
			for (int i = 0; i < chunk; i++){
				for (int j = 0; j < chunk; j++){
					C[i*chunk+j] += A[i*chunk + k] * B[k*chunk +j];
				}	 
			}
		}
	}


int main(int argc, char **argv)
{
// Initilize values
	double *my_block_A;
	double *my_block_B;
	double *blocks_array_A;
	double *blocks_array_B;
	double *C;
	double* big_C;
	double d;
  	double time, start_time, end_time, time_big = 0;
	int n;
	int k, l;
	int rank, size;
	int store;
	int tagA, tagB, tagAA, tagBB, destA, destB, sourceA, sourceB, source, dest, sourceBB, sourceAA, destAA, destBB;
	char* input = argv[1];
	char* output = argv[2];

// Start upp the MPI procs definitions

	MPI_Status status;
	MPI_Request request;
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int size_root = (int) sqrt(size);

	FILE *fp;
	// Send out n matrix size to all processors
	if(rank==0){
		fp=fopen(input,"r");
			store = fscanf(fp,"%d",&n); 
	  	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int chunk = n/size_root;
  	int chch = chunk*chunk;

	// Allocate memories
	my_block_A = (double*) malloc(chch * sizeof(double));
	my_block_B = (double*) malloc(chch * sizeof(double));
	big_C = (double*) malloc(n*n * sizeof(double));
	C = (double*) malloc(chch * sizeof(double));
	for(int i; i<chch; i++){
		C[i] = 0;}

	// Fill an array with the matrix elements from A and B in good order
	if(rank==0){
		assert(n % size_root == 0);
		blocks_array_A = (double*) malloc(size *chch* sizeof(double*));
		blocks_array_B = (double*) malloc(size *chch* sizeof(double*));

		int s=0;
	    //Fill A
		for(int k=0; k< size_root; k++){
			for(int l=0;l<chunk;l++){ 
		    	for(int i = 0; i < size_root; i++){
		    		s = i + k*size_root;
		    		for(int j=0; j<chunk; j++){
				    	store = fscanf(fp,"%lf",&d); 
		    			blocks_array_A[s*chch + l * chunk + j] = d;
	    		}
	    	}
		}
	}
	    //Fill B
		for(int k=0; k< size_root; k++){
			for(int l=0;l<chunk;l++){ 
		    	for(int i = 0; i < size_root; i++){
		    		s = i + k*size_root;
		    		for(int j=0; j<chunk; j++){
				    	store = fscanf(fp,"%lf",&d); 
		    			blocks_array_B[s*chch + l * chunk +j] = d;
		    	}
			}
	  	}
	}
  		fclose(fp);
}

	// Scatter the values for each processors block
  	MPI_Scatter(blocks_array_A, chch, MPI_DOUBLE, my_block_A , chch, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  	MPI_Scatter(blocks_array_B, chch, MPI_DOUBLE, my_block_B , chch, MPI_DOUBLE, 0, MPI_COMM_WORLD);


/////////////////////// Start the Cannon algo ///////////////////

		// Initialize rows and columns of processors
		int my_color_row = rank/size_root;  
		int my_color_kol = rank % size_root;

		MPI_Comm row_comm;
		MPI_Comm kol_comm;

		MPI_Comm_split(MPI_COMM_WORLD, my_color_row, rank, &row_comm);
		MPI_Comm_split(MPI_COMM_WORLD, my_color_kol, rank, &kol_comm);

		int row_rank, row_size;
		MPI_Comm_rank(row_comm, &row_rank);
		MPI_Comm_size(row_comm, &row_size);

		int kol_rank, kol_size;
		MPI_Comm_rank(kol_comm, &kol_rank);
		MPI_Comm_size(kol_comm, &kol_size);


  start_time = MPI_Wtime(); // start time measurement

// Initial shifts where each row/kolumn flips different amount of times (Cannon algo)
		for(int i = 1; i<size_root; i++){
			if(kol_rank == i){	
					destA = (row_rank + 2*row_size -i) % row_size;
					sourceA = (row_rank + i) % row_size;
					MPI_Sendrecv_replace(my_block_A, chch, MPI_DOUBLE, destA, 5, sourceA, 5, row_comm, &status);
				}
			if(row_rank == i){
					destB = (kol_rank + 2*kol_size -i) % kol_size;
					sourceB = (kol_rank + i) % kol_size;
					MPI_Sendrecv_replace(my_block_B, chch, MPI_DOUBLE, destB, 6, sourceB, 6, kol_comm, &status);
				}
			}
//// Time to do the multiplications 

		matrixMultiply(my_block_A, my_block_B, chunk, C);
	    ////Now do the whole movement. 
		for(int i=0; i<size_root-1; i++){
				destAA = (row_rank + row_size -1) % row_size;
				sourceAA = (row_rank + 1) % row_size;
				MPI_Sendrecv_replace(my_block_A, chch, MPI_DOUBLE, destAA, 7, sourceAA, 7, row_comm, &status);

				destBB = (kol_rank + kol_size -1) % kol_size;
				sourceBB = (kol_rank + 1) % kol_size;
				MPI_Sendrecv_replace(my_block_B, chch, MPI_DOUBLE, destBB, 8, sourceBB, 8, kol_comm, &status);
			
			matrixMultiply(my_block_A, my_block_B, chunk, C);
		}
//Gather the results
		MPI_Gather(C, chch, MPI_DOUBLE, big_C, chch, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	end_time = MPI_Wtime(); /// Now that we have a solution for a C stop timer
  	time = end_time - start_time;
	// Save largest time
  	MPI_Reduce(&time, &time_big, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    if(rank==0){
    	// Allocate memory for final C
    		double** C_final = (double**) malloc(n*n * sizeof(double*));
		for(int i=0; i<n; i++){
	    	C_final[i] = (double*)malloc(n * sizeof(double)); }

	    // Fill the C final with the elements in correct order
	  	for(int p=0; p<size; p++){
	  		for(int i=0; i<chunk; i++){
	  			for(int j=0; j<chunk; j++){
	  				l = p % size_root;
	  				k = p / size_root;
	  				C_final[k*chunk+i][l*chunk+j] = big_C[p*chch+i*chunk+j];
	  			}
	  	 	}
	  	}
	  
	  	// Print results to file
	  	
	     FILE *file=fopen(output,"w");
	     for (int i = 0; i < n; i++){
	     	for (int j = 0; j < n; ++j){
	         	fprintf(file, "%.6lf ", C_final[i][j]);
				}
			}
	     fclose(file);
	     

	     // print largest time
    //printf("For #proc = %d and matrix n = %d the time was %lf \n", size, n, time_big);
    printf("%lf", time_big);

//// Free all allocated mem
	free(blocks_array_A);
	free(blocks_array_B);

		for(int l=0; l<n; l++){
			free(C_final[l]);
		}
		free(C_final);
		free(big_C);

	} // End of rank 0

	free(my_block_A);
	free(my_block_B);

	MPI_Comm_free(&row_comm);
	MPI_Comm_free(&kol_comm);

	MPI_Finalize(); 

	return 0;

}