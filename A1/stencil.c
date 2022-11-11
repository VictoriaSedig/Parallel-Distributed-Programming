#define PI 3.14159265358979323846
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ROOT 0

int main(int argc, char **argv) {
	if (3 != argc) {
		printf("Usage: stencil num_values num_steps\n");
		return 0;
	}
    int num_values = atoi(argv[1]);
	int num_steps = atoi(argv[2]);
    int dest_up, dest_down, source_up, source_down;
	MPI_Status status;
	MPI_Request request;

    MPI_Init(&argc, &argv);
    int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int num_values_each = num_values/size;

	// Generate values for stencil operation
	double *input=(double*)malloc(num_values*sizeof(double));
	double *input_each =(double*)malloc(num_values_each*sizeof(double));

    double h = 2.0*PI/num_values;
    if (rank==0){
    	for (int i=0; i<num_values; i++) input[i]=sin(h*i);
    }

	// Stencil values
	const int STENCIL_WIDTH = 5;
	const int EXTENT = STENCIL_WIDTH/2;
	const double STENCIL[] = {1.0/(12*h), -8.0/(12*h), 0.0, 8.0/(12*h), -1.0/(12*h)};

    // Allocate data for result
    double *output =(double*)malloc(num_values*sizeof(double));
    double *output_each =(double*)malloc(num_values_each*sizeof(double));
	double *times = (double*)malloc(size*sizeof(double));

	double *to_send_up = (double*)malloc(2*sizeof(double));
	double *to_send_down = (double*)malloc(2*sizeof(double));
	double *to_recv_up = (double*)malloc(2*sizeof(double));
	double *to_recv_down = (double*)malloc(2*sizeof(double));

	// start parralelisation, scatter
	double start = MPI_Wtime();
	int index;
	MPI_Scatter(input, num_values_each, MPI_DOUBLE, input_each, num_values_each, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
	// Prepair where processors shall send their messages
    dest_up = (rank+1) % size;
    dest_down = (rank+size-1) % size;
    source_up = (rank+1) % size;
    source_down = (rank+size-1) % size;



	// Repeatedly apply stencil
	for (int s=0; s<num_steps; s++) {
		// Create the arrays containing the courner values that will be sent 
	 	to_send_down[0] = input_each[0];
	 	to_send_down[1] = input_each[1];
	 	to_send_up[0] = input_each[num_values_each-2];
	 	to_send_up[1] = input_each[num_values_each-1];

	    MPI_Send(to_send_up, 2, MPI_DOUBLE, dest_up, 0, MPI_COMM_WORLD);
	    MPI_Send(to_send_down, 2, MPI_DOUBLE, dest_down, 0, MPI_COMM_WORLD);
	    
	    MPI_Recv(to_recv_down, 2, MPI_DOUBLE, source_down , 0, MPI_COMM_WORLD, &status);
	    MPI_Recv(to_recv_up, 2, MPI_DOUBLE, source_up, 0, MPI_COMM_WORLD, &status);


        // Apply stencil on left boundary with periodic cond
		for (int i=0; i<num_values_each; i++) {

			// use received value
			double result = 0;
			for (int j=0; j<STENCIL_WIDTH; j++) {
				index = i + j;
				if (index == 0){
					result += STENCIL[j] * to_recv_down[0];
				}
				else if (index == 1){
					result += STENCIL[j] * to_recv_down[1];

				}
				else if (index == num_values_each+2 ){
					result +=STENCIL[j] * to_recv_up[0];
				}
				else if (index == num_values_each+3 ){
					result +=STENCIL[j] * to_recv_up[1];
				}
				else{

				result += STENCIL[j] * input_each[index-2];
				}
			}
			output_each[i] = result;
		}
        
       
		// Swap input and output
		if (s < num_steps-1) {
			double *tmp = input_each;
			input_each = output_each;
			output_each = tmp;

		}
	}

	MPI_Gather(output_each, num_values_each, MPI_DOUBLE, output, num_values_each, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Stop timer
	double my_execution_time = MPI_Wtime() - start;
		

    MPI_Gather(&my_execution_time, 1, MPI_DOUBLE, times, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	if(rank==0){
		double big = times[0];
		for (int i = 1; i < size; ++i)
		{
			if(times[i]>big)
			{
				big = times[i];
			}
		}
		printf("Largest time is %f s\n", big );
	}

    // Write to file
    if(rank==0){
     FILE *file=fopen("output.txt","w");
     for (int i = 0; i < num_values; i++){
         fprintf(file, "%.4f \n", output[i]);
		}
     fclose(file);
     }
     
	// Clean up
    free(input);
	free(output);
	free(output_each);
	free(input_each);
	free(to_send_up);
	free(to_send_down);
	free(to_recv_up);
	free(to_recv_down);
	free(times);
    MPI_Finalize();
	return 0;
}
