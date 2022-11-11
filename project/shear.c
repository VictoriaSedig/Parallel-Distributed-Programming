#include <mpi.h>
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ROOT 0

void *printlist(int *data, int len){
    for(int i=0; i<len; i++){
        printf("%d ",data[i] );}
    printf("\n");
}

int partition(int *data, int left, int right, int pivotIndex, char dir){
    int pivotValue,temp;
    int storeIndex,i;
    pivotValue = data[pivotIndex];
    temp=data[pivotIndex]; data[pivotIndex]=data[right]; data[right]=temp;
    storeIndex = left;
    for (i=left;i<right;i++)
    if(dir=='d'){	// Do for decending order
	    if (data[i] > pivotValue){
	        temp=data[i];data[i]=data[storeIndex];data[storeIndex]=temp;
	        storeIndex = storeIndex + 1;
    	}
	}
	else if(dir=='a'){	 // Do for acending order
	    if (data[i] <= pivotValue){
	        temp=data[i];data[i]=data[storeIndex];data[storeIndex]=temp;
	        storeIndex = storeIndex + 1;
    	}
	}
    temp=data[storeIndex];data[storeIndex]=data[right]; data[right]=temp;
    return storeIndex;
}

void quicksort(int *data, int left, int right, char dir){
    int pivotIndex, pivotNewIndex;
    if (right > left){
        pivotIndex = left+(right-left)/2;
        pivotNewIndex = partition(data, left, right, pivotIndex, dir);
        quicksort(data,left, pivotNewIndex - 1, dir);
        quicksort(data,pivotNewIndex + 1, right, dir);
    }
}



int main(int argc, char **argv)
{

// Assume that n is divisable of the number of processors
	// Initialize MPI 
    int rank, size;
	MPI_Status status;
	MPI_Request request;
    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	// Initialize values
	double time, start_time, end_time, time_big = 0;
	//int n = atoi(argv[1]);
	char* inputFile = argv[1];
	int store;
	int n;
	int *data, *my_data;
	FILE *fpInput;
	//// Let procs zero read the data

	if(rank==0){
		fpInput = fopen(inputFile, "r");
		store = fscanf(fpInput, "%d", &n);
		data = (int*) malloc(n*n*sizeof(int));
		for (int i=0; i<n*n; i++){
			store = fscanf(fpInput, "%d", &data[i]);

		}
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	int N = n*n;
	int d = ceil(log2(n));
	int num_steps = d+1;
	int part = n/size;
	int num_each = N/size;
	char ascending = 'a';
	char descending = 'd';
	int sendcount, recvcount;

	if(n%size != 0){
		printf("Error n is not divisable with number of procs\n");
		return 0;
	}

	// Do memory alllocations
    my_data=(int*)malloc(num_each*sizeof(int));
    int *to_send=(int*)malloc(part*part*sizeof(int));
    int *to_recv=(int*)malloc(part*part*sizeof(int));
    int *temp=(int*)malloc(part*part*sizeof(int));



  start_time = MPI_Wtime(); // start time measurement
	/// Proc 0 scatters the data to all procs
  	MPI_Scatter(data, num_each, MPI_INT, my_data , num_each, MPI_INT, 0, MPI_COMM_WORLD);

  	// The shear sort algorithm begins
  	for (int s = 0; s < d+1; s++){  // big loop


  		// Sort the arrays snakewise in rows
	  	if(rank%2==0 || part%2==0){
			for(int i=0; i<part; i+=2){
				quicksort(my_data,i*n, i*n+n-1, ascending);}
			for(int i=1; i<part; i+=2){
				quicksort(my_data,i*n, i*n+n-1, descending);}
		}
		else{
			for(int i=0; i<part; i+=2){
				quicksort(my_data,i*n, i*n+n-1, descending);}
			for(int i=1; i<part; i+=2){
				quicksort(my_data,i*n, i*n+n-1, ascending);}
		}


		// Send values to enable columnwise sort
		for (int i=0; i<size; i++){
		sendcount =0;
			for (int k = 0; k < part; k++){
				for(int j=0; j < part; j++){	
				to_send[sendcount]= my_data[k+j*n+i*part];
				sendcount++;}
	 		}
	 		if (rank != i){
	 			MPI_Sendrecv(to_send, part*part, MPI_INT, i, i, 
				to_recv, part*part, MPI_INT, i, rank, MPI_COMM_WORLD , &status);}
			else{
				for(int t=0; t<part*part; t++){
					to_recv[t]=to_send[t];}
			}
	 		recvcount =0;
			for (int k = 0; k < part; k++){
				for(int j=0; j<part; j++){
				my_data[i*part+j+k*part*size]=to_recv[recvcount];
				recvcount++;}
	 		}
		}


		// Sort the columns acsending order
	for(int i=0; i<part; i+=1){
		quicksort(my_data,i*n, i*n+n-1, ascending);
	}

	// Send back to row order in procs to be able to restart the loop
		for (int i=0; i<size; i++){
		sendcount =0;
			for (int k = 0; k < part; k++){
				for(int j=0; j < part; j++){	
				to_send[sendcount]= my_data[k+j*n+i*part];
				sendcount++;}
	 		}
	 		if (rank != i){
	 			MPI_Sendrecv(to_send, part*part, MPI_INT, i, i, to_recv, part*part, MPI_INT, i, rank, MPI_COMM_WORLD , &status);}
			else{
				for(int t=0; t<part*part; t++){
					to_recv[t]=to_send[t];}
			}
	 		recvcount =0;
			for (int k = 0; k < part; k++){
				for(int j=0; j<part; j++){
				my_data[i*part+j+k*part*size]=to_recv[recvcount];
				recvcount++;}
	 		}
		}

	}// End of big loop

	//Procs zero will gather the results
  	MPI_Gather(my_data, num_each, MPI_INT, data , num_each, MPI_INT, 0, MPI_COMM_WORLD);

  end_time = MPI_Wtime(); /// Now that we have a final result stop timer
  time = end_time - start_time;
  // Save largest time
  MPI_Reduce(&time, &time_big, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

  if(rank==0){
    printf("%d     %d     %lf\n", size, n ,time_big);
  }
   //Uncomment to write results to outputfile:
  
  if(rank==0){
  	FILE *fp=fopen("output.txt","w");
		for(int i=0; i<N; i++){
			if(i%n==0){
				fprintf(fp, "\n");
			}
			fprintf(fp, "%d ", data[i]);
		}
		fclose(fp);
	}
	
	
  	// Free allocations and finalize MPI 
  	if(rank==0){
		free(data);
	}
	free(to_send);
	free(to_recv);
	free(temp);
	free(my_data);

    MPI_Finalize(); 

	return 0;
}