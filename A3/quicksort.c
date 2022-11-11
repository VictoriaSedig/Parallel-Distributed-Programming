/**********************************************************************
 * Quick sort
 * Usage: ./a.out sequence length
 *
 **********************************************************************/
#define PI 3.14159265358979323846
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#define ROOT 0

int partition(double *data, int left, int right, int pivotIndex);
void quicksort(double *data, int left, int right);
void *merge(double *v1, int n1, double *v2, int n2, double* result);
void *printlist(double *data, int len);
double oneprocmedian(double *my_data, int mynum);
double oneprocmean(double *my_data, int mynum);
int find_split_index(double *data, int len, double pivot);
double get_wall_seconds(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    double seconds = tv.tv_sec + (double)tv.tv_usec/1000000;
    return seconds;
}
struct stuff {
    double *arr_smaller;
    int length;
};
typedef struct stuff Struct;
Struct paraquicksort(double* my_data, int mynum, int splitcount, MPI_Comm comm, int pivotmethod, int k);


int main(int argc, char *argv[]) {
    
// Start upp the MPI procs definitions
    int rank, size;
    MPI_Status status;
    MPI_Request request;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    /// Var definitions
    double *data_ref;
    double *data;
    int size_root = (int) sqrt(size);
    int sendcounts[size],displs[size];
    int splitcount = 0;
    int my_final_len;
    double time, start_time, end_time, time_big = 0, sec1, sec2;
    double drand48(void);
    int i,len,seq, pivotmethod, chunk;
    seq=atoi(argv[1]);
    len=atoi(argv[2]);
    pivotmethod=atoi(argv[3]);
    int rest = len % size;
    chunk = (len-rest)/size;
    /// There should be 2^k number of procs
    int k = (int) log2(size);
    /// Prepair sizes for the scatterv
    int count=0;
    for(int i=0; i<size; i++){
        if(i<(size-rest)){
            sendcounts[i]=chunk;
            displs[i]=chunk*i;
            }

        else{
            sendcounts[i] = chunk+1;
            displs[i] = chunk*(size-rest)+count*(chunk+1);
            count++;
        }
    }

    int mynum = sendcounts[rank];
    double *my_data =(double *)malloc(mynum*sizeof(double));

    // Fill data with random double values
    if(rank==0){
        data=(double *)malloc(len*sizeof(double));
        data_ref=(double *)malloc(len*sizeof(double));

        if (seq==0) {
            
            // Uniform random numbers
            for (i=0;i<len;i++){
            data[i]=drand48();
            }
        }
        else if (seq==1) {
            // Exponential distribution
            double lambda=10;
            for (i=0;i<len;i++){
            data[i]=-lambda*log(1-drand48());
            }
        }
        
        else if (seq==2) {
            // Normal distribution
            double x,y;
            for (i=0;i<len;i++){
                x=drand48(); y=drand48();
                data[i]=sqrt(-2*log(x))*cos(2*PI*y);

            }
        }

        else if (seq==3) {
            for (i=0;i<len;i++){
                data[len-i-1]=i*0.3;
            }
        }
        // Make ref data to verify correctness
        for (int i = 0; i < len; ++i)
            {
                data_ref[i]=data[i];
            }
    }

    // Give each procs its list
    MPI_Scatterv(data, sendcounts, displs, MPI_DOUBLE, my_data , mynum, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    

//// Start the parrallel quicksort
  start_time = MPI_Wtime(); // start time measurement
    quicksort(my_data,0,mynum-1);
    Struct result; // Making a result struckt to put results in
    result = paraquicksort(my_data, mynum, splitcount, MPI_COMM_WORLD, pivotmethod, k);
    double *final_arr;
    my_final_len = result.length;
    final_arr =(double *)malloc(my_final_len*sizeof(double));
    final_arr = result.arr_smaller;
    MPI_Gather(&my_final_len, 1, MPI_INT, sendcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int move=0;
    for(int i=0; i<size; i++){
        displs[i]=move;
        move += sendcounts[i];
    }

    MPI_Gatherv(final_arr, my_final_len, MPI_DOUBLE, data, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    end_time = MPI_Wtime(); /// Now that we have a solution for a C stop timer
    time = end_time - start_time; 

    MPI_Reduce(&time, &time_big, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    /// End of parrallell quicksort


    if(rank==0){
    // Check results
    int OK=1;
    for (i=0; i<len-1; i++) {
        if(data[i] > data[i+1]) {
            printf("Wrong result: data[%d] = %f, data[%d] = %f\n", i, data[i], i+1, data[i+1]);
            OK=0;
        }
    }
    
    if (OK) printf("Data sorted correctly!\n");
    
        /// Print time measurement
    printf("Para time for #procs=%d with seq= %d, len=%d, pivot=%d, was %lf \n",size, seq, len, pivotmethod,time_big);
    }

//Free memory allocations
    free(final_arr);
    if(rank==0){
        free(data);
        free(data_ref);
    }
    MPI_Finalize(); 

    return 0;
}

int partition(double *data, int left, int right, int pivotIndex){
    double pivotValue,temp;
    int storeIndex,i;
    pivotValue = data[pivotIndex];
    temp=data[pivotIndex]; data[pivotIndex]=data[right]; data[right]=temp;
    storeIndex = left;
    for (i=left;i<right;i++)
    if (data[i] <= pivotValue){
        temp=data[i];data[i]=data[storeIndex];data[storeIndex]=temp;
        storeIndex = storeIndex + 1;
    }
    temp=data[storeIndex];data[storeIndex]=data[right]; data[right]=temp;
    return storeIndex;
}

void quicksort(double *data, int left, int right){
    int pivotIndex, pivotNewIndex;
    
    if (right > left){
        pivotIndex = left+(right-left)/2;
        pivotNewIndex = partition(data, left, right, pivotIndex);
        quicksort(data,left, pivotNewIndex - 1);
        quicksort(data,pivotNewIndex + 1, right);
    }
}

void *merge(double *v1, int n1, double *v2, int n2, double *result)
{
    int i,j,k;
        
    i=0; j=0; k=0;
    while(i<n1 && j<n2)
        if(v1[i]<v2[j])
        {
            result[k] = v1[i];
            i++; k++;
        }
        else
        {
            result[k] = v2[j];
            j++; k++;
        }
    if(i==n1)
        while(j<n2)
        {
            result[k] = v2[j];
            j++; k++;
        }
    else
        while(i<n1)
        {
            result[k] = v1[i];
            i++; k++;
        }

}
void *printlist(double *data, int len){

    for(int i=0; i<len; i++){
        printf("%lf ",data[i] );
    }
    printf("\n");

}

double oneprocmedian(double *my_data, int mynum){
    

    if(mynum%2 == 0){
        return (my_data[mynum/2]+my_data[(mynum/2) -1])/2;
    }
    else{
        return my_data[(mynum-1)/2];
    }
}


double oneprocmean(double *my_data, int mynum){

    double mean=0;
    for(int i=0; i<mynum; i++){
        mean += my_data[i];
    }
    return mean/mynum;
}


int find_split_index(double *data, int len, double pivot){
 // Starting to search for the split from the middle of the array
   int half = (int) round((len)/2)+1;

   if(len==0){
    return -1;
   }
    else if(data[half]<pivot){
        for (int i = half; i < len; i++)
        {
            if(data[i]>pivot){
                return i-1;
            }
        }
    return len-1;
    }
    else if(data[half]>=pivot){
        for (int j = half; j >= 0; j-=1)
        {
            if(data[j]<=pivot){
                return j;
            }
        }
        return -1;
        }
    else{
        return half;
    }

}


Struct paraquicksort(double* my_data, int mynum, int splitcount, MPI_Comm comm, int pivotmethod, int k){
    int newsize, size, rank;
    double my_pivot, pivot;
    int midindex;
    double *new_list, *to_recv;

    MPI_Status status;
    MPI_Comm new_comm;
    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);


  if (splitcount == k) { /* We are done, do not have to proceed */
    Struct s;
    s.arr_smaller = my_data;
    s.length = mynum;
    return s;
  }
/// Calculating pivot depending on method
    if(pivotmethod==0){
        if (rank == 0){
        pivot = oneprocmedian(my_data, mynum);
        }
        MPI_Bcast(&pivot, 1, MPI_DOUBLE, 0, comm);   
    }
    else if(pivotmethod==1){
        my_pivot = oneprocmedian(my_data, mynum);
        double allmeds[size];
        MPI_Gather(&my_pivot, 1, MPI_DOUBLE, &allmeds, 1, MPI_DOUBLE, 0, comm);

        if(rank==0){
            quicksort(allmeds, 0, size-1);
            pivot = oneprocmedian(allmeds, size);
        }
        MPI_Bcast(&pivot, 1, MPI_DOUBLE, 0, comm);   
        // Reduce
    }
    else{
        my_pivot = oneprocmedian(my_data, mynum);
        double sum_pivot = 0;
        MPI_Reduce(&my_pivot, &sum_pivot, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if(rank==0){
            pivot = sum_pivot/size;
        }
        MPI_Bcast(&pivot, 1, MPI_DOUBLE, 0, comm);   
    }
    

    midindex = find_split_index(my_data, mynum, pivot);


    int size_to_send;
    int size_to_recv;
    int new_mynum;

    /// Perform the send and receivs from one group to another 
    if(rank<size/2){
        size_to_send = mynum-midindex-1;
        MPI_Send(my_data + midindex+1, size_to_send, MPI_DOUBLE, rank+size/2,1, comm);
        MPI_Probe(rank + size/2, 3, comm, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &size_to_recv);
        to_recv = (double*)malloc(size_to_recv*sizeof(double));
        MPI_Recv(to_recv, size_to_recv, MPI_DOUBLE, rank+size/2, 3, comm, &status);
        new_mynum = midindex+1 + size_to_recv;
        new_list= (double *) malloc(new_mynum*sizeof(double));
        merge(my_data, midindex+1, to_recv, size_to_recv, new_list);
        my_data = (double *) realloc(my_data,new_mynum*sizeof(double));

    }
    else{
        MPI_Probe(rank - size/2, 1, comm, &status);
        MPI_Get_count(&status, MPI_DOUBLE, &size_to_recv);
        to_recv = (double*)malloc(size_to_recv*sizeof(double));
        MPI_Recv(to_recv, size_to_recv, MPI_DOUBLE, rank-size/2, 1, comm, &status);
        size_to_send = midindex+1;
        MPI_Send(my_data, size_to_send, MPI_DOUBLE, rank-size/2,3, comm);
        new_mynum = mynum-midindex-1 + size_to_recv;
        new_list= (double *) malloc(new_mynum*sizeof(double));
        merge(my_data+midindex+1, mynum-midindex-1, to_recv, size_to_recv, new_list); 
        my_data = (double *) realloc(my_data,new_mynum*sizeof(double));

    }
    // put result in mydata
   for (int i = 0; i < new_mynum; ++i){
        my_data[i]=new_list[i];    
     }
     
    ///Split into two new groups
    int color = rank/(size/2);
    MPI_Comm_split(comm, color, rank, &new_comm);
    splitcount++;
    // Free before the recursive call
    free(new_list);
    free(to_recv);
    //Recursive call
    return paraquicksort(my_data, new_mynum, splitcount, new_comm, pivotmethod, k);
}