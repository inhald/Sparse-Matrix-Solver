#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "functions.h"




int main(int argc, char *argv[]) {


    if(argc != 2){
        printf("The correct usage of this program is: ./main filename.mtx");
        exit(1);
    }

    // <Handle the inputs here>x    
    const char *filename = argv[1];

    //reading CSRMatrix A
    CSRMatrix A;
    ReadMMtoCSR(filename, &A);



    //Initializing all the vector b (in Ax=b)
    double *b = (double *)malloc(A.num_rows * sizeof(double));
    // Set all elements of b to 1
    for (int i = 0; i < A.num_rows; ++i)
    {
        b[i] = 1.0;
    }

    //Initializng x and setting guess
    double *x = (double *)malloc(A.num_rows * sizeof(double));
    for(int i =0; i < A.num_rows; ++i){
        x[i] = 1.0;
    }

    //defining omega
    double omega = 1.1;

    //converting A to full matrix and calling solver

    fullMatrix(&A);

    //intializing clock and start+end times
    clock_t start_time, end_time;
    double cpu_time_used;
    start_time = clock();


    solver(&A, x, b, omega, 100000, 1e-10);

    end_time = clock();
    cpu_time_used = ((double)(end_time - start_time)) / CLOCKS_PER_SEC;


    //residual and result of Ax
    double *r = (double *)malloc(A.num_rows * sizeof(double));
    double *result = (double*)malloc(A.num_rows * sizeof(double));


    spmv_csr(&A, x, result);

    //calculating residual vector
    for(int i =0; i < A.num_rows; ++i){
        r[i] = result[i] - b[i];
    
    }
    

    //calculating norm
    double norm;
    norm = compute_norm(r, A.num_rows);


    //handling output
    printf("The matrix name: %s\n", filename);
    printf("The dimension of the matrix: %d by %d\n", A.num_rows, A.num_cols);
    printf("The number of non-zeros: %d\n", A.num_non_zeros);
    printf("CPU time: %f seconds\n", cpu_time_used);
    printf("Residual Norm: %.14lf\n", norm);

    

    //freeing memory
    free(b);
    free(x);
    free(r);
    free(result);
    free(A.csr_data);
    free(A.row_ptr);
    free(A.col_ind);

}
