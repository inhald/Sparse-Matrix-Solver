// Libraries
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

// Header
#include "functions.h"

void ReadMMtoCSR(const char *filename, CSRMatrix *matrix)
{

  // reading file
  FILE *file = fopen(filename, "r");
  char buffer[255];

  // catching errors
  if (file == NULL)
  {
    printf("Error Opening File.\n");
  }

  else
  {

    while ((fgets(buffer, 255, file)) != NULL)
    {
      char ch;
      ch = buffer[0];

      // ignores comments
      if (ch == '%')
      {
        continue;
      }
      // first line gives dimensions
      else
      {
        sscanf(buffer, "%d %d %d", &matrix->num_rows, &matrix->num_cols, &matrix->num_non_zeros);
        break;
      }
    }

    // building csr_data, col_indices and row_ptr
    matrix->csr_data = (double *)malloc(matrix->num_non_zeros * sizeof(double));
    matrix->col_ind = (int *)malloc(matrix->num_non_zeros * sizeof(int));
    matrix->row_ptr = (int *)malloc((matrix->num_rows + 1) * sizeof(int));
    int *row_ind = (int *)malloc((matrix->num_non_zeros) * sizeof(int));

    int row_index;

    // each non-zero value corresponds to an entry of the matrix
    for (int i = 0; i < matrix->num_non_zeros; ++i)
    {
      // scanning line and assigning it
      fgets(buffer, 255, file);
      sscanf(buffer, "%d %d %lf", &row_ind[i], &matrix->col_ind[i], &matrix->csr_data[i]);

      // row_ptr can be constructed as follows
      for (int j = row_ind[i]; j < matrix->num_rows + 1; ++j)
      {
        matrix->row_ptr[j]++;
      }
    }

    // selection sort to get the row indices in order
    int min;
    for (int i = 0; i < matrix->num_non_zeros - 1; ++i)
    {

      min = i;

      for (int j = i + 1; j < matrix->num_non_zeros; ++j)
      {

        if (row_ind[j] < row_ind[min])
        {
          min = j;
        }
      }

      if (min != i)
      {
        int temp = row_ind[i];
        row_ind[i] = row_ind[min];
        row_ind[min] = temp;

        double data = matrix->csr_data[i];
        matrix->csr_data[i] = matrix->csr_data[min];
        matrix->csr_data[min] = data;

        int col = matrix->col_ind[i];
        matrix->col_ind[i] = matrix->col_ind[min];
        matrix->col_ind[min] = col;
      }
    }

    // normalizing (going from matrix -> array type indicing)
    for (int i = 0; i < matrix->num_non_zeros; ++i)
    {
      matrix->col_ind[i]--;
    }

    //freeing memory 
    free(row_ind);
  }

  fclose(file);
}


//function for printing CSRMatrix
void printCSRMatrix(const CSRMatrix *matrix){

    //printing num_rows, num_cols and num_non_zeros
    printf("Number of rows: %d\n", matrix -> num_rows);
    printf("Number of columns: %d\n", matrix -> num_cols);
    printf("Number of non-zero elements %d\n", matrix->num_non_zeros);

    printf("CSR data\n");

    //printing data
    for(int i=0; i < matrix->num_non_zeros; ++i){
        printf("%lf ", matrix->csr_data[i]);
    }
    printf("\n");
    

    //printing column indices
    printf("col indices: ");

    for(int i=0; i < matrix->num_non_zeros; ++i){
        printf("%d ", matrix->col_ind[i]);
    }
    printf("\n");
    
    //printing row_ptr
    printf("row pointer: ");

    for(int i=0; i <=matrix->num_rows; ++i){
        printf("%d ", matrix->row_ptr[i]);
    }

    printf("\n");

}






void spmv_csr(const CSRMatrix *A, const double *x, double *y)
{


  //iterating over rows
  for (int i = 0; i < A->num_rows; ++i)
  {
    //iterating over each non-zero element in each row
    for (int j = A->row_ptr[i]; j < A->row_ptr[i + 1]; ++j)
    {
      y[i] += A->csr_data[j] * x[A->col_ind[j]];

    }
  }

}

double compute_norm(const double *r, int size)
{

  // computing dot product

  double sum = 0.0;

  for (int i = 0; i < size; ++i)
  {
    sum += r[i] * r[i];
  }

  // sqrt of dot product is norm

  return sqrt(sum);
}


void fullMatrix(CSRMatrix *A){

  //this function converts the lower triangular matrix to a symmetric matrix

  //calculating new num_non_zeros
  int num_non_zeros = A->num_non_zeros*2 - A->num_rows;

  //initializing
  double *csr_data = (double *)malloc(num_non_zeros*sizeof(double));
  int *col_ind = (int *)malloc(num_non_zeros*sizeof(int));
  int *row_ptr = (int *)malloc((A->num_rows+1)*sizeof(int));


  //calculating # of non-diag elements in cols
  int *colCount = (int *)calloc(A->num_rows, sizeof(int));

  for(int i =0; i < A->num_non_zeros; ++i){
    colCount[A->col_ind[i]]++;
  }

  row_ptr[0] = 0;

  for(int i =0; i < A->num_rows; ++i){
    colCount[i]--;
  }


  //changing it to be a cumulative sum like row_ptr
  for(int i=1; i < A->num_rows; ++i){
    colCount[i] = colCount[i] + colCount[i-1];
  } 

  //adding this to row_ptr adjusts it for the full matrix
  for(int i =1; i < A->num_rows+1; ++i){
    row_ptr[i] = A->row_ptr[i] + colCount[i-1];
  }


  //filling csr_data and col_ind
  for(int i =0; i < A->num_rows; ++i){
    //iterating over each non-zero element in lower triangular matrix
    for(int j =A->row_ptr[i]; j < A->row_ptr[i+1]; ++j){
      int col = A->col_ind[j];
      double value = A->csr_data[j];

      //two cases: 
      //1 for diagonal elems
      if(i != col){
        csr_data[row_ptr[i]] = value;
        col_ind[row_ptr[i]] = col;
        row_ptr[i]++;

        csr_data[row_ptr[col]]=value;
        col_ind[row_ptr[col]]=i;
        row_ptr[col]++;

      }
      //1 for non-diag
      else{
        csr_data[row_ptr[i]] =value;
        col_ind[row_ptr[i]] = col;
        row_ptr[i]++;
      }

    }
  }

  //correcting row_ptr
  for(int i = A->num_rows; i > 0; --i){
    row_ptr[i] = row_ptr[i-1];
  }

   row_ptr[0] = 0;

  //freeing old memory
  free(A->row_ptr);
  free(A->csr_data);
  free(A->col_ind);


  //replacing it 
  A->num_non_zeros = num_non_zeros;
  A->csr_data = csr_data;
  A->row_ptr = row_ptr;
  A->col_ind = col_ind;

}


void sorIteration(const CSRMatrix *A, double *x, double *b, double omega){
  
  int n= A->num_rows;
  //defining diagonal element and sigma
  double aii;
  double sigma;


  for(int i =0; i <n; ++i){
    
    sigma =0.0;
    
    for(int j=A->row_ptr[i]; j < A->row_ptr[i+1]; ++j){
      int col_ind = A->col_ind[j];

      //initializing aii
      if(i == col_ind){
        aii = A->csr_data[j];
      }
      //sigma is result of matrix multiplication 
      else{
        sigma +=  A->csr_data[j]*x[col_ind];
      }

    }
    //updating x
    x[i] = (1-omega)*x[i] + omega*(b[i] - sigma)/aii;



  }

}



void solver(CSRMatrix *A, double *x, double *b, double omega, double max_iter, double tolerance){


  //defining and initializing residual 
  double *residual = (double *)calloc(A->num_rows, sizeof(double));

  for(int i =0; i < A->num_rows; ++i){
    residual[i] = 1.0;
  }
  double resi_norm;
  int n = A->num_rows;
  int k=0;

  //loop will break after max_iter or if tolerance is acheived (1e-14)

  while((compute_norm(residual, n) > tolerance) && (k < max_iter)){
    spmv_csr(A, x, residual);

    //calc residual
    for(int i =0; i < n; ++i){
      residual[i]= b[i] - residual[i];
    }

    resi_norm = compute_norm(residual, n);

    if(resi_norm < tolerance){
      break;
    }

    //perform iteration of SOR
    sorIteration(A,x,b, omega);
    k++;
  }

  //freeing residual

  free(residual);


}



