#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define TOL 1.0e-2

void communicate(double *matrix, int n_loc, int dim, int npes, int me){
  // Send up and send down
  if(me != 0){
    // send up the 2nd row (idx = 1) of previous processor (tag = 1)
    MPI_Send(matrix + dim, dim, MPI_DOUBLE, me-1, 1, MPI_COMM_WORLD);
  }
  if(me!= npes-1){
    // send down the second-to-the-last row (idx = n_loc-2) of next processor (tag = 2)
    MPI_Send(matrix + (n_loc-2)*dim, dim, MPI_DOUBLE, me+1, 2, MPI_COMM_WORLD);
  }

  // Receive from above and and from below
  if(me!= 0){
    // receive from above and update 1st row (tag = 2)
    MPI_Recv(matrix, dim, MPI_DOUBLE, me-1, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  if(me!= npes-1){
    // receive from below and update last row (tag = 1)
    MPI_Recv(matrix + (n_loc-1)*dim, dim, MPI_DOUBLE, me+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
}

void evolve(double* matrix, int n_loc, int dim, int npes, int me, int offset, int max_iter){
	double *new_mat = (double*)malloc(sizeof(double) * dim * n_loc);
  if (new_mat == NULL){
    printf("[ERROR] Check memory allocation for new matrix.\n");
    exit(-1);
  }

  FILE *file;
  file = fopen( "error.dat", "w");

  double west, east, h = 0.1;
  double increment = 100.0 / (dim+1);
  for (size_t it = 0; it < max_iter; it++) {
		double err_chunk = 0;
		for (int i = 1; i < n_loc-1; i++){
			for(int j = 0; j < dim; j++){
        // stencil west
        if (j == 0) west = (me*(n_loc-2) + i + offset) * increment;
        else        west = matrix[i*dim + j-1];
  			// stencil east
        if (j == dim-1) east = 0;
        else            east = matrix[i*dim + j+1];

        new_mat[i*dim + j] = 0.25*( matrix[(i-1)*dim + j] + matrix[(i+1)*dim + j] + west + east );

				err_chunk += pow(new_mat[i*dim +j] - matrix[i*dim +j], 2);
			}
		}

    // Replace matrix values by new_mat values
    for (int i=1; i < n_loc-1; i++){
      for(int j=0; j < dim; j++){
        matrix[i*(dim) + j] = new_mat[i*(dim) + j ];
      }
    }
    communicate(matrix, n_loc, dim, npes, me);

		double error = 0;
		MPI_Allreduce(&err_chunk, &error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if (me == 0) {
      printf("it %zu:\t %.7f\n", it, sqrt(error) );
      fprintf(file, "%.7f\n", sqrt(error) );
    }
    if (sqrt(error) < TOL) break;
	}
  fclose(file);
	free(new_mat);
}


void export_solution(double* matrix, int me, int dim, int n_loc){
  FILE *file;
  file = fopen( "solution.dat", "w");

  double *sol_mat = (double*)malloc(sizeof(double) * dim * dim);
  if (sol_mat == NULL){
    printf("[ERROR] Check memory allocation for the solution matrix.\n");
    exit(-1);
  }

  MPI_Gather(matrix + dim, dim*(n_loc-2), MPI_DOUBLE, sol_mat, dim*(n_loc-2), MPI_DOUBLE, 0, MPI_COMM_WORLD);

  const double h = 0.1;
  if (me == 0){
    for (int i = 0; i < dim; i++){
      for (int j=0; j < dim; j++){
          fprintf(file, "%f %f %f\n", (j+1)*h, -(i+1)*h, sol_mat[i*dim + j]);
      }
    }
    fclose(file);
  }
  free(sol_mat);
}


int main( int argc, char *argv[] )
{
  if (argc < 3){
    printf("Argument missing.\n");
    exit(-1);
  }
  int dim = atoi(argv[1]);
  int max_iter = atoi(argv[2]);

  int me, npes, n_loc, rest, offset=0;
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &me );
  MPI_Comm_size( MPI_COMM_WORLD, &npes );
  MPI_Barrier(MPI_COMM_WORLD);  //synchronize all processes

  double t_start, t_end, t_elapsed;
  t_start = MPI_Wtime();

  n_loc = dim/npes + 2;
  rest = dim % npes;
  if( me < rest ) n_loc += 1;
  else offset = rest;

  double* matrix = NULL;
  matrix = (double *)  malloc(sizeof(double) * dim * n_loc);
  if (matrix == NULL){
    printf("[ERROR] Check memory allocation for the matrix.\n");
    exit(-1);
  }
  memset(matrix, 0.5, sizeof(double) * dim * n_loc);


  double increment = 100.0 / (dim+1);
  if (me==0) {  // top
    for (size_t j = 0; j < dim; j++) *(matrix + j) = 0;
  }
  if (me == npes-1){  // bottom
    for (int j = 0; j < dim; j++){
      matrix[(n_loc-1)*dim + j] = (dim - j) * increment;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  communicate(matrix, n_loc, dim, npes, me);
  evolve(matrix, n_loc, dim, npes, me, offset, max_iter);
  t_end = MPI_Wtime();
  t_elapsed = t_end - t_start;

  double total_time = 0.0;
  MPI_Reduce(&t_elapsed, &total_time, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if (me==0) printf("Time elapsed: %.7f\n", total_time);

  export_solution(matrix, me, dim, n_loc);
  MPI_Finalize();
  free(matrix);
  return 0;
}
