#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "MPI_Additions.h"

int main(void){

  MPI_Init(NULL, NULL);
  init_MPI_Additions();
  pingpong(1,10000,20);
  MPI_Finalize();

  return 0;
}
