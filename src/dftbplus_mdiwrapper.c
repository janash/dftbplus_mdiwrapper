#include <mpi.h>
#include <string.h>
#include "mdi.h"
#include "dftbplus.h"
#include <stdlib.h>

int main(int argc, char **argv) {

  // Initialize the MPI environment
  MPI_Comm world_comm = MPI_COMM_WORLD;
  MPI_Init(&argc, &argv);

  DftbPlus *instance = malloc(sizeof(DftbPlus));
  dftbp_init(instance, NULL);

  // Synchronize all MPI ranks
  MPI_Barrier(world_comm);
  MPI_Finalize();

  return 0;
}
