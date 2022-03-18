#include <iostream>
#include <mpi.h>
#include <stdexcept>
#include <string.h>
#include "mdi.h"
#include "dftbplus.h"

using namespace std;

int main(int argc, char **argv) {

  // Initialize the MPI environment
  MPI_Comm world_comm;
  MPI_Init(&argc, &argv);

  DftbPlus *instance;
  dftbp_init(instance, NULL);

  // Synchronize all MPI ranks
  MPI_Barrier(world_comm);
  MPI_Finalize();

  return 0;
}
