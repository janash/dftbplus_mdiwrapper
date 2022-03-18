#include <iostream>
#include <mpi.h>
#include <stdexcept>
#include <string.h>
#include "mdi.h"

using namespace std;

int main(int argc, char **argv) {

  // Initialize the MPI environment
  MPI_Comm world_comm;
  MPI_Init(&argc, &argv);

  // Initialize MDI
  if ( MDI_Init(&argc, &argv) ) {
    throw std::runtime_error("The MDI library was not initialized correctly.");
  }

  // Confirm that MDI was initialized successfully
  int initialized_mdi;
  if ( MDI_Initialized(&initialized_mdi) ) {
    throw std::runtime_error("MDI_Initialized failed.");
  }
  if ( ! initialized_mdi ) {
    throw std::runtime_error("MDI not initialized: did you provide the -mdi option?.");
  }

  // Get the correct MPI intra-communicator for this code
  if ( MDI_MPI_get_world_comm(&world_comm) ) {
    throw std::runtime_error("MDI_MPI_get_world_comm failed.");
  }

  // Connect to the engines
  MDI_Comm mm_comm = MDI_COMM_NULL;
  MDI_Comm qm_comm = MDI_COMM_NULL;
  int nengines = 2;
  for (int iengine=0; iengine < nengines; iengine++) {
    MDI_Comm comm;
    MDI_Accept_communicator(&comm);
 
    // Determine the name of this engine
    char* engine_name = new char[MDI_NAME_LENGTH];
    MDI_Send_command("<NAME", comm);
    MDI_Recv(engine_name, MDI_NAME_LENGTH, MDI_CHAR, comm);
 
    cout << "Engine name: " << engine_name << endl;
 
    if ( strcmp(engine_name, "MM") == 0 ) {
      if ( mm_comm != MDI_COMM_NULL ) {
	throw runtime_error("Accepted a communicator from a second MM engine.");
      }
      mm_comm = comm;
    }
    else if ( strcmp(engine_name, "QM") == 0 ) {
      if ( qm_comm != MDI_COMM_NULL ) {
	throw runtime_error("Accepted a communicator from a second QM engine.");
      }
      qm_comm = comm;
    }
    else {
      throw runtime_error("Unrecognized engine name.");
    }
 
    delete[] engine_name;
  }

  // Perform the simulation
  int niterations = 10;  // Number of MD iterations
  int natoms;
  double qm_energy;
  double mm_energy;
 
  // Receive the number of atoms from the MM engine
  MDI_Send_command("<NATOMS", mm_comm);
  MDI_Recv(&natoms, 1, MDI_INT, mm_comm);
 
  // Allocate the arrays for the coordinates and forces
  double coords[3*natoms];
  double forces[3*natoms];

  // Have the MM engine initialize a new MD simulation
  MDI_Send_command("@INIT_MD", mm_comm);
 
  // Perform each iteration of the simulation
  for (int iiteration = 0; iiteration < niterations; iiteration++) {

    // Have the MM engine proceed to the @FORCES node
    MDI_Send_command("@FORCES", mm_comm);

    // Receive the coordinates from the MM engine
    MDI_Send_command("<COORDS", mm_comm);
    MDI_Recv(&coords, 3*natoms, MDI_DOUBLE, mm_comm);
 
    // Send the coordinates to the QM engine
    MDI_Send_command(">COORDS", qm_comm);
    MDI_Send(&coords, 3*natoms, MDI_DOUBLE, qm_comm);
  
    // Get the QM energy
    MDI_Send_command("<ENERGY", qm_comm);
    MDI_Recv(&qm_energy, 1, MDI_DOUBLE, qm_comm);
 
    // Get the MM energy
    MDI_Send_command("<ENERGY", mm_comm);
    MDI_Recv(&mm_energy, 1, MDI_DOUBLE, mm_comm);
 
    // Receive the forces from the QM engine
    MDI_Send_command("<FORCES", qm_comm);
    MDI_Recv(&forces, 3*natoms, MDI_DOUBLE, qm_comm);
 
    // Send the forces to the MM engine
    MDI_Send_command(">FORCES", mm_comm);
    MDI_Send(&forces, 3*natoms, MDI_DOUBLE, mm_comm);
 
    cout << "timestep: " << iiteration << " " << mm_energy << " " << qm_energy << endl;
  }

  // Send the "EXIT" command to each of the engines
  MDI_Send_command("EXIT", mm_comm);
  MDI_Send_command("EXIT", qm_comm);

  // Synchronize all MPI ranks
  MPI_Barrier(world_comm);
  MPI_Finalize();

  return 0;
}
