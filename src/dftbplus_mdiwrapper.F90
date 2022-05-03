PROGRAM DRIVER_F90

USE mpi
USE mdi
USE dftbp_mmapi

IMPLICIT NONE

   INTEGER :: world_comm, ierr
   type(TDftbPlus) :: instance

   ! Initialize the MPI environment
   world_comm = MPI_COMM_WORLD;
   call MPI_Init(ierr)

   call TDftbPlus_init(instance, 6)

   ! Synchronize all MPI ranks
   call MPI_Barrier( world_comm, ierr )
   call MPI_Finalize( ierr )

END PROGRAM DRIVER_F90
