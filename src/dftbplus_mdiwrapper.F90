PROGRAM DRIVER_F90

USE mpi
USE mdi
USE dftbp_mmapi
USE dftbp_dftbplus_main, only : runDftbPlus

IMPLICIT NONE

   INTEGER :: world_comm, ierr
   type(TDftbPlus) :: instance
   type(TDftbPlusInput) :: input
   character(len=255) :: cwd
   character(len=255) :: hsd_file
   real*8 :: energy

   ! Initialize the MPI environment
   world_comm = MPI_COMM_WORLD;
   call MPI_Init(ierr)

   call TDftbPlus_init(instance, 6)

   ! Set up file iinformation
   hsd_file = "/dftb_in.hsd"
   call getcwd(cwd)

   hsd_file = TRIM(cwd) // TRIM(hsd_file)

   ! read file
   call instance%getInputFromFile(hsd_file, input)
   call instance%setupCalculator(input)

   ! get the energy
   call instance%getEnergy(energy)
   !write(*,*) energy

   ! This is a function that Jessica wrote.
   call instance%run()

   ! run the calculation
   ! call runDftbPlus(instance%main, instance%env)
   !write(*,*) energy


   ! Synchronize all MPI ranks
   call MPI_Barrier( world_comm, ierr )
   call MPI_Finalize( ierr )

END PROGRAM DRIVER_F90
