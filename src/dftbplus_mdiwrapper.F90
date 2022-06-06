PROGRAM DRIVER_F90

USE mpi
USE mdi

! From dftbplus.F90
use dftbp_common_environment, only : TEnvironment, TEnvironment_init
use dftbp_common_globalenv, only : initGlobalEnv, destructGlobalEnv
use dftbp_dftbplus_hsdhelpers, only : parseHsdInput
use dftbp_dftbplus_initprogram, only : TDftbPlusMain
use dftbp_dftbplus_inputdata, only : TInputData
use dftbp_dftbplus_main, only : runDftbPlus
use dftbp_io_formatout, only : printDftbHeader

! main api
use dftbp_dftbplus_mainapi, only : getEnergy, nrOfAtoms

IMPLICIT NONE

   INTEGER :: world_comm, ierr
   character(len=255) :: cwd
   character(len=255) :: hsdFileName
   real*8 :: initialEnergy
   real*8 :: finalEnergy


   type(TEnvironment) :: env
   type(TInputData), allocatable :: input

   ! This has all the system info
   type(TDftbPlusMain) :: main

   ! Initialize the MPI environment
   world_comm = MPI_COMM_WORLD;
   call MPI_Init(ierr)

   ! Initialize DFTB+ environment
   call initGlobalEnv()

   ! Set up file iinformation
   ! hsdFileName = "/dftb_in_.hsd"
   call getcwd(cwd)

   ! hsdFileName = TRIM(cwd) // TRIM(hsdFileName)

   ! Allocate the input
   ! Unsure if it needs to be done this way
   allocate(input)

   !! Read HSD and put into input object.

   ! This assumes the filename is dftb_in.hsd.
   ! This is set in hsdhelpers.F90, where parseHsdInput is defined.
   ! Unsure if there is a way to change the file input name if this function is used.
   ! File input name could be different if this function is not used (could follow example)
   ! from this routine definition. Involves several more steps and directly dealing with
   ! hsd tree. 
   call parseHsdInput(input)

   ! Initialize the environment.
   call TEnvironment_init(env)

   ! I'm guessing this populates everything
   call main%initProgramVariables(input, env)

   ! Get the energy using the main api
   call getEnergy(env, main, initialEnergy)
   
   ! run the calculation
   call runDftbPlus(main, env)
   

   ! Get the energy using the main api
   call getEnergy(env, main, finalEnergy)

   ! write the energy
   write(*,*) "Initial Energy:", initialEnergy
   write(*,*) "Final Energy:", finalEnergy

   ! Access number of atoms
   write(*,*) "Number of Atoms (from main object):", main%nAtom
   write(*,*) "Number of Atoms (using main api):", nrOfAtoms(main)


   ! Synchronize all MPI ranks
   call MPI_Barrier( world_comm, ierr )
   call MPI_Finalize( ierr )

   ! Destruct
   call env%destruct()
   call destructGlobalEnv()


END PROGRAM DRIVER_F90
