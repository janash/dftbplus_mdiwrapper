MODULE DFTBPLUS_ENGINE

  USE mpi
  USE ISO_C_binding
  USE mdi,              ONLY : MDI_Init, MDI_Send, MDI_INT, MDI_CHAR, MDI_NAME_LENGTH, &
       MDI_Accept_communicator, MDI_Recv_command, MDI_Recv, &
       MDI_Set_execute_command_func, MDI_MPI_get_world_comm, MDI_DOUBLE, MDI_BYTE, &
       MDI_ENGINE, MDI_Get_role, MDI_Register_command, MDI_Register_node, &
       MDI_Register_callback, MDI_COMMAND_LENGTH, MDI_MPI_get_world_comm, &
       MDI_Plugin_get_argc, MDI_Plugin_get_arg
  
  
  USE dftbp_common_environment, ONLY : TEnvironment, TEnvironment_init
  USE dftbp_dftbplus_inputdata, ONLY : TInputData
  USE dftbp_dftbplus_initprogram, ONLY : TDftbPlusMain

  IMPLICIT NONE

  INTEGER, PARAMETER :: dp = selected_real_kind(15, 307)

  ! MDI Communicator to the driver
  INTEGER :: comm

  ! MPI intra-communicator for this code
  INTEGER :: world_comm

  ! Flag to terminate MDI response function
  LOGICAL :: terminate_flag = .false.

  ! Create DFTB+ instance and environment
  type(TEnvironment) :: env
  type(TDftbPlusMain) :: main
  
  ! Create DFTB+ input
  type(TInputData), allocatable :: input

CONTAINS

  FUNCTION MDI_Plugin_init_dftbplus() bind ( C, name="MDI_Plugin_init_dftbplus" )

    INTEGER :: MDI_Plugin_init_dftbplus
    INTEGER :: ierr
    INTEGER :: argc
    INTEGER :: iarg
    CHARACTER(LEN=1024) :: option
    CHARACTER(LEN=1024) :: mdi_options
    LOGICAL :: mdi_options_found

    real*8 :: initialEnergy

    ! Get the command-line options from the driver
    mdi_options_found = .false.
    CALL MDI_Plugin_get_argc(argc, ierr)
    DO iarg=0, argc-1
       CALL MDI_Plugin_get_arg(iarg, option, ierr)
       IF ( (TRIM(option) .eq. "-mdi") .or. (TRIM(option) .eq. "--mdi") ) THEN
          IF ( argc .gt. (iarg+1) ) THEN
             CALL MDI_Plugin_get_arg(iarg+1, mdi_options, ierr)
             mdi_options_found = .true.
          ELSE
             WRITE(6,*)'ERROR: argument to -mdi option not provided'
             MDI_Plugin_init_dftbplus = 1
             RETURN
          END IF
       END IF
    END DO
    IF ( .not. mdi_options_found ) THEN
       WRITE(6,*)'ERROR: -mdi option not provided'
       MDI_Plugin_init_dftbplus = 1
       RETURN
    END IF

    ! Call MDI_Init
    world_comm = MPI_COMM_WORLD
    CALL MDI_Init(mdi_options, ierr)

    ! Get the MPI intra-communicator over which this plugin will run
    CALL MDI_MPI_get_world_comm(world_comm, ierr);

    ! Perform one-time operations required to establish a connection with the driver
    CALL initialize_mdi()

    ! Respond to commands from the driver
    CALL respond_to_commands()

    MDI_Plugin_init_dftbplus = 0
  END FUNCTION MDI_Plugin_init_dftbplus


  SUBROUTINE initialize_mdi()
    USE dftbp_dftbplus_hsdhelpers, only : parseHsdInput 

    INTEGER :: ierr, role

    PROCEDURE(execute_command), POINTER :: generic_command => null()
    TYPE(C_PTR)                         :: class_obj
    generic_command => execute_command

    ! Confirm that the code is being run as an ENGINE
    call MDI_Get_role(role, ierr)
    IF ( role .ne. MDI_ENGINE ) THEN
       WRITE(6,*)'ERROR: Must run dftbplus as an ENGINE'
    END IF

    ! Register the commands
    CALL MDI_Register_node("@DEFAULT", ierr)
    CALL MDI_Register_command("@DEFAULT", "EXIT", ierr)
    CALL MDI_Register_command("@DEFAULT","<CELL", ierr)
    CALL MDI_Register_command("@DEFAULT", "<CELL_DISPL", ierr)
    CALL MDI_Register_command("@DEFAULT", "<CHARGES", ierr)
    CALL MDI_Register_command("@DEFAULT", "<COORDS", ierr)
    CALL MDI_Register_command("@DEFAULT","<ENERGY", ierr)
    CALL MDI_Register_command("@DEFAULT","<NATOMS", ierr)

    ! Connct to the driver
    CALL MDI_Accept_communicator(comm, ierr)

    ! Set the generic execute_command function
    CALL MDI_Set_execute_command_func(generic_command, class_obj, ierr)

    ! Initialize DFTB+
    ! Initialize DFTB+ Environment
    CALL TEnvironment_init(env)

    ! Parse DFTB+ input
    allocate(input)
    CALL parseHsdInput(input)

    ! Initialize variables according to input
    call main%initProgramVariables(input, env)

  END SUBROUTINE initialize_mdi


  SUBROUTINE respond_to_commands()
    CHARACTER(len=:), ALLOCATABLE :: command
    INTEGER :: ierr

    ALLOCATE( character(MDI_COMMAND_LENGTH) :: command )

    ! Respond to the driver's commands
    response_loop: DO

       ! Receive a command from the driver and broadcast it to all ranks
       CALL MDI_Recv_command(command, comm, ierr)
       CALL MPI_Bcast(command, MDI_COMMAND_LENGTH, MPI_CHAR, 0, world_comm, ierr)

       CALL execute_command(command, comm, ierr)

       IF ( terminate_flag ) EXIT

    END DO response_loop

    DEALLOCATE( command )

  END SUBROUTINE respond_to_commands


  SUBROUTINE execute_command(command, comm, ierr)
    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN)  :: command
    INTEGER, INTENT(IN)           :: comm
    INTEGER, INTENT(OUT)          :: ierr
 
    SELECT CASE( TRIM(command) )
    CASE( "EXIT" )
       terminate_flag = .true.
   case( "<CELL" )
      call send_cell(comm)
   case( "<CELL_DISPL" )
      call send_cell_displ(comm)
   case( "<CHARGES" )
      call send_charges(comm)
   case( "<COORDS" )
      call send_coords(comm)
   case( "<ENERGY" )
      call send_energy(comm)
   case( "<NATOMS" )
      call send_natoms(comm)
    CASE DEFAULT
       WRITE(6,*)'Error: command not recognized'
       STOP 1
    END SELECT

    ierr = 0
  END SUBROUTINE execute_command


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Corresponds to <ENERGY
   ! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE send_energy(comm)
   USE mdi , ONLY    : MDI_DOUBLE, MDI_Send
   USE dftbp_dftbplus_mainapi, ONLY : getEnergy

   implicit none
   integer, intent(in)          :: comm
   integer                      :: ierr
   real*8                       :: etotal

   ! DFTB+ uses Hartrees, so no need to convert.
   call getEnergy(env, main, etotal)
   call MDI_Send(etotal, 1, MDI_DOUBLE, comm, ierr)

  END SUBROUTINE send_energy

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Corresponds to <CELL
   ! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE send_cell(comm)
   USE mdi, ONLY     : MDI_DOUBLE, MDI_Send

   implicit none
   integer, intent(in)     :: comm
   integer                 :: ierr
   real*8, dimension(9)    :: latVecs

   ! DFTB+ uses atomic length units, so no need to convert.
   latVecs = RESHAPE(main%latVec, (/9/))
   call MDI_Send(latVecs, 9, MDI_Double, comm, ierr)

  END SUBROUTINE send_cell


   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Corresponds to <CELL
   ! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE send_cell_displ(comm)
   USE mdi, only     : MDI_DOUBLE, MDI_Send

   implicit none
   integer, intent(in)  :: comm
   integer              :: ierr
   real*8, dimension(3) :: origin

   ! retrieve the origin of the cell
   call MDI_Send(origin, 3, MDI_Double, comm, ierr)

  END SUBROUTINE send_cell_displ

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Corresponds to <CHARGES
   ! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE send_charges(comm)
   USE mdi, ONLY       : MDI_DOUBLE, MDI_Send
   USE dftbp_dftbplus_mainapi, ONLY : getGrossCharges, nrOfAtoms

   implicit none
   integer, intent(in)        :: comm
   integer                    :: ierr
   integer                    :: nAtoms
   real*8, allocatable        :: atomCharges(:)

   nAtoms = nrOfAtoms(main)

   allocate( atomCharges(nAtoms) )

   call getGrossCharges(env, main, atomCharges)
   call MDI_Send(atomCharges, nAtoms, MDI_DOUBLE, comm, ierr)

   END SUBROUTINE send_charges

   SUBROUTINE send_coords(comm)
      USE mdi, ONLY          : MDI_DOUBLE, MDI_Send
      USE dftbp_dftbplus_mainapi, ONLY : nrOfAtoms

      implicit none
      integer, intent(in)        :: comm
      integer                    :: ierr
      integer                    :: nAtoms
      integer                    :: nCoords
      real*8, allocatable        :: atomCoords(:)

      nAtoms = nrOfAtoms(main)
      nCoords = 3 * nAtoms

      allocate( atomCoords(nCoords) )

      atomCoords = RESHAPE( main%coord0, (/nCoords/) )

      call MDI_Send(atomCoords, nCoords, MDI_DOUBLE, comm, ierr)

   END SUBROUTINE send_coords

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !
   ! Corresponds to <NATOMS
   ! 
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   SUBROUTINE send_natoms(comm)
      USE mdi, ONLY          : MDI_INT, MDI_Send
      USE dftbp_dftbplus_mainapi, ONLY : nrOfAtoms

      implicit none
      integer, intent(in)        :: comm
      integer                    :: ierr
      integer                    :: nAtoms

      nAtoms = nrOfAtoms(main)
      
      call MDI_Send(nAtoms, 1, MDI_INT, comm, ierr)

   END SUBROUTINE send_natoms



END MODULE DFTBPLUS_ENGINE
