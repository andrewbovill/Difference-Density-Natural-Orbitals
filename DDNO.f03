      Program DDNO
!
!     This program carries out the Density Difference Natural Orbital analysis
!     Computed the promotion number, excitation number, and the DDNO's for
!     an excited state calculation
!
!     -A. J. Bovill, 2022
!
!     Excitation number was proposed by Peter M. W. Gill to 
!     represent the number of electrons in the excited state that lie 
!     in the unoccupied space of the ground state
!     "Excitation Number: Characterizing Mutiply Excited States"
!     DOI: 10.1021/acs.jctc.7B00963
!
!     USE Connections

      use DDNO_Mod
!
!     Variable Declarations
!
      implicit none 
!     First file  filename_GS is Ground state
!     Second file filename_EX is Excited state
      character(len=512) :: filename_GS, filename_EX
      integer :: N_Electrons_GS, N_Electrons_GS_Alpha, N_Electrons_GS_Beta 
      integer :: N_Electrons_EX, N_Electrons_EX_Alpha, N_Electrons_EX_Beta 
      integer :: N_Basis_GS, N_Basis_EX
      integer :: i, j
      real :: Excitation_Number, Alpha_Excitation, Beta_Excitation

!     Arrays
      real, allocatable, dimension(:) :: Alpha_MO_Array_GS, Beta_MO_Array_GS
      real, allocatable, dimension(:) :: Alpha_MO_Array_EX, Beta_MO_Array_EX

!     Matrices
      real, allocatable, dimension(:,:) :: Alpha_MO_Matrix_GS, Beta_MO_Matrix_GS
      real, allocatable, dimension(:,:) :: Alpha_MO_Matrix_EX, Beta_MO_Matrix_EX

!     Matrices of occupied orbitals.
      real, allocatable, dimension(:,:) :: Alpha_OCC_AB, Beta_OCC_AB

!     Overlap matrix (identity matrix is just test)
      real, allocatable, dimension(:,:) :: identity_matrix
      real, allocatable, dimension(:,:) :: Overlap_Alpha
      real, allocatable, dimension(:,:) :: Overlap_Beta

!     Integers for write files.
      integer,parameter :: IUnit_GS=10
      integer,parameter :: IUnit_EX=20


!
!     Format Statements
!

 1000 Format(1x,'Enter Program DDNO.')
 1010 Format(1x,'Matrix File 1: ',A,/,  &
             1x,'Matrix File 2: ',A,/)

!     Read in Electrons and Basis number for Ground State
      write(*,1000)
      write(*,*) "What's the Ground State fchk file?" 
      call Get_Command_Argument(1,filename_GS)
      write(*,*) filename_GS
      open(File=trim(filename_GS), Unit = IUnit_GS)
      call Number_of_Electrons(filename_GS,N_Electrons_GS_Alpha,IUnit_GS,1)
      call Number_of_Electrons(filename_GS,N_Electrons_GS_Beta,IUnit_GS,2)
      call Number_of_Basis_Functions(filename_GS,N_Basis_GS,IUnit_GS)

!
      N_Electrons_GS = N_Electrons_GS_Alpha + N_Electrons_GS_Beta
      write(*,*) 'Number of Alpha Electrons in Ground State = ', N_Electrons_GS_Alpha  
      write(*,*) 'Number of Beta Electrons in Ground State = ', N_Electrons_GS_Beta  
      write(*,*) 'Number of Electrons in Ground State = ', N_Electrons_GS  
      write(*,*) 'Number of Basis in Ground State = ', N_Basis_GS  

!     Read into Arrays for both Alpha and Beta Ground States
      Allocate (Alpha_MO_Array_GS(N_Basis_GS*N_Basis_GS))
      call GET_MO_Array(IUnit_GS,N_Basis_GS,Alpha_MO_Array_GS,1)
      Alpha_MO_Matrix_GS = reshape(Alpha_MO_Array_GS, (/ N_Basis_GS,N_Basis_GS /))

      Allocate (Beta_MO_Array_GS(N_Basis_GS*N_Basis_GS))
      call GET_MO_Array(IUnit_GS,N_Basis_GS,Beta_MO_Array_GS,2)
      Beta_MO_Matrix_GS = reshape(Beta_MO_Array_GS, (/ N_Basis_GS,N_Basis_GS /))

      close (IUnit_GS, status='keep') 

!     Read in Electrons and Basis number for Excited State
      write(*,*)
      write(*,*) "What's the Excited State fchk file?" 
      call Get_Command_Argument(2,filename_EX)
      write(*,*) filename_EX
      open(File=trim(filename_EX), Unit = IUnit_EX)
      call Number_of_Electrons(filename_EX,N_Electrons_EX_Alpha,IUnit_EX,1)
      call Number_of_Electrons(filename_EX,N_Electrons_EX_Beta,IUnit_EX,2)
      call Number_of_Basis_Functions(filename_EX,N_Basis_EX,IUnit_EX)

      N_Electrons_EX = N_Electrons_EX_Alpha + N_Electrons_EX_Beta
      write(*,*) 'Number of Alpha Electrons in Excited State = ', N_Electrons_EX_Alpha  
      write(*,*) 'Number of Beta Electrons in Excited State = ', N_Electrons_EX_Beta  
      write(*,*) 'Number of Electrons in Excited State = ', N_Electrons_EX  
      write(*,*) 'Number of Basis in Excited State = ', N_Basis_EX  

!     Read into Arrays for both Alpha and Beta Excited States
      Allocate (Alpha_MO_Array_EX(N_Basis_EX*N_Basis_EX))
      call GET_MO_Array(IUnit_EX,N_Basis_EX,Alpha_MO_Array_EX,1)
      Alpha_MO_Matrix_EX = reshape(Alpha_MO_Array_EX, (/ N_Basis_EX,N_Basis_EX /))

      Allocate (Beta_MO_Array_EX(N_Basis_EX*N_Basis_EX))
      call GET_MO_Array(IUnit_EX,N_Basis_EX,Beta_MO_Array_EX,2)
      Beta_MO_Matrix_EX = reshape(Beta_MO_Array_EX, (/ N_Basis_EX,N_Basis_EX /))

      close (IUnit_EX, status='keep') 
     
      if(N_Electrons_GS.ne.N_Electrons_EX) then
        write(*,*) "The number of the electrons in both files do not match!"
        STOP
      end if

      if(N_Basis_GS.ne.N_Basis_EX) then
        write(*,*) "The number of the basis in both files do not match!"
        STOP
      end if
!
!     Generate Overlap for Alpha & Beta matrices
!     

      Overlap_Alpha = matmul(inv(transpose(Alpha_MO_Matrix_GS)),inv(Alpha_MO_Matrix_GS))
      Overlap_Beta = matmul(inv(transpose(Beta_MO_Matrix_GS)),inv(Beta_MO_Matrix_GS))

!
!     Compute the Alpha & Beta occupied matrix from the generated Overlap_Alpha &
!     Overlap_Beta
!
      ALLOCATE(Alpha_OCC_AB(N_Electrons_GS_Alpha,N_Electrons_EX_Alpha))
      ALLOCATE(Beta_OCC_AB(N_Electrons_GS_Beta,N_Electrons_EX_Beta))

      Alpha_OCC_AB = MatMul(Transpose(Alpha_MO_Matrix_GS(:,1:N_Electrons_GS_Alpha)),  &
         MatMul(Overlap_Alpha,Alpha_MO_Matrix_EX(:,1:N_Electrons_EX_Alpha)))

      Beta_OCC_AB = MatMul(Transpose(Beta_MO_Matrix_GS(:,1:N_Electrons_GS_Beta)),  &
         MatMul(Overlap_Beta,Beta_MO_Matrix_EX(:,1:N_Electrons_EX_Beta)))

      End Program DDNO
