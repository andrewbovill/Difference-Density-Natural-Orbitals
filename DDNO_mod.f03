      Module DDNO_mod
!
!
!     Program uses the function "function inv(a)" from "LAPACK_Helper" to invert a GE Matrix 
!     You can find the source code here.
!     https://github.com/b-fg/LAPACK_helper/blob/master/lapack_helper.f90
!
!


      Contains

      function inv(A) result(Ainv)
        implicit none
        real,intent(in) :: A(:,:)
        real            :: Ainv(size(A,1),size(A,2))
        real            :: work(size(A,1))            ! work array for LAPACK
        integer         :: n,info,ipiv(size(A,1))     ! pivot indices

      ! Store A in Ainv to prevent it from being overwritten by LAPACK
        Ainv = A
        n = size(A,1)
      ! SGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
        call SGETRF(n,n,Ainv,n,ipiv,info)
        if (info.ne.0) stop 'Matrix is numerically singular!'
        ! SGETRI computes the inverse of a matrix using the LU factorization
      ! computed by SGETRF.
        call SGETRI(n,Ainv,n,ipiv,work,n,info)
        if (info.ne.0) stop 'Matrix inversion failed!'
      end function inv

      Subroutine Number_of_Electrons(filename,NElectrons,IUnit,Spin)

!     This subroutine reads the input file line by line for the number of 
!     Alpha and Beta electrons

      implicit none 
      character(len = 512) :: str1, str2, str3
      character,intent(in) :: filename
      integer,intent(in) :: Spin, IUnit
      integer,intent(out) :: NElectrons
      logical :: found = .false.
      integer :: io

1000  format (A61)

!     str1 is the whole line by line for the first 61 columns
!     You can change the format statement above depending on how much
!     you want to read into the line.

!     str2 is the check string against str1. Change to whatever
!     phrase you want. Make sure to index through correct values!

      str2 = ('Number of alpha electrons')
      str3 = ('Number of beta electrons')

!     Do loop through line by line of input file
      do while(.not.found)
        ! Read in each line for str1 and check against str2
        read (IUnit,1000,iostat=io) str1
        ! If yes change logical to true
        ! Change index to match length of string.
        if(str1(1:25) == str2.and.Spin == 1) then 
          found = .true.
          write(*,*) 'Found Alpha Electrons'
        end if

        if(str1(1:24) == str3.and.Spin == 2) then 
          found = .true.
          write(*,*) 'Found Beta Electrons'
        end if
        
        ! Check if End of file. If yes exit 
        if (io.ne.0) then 
          write(*,*) 'Never Found Electrons'
          stop
        end if

      end do
!     Convert end of str1 into integer value for Nelectrons
!     Change index above if your column in the file of the line
!     is different in your program.
      read(str1(59:61),'(I3)') NElectrons

      found =.false.
      write(*,*)
      End Subroutine Number_of_Electrons


      Subroutine Number_of_Basis_Functions(filename1,NBasis,IUnit)

!     This subroutine calls the string in question (in this case number
!     of basis sets) You can use this subroutine and modify it for any 
!     string and read the values at the end of the string.     


      implicit none 
      character(len = 512) :: str1, str2
      character,intent(in) :: filename1
      logical :: found = .false.
      integer,intent(out):: NBasis
      integer:: io, IUnit

1000  format (A61)

!     str1 is the whole line by line for the first 61 columns
!     You can change the format statement above depending on how much
!     you want to read into the line.

!     str2 is the check string against str1. Change to whatever
!     phrase you want. Make sure to index through correct values!

      str2 = ('Number of basis functions')

!     Do loop through line by line of input file
      do while(.not.found)
        ! Read in each line for str1 and check against str2
        read (IUnit,1000,iostat=io) str1
        ! If yes change logical to true
        ! Change index to match length of string.
        if(str1(1:25) == str2) then 
          found = .true.
        end if
        
        ! Check if End of file. If yes exit 
        if (io.ne.0) then 
          write(*,*) 'Never Found # of Basis Functions.'
          stop
        end if

      end do
!     Convert end of str1 into integer value for Nelectrons
!     Change index above if your column in the file of the line
!     is different in your program.
      read(str1(55:61),'(I7)') NBasis

      found =.false.
      
      End Subroutine Number_of_Basis_Functions
      
      Subroutine GET_MO_Array(IUnit,NBasis,MO_Array,Spin)

!     This subroutine does the following:
!     1. Read in string to find MO coefficients for Alpha or Beta
!     2. Read in MO coefficients as an array of NBasis*NBasis length
!     3. Sends MO coefficient matrix back as array

      
      implicit none
      character(len = 512) :: str1, str2, str3
      logical :: found = .false.
      integer, intent(in) :: Spin,NBasis
      real,dimension(NBasis*NBasis), intent(out) :: MO_Array  
      integer:: io, IUnit

1000  format (A61)

      str2 = ('Alpha MO coefficients')
      str3 = ('Beta MO coefficients')

      write(*,*)
!     Resets logical for multiple subroutine runs
      found = .false.

!     Do loop through line by line of input file
      do while(.not.found)
        ! Read in each line for str1 and check against str2
        read (IUnit,1000,iostat=io) str1
        ! If yes change logical to true
        ! Change index to match length of string.
        if(str1(1:21) == str2.and.Spin == 1) then 
          found = .true.
          write(*,*) 'Found Alpha MO Coefficients'
          !Read into Alpha Array
          read (IUnit,*) MO_Array
        end if

        if(str1(1:20) == str3.and.Spin == 2) then 
          found = .true.
          write(*,*) 'Found Beta MO Coefficients'
          !Read into Beta Array
          read (IUnit,*) MO_Array
        end if
        
        ! Check if End of file. If yes exit 
        if (io.ne.0) then 
          write(*,*) 'Never Found MO Coefficients'
          stop
        end if

      end do

      End Subroutine GET_MO_Array


      End Module DDNO_mod
