!----------------------------------------------------------------------
!
!   module name        - matrices
!
!----------------------------------------------------------------------
!
!   latest revision    - May 2
!
!   purpose            - store element matrices for the first layer
!                        elements
!
!----------------------------------------------------------------------
!
module matrices
!
  use parametersDPG
  implicit none
#if C_MODE
#define VTYPE  complex*16
#else
#define VTYPE double precision
#endif
!
! max # of elements in the first layer
  integer, parameter :: MAXNRFL=8
!
! # of stored elements in the first layer
  integer :: NRFL
!!$OMP THREADPRIVATE (NRFL)
!
! order of elements
  integer, parameter :: MYP=6
!
! matrix dimensions
  integer, parameter :: MYE = 3*MYP*(MYP+1)**2*2
!
! stiffness matrices to store
  VTYPE, dimension(MYE,MYE,MAXNRFL) :: ZFL_EE
!!$OMP THREADPRIVATE (ZFL_EE)
  VTYPE, dimension(MYE,MAXNRFL) :: ZLOADFL_E
!
! xy coordinates of the first vertex node
  real*8, dimension(2,MAXNRFL) :: XYVERT
!!$OMP THREADPRIVATE (XYVERT)
!
end module matrices



