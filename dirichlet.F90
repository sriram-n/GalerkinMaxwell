!----------------------------------------------------------------------------------
!> Purpose : calculate dirichlet boundary condition
!!
!! @param[in]  Mdle  - an element (middle node) number
!! @param[in]  X     - physical coordinates of a point
!! @param[in]  Icase - node case
!!
!! @param[out] ValH  - value of the H1 solution
!! @param[out] DvalH - H1 corresponding first derivatives
!! @param[out] ValE  - value of the H(curl) solution
!! @param[out] DvalE - H(curl) corresponding first derivatives
!! @param[out] ValV  - value of the H(div) solution
!! @param[out] DvalV - H(div) corresponding first derivatives
!----------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine dirichlet(Mdle,X,Icase, ValH,DvalH, ValE,DvalE, ValV,DvalV)
!
      use control    , only : NEXACT, GEOM_TOL
      use parameters , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ, ZERO
!
      implicit none
      integer,                       intent(in)  :: Mdle
      real*8, dimension(3),          intent(in)  :: X
      integer,                       intent(in)  :: Icase
      VTYPE,dimension(  MAXEQNH    ),intent(out) ::   ValH
      VTYPE,dimension(  MAXEQNH,3  ),intent(out) ::  DvalH
      VTYPE,dimension(  MAXEQNH,3,3)             :: d2valH
      VTYPE,dimension(3,MAXEQNE    ),intent(out) ::   ValE
      VTYPE,dimension(3,MAXEQNE,3  ),intent(out) ::  DvalE
      VTYPE,dimension(3,MAXEQNE,3,3)             :: d2valE
      VTYPE,dimension(3,MAXEQNV    ),intent(out) ::   ValV
      VTYPE,dimension(3,MAXEQNV,3  ),intent(out) ::  DvalV
      VTYPE,dimension(3,MAXEQNV,3,3)             :: d2valV
      VTYPE,dimension(  MAXEQNQ    )             ::   valQ
      VTYPE,dimension(  MAXEQNQ,3  )             ::  dvalQ
      VTYPE,dimension(  MAXEQNQ,3,3)             :: d2valQ
!
!     printing flag
      integer :: iprint
!----------------------------------------------------------------------------------
!
!     printing flag : 0 - silent ; 1 - verbose
      iprint=0
!
!     initialize
      ValH=ZERO ; DvalH=ZERO ; ValE=ZERO ; DvalE=ZERO ; ValV=ZERO ; DvalV=ZERO
!
      select case(NEXACT)
!     exact solution UNKNOWN
      case(0)
      stop 1
!
!     exact solution KNOWN
      case(1,2)
!       use the exact solution to determine Dirichlet data
        call exact(X,Mdle, ValH,DvalH,d2valH, ValE,DvalE,d2valE,  &
                           ValV,DvalV,d2valV, valQ,dvalQ,d2valQ)
!
      case default
        write(*,1000) NEXACT
 1000   format(' dirichlet: NEXACT = ',i2)
        stop
      endselect
!
      if (iprint == 1) then
        write(*,1001) X(1:3),ValE(1:3,1)
 1001   format(' dirichlet: X    = ',3(f8.3,2x),/, &
               '            ValE = ',2(3(2e12.5,2x),4x))
      endif
!
!
endsubroutine dirichlet
