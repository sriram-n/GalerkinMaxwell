!------------------------------------------------------------------------------
!> Purpose : source term
!!
!! @param[in]  Mdle  - element (middle node) number
!! @param[in]  X     - physical coordinates
!! @param[out] Fval  - rhs
!------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine getf(Mdle,X, Jval)
!
      use control          , only : NEXACT
      use assembly         , only : NR_RHS
      use data_structure3D , only : NR_COMP,ADRES,NRINDEX
      use parameters       , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
      use problem          , only : ZI,MU,EPSILON,SIGMA,OMEGA
!
      implicit none
      integer,                  intent(in)  :: Mdle
      real*8,dimension(3),      intent(in)  :: X
      VTYPE,dimension(3),       intent(out) :: Jval
!------------------------------------------------------------------------------
!
!     exact solution
      VTYPE,dimension(  MAXEQNH    ) ::   valH
      VTYPE,dimension(  MAXEQNH,3  ) ::  dvalH
      VTYPE,dimension(  MAXEQNH,3,3) :: d2valH
      VTYPE,dimension(3,MAXEQNE    ) ::   valE
      VTYPE,dimension(3,MAXEQNE,3  ) ::  dvalE
      VTYPE,dimension(3,MAXEQNE,3,3) :: d2valE
      VTYPE,dimension(3,MAXEQNV    ) ::   valV
      VTYPE,dimension(3,MAXEQNV,3  ) ::  dvalV
      VTYPE,dimension(3,MAXEQNV,3,3) :: d2valV
      VTYPE,dimension(  MAXEQNQ    ) ::   valQ
      VTYPE,dimension(  MAXEQNQ,3  ) ::  dvalQ
      VTYPE,dimension(  MAXEQNQ,3,3) :: d2valQ
!
!     miscellaneus
      integer :: iload,ivar,ibeg,icomp,jcomp,k,l
      complex*16 :: zk2
!
!     printing flag
      integer :: iprint
!------------------------------------------------------------------------------
!
!     printing flag : 0 - silent ; 1 - verbose
      iprint=0
!
      if (iprint.eq.1) then
        write(*,7001) Mdle,X
 7001   format(' getf: Mdle,X = ',i8,2x,3(f8.3,2x))
      endif
!
!     initialize source terms
      Jval = ZERO
!
      select case(NEXACT)
!==============================================================================
!  UNKOWN EXACT SOLUTION                                                      |
!==============================================================================
      case(0)
!!!        Jval(1) = 1.d0
!
!==============================================================================
!  KNOWN EXACT SOLUTION, NON-ZERO RHS
!==============================================================================
      case(1)
!
!     compute exact solution
      call exact(X,Mdle, valH,dvalH,d2valH, valE,dvalE,d2valE, &
                         valV,dvalV,d2valV, valQ,dvalQ,d2valQ)
!
      Jval(1) = d2valE(2,1,1,2) - d2valE(1,1,2,2) - d2valE(1,1,3,3) + d2valE(3,1,1,3)
      Jval(2) = d2valE(3,1,2,3) - d2valE(2,1,3,3) - d2valE(2,1,1,1) + d2valE(1,1,2,1)
      Jval(3) = d2valE(1,1,3,1) - d2valE(3,1,1,1) - d2valE(3,1,2,2) + d2valE(2,1,3,2)
      Jval(1:3) = Jval(1:3)/MU
!
      zk2 = OMEGA**2*EPSILON - ZI*OMEGA*SIGMA
      Jval(1:3) = Jval(1:3) - zk2*valE(1:3,1)
      Jval(1:3) = Jval(1:3)/(-ZI*OMEGA)
      !Jval(1:3) = Jval(1:3)

!
!==============================================================================
!  KNOWN EXACT SOLUTION , HOMOGENEOUS RHS                                     |
!==============================================================================
      case(2)
!
      endselect
!
      if (iprint.eq.1) then
        write(*,7010) Jval
 7010   format(' getf: Jval = ',e12.5)
      endif
!
!
endsubroutine getf


subroutine get_bdSource(Mdle,X,Rn, Imp_val)
!
      use control          , only : NEXACT
      use assembly         , only : NR_RHS
      use data_structure3D , only : NR_COMP,ADRES,NRINDEX
      use parameters       , only : MAXEQNH,MAXEQNE,MAXEQNV,MAXEQNQ,ZERO
      use problem
!
      implicit none
      integer,                  intent(in)  :: Mdle
      real*8,dimension(3),      intent(in)  :: X
      real*8,dimension(3),      intent(in)  :: Rn
      VTYPE,dimension(3),       intent(out) :: Imp_val
!------------------------------------------------------------------------------
!
!     exact solution
      VTYPE,dimension(  MAXEQNH    ) ::   zvalH
      VTYPE,dimension(  MAXEQNH,3  ) ::  zdvalH
      VTYPE,dimension(  MAXEQNH,3,3) :: zd2valH
      VTYPE,dimension(3,MAXEQNE    ) ::   zvalE
      VTYPE,dimension(3,MAXEQNE,3  ) ::  zdvalE
      VTYPE,dimension(3,MAXEQNE,3,3) :: zd2valE
      VTYPE,dimension(3,MAXEQNV    ) ::   zvalV
      VTYPE,dimension(3,MAXEQNV,3  ) ::  zdvalV
      VTYPE,dimension(3,MAXEQNV,3,3) :: zd2valV
      VTYPE,dimension(  MAXEQNQ    ) ::   zvalQ
      VTYPE,dimension(  MAXEQNQ,3  ) ::  zdvalQ
      VTYPE,dimension(  MAXEQNQ,3,3) :: zd2valQ
!
!     miscellaneus
      integer :: iload,ivar,ibeg,icomp,jcomp,k,l
      complex*16 :: zaux
      VTYPE,dimension(3) ::   rntimesE,rn2timesE
      VTYPE,dimension(3) ::   curlE,rntimescurlE
      real*8                    :: impedanceConstant
!      real*8                    :: E   ! vector field
!      real*8, dimension(3)      :: dE  ! 1st derivative
!      real*8, dimension(3,3)    :: d2E ! 2nd derivative
!
!     printing flag
      integer :: iprint
!------------------------------------------------------------------------------
!
!     printing flag : 0 - silent ; 1 - verbose
      iprint=0
!
      if (iprint.eq.1) then
        write(*,7001) Mdle,X
 7001   format(' get_bdSource: Mdle,X = ',i8,2x,3(f8.3,2x))
      endif
!
!     initialize source terms
      Imp_val = ZERO
!
      select case(NEXACT)
!  ... known HOMOGENEOUS solution: do nothing
      case(0,2)

      case(1)
!
!  .....impedance BC
          call exact(X,Mdle, zvalH,zdvalH,zd2valH, &
                        zvalE,zdvalE,zd2valE, &
                        zvalV,zdvalV,zd2valV, &
                        zvalQ,zdvalQ,zd2valQ)

          call my_cross_product(Rn,zvalE(1:3,1), rntimesE)
          !  ...value
          curlE(1) = zDvalE(3,1,2) - zDvalE(2,1,3)
          curlE(2) = zDvalE(1,1,3) - zDvalE(3,1,1)
          curlE(3) = zDvalE(2,1,1) - zDvalE(1,1,2)

          call my_cross_product(Rn,rntimesE, rn2timesE)
          call my_cross_product(Rn,curlE, rntimescurlE)
!
          Imp_val = rntimescurlE/(-ZI*OMEGA*MU) -(GAMMA)*rn2timesE
!
!  ...exact solution unknown
      case default
        write(*,*) 'get_bdSource: UNSPECIFIED NEXACT';stop
      end select
      if (iprint.eq.1) then
        write(*,7002) Imp_val
 7002   format('get_bsource: Imp_val = ',2e12.5)
      endif
!
end subroutine get_bdSource

subroutine my_cross_product(A,Zb,Zc)
!
      use problem
!
      implicit none
      real*8,dimension(3),      intent(in)  :: A
      VTYPE,dimension(3),       intent(in)  :: Zb
      VTYPE,dimension(3),       intent(out) :: Zc
!
      Zc(1) =   A(2)*Zb(3) - A(3)*Zb(2)
      Zc(2) = - A(1)*Zb(3) + A(3)*Zb(1)
      Zc(3) =   A(1)*Zb(2) - A(2)*Zb(1)
!
end subroutine my_cross_product
