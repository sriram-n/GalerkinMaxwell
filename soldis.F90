!---------------------------------------------------------------------------------------
!> Purpose : display user-defined quantity
!! @param[in] Mdle   - element (middle node) number 
!! @param[in] Xi     - master element coordinates
!! @param[in] X      - physical coordinates
!! @param[in] Rn     - outward normal unit vector
!! @param[in] ZsolH  - H1    sol
!! @param[in] ZgradH - H1    grad
!! @param[in] ZsolE  - Hcurl sol
!! @param[in] ZcurlE - Hcurl curl
!! @param[in] ZsolV  - Hdiv  sol
!! @param[in] ZdivV  - Hdiv  div
!! @param[in] ZsolQ  - L2    sol
!! 
!! @param[out] val   - quantity to display
!---------------------------------------------------------------------------------------
!
      subroutine soldis(Mdle,Xi,X,Rn,ZsolH,ZgradH,ZsolE,ZcurlE,ZsolV,ZdivV,ZsolQ, Val)
!
      use data_structure3D
      use problem
!
#include "implicit_none.h"
      integer,                     intent(in)  :: Mdle
      real*8,dimension(3),         intent(in)  :: Xi,X,Rn
      VTYPE,dimension(  MAXEQNH  ),intent(in)  :: ZsolH
      VTYPE,dimension(  MAXEQNH,3),intent(in)  :: ZgradH
      VTYPE,dimension(3,MAXEQNE  ),intent(in)  :: ZsolE
      VTYPE,dimension(3,MAXEQNE  ),intent(in)  :: ZcurlE
      VTYPE,dimension(3,MAXEQNV  ),intent(in)  :: ZsolV
      VTYPE,dimension(  MAXEQNV  ),intent(in)  :: ZdivV
      VTYPE,dimension(  MAXEQNQ  ),intent(in)  :: ZsolQ
      real*8,                      intent(out) :: Val
!
      VTYPE :: &
         zvalH(MAXEQNH), &
         zdvalH(MAXEQNH,3),zd2valH(MAXEQNH,3,3), &
         zvalE(3,MAXEQNE), &
         zdvalE(3,MAXEQNE,3),zd2valE(3,MAXEQNE,3,3), &
         zvalV(3,MAXEQNV), &
         zdvalV(3,MAXEQNV,3),zd2valV(3,MAXEQNV,3,3), &
         zvalQ(MAXEQNQ), &
         zdvalQ(MAXEQNQ,3),zd2valQ(MAXEQNQ,3,3)
      integer :: icase = 1, iprint
      VTYPE :: zval
!---------------------------------------------------------------------------------------
!
      iprint=10
!
      select case (IEXACT_DISP)
!
!  ...compute the exact solution if needed
      case(1,3)
        call exact(X, icase, &
                   zvalH,zdvalH,zd2valH, &
                   zvalE,zdvalE,zd2valE, &
                   zvalV,zdvalV,zd2valV, &
                   zvalQ,zdvalQ,zd2valQ)
      end select
!
      select case (IEXACT_DISP)
!
!  ...display exact solution
      case(1)
        zval = ZvalE(iabs(ICHOOSE_DISP),1)
!
!  ...display approximate solution
      case(2)
        zval = ZsolE(iabs(ICHOOSE_DISP),1)
!
!  ...display error
      case(3)
        zval = ZvalE(iabs(ICHOOSE_DISP),1) - ZsolE(iabs(ICHOOSE_DISP),1)
      end select
!
!  ...choose between real and imaginary parts
      if (ICHOOSE_DISP.gt.0) then
        Val = dreal(zval)
      else
        Val = dimag(zval)
      endif
      if (iprint.eq.1) then
        write(*,1010) IEXACT_DISP,ICHOOSE_DISP,Val
 1010   format('soldis: IEXACT_DISP,ICHOOSE_DISP = ',2i3,' Val = ',e12.5)
        call pause
      endif
!  
      end subroutine soldis
!
!
!
!---------------------------------------------------------------------------------------
!> Purpose : show the quantities to display 
!---------------------------------------------------------------------------------------
      subroutine soldis_select
      use problem
      use control,  only : NEXACT

!
      implicit none
      integer, save :: ivis=0 ! visitation flag
      integer :: iprev
!---------------------------------------------------------------------------------------
!
      if (ivis.eq.1) then
   10   write(*,*) 'USE PREVIOUS SELECTION? 1-Y,0-N'
        read(*,*) iprev
        select case(iprev)
        case(0)
        case(1)
          call disp_soldis(NSTD_OUT) ; return
        case default
          go to 10
        end select
      else
        ivis=1
      endif
!
      if ((NEXACT.eq.1).or.(NEXACT.eq.2)) then
   20   write(*,*) 'SET VARIABLE     1) EXACT 2) APPROX 3) ERROR'
        read(*,*) IEXACT_DISP
        select case(IEXACT_DISP)
        case(1,2,3)
        case default
          go to 20
        end select
      else
        IEXACT_DISP=2
      endif  
!
   30 write(*,*) 'SET VARIABLE: 1) E1, 2) E2, 3) E3'
      write(*,*) 'USE NEGATIVE FLAGS FOR IMAGINARY PARTS'
      read(*,*) ICHOOSE_DISP
      select case(ICHOOSE_DISP)
      case(1,2,3)
      case(-1,-2,-3)
      case default
        go to 30
      end select
!
      call disp_soldis(NSTD_OUT)
!  
      end subroutine soldis_select
!
      subroutine disp_soldis(Nstream)
      use problem
      implicit none
      integer, intent(in) :: Nstream
!
      write(Nstream,1005)
      write(Nstream,1000)
      select case (IEXACT_DISP)
      case(1); write(Nstream,1010)
      case(2); write(Nstream,1020)
      case(3); write(Nstream,1030)
      end select
      if (ICHOOSE_DISP.gt.0) write(Nstream,1040)  ICHOOSE_DISP
      if (ICHOOSE_DISP.lt.0) write(Nstream,1050) -ICHOOSE_DISP
      write(Nstream,1000)
!
 1000 format('-----------------------')
 1005 format('DISPLAY SETUP')
 1010 format('EXACT  SOLUTION SELECTED')
 1020 format('APPROX SOLUTION SELECTED')
 1030 format('ERROR SELECTED')
 1040 format('REAL PART OF E',i1,' SELECTED')
 1050 format('IMAG PART OF E',i1,' SELECTED')
!
      end subroutine disp_soldis
