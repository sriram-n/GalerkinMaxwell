c--------------------------------------------------------------------
c                                                                     
c     routine name      - my_own_tests
c                                                                     
c-------------------------------------------------------------------- 
c                                                                     
c     latest revision:  - Apr 14
c                                                                     
c     purpose:          - driver for quick testing routines
c                                                                    
c-----------------------------------------------------------------------
c
      subroutine my_own_tests
c
      use control
      use data_structure3D
      use assembly
#include "syscom.blk"
c
c  ...work space for routine exact
      real*8,dimension(3)            :: xp
      integer                        :: mdle
      real*8,dimension(  MAXEQNH    )::   valH
      real*8,dimension(  MAXEQNH,3  )::  dvalH
      real*8,dimension(  MAXEQNH,3,3):: d2valH
      real*8,dimension(3,MAXEQNE    )::   valE
      real*8,dimension(3,MAXEQNE,3  )::  dvalE
      real*8,dimension(3,MAXEQNE,3,3):: d2valE
      real*8,dimension(3,MAXEQNV    )::   valV
      real*8,dimension(3,MAXEQNV,3  )::  dvalV
      real*8,dimension(3,MAXEQNV,3,3):: d2valV
      real*8,dimension(  MAXEQNQ    )::   valQ
      real*8,dimension(  MAXEQNQ,3  )::  dvalQ
      real*8,dimension(  MAXEQNQ,3,3):: d2valQ
c
 10   write(*,*) 'SELECT'
      write(*,*) '  '
      write(*,*) 'Exit...........................................0'
      write(*,*) 'Test exact solution routine....................1'
c
      read(*,*) iselect
      select case(iselect)
      case(0)
        return
      case(1)
        write(*,*) 'SET mdle,xp'
        read(*,*) mdle,xp
        write(*,*) 'mdle,xp = ',mdle,xp(1:3)
        call  exact(xp,mdle, valH,dvalH,d2valH, valE,dvalE,d2valE, 
     .                       valV,dvalV,d2valV, valQ,dvalQ,d2valQ)
        write(*,7010) valH(1)
 7010   format('valH(1) = ',e12.5)
      end select
      go to 10
c
      end subroutine my_own_tests
