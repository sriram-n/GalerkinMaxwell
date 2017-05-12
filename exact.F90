!------------------------------------------------------------------------------
!> Purpose : exact (manufactured, no miracle!) solution
!!
!> @param[in]  X      - a point in physical space
!> @param[in]  Mdle   - element (middle node) number
!> @param[out] ValH   - value of the H1 solution
!> @param[out] DvalH  - corresponding first derivatives
!> @param[out] D2valH - corresponding second derivatives
!> @param[out] ValE   - value of the H(curl) solution
!> @param[out] DvalE  - corresponding first derivatives
!> @param[out] D2valE - corresponding second derivatives
!> @param[out] ValV   - value of the H(div) solution
!> @param[out] DvalV  - corresponding first derivatives
!> @param[out] D2valV - corresponding second derivatives
!> @param[out] ValQ   - value of the H(div) solution
!> @param[out] DvalQ  - corresponding first derivatives
!> @param[out] D2valQ - corresponding second derivatives
!------------------------------------------------------------------------------
!
#include "implicit_none.h"
!
subroutine exact(Xp,Mdle, ValH,DvalH,D2valH, ValE,DvalE,D2valE, &
                          ValV,DvalV,D2valV, ValQ,DvalQ,D2valQ)
!
      use data_structure3D
      use problem
!
      implicit none
      real*8,dimension(3),            intent(in)  :: Xp
      integer                       , intent(in)  :: Mdle
      VTYPE,dimension(  MAXEQNH    ), intent(out) ::   ValH
      VTYPE,dimension(  MAXEQNH,3  ), intent(out) ::  DvalH
      VTYPE,dimension(  MAXEQNH,3,3), intent(out) :: D2valH
      VTYPE,dimension(3,MAXEQNE    ), intent(out) ::   ValE
      VTYPE,dimension(3,MAXEQNE,3  ), intent(out) ::  DvalE
      VTYPE,dimension(3,MAXEQNE,3,3), intent(out) :: D2valE
      VTYPE,dimension(3,MAXEQNV    ), intent(out) ::   ValV
      VTYPE,dimension(3,MAXEQNV,3  ), intent(out) ::  DvalV
      VTYPE,dimension(3,MAXEQNV,3,3), intent(out) :: D2valV
      VTYPE,dimension(  MAXEQNQ    ), intent(out) ::   ValQ
      VTYPE,dimension(  MAXEQNQ,3  ), intent(out) ::  DvalQ
      VTYPE,dimension(  MAXEQNQ,3,3), intent(out) :: D2valQ
!
!------------------------------------------------------------------------------
!     Space for temporary solutions
!

      VTYPE                    :: E
      VTYPE,dimension(3)       :: dE
      VTYPE, dimension(3,3)    :: d2E
      integer                  :: icomp

!
!     initialize exact solution
      ValH=ZERO ; DvalH=ZERO ; D2valH=ZERO
      ValE=ZERO ; DvalE=ZERO ; D2valE=ZERO
      ValV=ZERO ; DvalV=ZERO ; D2valV=ZERO
      ValQ=ZERO ; DvalQ=ZERO ; D2valQ=ZERO

      icomp = ICOMP_EXACT

!     initialize variables
      E = ZERO
      dE = ZERO
      d2E = ZERO

      call hcurl_solution(Xp, E,dE,d2E)


 !  .....value
        ValE(icomp,1    ) =   E
!
!  .....1st order derivatives
        DvalE(icomp,1,1  ) =  dE(1)
        DvalE(icomp,1,2  ) =  dE(2)
        DvalE(icomp,1,3  ) =  dE(3)
!
!  .....2nd order derivatives
        D2valE(icomp,1,1,1) = d2E(1,1)
        D2valE(icomp,1,1,2) = d2E(1,2)
        D2valE(icomp,1,1,3) = d2E(1,3)
        D2valE(icomp,1,2,1) = d2E(2,1)
        D2valE(icomp,1,2,2) = d2E(2,2)
        D2valE(icomp,1,2,3) = d2E(2,3)
        D2valE(icomp,1,3,1) = d2E(3,1)
        D2valE(icomp,1,3,2) = d2E(3,2)
        D2valE(icomp,1,3,3) = d2E(3,3)


  end subroutine exact


    subroutine hcurl_solution(Xp, E,dE,d2E)
    use data_structure3D
    use problem
#if C_MODE
#define V_TYPE  complex*16
#else
#define V_TYPE double precision
#endif
    implicit none
    !-----------------------------------------------------------------------------------
    real*8, dimension(3),     intent(in)  :: Xp
!
    V_TYPE,                   intent(out) :: E   ! vector field
    V_TYPE, dimension(3),     intent(out) :: dE  ! 1st derivative
    V_TYPE, dimension(3,3),   intent(out) :: d2E ! 2nd derivative
    ! !-----------------------------------------------------------------------------------

    real*8  :: x1,x2,x3
    real*8 :: nn        ! for EXPONENTIAL solution
    real*8 :: cn,dn     ! for singular solution
    real*8 :: np_x,np_y,np_z,r0,k0,w0,phase,amplitude
    real*8 :: impedanceConstant,om
    V_TYPE :: f_x,f_y,f_z,df_x,df_y,df_z,ddf_x,ddf_y,ddf_z
    integer :: icomp
    V_TYPE  :: c2z,uz,uz_x,uz_y,uz_z,uz_xx,uz_xy,uz_xz,uz_yy,uz_yx,uz_yz
    V_TYPE  :: uz_zz,uz_zy,uz_zx
    V_TYPE  :: pz,pz_x,pz_y,pz_z,pz_xx,pz_xy,pz_xz,pz_yy,pz_yx,pz_yz
    V_TYPE  :: pz_zz,pz_zy,pz_zx
    ! !-----------------------------------------------------------------------------------

    ! initialize variables
    E = ZERO; dE = ZERO; d2E = ZERO;
    icomp = ICOMP_EXACT
    f_x = ZERO; f_y = ZERO; f_z = ZERO
    df_x = ZERO; df_y = ZERO; df_z = ZERO
    ddf_x = ZERO; ddf_y = ZERO; ddf_z = ZERO


    ! !
    ! !-----------------------------------------------------------------------------------
    ! !      D E C L A R E    S O L U T I O N    V A R I A B L E S                       |
    ! !-----------------------------------------------------------------------------------
    ! !
    x1 = Xp(1); x2 = Xp(2); x3 = Xp(3)


    !--------------- 1st prob -------------------------------------------------------
    if (ISOL .eq. 1) then

         np_x=real(NPX,8); np_y=real(NPY,8); np_z=real(NPZ,8)
!
!       value
        f_x = x1**np_x
        f_y = x2**np_y
        f_z = x3**np_z
!
!       derivatives
        select case(int(np_x))
          case(0); df_x = 0.d0; ddf_x = 0.d0
          case(1); df_x = 1.d0; ddf_x = 0.d0
          case default
            df_x = np_x * x1**(np_x-1.d0)
            ddf_x = np_x * (np_x-1.d0) * x1**(np_x-2.d0)
        end select
        select case(int(np_y))
          case(0); df_y = 0.d0; ddf_y = 0.d0
          case(1); df_y = 1.d0; ddf_y = 0.d0
          case default
            df_y = np_y * x2**(np_y-1.d0)
            ddf_y = np_y * (np_y-1.d0) * x2**(np_y-2.d0)
        end select
        select case(int(np_z))
          case(0); df_z = 0.d0; ddf_z = 0.d0
          case(1); df_z = 1.d0; ddf_z = 0.d0
          case default
            df_z = np_z * x3**(np_z-1.d0)
            ddf_z = np_z * (np_z-1.d0) * x3**(np_z-2.d0)
        end select



!----------- 2nd prob -------------------------------------------------------
!  ...a smooth solution
      elseif (ISOL .eq. 2) then
      om = OMEGA

      f_x= sin(om*Xp(1))!dexp(-Xp(1)**2)/20.d0!sin(OMEGA*Xp(1))
      f_y= sin(om*Xp(2))!dexp(-Xp(2)**2)/20.d0!sin(OMEGA*Xp(2))
      f_z= sin(om*Xp(3))
!
!     1st order derivatives
      df_x=(om)*cos(om*Xp(1))!-Xp(1)*dexp(-Xp(1)**2)/10.d0!(OMEGA)*cos(OMEGA*Xp(1))
      df_y=(om)*cos(om*Xp(2))!-Xp(2)*dexp(-Xp(2)**2)/10.d0!(OMEGA)*cos(OMEGA*Xp(2))
      df_z=(om)*cos(om*Xp(3))
!
!     2nd order derivatives
      ddf_x=-om**2*f_x!dexp(-Xp(1)**2)*(2.d0*Xp(1)**2-1.d0)/10.d0!-OMEGA**2*f_x
      ddf_y=-om**2*f_y!dexp(-Xp(2)**2)*(2.d0*Xp(2)**2-1.d0)/10.d0!-OMEGA**2*f_y
      ddf_z=-om**2*f_z

        !     an exponential
      elseif (ISOL .eq. 3) then
!
!       value
        f_x=exp(Xp(1))
        f_y=1.d0!exp(Xp(2))
        f_z=1.d0!exp(Xp(3))
!
!       1st order derivatives
        df_x=exp(Xp(1))
        df_y=0.d0!exp(Xp(2))
        df_z=0.d0!exp(Xp(3))
!
!       2nd order derivatives
        ddf_x=exp(Xp(1))
        ddf_y=0.d0!exp(Xp(2))
        ddf_z=0.d0!exp(Xp(3))
!

!
      elseif (ISOL .eq. 4) then
!
!     value
      f_x=1.d0
      f_y=1.d0
      f_z=Xp(3)**2 * (1.d0-Xp(3))
!
!     1st order derivatives
      df_x=0.d0
      df_y=0.d0
      df_z=2.d0*Xp(3)-3.d0*Xp(3)**2
!
!     2nd order derivatives
      ddf_x=0.d0
      ddf_y=0.d0
      ddf_z=2.d0-6.d0*Xp(3)
!
      elseif (ISOL .eq. 5) then
!
!     value
      f_x=1.d0
      f_y=1.d0
      f_z=Xp(3)**2
!
!     1st order derivatives
      df_x=0.d0
      df_y=0.d0
      df_z=2.d0*Xp(3)
!
!     2nd order derivatives
      ddf_x=0.d0
      ddf_y=0.d0
      ddf_z=2.

      !  ...fundamental TE10 mode for rectangular waveguide
      elseif (ISOL .eq. 9) then

      f_x=-ZI*(OMEGA/PI)*sin(PI*Xp(1))
      f_y= 1.d0
      f_z=cdexp(-ZI*OMEGA*Xp(3)*GAMMA)
!
!     1st order derivatives
      df_x=-ZI*(OMEGA/PI)*PI*cos(PI*Xp(1))
      df_y=0.d0
      df_z=(-ZI*OMEGA*GAMMA)*f_z
!
!     2nd order derivatives
      ddf_x=-PI**2*f_x
      ddf_y=0.d0
      ddf_z=(-ZI*OMEGA*GAMMA)*df_z

      !----------- prob -------------------------------------------------------
!  ...a smooth solution
      elseif (ISOL .eq. 50) then

      f_x= Xp(1)*(1.d0-Xp(1))!dexp(-Xp(1)**2)/20.d0!sin(OMEGA*Xp(1))
      f_y= Xp(2)*(1.d0-Xp(2))!dexp(-Xp(2)**2)/20.d0!sin(OMEGA*Xp(2))
      f_z= Xp(3)*(128.d0-Xp(3))!/1000.0d0
!
!     1st order derivatives
      df_x= 1.d0-2.d0*Xp(1)!-Xp(1)*dexp(-Xp(1)**2)/10.d0!(OMEGA)*cos(OMEGA*Xp(1))
      df_y=1.d0-2.d0*Xp(2)!-Xp(2)*dexp(-Xp(2)**2)/10.d0!(OMEGA)*cos(OMEGA*Xp(2))
      df_z=(128.d0-2.d0*Xp(3))!/1000.0d0
!
!     2nd order derivatives
      ddf_x=-2.d0!dexp(-Xp(1)**2)*(2.d0*Xp(1)**2-1.d0)/10.d0!-OMEGA**2*f_x
      ddf_y=-2.d0!dexp(-Xp(2)**2)*(2.d0*Xp(2)**2-1.d0)/10.d0!-OMEGA**2*f_y
      ddf_z=-2.d0!/1000.0d0

        !     an exponential
    endif

!     !  .....1st order derivatives
         dE(1)=  df_x *   f_y *   f_z
         dE(2) =   f_x *  df_y *   f_z
         dE(3)=   f_x *   f_y *  df_z
! !
! !  .....2nd order derivatives
         d2E(1,1) = ddf_x *   f_y *   f_z
         d2E(1,2) =  df_x *  df_y *   f_z
         d2E(1,3) =  df_x *   f_y *  df_z
         d2E(2,1) =  d2E(1,2)
         d2E(2,2) =   f_x * ddf_y *   f_z
         d2E(2,3) =   f_x *  df_y *  df_z
         d2E(3,1) =  d2E(1,3)
         d2E(3,2) =  d2E(2,3)
         d2E(3,3) =   f_x *   f_y * ddf_z

         E=f_x*f_y*f_z
     end subroutine hcurl_solution
