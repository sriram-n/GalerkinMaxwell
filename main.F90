!---------------------------------------------------------
!  The code is run through a bash script
!---------------------------------------------------------
program main
!
      use environment
      use paraview
      use control
      use parametersDPG
      use GMP
      use data_structure3D
      use physics
      use uhm2
      use problem
      use matrices

!
      implicit none
      character(len=1024) :: argv
      real*8  :: err,rnorm,rvoid,factor,t
      integer :: mdle,i,kref,idec,nvoid,niter,ibeg,iend,nreflag,istep,nrdof_old,nrdof_new,nstop,idec_solve
      integer :: iso_ans, iii,info,info1,ref_xyz,numRef,iref
      integer, dimension(1) :: flag
!
!----------------------------------------------------------------------
!
!  ...initialize environment
      call begin_environment
!
!  ...read in HP3D input files location (if option is not present, the default value is used)
!
!                             option label      // explanation                // default value     // parameter
      call get_option_string( '-file-control'    , 'Control file'              , './files/control'  , FILE_CONTROL)
      call get_option_string( '-file-geometry'   , 'Geometry file'             , './files/cube_waveguide512', FILE_GEOM   )
      call get_option_string( '-file-phys'       , 'Physics file'              , './files/physics'  , FILE_PHYS   )
      call get_option_string( '-file-refinement' , 'Refinement files location' , '../../files/ref'  , FILE_REFINE )
      call get_option_string( '-file-history'    , 'History file'              , './files/history'  , FILE_HISTORY)
      call get_option_string( '-file-err'        , 'Error file'                , './files/dump_err' , FILE_ERR    )
!
!  ...read in problem dependent parameters (defined in module parametersDPG,DPGH1)
!
!                              option label      // explanation                // default value     // parameter
      call get_option_int(    '-order-approx'       , 'ORDER_APPROX'              , 3                  , ORDER_APPROX)
      call get_option_int(    '-orderx'             , 'NPX'                       , 2                  , NPX         )
      call get_option_int(    '-ordery'             , 'NPY'                       , 1                  , NPY         )
      call get_option_int(    '-orderz'             , 'NPZ'                       , 0                  , NPZ         )
      call get_option_int(    '-comp'               , 'ICOMP_EXACT'               , 2                  , ICOMP_EXACT )
      call get_option_int(    '-isol'               , 'ISOL'                      , 9                  , ISOL        )
!
      call get_option_real(   '-mu'                 , 'MU'                        , 1.d0               , MU          )
      call get_option_real(   '-epsilon'            , 'EPSILON'                   , 1.d0               , EPSILON     )
      call get_option_real(   '-sigma'              , 'SIGMA'                     , 0.0d0               , SIGMA       )
      PI = dacos(-1.d0)
      call get_option_real(   '-omega'              , 'OMEGA'                     , 1.5d0*PI           , OMEGA       )
      GAMMA = sqrt(1.d0-(PI**2)/(OMEGA**2))
      !call get_option_real(   '-omega'              , 'OMEGA'                     , 1.0d0           , OMEGA       )
      !GAMMA = 1.d0
!      call get_option_real(   '-gamma'              , 'GAMMA'                     , 1.d0               , GAMMA       )
!     -- Parview Interface --
      ! Variables relevant to src/modules/paraview
!                        option label     // explanation                        // default value     // parameter
      call get_option_string( '-prefix'          , 'Prefix paraview file'      ,'galerkinMaxwell'            , PREFIX      )
      call get_option_string('-file-vis-upscale','Visualization upscale file location','../../files/vis',FILE_VIS          )
      call get_option_string('-vis-level'       ,'Visualization upscale level (0-3)'  ,'2'                 ,VLEVEL            )
      call get_option_string('-dir-paraview'    ,'Paraview root directory'            ,'./output/figures/test512_p4/' ,PARAVIEW_DIR      )
      call get_option_bool(  '-paraview-geom'   ,'Dump geom at every Paraview call'   ,.TRUE.              ,PARAVIEW_DUMP_GEOM)
      call get_option_bool(  '-paraview-attr'   ,'Dump solution to Paraview'          ,.TRUE.              ,PARAVIEW_DUMP_ATTR)
!
!  ...finalize
      call end_environment
      !call set_problem
!
!  ...print fancy header
      write(*,*)'                      '
      write(*,*)'// --  GALERKIN  FOR MAXWELL EQUATION -- //'
      write(*,*)'                      '
!
!  ...initialize problem
      call initialize

!     Kyungjoo's magic...
      UHM_VERBOSE            = .FALSE.
      UHM_SOLVER_REPORT_SHOW = .FALSE.
!
      call get_command(argv)
      call uhm_initialize(argv)
!
      call uhm_option_begin
      call uhm_option_end
!
      call uhm_time_begin(UHM_TIME_LEVEL_FRONT)
      call uhm_direct_lu_piv_initialize( &
!              UHM_DOUBLE, NR_RHS_PROB, 256, UHM_SOLVER_PTR)
              UHM_DOUBLE, MY_NR_RHS, 256, UHM_SOLVER_PTR)
!
!  ...display menu in infinite loop
 10   continue
      write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
      write(*,*) 'Quit ...................................0'
      write(*,*) '                                         '
      write(*,*) 'Geometry graphics (X11) ................1'
      write(*,*) 'HP3D graphics (X11) ....................3'
      write(*,*) 'Paraview ...............................4'
      write(*,*) 'Print Data Structure arrays ............5'
      write(*,*) 'Dumpout Data Structure .................7'
      write(*,*) '                                         '
      write(*,*) ' --  Geometry & Refinements  --          '
      write(*,*) 'Geometry error ........................12'
      write(*,*) 'Interactive H-refinements .............31'
      write(*,*) 'Uniform     H-refinements .............32'
      write(*,*) 'anisotropic H-refinements .............33'
      write(*,*) '                                         '
      write(*,*) 'Solve (frontal) .......................40'
      write(*,*) 'Solve (MUMPS) .........................50'
      write(*,*) 'Solve (UHM) ...........................60'
      write(*,*) '                                         '
      write(*,*) 'Compute Hcurl error ..................100'
      write(*,*) 'Compute Hcurl error rates ............110'
      write(*,*) '2 uniform n anisotropic refinement....120'
      write(*,*) 'Compute BC data interpolation error...130'
      write(*,*) '                                         '
      write(*,*) 'My tests..... ........................200'
      write(*,*) 'My own tests..... ....................210'
      write(*,*) '=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-='
!
      read(*,*) idec
!
!----------------------------------------------------------------------
!
      select case(idec)
!  ...quit
      case( 0) ; call finalize ; stop
!
!  ...GMP graphics
      case( 1) ; call graphg
!
!  ...hp3d graphics
      case( 3) ; call graphb
!
!  Paraview graphics
        case(4) ; call paraview_driver(iParAttr)
!
!  ...print data structure
      case( 5) ; call result
!
!  ...dump out data structure
      case( 7)
        write(*,*)'dumping out GMP...'
        call dumpout_GMP
        write(*,*)'dumping out physics...'
        call dumpout_physics
        write(*,*)'dumping out HP3D...'
        call dumpout_hp3d('./files/dumpc3Dhp')
!
!  ...geometry error
      case(12) ; call geom_error(err,rnorm)
!
!  ...interactive refinements
      case(31)
        write(*,*)'Active elements:'
        mdle=0
        do i=1,NRELES
          call nelcon(mdle, mdle)
!
          select case(NODES(mdle)%type)
          case('mdlp') ; write(*,7000) mdle
          case('mdlb') ; write(*,7001) mdle
          case('mdln') ; write(*,7002) mdle
          case('mdld') ; write(*,7003) mdle
          endselect
 7000     format(' mdle = ',i6,' ; PRISM')
 7001     format(' mdle = ',i6,' ; BRICK')
 7002     format(' mdle = ',i6,' ; TET')
 7003     format(' mdle = ',i6,' ; PYRAMID')
!
        enddo
        call display_ref_kind
        write(*,7010)
 7010   format(' mdle,kref =')
        read(*,*) mdle,kref
!
!       refine element
        call refine(mdle,kref)
!       recover 1-irregular mesh, update geometry and Dirichlet dof's
        call close_mesh ; call update_gdof ; call update_ddof
!
!     uniform global H-refinements
      case(32)
!       refine elements
        call global_href
!       recover 1-irregular mesh, update geometry and Dirichlet dof's
        call close_mesh ; call update_gdof ; call update_Ddof

!  ..... anisotropic h-refinement of cylinder
        case(33)
        write(*,*)'set number of anisotropic refinements'
        write(*,*)'-1 for infinite loop'
        read(*,*), iso_ans
        info=0
        if(iso_ans.ne.-1) then
        !ref_xyz=1
        !call setAnisoRef(info,ref_xyz, info)
        ref_xyz=2
        do iii=1,iso_ans
          call setAnisoRef(info,ref_xyz, info)
        enddo
        call update_gdof
        call update_Ddof
        else
        !ref_xyz=1
        !call setAnisoRef(info,ref_xyz, info)
        ref_xyz=2
        do while(1.gt.0)
          call setAnisoRef(info,ref_xyz, info)
        enddo
        call update_gdof
        call update_Ddof
        endif
!     adaptive H-refinements
!!      case(33) ; call adaptivity_geom(0.3d0, nvoid,rvoid)
!
!     frontal solve
      case(40)
        call solve1(MY_NR_RHS)
!
!     MUMPS solve
      case(50)
        NRFL = 0
        call uhm_time_in
        !call hack_impedance
        call mumps_interf(MY_NR_RHS)
        !call mumps_sc_3D
        !call unhack_impedance
        call uhm_time_out(t)
        write(*,*) 'time mumps = ', t
!
!     UHM solve
      case(60)
        call uhm_solve
        call uhm_solver_flush(UHM_SOLVER_PTR)
!
!     compute Hcurl error for the E-field only
      case(100)
        flag(1)=1
        call compute_error(flag,1)
!
! ......... rate test
        case(110)
        write(*,*) 'Testing h-convergence rates'
!  ........ Get #refinements to be done
        write(*,*) 'Enter number of uniform H-refinements:'
        read(*,*) numRef
          iref = 0
          do while(iref.lt.numRef)
!  ........ Solve the problem
            !call mumps_sc_3D
            call mumps_interf(MY_NR_RHS)
!  ........ Set error flags
            flag(1) = 1
!  ........ Compute error
            call compute_error(flag,1)
!  ........ Do uniform h-refinements
            call global_href
            call close_mesh
            call update_gdof
            call update_Ddof
            iref = iref+1
          enddo

! ......... 1 uniform n anisotropic refinements and solve
        case(120)
        write(*,*) '2 uniform n anisotropic refinements and solve'
!  ........ Get #refinements to be done
        write(*,*) 'Enter number of anisotropic H-refinements:'
        read(*,*) numRef
        info = 0
        call global_href
        !call close_mesh
        call update_gdof
        call update_Ddof
        call mumps_interf(MY_NR_RHS)
        flag(1) = 1
        write(*,*) 'after one uniform refinement'
        call compute_error(flag,1)

        call global_href
        !call close_mesh
        call update_gdof
        call update_Ddof
        call mumps_interf(MY_NR_RHS)
        flag(1) = 1
        write(*,*) 'after 2 uniform refinement'
        call compute_error(flag,1)

          iref = 0
          do while(iref.lt.numRef)
!  ........ Do anisotropic h-refinements
            ref_xyz=2
            call setAnisoRef(info,ref_xyz, info)
            call update_gdof
            call update_Ddof
! ........... solve
            call mumps_interf(MY_NR_RHS)
!  ........ Set error flags
            flag(1) = 1
!  ........ Compute error
            write(*,*) 'after anisotropic refinement'
            call compute_error(flag,1)
            iref = iref+1
          enddo
!     compute BC data interpolation error
      case(130)
        call compute_BC_interp_error
!
      case(200) ; call my_tests
!
      case(210) ; call my_own_tests
      endselect
!
!  ...go back to menu
      goto 10
!
!
endprogram main

      subroutine my_tests
      end subroutine my_tests



subroutine setAnisoRef(info,ref_xyz, info1)
! This routine anisotropically refines
! prism and hexa elements
! info = 0 if success
! ref_xyz: refinement flag to set kref
!          1 for xy plane direction
!          2 for z direction
! info = 1 otherwise
#include "implicit_none.h"
    !
    use control
    use data_structure3D
    use uhm2
    use parameters, only : ZERO,ZONE,ZIMG
      implicit none
      integer,                       intent(in)  :: info
      integer,                       intent(in)  :: ref_xyz
      integer,                       intent(out)  :: info1
      integer, allocatable, dimension(:) :: list_elem
      integer :: i,iprint,ic,mdle,iel,kref,nr_elem_to_refine
      write(*,*) 'From IANISOREF before anything'
        !call pause

        allocate (list_elem(NRELES))
        write(*,*) 'NRELES is: ', NRELES
        !call pause
        ic=0
        mdle=0
        do iel=1,NRELES
            call nelcon(mdle, mdle)
            ic=ic+1
            list_elem(ic) = mdle
        enddo
        nr_elem_to_refine = NRELES
        if (nr_elem_to_refine.gt.0) then
            !      ...refine the elements from the list
            do iel=1,nr_elem_to_refine
                mdle = list_elem(iel)
                select case(NODES(mdle)%type)
!               PRISM
                case('mdlp')
                  if(ref_xyz.eq.2) then
                    kref=01
                  elseif(ref_xyz.eq.1) then
                    kref=10
                  else
                    write(*,*) 'error from IANISOREF: invalid ref_xyz in prism'
                    stop
                  endif
!
!               BRICK
                case('mdlb')
                  if(ref_xyz.eq.2) then
                    kref=001
                  elseif(ref_xyz.eq.1) then
                    kref=010
                  else
                    write(*,*) 'error from IANISOREF: invalid ref_xyz in bric'
                    stop
                  endif
                end select
                call refine(mdle,kref)
            enddo
            !      ...close the mesh
            call close
            !      ...update geometry and Dirichlet flux dof after the refinements
            !call update_gdof
            !call update_Ddof
        endif
        info1 = 1
!
!
end subroutine setAnisoRef

subroutine hack_impedance

end subroutine hack_impedance

