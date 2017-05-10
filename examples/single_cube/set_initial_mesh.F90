!
!  ...single cube example
!
      subroutine set_initial_mesh(Nelem_order)
!
      use parameters       , only : MAXP
      use data_structure3D 
      use DPGHcurl         , only : ORDER_APPROX
!
      implicit none
      integer,dimension(NRELIS),intent(out) :: Nelem_order
!
      integer :: iel,ndom,if,nrf
      integer,dimension(6,2) :: ibc
      integer,save :: ivis=0   ! visitation flag   
!
!----------------------------------------------------------------------
!
!     check if have not exceeded the maximum order
      if (ORDER_APPROX.gt.MAXP) then
        write(*,*) 'set_initial_mesh: ORDER_APPROX, MAXP = ', ORDER_APPROX,MAXP
        stop
      endif
      write(*,*) 'set_initial_mesh: SINGLE CUBE EXAMPLE'; call pause
      if (NRELIS.ne.1) then
        write(*,*) 'set_initial_mesh: NRELIS = ',NRELIS; stop1
      endif
! 
!     loop through initial mesh elements
      do iel=1,NRELIS
!      
!  STEP 1 : set up order of approximation (elements can have different
!           orders in each direction)
        select case(ELEMS(iel)%Type)
        case('pris') ; Nelem_order(iel)=ORDER_APPROX*11
        case('bric') ; Nelem_order(iel)=ORDER_APPROX*111
        case('tetr') ; Nelem_order(iel)=ORDER_APPROX*1
        case('pyra') ; Nelem_order(iel)=ORDER_APPROX*1
        endselect
!
!  STEP 2 : set up physics
!
!       set up the number of physical attributes supported by the element
        ELEMS(iel)%nrphysics=2
        allocate(ELEMS(iel)%physics(2))
        ELEMS(iel)%physics(1)='Elfld'
        ELEMS(iel)%physics(2)='curnt'
!
!       initialize BC flags
!
!  .....remember that you can use only the GMP attributes to set up the BC!
!!        call domain_number(iel, ndom)
!
!  .....Dirichlet BC
        ibc(1:6,1) = 1 ; ibc(1:6,2) = 0
!
!       allocate BC flags (1 per attribute)
        allocate(ELEMS(iel)%bcond(2)) 
!
!       for each attribute, encode face BC into a single BC flag
        call encodg(ibc(1:6,1),10,6, ELEMS(iel)%bcond(1))
        call encodg(ibc(1:6,2),10,6, ELEMS(iel)%bcond(2))
!
!     loop through initial mesh elements
      enddo
!
!
      end subroutine set_initial_mesh
