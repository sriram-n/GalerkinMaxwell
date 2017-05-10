
!  ...user defined routine to define the initial mesh
!
      subroutine set_initial_mesh(Nelem_order)
!
      use GMP
      use parameters
      use data_structure3D
      use problem
!
      implicit none
      integer,dimension(NRELIS),intent(out) :: Nelem_order ! polynomial order for initial mesh elements
!
      integer :: iel,ndom,i,ifc,neig
      integer,dimension(6,4) :: ibc
!
!----------------------------------------------------------------------
!
!  ...check if have not exceeded the maximum order
      if (ORDER_APPROX.gt.MAXP) then
        write(*,*) 'set_initial_mesh: ORDER_APPROX, MAXP = ', ORDER_APPROX,MAXP
        stop
      endif
!
!  ...loop through initial mesh elements
      do iel=1,NRELIS
!
!  STEP 1 : set up order of approximation (elements may have different
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
        ELEMS(iel)%nrphysics=1
        allocate(ELEMS(iel)%physics(1))
        ELEMS(iel)%physics(1)='Elfld'
!
!       initialize BC flags
!
! .... single cube: impedance
        ibc(1:6,1) = (/1,9,1,1,1,1/)
! !  .....single cube with Dirichlet, PEC BC on 5 out of 6 faces
! !  .....Impedance one one face
!         ibc(1:6,1) = (/1,9,1,1,1,1/)
! !  .......prism core geometry
!          ibc(1:6,1) = (/1,9,1,1,1,1/)
!
!       allocate BC flags (one per attribute)
        allocate(ELEMS(iel)%bcond(1))
!
!       for each attribute, encode face BC into a single BC flag
        do i=1,1
          call encodg(ibc(1:6,i),10,6, ELEMS(iel)%bcond(i))
        enddo
!
!  ...end of loop through initial mesh elements
      enddo
!
!
      end subroutine set_initial_mesh
