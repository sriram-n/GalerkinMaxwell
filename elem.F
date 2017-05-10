c--------------------------------------------------------------------
c
c     routine name      - elem
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns unconstrained (ordinary)
c                         stiffness matrix and load vector
c                         for the standard Galerkin method for Maxwell
c                         eqn
c
c
c     arguments:
c
c     in:
c             Mdle      - an element middle node number, identified
c                         with the element
c             Maxdof    - column length of Zaloc
c     out:
c             Nrdof     - number of dof for a single component
c             Itest_e,Itest_a,Itrial_e,Itrial_a - flags indicating
c                         presence of corresponding load vectors
c                         and stiffness matrices
c             Zbloc     - load vector
c             Zaloc     - stiffness matrix
c
c-----------------------------------------------------------------------
c
c  ...this is a system routine, the head cannot be changed
      subroutine elem(Mdle, Itest,Itrial)
c
      use control
      use data_structure3D
      use assembly
#include "syscom.blk"
c
      dimension Itest(NR_PHYSA),Itrial(NR_PHYSA)
c
      Itest(1:NR_PHYSA)=0; Itrial(1:NR_PHYSA)=0
c
      select case(NODES(Mdle)%case)
c
      case(1)
c
c  .....this is user-defined routine, you are in charge of it
        Itest(1)=1; Itrial(1)=1
        call elem_Hcurl(Mdle,BLOC(1)%nrow,
     .                  BLOC(1)%array,ALOC(1,1)%array)
c
      case default
        write(*,*) 'elem: Mdle,NODES(Mdle)%case = ',
     .              Mdle,NODES(Mdle)%case
        stop1
      end select

c
c
      end subroutine elem

c--------------------------------------------------------------------
c
c     routine name      - elem_Hcurl
c
c--------------------------------------------------------------------
c
c     latest revision:  - Apr 17
c
c     purpose:          - routine returns unconstrained (ordinary)
c                         stiffness matrix and load vector
c                         for the Galerkin formulation for time-
c                         harmonic Maxwell equations
c
c     arguments:
c
c     in:
c             Mdle      - an element middle node number, identified
c                         with the element
c             MdE       - column length of ZalocEE
c     out:
c             ZblocE    - load vectors
c             ZalocEE   - stiffness matrix
c
c---------------------------------------------------------------------
c
      subroutine elem_Hcurl(Mdle,MdE,
     .                      ZblocE,ZalocEE)
c
      use control
      use parameters
      use element_data
      use data_structure3D
      use problem
#include "syscom.blk"
c
      dimension ZblocE(MdE),ZalocEE(MdE,MdE)
c
      character(len=4) :: etype,ftype
c
c  ...element order, orientation for edges and faces
      dimension norder(19),norient_edge(12),norient_face(6)
c
c  ...face order
      dimension norderf(5)
c
c  ...geometry dof
      dimension xnod(3,MAXbrickH)
c
c  ...geometry
      dimension xi(3),dxidt(3,2),x(3),dxdxi(3,3),dxidx(3,3),
     .          dxdt(3,2),rt(3,2),rn(3),t(2)
c
c  ...H1 shape functions (for geometry)
      dimension shapH(MAXbrickH),gradH(3,MAXbrickH)
c
c  ...Hcurl shape functions
      dimension shapE(3,MAXbrickE),curlE(3,MAXbrickE)
c
c  ...3D quadrature data
      dimension xiloc(3,MAX_NINT3),waloc(MAX_NINT3)
c
c  ...2D quadrature data
      dimension tloc(2,MAX_NINT2),wtloc(MAX_NINT2)
c
c  ...BC's flags
      dimension ibc(6,NR_PHYSA)
c
c  ...impressed current
      dimension zJ(3)
      integer, dimension(27) :: nodesl, norientl
c
c  ...shape functions and their curl after Piola transform, rotated
c     tangential flux
      dimension E1(3),curlE1(3),E2(3),curlE2(3),rn2timesE(3),rntimesE(3)
      dimension zImp(3)
c
c  ...for debug printing
      dimension aux(10)
c
c---------------------------------------------------------------------
c
      iprint=0
c      write(*,*) 'elem: Mdle = ',Mdle
c
c  ...element type
      etype = NODES(Mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
c
c  ...determine order of approximation
      call find_order(Mdle, norder)
c
c  ...determine edge and face orientations
      call find_orient(Mdle, norient_edge,norient_face)
c
c  ...determine nodes coordinates
      call nodcor(Mdle, xnod)
c
c  ...get the element boundary conditions flags
      call find_bc(Mdle, ibc)
C       call elem_nodes(Mdle, nodesl,norientl)

C       do i=nrv+1,nrv+nre
C         write(*,*) 'Mdle, nodesl(i), bcond =', Mdle, nodesl(i),
C      .   NODES(nodesl(i))%bcond
C         call pause
C       enddo

c      iprint = 2
      if (iprint.ge.1) then
        write(*,7001) Mdle
 7001   format('elem_Hcurl: BC FLAGS FOR Mdle = ',i5)
        do i=1,NR_PHYSA
          write(*,7002) PHYSA(i), ibc(1:nrf,i)
 7002     format('          ATTRIBUTE = ',a6,' FLAGS = ',6i2)
        enddo
      endif
c      iprint = 0
c
C C c..............FOR LONG WAVEGUIDE: RE-USING STIFFNESS MATRICES
C c ...check if one needs to recompute element matrices
C       idec=0
C       call copy_element_matrices(MdE,xnod,ibc, idec,
C      .                        ZalocEE,ZblocE)
C       if (idec.eq.1) return
C C c..............FOR LONG WAVEGUIDE: RE-USING STIFFNESS MATRICES
c  ...clear space for stiffness matrix and rhsv:
      ZblocE = ZERO; ZalocEE = ZERO
c
c  ...complex wave number
      zk2 = OMEGA**2*EPSILON - ZI*OMEGA*SIGMA
c
c-----------------------------------------------------------------------
c
c  ...element integrals...
c
c  ...set the quadrature
      call set_3Dint(etype,norder, nint,xiloc,waloc)
      do l=1,nint
        xi(1:3) = xiloc(1:3,l)
        wa = waloc(l)
c
c  .....determine element H1 shape functions (for geometry)
        call shape3H(etype,xi,norder,norient_edge,norient_face,
     .               nrdofH,shapH,gradH)
c
c  .....determine element H(curl) shape functions
        call shape3E(etype,xi,norder,norient_edge,norient_face,
     .               nrdofE,shapE,curlE)
c
c  .....geometry
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,
     .              x,dxdxi,dxidx,rjac,iflag)
c
c  .....integration weight
        weight = rjac*wa
c
c  .....get the RHS
        call getf(Mdle,x, zJ)
c
c  .....loop through H(curl) test functions
        do k1=1,nrdofE
          E1(1:3) = shapE(1,k1)*dxidx(1,1:3)
     .            + shapE(2,k1)*dxidx(2,1:3)
     .            + shapE(3,k1)*dxidx(3,1:3)
          curlE1(1:3) = dxdxi(1:3,1)*curlE(1,k1)
     .                + dxdxi(1:3,2)*curlE(2,k1)
     .                + dxdxi(1:3,3)*curlE(3,k1)
          curlE1(1:3) = curlE1(1:3)/rjac
c
c  .......compute the RHS
          ZblocE(k1) = ZblocE(k1) +
     .      - ZI*OMEGA*(zJ(1)*E1(1)+zJ(2)*E1(2)+zJ(3)*E1(3))*weight
c
c  .......loop through H(curl) trial functions
          do k2=1,nrdofE
            E2(1:3) = shapE(1,k2)*dxidx(1,1:3)
     .              + shapE(2,k2)*dxidx(2,1:3)
     .              + shapE(3,k2)*dxidx(3,1:3)
            curlE2(1:3) = dxdxi(1:3,1)*curlE(1,k2)
     .                  + dxdxi(1:3,2)*curlE(2,k2)
     .                  + dxdxi(1:3,3)*curlE(3,k2)
            curlE2(1:3) = curlE2(1:3)/rjac
c
c  .........accumulate for the stiffness matrix
            ZalocEE(k1,k2) = ZalocEE(k1,k2)
     .      + ((curlE1(1)*curlE2(1)+curlE1(2)*curlE2(2)
     .                   +curlE1(3)*curlE2(3))/MU
     .        -zk2*(E1(1)*E2(1)+E1(2)*E2(2)+E1(3)*E2(3)))*weight
           enddo
        enddo
      enddo
      if (iprint.eq.1) then
        write(*,*) 'elem_Hcurl: ZalocEE = '
        do i=1,10
          write(*,7011) ZalocEE(i,1:10)
        enddo
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  ...boundary integrals
c
c  ...loop through element faces
      do if=1,nrf
        if (ibc(if,1).ne.9) cycle
c         if (IBCFLAG.ne.) cycle
c
c      write(*,*) 'in the ibc boundary'
c      call pause
c  .....sign factor to determine the OUTWARD normal unit vector
        nsign = nsign_param(etype,if)
c
c  .....face type
        ftype = face_type(etype,if)
c
c  .....face order of approximation
        call face_order(etype,if,norder, norderf)
c
c  .....set 2D quadrature
        call set_2Dint(ftype,norderf, nint,tloc,wtloc)
c
c  .....loop through integration points
        do l=1,nint
c
c  .......face coordinates
          t(1:2) = tloc(1:2,l)
c
c  .......face parametrization
          call face_param(etype,if,t, xi,dxidt)
c
c  .......determine element H1 shape functions (for geometry)
          call shape3H(etype,xi,norder,norient_edge,norient_face,
     .                 nrdofH,shapH,gradH)
c
c  .......determine element H(curl) shape functions (for fluxes)
          call shape3E(etype,xi,norder,norient_edge,norient_face,
     .                 nrdofE,shapE,curlE)
c
c  .......geometry
          call bgeom3D(Mdle,xi,xnod,shapH,gradH,nrdofH,dxidt,nsign,
     .                 x,dxdxi,dxidx,rjac,dxdt,rn,bjac)
          weight = bjac*wtloc(l)

          call  get_bdSource(Mdle,x,rn, zImp)
c
c  .......loop through the H(curl) test functions
          do k1=1,nrdofE
            E1(1:3) = shapE(1,k1)*dxidx(1,1:3)
     .              + shapE(2,k1)*dxidx(2,1:3)
     .              + shapE(3,k1)*dxidx(3,1:3)
c
c  .........accumulate for IBC load
c  .......compute the RHS
            ZblocE(k1) = ZblocE(k1) +
c     .       (zImp(1)*E1(1)+zImp(2)*E1(2)+zImp(3)*E1(3))*weight
     .      ZI*OMEGA*(zImp(1)*E1(1)+zImp(2)*E1(2)+zImp(3)*E1(3))*weight
c  .........loop through H(curl) trial functions
            do k2=1,nrdofE
              E2(1:3) = shapE(1,k2)*dxidx(1,1:3)
     .                + shapE(2,k2)*dxidx(2,1:3)
     .                + shapE(3,k2)*dxidx(3,1:3)
c              call dot_product(rn,E2, rndotE2)
              call cross_product(rn,E2, rntimesE)
              call cross_product(rn,rntimesE, rn2timesE)

c  ...........accumulate for the extended stiffness matrix
              ZalocEE(k1,k2) = ZalocEE(k1,k2)
     .   - ZI*OMEGA*(GAMMA)*(E1(1)*rn2timesE(1)+E1(2)*rn2timesE(2)
     .                          +E1(3)*rn2timesE(3))*weight
            enddo
          enddo
        enddo
      enddo
c
C       if (idec.ne.2) then
C         write(*,*) 'elem_dpgMaxwell: idec = ',idec
C         stop 1
C       endif
C       call copy_element_matrices(MdE,xnod,ibc, idec,
C      .                        ZalocEE,ZblocE)
c
      if (iprint.ge.1) then
        write(*,7010)
 7010   format('elem_Hcurl: ZblocE = ')
        write(*,7011) ZblocE(1:NrdofE)
 7011   format(10e12.5)
        write(*,7012)
 7012   format('elem_Hcurl: ZalocEE = ')
        do i=1,NrdofE
          write(*,7013) i,ZalocEE(i,1:NrdofE)
 7013     format('i = ',i3,10(/,5(2e12.5,2x)))
        enddo
        call pause
      endif
c
      end subroutine



      subroutine copy_element_matrices(MdE,Xnod,ibc, Idec,
     .                            ZalocEE,ZblocE)
c
      use control
      use parameters
      use matrices
      implicit none
c
#if C_MODE
#define VTYPE  complex*16
#else
#define VTYPE double precision
#endif
c.......declare input/output variables
      integer,                     intent(in)    :: MdE
      integer, dimension(6,1), intent(in) :: ibc
      VTYPE, dimension(MdE,MdE) :: ZalocEE
      VTYPE, dimension(MdE) :: ZblocE
c
c ...element geometry dof
      real*8, dimension(3,MAXbrickH),intent(in) :: Xnod
      integer  :: Idec
c
      integer :: iprint,i,checkIBC,j
      !checkIBC = 0
c
      iprint=1
c
      select case(Idec)
c
c ..check if one should recompute the element matrices
      case(0)
c ... check if IBC is imposed, if so, always compute matrices
      do j=1,6
        if(ibc(j,1).eq.9) then
         ! checkIBC = checkIBC+1
        !endif
        !if(checkIBC.eq.5) then
          write(*,*) 'ibc is imposed'
          Idec = 2
          return
        endif
      enddo
c
c ....find the corresponding element in the first layer
      do i=1,NRFL
        if ((abs(Xnod(1,1)-XYVERT(1,i)).lt.GEOM_TOL).and.
     .       (abs(Xnod(2,1)-XYVERT(2,i)).lt.GEOM_TOL)) go to 10
      enddo
c      cwrite(*,*) 'copy_element_matrices: NRFL is: ', NRFL
c      cwrite(*,*) 'copy_element_matrices: Xnod(1:3,1) is: ', Xnod(1:2,1)
c      cwrite(*,*) 'i is ', i
c      cwrite(*,*) 'copy_element_matrices: XYVERT(1:2,i) is: ', XYVERT(1:2,i)
c
c ....no element with this position has been computed yet
      Idec=2; return
c
c ....element matrices have been computed
10    Idec=1
c       if (checkIBC.eq.1) then
c         Idec=2
c         return
c       endif
c
      write(*,*) 'not computing the stiffness matrices'
c      cwrite(*,*) 'copy_element_matrices: NRFL is: ', NRFL
c      cwrite(*,*) 'copy_element_matrices: Xnod(1:3,1) is: ', Xnod(1:2,1)
c      cwrite(*,*) 'i is ', i
c      cwrite(*,*) 'copy_element_matrices: XYVERT(1:2,i) is: ', XYVERT(1:2,i)
c ....copy the element matrices from the temporary data structure
      ZalocEE(1:MdE,1:MdE) = ZFL_EE(1:MdE,1:MdE,i)
      ZblocE(1:MdE) = ZLOADFL_E(1:MdE,i)
c
c ....finish generating the matrices
c      do i=1,MdQ
c        ZalocQE(i,1:MdE) = conjg(ZalocEQ(1:MdE,i))
c      enddo
c      return
c
c ..store the element matrices
      case(2)
c
c ....store first vertex node xy coordinates
      NRFL=NRFL+1
      if (NRFL.gt.MAXNRFL) then
       write(*,*) 'copy_element_matrices: NRFL is: ', NRFL
c        write(*,*) 'copy_element_matrices: Xnod(1:3,1) is: ', Xnod(1:3,1)
        write(*,*) 'copy_element_matrices: INCREASE MAXNRFL'
        stop 1
      endif
      XYVERT(1:2,NRFL) = Xnod(1:2,1)
c
c ....copy the element matrices from the temporary data structure
      ZFL_EE(1:MdE,1:MdE,NRFL) = ZalocEE(1:MdE,1:MdE)
      ZLOADFL_E(1:MdE,NRFL) =  ZblocE(1:MdE)
c
      end select
c
c
      end subroutine copy_element_matrices







