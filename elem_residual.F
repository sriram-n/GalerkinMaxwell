c--------------------------------------------------------------------
c                                                                     
c     routine name      - elem_residual
c                                                                     
c-------------------------------------------------------------------- 
c                                                                     
c     latest revision:  - May 14
c                                                                     
c     purpose:          - routine returns element residual (squared)
c                         for the DPGHcurl formulation for Maxwell
c                         equations
c                                                                    
c     arguments:                                                     
c                                                                     
c     in:              
c             Mdle      - an element middle node number, identified
c                         with the element
c     out:              
c             Resid     - element residual (squared)
c             Nref_flag - suggested h-refinement flag
c
c---------------------------------------------------------------------
c
      subroutine elem_residual(Mdle, Resid,Nref_flag)
c
      use control
      use parametersDPG
      use element_data
      use data_structure3D
      use DPGHcurl
      use problem
#include "syscom.blk"
c
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
c  ...solution dof
      dimension zdofH(MAXEQNH,MAXbrickH),
     .          zdofE(MAXEQNE,MAXbrickE),
     .          zdofV(MAXEQNV,MAXbrickV),
     .          zdofQ(MAXEQNQ,MAXbrickQ)
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
      dimension xiloc(3,MAXNINT3ADD),waloc(MAXNINT3ADD)
c
c  ...2D quadrature data
      dimension tloc(2,MAXNINT2ADD),wtloc(MAXNINT2ADD)
c
c  ...impressed current
      dimension zJ(3)
c
c  ...shape functions and their curl after Piola transform, rotated
c     tangential flux
      dimension E1(3),curlE1(3),E2(3),curlE2(3),rntimesE(3)
c
c  ...approximate solution
      dimension zsolExi(3),zcurlExi(3),zsolE(3),zcurlE(3),
     .          zsolFxi(3),zsolF(3),zflux(3),
     .          zcurl_exact(3),zflux_exact(3)
c
c  ...directional contributions to element residual
      dimension residd(0:4),zresidd(0:4),nref(3)
c
c  ...error representation function
      dimension zpsi(3),zcurl_xi_psi(3)
c
      dimension aux(10)
c
c  ...exact solution (for debugging)
      dimension  zvalH(  MAXEQNH    ),
     .          zdvalH(  MAXEQNH,3  ),
     .         zd2valH(  MAXEQNH,3,3),
     .           zvalE(3,MAXEQNE    ),
     .          zdvalE(3,MAXEQNE,3  ),
     .         zd2valE(3,MAXEQNE,3,3),
     .           zvalV(3,MAXEQNV    ),
     .          zdvalV(3,MAXEQNV,3  ),
     .         zd2valV(3,MAXEQNV,3,3),
     .           zvalQ(  MAXEQNQ    ),
     .          zdvalQ(  MAXEQNQ,3  ),
     .         zd2valQ(  MAXEQNQ,3,3)
c
      character uplo
      nk(k1,k2) = (k2-1)*k2/2+k1
c
c---------------------------------------------------------------------
c
      select case(Mdle)
      case(1)
        iprint=0
      case default
        iprint=0
      end select
c
c  ...element type
      etype = NODES(Mdle)%type
      nrv = nvert(etype); nre = nedge(etype); nrf = nface(etype)
c
c  ...determine order of approximation
      call find_order(Mdle, norder)
c
c  ...set the enriched order of appoximation
      select case(etype)
      case('mdlb'); nordP = NODES(Mdle)%order+NORD_ADD*111
      case('mdln','mdld'); nordP = NODES(Mdle)%order+NORD_ADD
      case('mdlp'); nordP = NODES(Mdle)%order+NORD_ADD*11
      end select
c
c  ...determine edge and face orientations
      call find_orient( Mdle, norient_edge,norient_face)
c                                                                     
c  ...determine nodes coordinates 
      call nodcor(Mdle, xnod)
c
c  ...determine solution dof
      call solelm(Mdle, zdofH,zdofE,zdofV,zdofQ)
      if (iprint.eq.1) then
        write(*,7020) xnod(1,1:nrv)
 7020   format('elem_residual: xnod  = ',8(f8.3,2x))
        write(*,7025) xnod(2,1:nrv)
        write(*,7025) xnod(3,1:nrv)
 7025   format('                       ',8(f8.3,2x))
        write(*,7030) 1,zdofE(1,1:nre)
        write(*,7030) 2,zdofE(2,1:nre)
 7030   format('elem_residual: zdofE(',i1',*) = ',2(/,6(2e12.5,2x)))
        call pause
      endif
c
c  ...clear space for auxiliary matrices
      BLOADE = ZERO; AP = ZERO
c
c  ...complex wave number
      zk2 = OMEGA**2*EPSILON - ZI*OMEGA*SIGMA
c
c-----------------------------------------------------------------------
c
c  ...element integrals...
c
c  ...use the enriched order to set the quadrature
      INTEGRATION = NORD_ADD
      call set_3Dint(etype,norder, nint3,xiloc,waloc)
      INTEGRATION = 0
      do l=1,nint3
        xi(1:3) = xiloc(1:3,l)
        wa = waloc(l)
c
c  .....determine element H1 shape functions
        call shape3H(etype,xi,norder,norient_edge,norient_face, 
     .               nrdofH,shapH,gradH)
c
c  .....determine element H(curl) shape functions
        call shape3E(etype,xi,norder,norient_edge,norient_face, 
     .               nrdofE,shapE,curlE)
c
c  .....determine discontinuous H(curl) shape functions
        call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
c
c  .....geometry
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, 
     .              x,dxdxi,dxidx,rjac,iflag)
c
c  .....integration weight 
        weight = rjac*wa
c
c  .....compute the approximate solution
        zsolExi = ZERO; zcurlExi(1:3) = ZERO
        do k=1,nrdofE
          zsolExi(1:3)  = zsolExi(1:3)  + zdofE(1,k)*shapE(1:3,k)
          zcurlExi(1:3) = zcurlExi(1:3) + zdofE(1,k)*curlE(1:3,k)
        enddo
        zsolE(1:3) = zsolExi(1)*dxidx(1,1:3) 
     .             + zsolExi(2)*dxidx(2,1:3) 
     .             + zsolExi(3)*dxidx(3,1:3) 
        zcurlE(1:3) = dxdxi(1:3,1)*zcurlExi(1)
     .              + dxdxi(1:3,2)*zcurlExi(2)
     .              + dxdxi(1:3,3)*zcurlExi(3)
        zcurlE(1:3) = zcurlE(1:3)/rjac
c
        call exact(x,Mdle, zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE, 
     .                     zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
        zcurl_exact(1) = zdvalE(3,1,2)-zdvalE(2,1,3)
        zcurl_exact(2) = zdvalE(1,1,3)-zdvalE(3,1,1)
        zcurl_exact(3) = zdvalE(2,1,1)-zdvalE(1,1,2)

        if (iprint.eq.1) then
          write(*,8010) Mdle,x(1:3)
 8010     format('elem_residual: Mdle = ',i6,',  x = ',3f8.3)
          write(*,8020) zsolE(1:3), zcurlE(1:3)
 8020     format('  approximate solution and curl = ',
     .           3(2e12.5,2x),3x,3(2e12.5,2x))
          write(*,8030) zvalE(1:3,1),zcurl_exact(1:3)
 8030     format('  exact solution and curl       = ',
     .           3(2e12.5,2x),3x,3(2e12.5,2x))

        endif
c
c  .....debugging with the exact solution...
        if (iprint.eq.10) then
          zsolE(1:3)  = zvalE(1:3,1)
          zcurlE(1:3) = zcurl_exact(1:3)
        endif
c
c  .....get the RHS
        call getf(Mdle,x, zJ)
c
c  .....loop through enriched H(curl) test functions
        do k1=1,nrdofEE
          E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) 
     .            + shapEE(2,k1)*dxidx(2,1:3) 
     .            + shapEE(3,k1)*dxidx(3,1:3) 
          curlE1(1:3) = dxdxi(1:3,1)*curlEE(1,k1)
     .                + dxdxi(1:3,2)*curlEE(2,k1)
     .                + dxdxi(1:3,3)*curlEE(3,k1)
          curlE1(1:3) = curlE1(1:3)/rjac
c
c  .......accumulate for the load
          BLOADE(k1) = BLOADE(k1)
     .      + ((curlE1(1)*zcurlE(1)
     .         +curlE1(2)*zcurlE(2)
     .         +curlE1(3)*zcurlE(3))/MU
     .        -zk2*(E1(1)*zsolE(1)+E1(2)*zsolE(2)+E1(3)*zsolE(3))
     .      + zi*OMEGA*(zJ(1)*E1(1)+zJ(2)*E1(2)+zJ(3)*E1(3)))*weight
ccc          write(*,*) 'l,k1,BLOADE(k1) = ', l,k1,BLOADE(k1)
ccc          if (k1.eq.nrdofEE) call pause
c
c  .......loop through enriched H(curl) trial functions
          do k2=k1,nrdofEE
            E2(1:3) = shapEE(1,k2)*dxidx(1,1:3) 
     .              + shapEE(2,k2)*dxidx(2,1:3) 
     .              + shapEE(3,k2)*dxidx(3,1:3) 
            curlE2(1:3) = dxdxi(1:3,1)*curlEE(1,k2)
     .                  + dxdxi(1:3,2)*curlEE(2,k2)
     .                  + dxdxi(1:3,3)*curlEE(3,k2)
            curlE2(1:3) = curlE2(1:3)/rjac
c
c  .........accumulate for the test stiffness matrix
            k = nk(k1,k2)
            select case(INNER_PRODUCT)
            case(1)
              AP(k) = AP(k) 
     .              + (E1(1)*E2(1)+E1(2)*E2(2)+E1(3)*E2(3)
     .                + curlE1(1)*curlE2(1)+curlE1(2)*curlE2(2)
     .                + curlE1(3)*curlE2(3))*weight
            end select
          enddo
        enddo
      enddo
      if (iprint.eq.1) call pause
      if (iprint.eq.2) then
        do i=1,10
          do j=1,i-1
            aux(j) = AP(nk(j,i))
          enddo
          do j=i,10
            aux(j) = AP(nk(i,j))
          enddo
          write(*,7011) aux
 7011   format(10e12.5)
        enddo
        call pause
      endif 
      if (iprint.ge.1) then
        write(*,7014) BLOADE(1:nrdofEE)
 7014   format('elem_residual: BLOADE AFTER VOL INT = ',
     .         10(/,6(2e12.5,2x)))
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  ...boundary integrals
c
c  ...loop through element faces
      do if=1,nrf
c
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
        INTEGRATION = NORD_ADD
        call set_2Dint(ftype,norderf, nint,tloc,wtloc)
        INTEGRATION = 0
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
c  .......determine discontinuous H(curl) shape functions
          call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
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
c
c  .......compute approximate flux at the point
          zsolFxi = ZERO
          do k=1,nrdofE
            zsolFxi(1:3)  = zsolFxi(1:3)  + zdofE(2,k)*shapE(1:3,k)
          enddo
          zsolF(1:3) = zsolFxi(1)*dxidx(1,1:3) 
     .               + zsolFxi(2)*dxidx(2,1:3) 
     .               + zsolFxi(3)*dxidx(3,1:3)
          call zcross_product(rn,zsolF, zflux) 
          call exact(x,Mdle, 
     .              zvalH,zdvalH,zd2valH, zvalE,zdvalE,zd2valE, 
     .              zvalV,zdvalV,zd2valV, zvalQ,zdvalQ,zd2valQ)
          call zcross_product(rn,zvalE(1:3,2), zflux_exact) 
          if (iprint.eq.1) then
            write(*,8040) if,x(1:3)
 8040       format('elem_residual: if = ',i6,' x = ',3f8.3)
            write(*,8050) zflux
 8050       format('  approximate flux = ',3(2e12.5,2x))
            write(*,8060) zflux_exact
 8060       format('  exact flux       = ',3(2e12.5,2x))
          endif
c
c  .......debugging with the exact solution...
          if (iprint.eq.10) then
            zflux(1:3) = zflux_exact(1:3)
          endif
c
c  .......loop through enriched test functions
          do k1=1,nrdofEE
            E1(1:3) = shapEE(1,k1)*dxidx(1,1:3) 
     .              + shapEE(2,k1)*dxidx(2,1:3) 
     .              + shapEE(3,k1)*dxidx(3,1:3) 
c
c  .........accumulate for the load vector
            BLOADE(k1) = BLOADE(k1) 
     .      + (E1(1)*zflux(1)+E1(2)*zflux(2)+E1(3)*zflux(3))
     .       *weight
          enddo
        enddo
        if (iprint.eq.1) call pause
      enddo
      if (iprint.ge.1) then
        write(*,7015) BLOADE(1:nrdofEE)
 7015   format('elem_residual: FINAL BLOADE = ',10(/,6(2e12.5,2x)))
        call pause
      endif
c
c-----------------------------------------------------------------------
c
c  ...factorize the test stiffness matrix
      uplo = 'U'
      call ZPPTRF(uplo, nrdofEE, AP, info) 
      if (info.ne.0) then
        write(*,*) 'elem_dpgH1: info = ',info
        stop1
      endif
c
c  ...save copies of the RHS to compute later the residual
      BLOADEc = BLOADE
c
c  ...compute the product of inverted test Gram matrix with RHS,
c     BLOADE is overwritten with the solution
      call ZPPTRS(uplo, nrdofEE, 1, AP, BLOADE, MAXbrickEE, info1 )
      if (info1.ne.0) then
        write(*,*) 'elem_dpgH1: info1 = ',info1
        stop1
      endif
c
c  ...compute the residual
      zresid = ZERO
      do k=1,NRDOFEE
        zresid = zresid + BLOADEc(k)*conjg(BLOADE(k))
      enddo
      Resid = zresid
c
c-----------------------------------------------------------------------
c
c  ...recompute the element residual through direct integration to
c     establish anisotropy flags
      zresidd(0:3) = ZERO
      do l=1,nint3
        xi(1:3) = xiloc(1:3,l)
        wa = waloc(l)
c
c  .....determine element H1 shape functions (for geometry)
        call shape3H(etype,xi,norder,norient_edge,norient_face, 
     .               nrdofH,shapH,gradH)
c
c  .....determine discontinuous H(curl) shape functions
        call shape3EE(etype,xi,nordP, nrdofEE,shapEE,curlEE)
c
c  .....geometry
        call geom3D(Mdle,xi,xnod,shapH,gradH,nrdofH, 
     .              x,dxdxi,dxidx,rjac,iflag)
c
c  .....integration weight 
        weight = rjac*wa
c
c  .....compute the error representation function
        zpsi = ZERO; zcurl_xi_psi(1:3) = ZERO
        do k=1,nrdofEE
          zpsi(1:3) = zpsi(1:3) + BLOADE(k)*shapEE(1:3,k)
          zcurl_xi_psi(1:3) = zcurl_xi_psi(1:3) 
     .                      + BLOADE(k)*curlEE(1:3,k)
        enddo
c
        select case(INNER_PRODUCT)
        case(1)
          do i=1,3; do j=1,3
            axi = 0.d0; bxi = 0.d0
            do k=1,3
              axi = axi + dxidx(i,k)*dxidx(j,k)
              bxi = bxi + dxdxi(k,i)*dxdxi(k,j)
            enddo
            bxi = bxi/rjac**2
c
c  .........only the vector contributes to directional residuals
            if (i.eq.j) then
              zresidd(i) = zresidd(i) + axi*abs(zpsi(i))**2*weight
            else
              zresidd(0) = zresidd(0) 
     .                   + axi*zpsi(i)*conjg(zpsi(j))*weight
            endif
            zresidd(0) = zresidd(0) 
     .      + bxi*zcurl_xi_psi(i)*conjg(zcurl_xi_psi(j))*weight
          enddo; enddo
        end select
c
c  ...end of loop through integration points
      enddo
      residd(0:3) = zresidd(0:3)
      diff = residd(0)+residd(1)+residd(2)+residd(3) - Resid
      if (abs(diff).gt.1.d-8*abs(Resid)) then
        write(*,*) 'Resid = ',Resid,
     .              residd(0)+residd(1)+residd(2)+residd(3)
        write(*,*) 'residd = ',residd(0:3)
        call pause
      endif
c
c  ...determine the refinement flag
      select case(etype)
      case('mdlb')
        if (residd(0).lt..1d0*Resid) then
          nref(1:3) = 1
          do i=1,3
            if (residd(i).lt..1d0*Resid) nref(i)=0
          enddo
          Nref_flag = nref(1)*100+nref(2)*10+nref(3)
        else
          Nref_flag = 111
        endif
      case('mdln','mdld')
        Nref_flag = 1
      case('mdlp')
        if (residd(0).lt..1d0*Resid) then
          nref(1:2) = 1
          if (residd(1)+residd(2).lt..2d0*Resid) nref(1)=0
          if (residd(3).lt..1d0*Resid) nref(2)=0
          Nref_flag = nref(1)*10+nref(2)
        else
          Nref_flag = 111
        endif
      end select
ccc      write(*,*) 'residd = ',residd(0:3)
ccc      write(*,*) 'Mdle,Nref_flag = ',Mdle,Nref_flag
ccc      call pause
c
c-----------------------------------------------------------------------
c
      if (iprint.eq.1) then
        write(*,7010) Mdle, Resid
 7010   format('elem_residual: Mdle, Resid = ',i5,3x,e12.5)
        call pause
      endif
c
      end subroutine elem_residual    

c-----------------------------------------------------------------------
c                                                                     
c     routine name      - compute_residual
c                                                                     
c-------------------------------------------------------------------- 
c                                                                     
c     latest revision:  - Apr 14
c                                                                     
c     purpose:          - routine returns global residual 
c                         for the DPGH1 formulation for Laplace
c                         equation
c                                                                    
c---------------------------------------------------------------------
c
      subroutine compute_residual
c
      use data_structure3D
      use environment      , only : QUIET_MODE
#include "syscom.blk"
c
c  ...visitation flag
      integer, save :: ivis = 0
c
c  ...total number of dof for the old mesh
      integer, save :: nrdof_total_old
c
c  ...residual for the old mesh
      real*8,  save :: residual_old
c
c  ...residuals and rates to display
      real*8 , dimension(2,10), save :: rwork
c
c  ...number of dof to display
      integer , dimension(10), save :: iwork
    
c
      iprint=1
c
c  ...compute total residual and number of dof
      nrdof_total = 2*NRDOFSE; residual = 0.d0
      mdle = 0
      do iel=1,NRELES
        call nelcon(mdle, mdle)
        call elem_residual(mdle, resid,nref_flag)
        if (iprint.eq.1) then
          write(*,7010) iel, mdle, resid
 7010     format('compute_residual: iel, mdle = ',2i8,
     .           ' element residual = ',e12.5)
        endif
        residual = residual + resid
        call find_ndof(mdle, ndofH,ndofE,ndofV,ndofQ)
c
c  .....subtract the middle node H(div) dof from the global count
        nrdof_total = nrdof_total - ndofE
      enddo
      residual = sqrt(residual)
c
c  ...compute rate
      rate = 0.d0 
      if (ivis.ne.0) then
        if (nrdof_total.gt.nrdof_total_old) then
          rate = log(residual_old/residual)
     .          /log(float(nrdof_total_old)/float(nrdof_total))
        endif
      endif
c
c  ...save current data
      ivis = ivis+1
      nrdof_total_old = nrdof_total
      residual_old = residual
c
c  ...store data to display
      iwork(ivis) = nrdof_total
      rwork(1,ivis) = residual; rwork(2,ivis) = rate
c
c  ...display the convergence history
      if (.NOT. QUIET_MODE) then
        write(*,*)''
        write(*,*)'         -- Error Report --'
        write(*,7100)
 7100   format(' Mesh  ','  Nrdof  ', ' Residual   ','     Rate ')
        do i=1,ivis     
          write(*,7110)i,iwork(i),rwork(1:2,i)
 7110     format(i3,4x,i7,2x,e12.5,2x,f8.3)
        enddo
        write(*,*)''
      endif
c
c
      end subroutine compute_residual






