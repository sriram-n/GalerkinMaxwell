subroutine set_mesh_refinement
!
      use data_structure3D , only : NRELIS,NODES
!
      implicit none
      integer :: i,ndom
!      
!-----------------------------------------------------------------      
!
!     loop over initial mesh elements
      do i=1,NRELIS
!      
!        EITHER set refinement filter based on element shape...
         select case(NODES(i)%type)
!        prism          
         case('mdlp') ; NODES(i)%ref_filter=10
!        hexa                  
         case('mdlb') ; NODES(i)%ref_filter=100
!        tet                  
         case('mdln') ; NODES(i)%ref_filter=24
!        pyramid                  
         case('mdld') ; NODES(i)%ref_filter=10
         endselect        
!
!        OR set filters based on subdomain number
         call domain_number(i, ndom)
         select case(ndom)
         case(1)      ; NODES(i)%ref_filter=10   
         case(2)      ; NODES(i)%ref_filter=24
         case(3)      ; NODES(i)%ref_filter=32   
         endselect
!         
!        OR set them some other way that is convenient 
!
!
      enddo        
!
!      
endsubroutine set_mesh_refinement
