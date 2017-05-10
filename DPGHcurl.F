c----------------------------------------------------------------------
c
c   module name        - DPGHcurl
c
c----------------------------------------------------------------------
c
c   latest revision    - May 14
c
c   purpose            - work space for the primal DPG method related 
c                        computations of optimal test functions and 
c                        energy error 
c
c----------------------------------------------------------------------
c
      module DPGHcurl
c
      use parametersDPG
#if C_MODE
#define V_TYPE  complex*16
#else     
#define V_TYPE double precision 
#endif
c
c  ...H1 discontinuous shape functions
      double precision shapEE(3,MAXbrickEE),curlEE(3,MAXbrickEE)
c
c  ...stiffnes matrix for the local Riesz Hcurl matrix in Lapack format
      V_TYPE AP(MAXbrickEE*(MAXbrickEE+1)/2)
c
c  ...load vector for the enriched space
      V_TYPE BLOADE(MAXbrickEE)
c
c  ...copy of load vector for the enriched space
      V_TYPE BLOADEc(MAXbrickEE)
c
c  ...stiffnes matrices for the enriched test space
      V_TYPE STIFFTE(MAXbrickEE,MAXbrickE)
      V_TYPE STIFFTF(MAXbrickEE,MAXbrickE)
c
c  ...copies of the same matrices
      V_TYPE STIFFTEc(MAXbrickEE,MAXbrickE)
      V_TYPE STIFFTFc(MAXbrickEE,MAXbrickE)
c
c  ...additional parameters including those required by the system
      integer :: ORDER_APPROX
      integer :: NPX, NPY, NPZ, ICOMP_EXACT
      integer,parameter :: MY_NR_RHS=1
      integer :: ICHOOSE_DISP
      integer :: IEXACT_DISP
c
c
      end module DPGHcurl


