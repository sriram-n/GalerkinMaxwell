c----------------------------------------------------------------------
c
c   module name        - matrices
c
c----------------------------------------------------------------------
c
c   latest revision    - May 2
c
c   purpose            - store element matrices for the first layer
c                        elements
c
c----------------------------------------------------------------------
c
        module matrices
c
        use parametersDPG
        implicit none
#if C_MODE
#define VTYPE  complex*16
#else
#define VTYPE double precision
#endif
c
c max # of elements in the first layer
        integer, parameter :: MAXNRFL=32
c
c # of stored elements in the first layer
        integer :: NRFL
cc$OMP THREADPRIVATE (NRFL)
c
c order of elements
        integer, parameter :: MYP=6
c
c matrix dimensions
        integer, parameter :: MYE = 3*MYP*(MYP+1)**2*2
c
c stiffness matrices to store
        VTYPE, dimension(MYE,MYE,MAXNRFL) :: ZFL_EE
cc$OMP THREADPRIVATE (ZFL_EE)
        VTYPE, dimension(MYE,MAXNRFL) :: ZLOADFL_E
c
c xy coordinates of the first vertex node
        real*8, dimension(2,MAXNRFL) :: XYVERT
cc$OMP THREADPRIVATE (XYVERT)
c
        end module matrices



