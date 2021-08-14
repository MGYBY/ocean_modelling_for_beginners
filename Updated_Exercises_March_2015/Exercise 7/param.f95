MODULE param

! DECLARATION SECTION
 
INTEGER, PARAMETER :: nx = 101 ! number of cells in x direction	
INTEGER, PARAMETER :: nz = 10 ! number of layers	

REAL :: htotal(0:nx+1) ! initial bathymetry
REAL :: hzero(nz,0:nx+1) ! initial layer thicknesses
REAL :: h(nz,0:nx+1) ! actual layer thicknesses
REAL :: dp(0:nz,0:nx+1) ! dynamic pressure 
REAL :: eta(nz+1,0:nx+1) ! actual interface displacements
REAL :: etan(nz,0:nx+1) ! interface displacements at time level n+1
REAL :: eta0(nz+1,0:nx+1) ! initial interface displacements 
REAL :: dhdt(nz,0:nx+1) ! actual layer-thickness change
REAL :: u(nz,0:nx+1) ! actual lateral velocity
REAL :: un(nz,0:nx+1) ! lateral velocity at time level n+1
REAL :: dt, dx, g ! time step, grid spacing, acceleration due to gravity
REAL :: eps ! parameter for Shapiro filter
REAL :: rho(0:nz) ! layer densities
REAL :: hmin ! minimum layer thickness
INTEGER :: i,k  ! layer index and cell index
INTEGER :: wet(nz,0:nx+1)  ! wet and dry pointer

END MODULE param