MODULE param

INTEGER, PARAMETER :: nx = 101	
INTEGER, PARAMETER :: ny = 51	

REAL :: hzero(0:ny+1,0:nx+1), h(0:ny+1,0:nx+1)
REAL :: eta(0:ny+1,0:nx+1),etan(0:ny+1,0:nx+1)
REAL :: u(0:ny+1,0:nx+1), un(0:ny+1,0:nx+1)
REAL :: v(0:ny+1,0:nx+1), vn(0:ny+1,0:nx+1)
REAL :: dt, dx, dy, g, r, taux, tauy, rho, ah
REAL :: slip
REAL :: eps ! parameter for Shapiro filter

INTEGER :: j,k

INTEGER :: wet(0:ny+1,0:nx+1)
REAL :: hmin

REAL :: CuP(0:ny+1,0:nx+1), CuN(0:ny+1,0:nx+1)
REAL :: CvP(0:ny+1,0:nx+1), CvN(0:ny+1,0:nx+1)
REAL :: Cu(0:ny+1,0:nx+1), Cv(0:ny+1,0:nx+1)
REAL :: B(0:ny+1,0:nx+1), BN(0:ny+1,0:nx+1)
REAL :: T(0:ny+1,0:nx+1), TN(0:ny+1,0:nx+1)

INTEGER :: mode

END MODULE param