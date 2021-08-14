MODULE param

! DECLARATION SECTION

INTEGER, PARAMETER :: nx = 101	

REAL :: hzero(0:nx+1), h(0:nx+1)
REAL :: eta(0:nx+1),etan(0:nx+1)
REAL :: u(0:nx+1), un(0:nx+1)
REAL :: dt,dx,g
REAL :: eps ! parameter for Shapiro filter

INTEGER :: k

!Parameters for wetting-drying algorithm
INTEGER :: wet(0:nx+1)
REAL :: hmin

END MODULE param