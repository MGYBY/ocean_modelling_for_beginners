MODULE param

! DECLARATION SECTION
! 
! Variables and parameters shared by 
! all components of the code.
!
INTEGER, PARAMETER :: nx = 101	

REAL :: hzero(0:nx+1), h(0:nx+1)
INTEGER :: wet(0:nx+1)
REAL :: eta(0:nx+1),etan(0:nx+1)
REAL :: u(0:nx+1), un(0:nx+1)
REAL :: w(0:nx+1), wu(0:nx+1)
REAL :: dt,dx,g
REAL :: eps ! parameter for Shapiro filter

INTEGER :: k  ! grid index

END MODULE param