MODULE sub
USE param

CONTAINS

!=======================
SUBROUTINE init

hmin = 0.1

! grid parameters
dx = 10.0
dt = 0.1

! physical parameters
g = 9.81

! initial conditions
DO k = 1,nx
  hzero(k) = -20.+10.*real(k)/real(nx)
END DO

DO k = nx/2-5,nx/2+5
  hzero(k) = -14.
END DO

hzero(0) = hzero(1)
hzero(nx+1) = hzero(nx)

DO k = 0,nx+1
  eta(k) = -MIN(0.0,hzero(k))
  etan(k) = eta(k)
END DO

DO k = 0,nx+1
  h(k) = hzero(k)+eta(k)
!XXXXX wet = 1 defines "wet" grid cells
!XXXXX wet = 0 defines "dry" grid cells
  wet(k) = 1
  if(h(k) < hmin)wet(k) = 0
!XXXXXXXXXXX
  u(k) = 0.
  un(k) = 0.
END DO

END SUBROUTINE init

!================
SUBROUTINE dyn

! local parameters
REAL :: pgradx 
REAL :: hue, huw, hwp,hwn,hen,hep

DO k = 1,nx

!XXXXXXXXXXXXX NEW XXXXXXXXXXXXXXXXXXX 
! velocity predictor for wet grid cells

pgradx = -g*(eta(k+1)-eta(k))/dx
un(k) = 0.0
IF(wet(k)==1) THEN
  IF((wet(k+1)==1).or.(pgradx>0)) un(k) = u(k)+dt*pgradx
ELSE
  IF((wet(k+1)==1).and.(pgradx<0)) un(k) = u(k)+dt*pgradx
END IF

END DO

! un conditions for boundaries 
un(0) = 0.
un(nx+1) = 0.

! sea level predictor
DO k = 1,nx
     hep = 0.5*(un(k)+abs(un(k)))*h(k)
     hen = 0.5*(un(k)-abs(un(k)))*h(k+1)
     hue = hep+hen
     hwp = 0.5*(un(k-1)+abs(un(k-1)))*h(k-1)
     hwn = 0.5*(un(k-1)-abs(un(k-1)))*h(k)
     huw = hwp+hwn
     etan(k) = eta(k)-dt*(hue-huw)/dx
END DO

END SUBROUTINE dyn

!======================
SUBROUTINE shapiro

!local parameters
REAL :: term1,term2

! 1-order Shapiro filter

DO k = 1,nx

IF(wet(k)==1)THEN
  term1 = (1.0-0.5*eps*(wet(k+1)+wet(k-1)))*etan(k)
  term2 = 0.5*eps*(wet(k+1)*etan(k+1)+wet(k-1)*etan(k-1))
  eta(k) = term1+term2
ELSE
  eta(k) = etan(k)
END IF

END DO

END SUBROUTINE shapiro


END MODULE sub