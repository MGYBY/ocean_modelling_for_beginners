MODULE sub
USE param

CONTAINS

!=======================
SUBROUTINE init

! grid parameters
dx = 10.0
dt = 0.1

! physical parameters
g = 9.81

! initial conditions
DO k = 0,nx+1
  hzero(k) = 10.0 
  eta(k) = 0.0
  etan(k) = 0.0
END DO

! specify closed channel
hzero(0) = 0.
hzero(nx+1) = 0.

DO k = 0,nx+1
  h(k) = hzero(k)+eta(k)
  wet(k) = 1
  if(h(k) <= 0.) wet(k) = 0
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

! velocity predictor for wet grid cells
IF(wet(k) == 1) THEN
  IF (wet(k+1) == 1) THEN
     pgradx = -g*(eta(k+1)-eta(k))/dx
     un(k) = u(k)+dt*pgradx
  ELSE
     un(k) = 0.
  END IF       
END IF

END DO

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