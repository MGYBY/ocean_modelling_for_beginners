MODULE sub
USE param

CONTAINS
!=======================
SUBROUTINE init
REAL :: htot, hini(nz)

hmin = 0.01 ! minimum layer thickness
dx = 10.0 ! grid spacing
g = 9.81 ! acceleration due to gravity

! bathymetry
DO k = 1,nx
htotal(k) = 100.0
END DO

! triangle-shaped island 
DO k = 31,51
 htotal(k) = 100.-95.*real(k-30)/21.
END DO

DO k = 52,71
 htotal(k) = 100.-95.*real(71-k+1)/20.
END DO

htotal(0) = -10.0
htotal(nx+1) = -10.0

! undisturbed layer thicknesses & interface displacements
DO i = 1,nz
hini(i) = 10.0
END DO
DO k = 0,nx+1
 htot = htotal(k)
 DO i = 1,nz
 hzero(i,k) = MAX(MIN(hini(i),htot),0.0)
 eta(i,k) = MAX(0.0,-htot)
 htot = htot - hini(1)
 END DO
END DO

! layer densities
rho(0) = 0.0 ! air density ignored
rho(1) = 1025.0
DO i = 2,nz
  rho(i) = 1026.0+REAL(i-2)/REAL(nz-2)*0.5
END DO

! boundary values for dp and eta
DO k = 0,nx+1
  dp(0,k) = 0.0 ! air pressure ignored
  eta(nz+1,k) = 0.0 ! sea floor is rigid
END DO

! store initial interface displacements
DO k = 0,nx+1
DO i = 1,nz+1
  eta0(i,k) = eta(i,k)
END DO
END DO

! layer thicknesses, wet\dry pointers and velocities
DO i = 1,nz
DO k = 0,nx+1
  h(i,k) = hzero(i,k)
  wet(i,k) = 1
  if(h(i,k) < hmin) wet(i,k) = 0
  u(i,k) = 0.
  un(i,k) = 0.
END DO
END DO

END SUBROUTINE init

!================
SUBROUTINE dyn

! local parameters
REAL :: pgradx 
REAL :: hue, huw, hwp, hwn, hen, hep
REAL :: deta

! calculate dynamic pressure
  DO k = 0,nx+1
  DO i = 1,nz
    dp(i,k) = dp(i-1,k)+(rho(i)-rho(i-1))*g*eta(i,k)  
  END DO
  END DO

DO k = 1,nx
DO i = 1,nz
 
! velocity predictor for wet grid cells
  pgradx = -(dp(i,k+1)-dp(i,k))/rho(i)/dx
  un(i,k) = 0.0
  IF(wet(i,k)==1) THEN
    IF((wet(i,k+1)==1).or.(pgradx>0.0)) un(i,k) = u(i,k)+dt*pgradx
  ELSE
    IF((wet(i,k+1)==1).and.(pgradx<0.0)) un(i,k) = u(i,k)+dt*pgradx
  END IF
END DO
END DO

! layer-thickness change predictor
DO k = 1,nx
DO i = 1,nz
  hep = 0.5*(un(i,k)+abs(un(i,k)))*h(i,k)
  hen = 0.5*(un(i,k)-abs(un(i,k)))*h(i,k+1)
  hue = hep+hen
  hwp = 0.5*(un(i,k-1)+abs(un(i,k-1)))*h(i,k-1)
  hwn = 0.5*(un(i,k-1)-abs(un(i,k-1)))*h(i,k)
  huw = hwp+hwn
  dhdt(i,k) = -(hue-huw)/dx
END DO
END DO

! update interface displacements
DO k = 1,nx
  deta = 0.0
  DO i = nz,1,-1
    deta = deta + dhdt(i,k)
    etan(i,k)= eta(i,k)+dt*deta
  END DO
END DO

! apply Shapiro filter
CALL shapiro

! update layer thicknesses, lateral velocities and wet/dry pointers
DO k = 1,nx
   DO i = 1,nz
    h(i,k) = hzero(i,k)+eta(i,k)-eta(i+1,k)-eta0(i,k)+eta0(i+1,k)
    u(i,k) = un(i,k)
    wet(i,k) = 1
    if(h(i,k)<hmin)wet(i,k) = 0
  END DO
END DO

END SUBROUTINE dyn

!======================
SUBROUTINE shapiro

!local parameters
REAL :: term1,term2

! 1-order Shapiro filter 

DO i = 1,nz
DO k = 1,nx

IF(wet(i,k)==1)THEN
  term1 = (1.0-0.5*eps*(wet(i,k+1)+wet(i,k-1)))*etan(i,k)
  term2 = 0.5*eps*(wet(i,k+1)*etan(i,k+1)+wet(i,k-1)*etan(i,k-1))
  eta(i,k) = term1+term2
ELSE
  eta(i,k) = etan(i,k)
END IF

END DO
END DO

END SUBROUTINE shapiro


END MODULE sub