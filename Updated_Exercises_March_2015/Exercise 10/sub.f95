MODULE sub
USE param

CONTAINS

!=======================
SUBROUTINE init

hmin = 0.05

! grid parameters
dx = 100.0
dy = 100.0
dt = 3.0

! physical parameters
g = 9.81 ! acceleration due to gravity 
r = 1.e-3 ! bottom-drag coefficient
rho = 1028.0 ! mean seawater density

! read bathymetry file
OPEN(10,file ='topo.dat',form='formatted')
DO j = 1,ny
  READ(10,'(51F12.6)')(hzero(j,k),k=1,nx)
END DO

! land boundaries 
DO k = 0,nx+1
 hzero(0,k) = -10.0
 hzero(ny+1,k) = -10.0
END DO

DO j = 0,ny+1
 hzero(j,0) = -10.0
 hzero(j,nx+1) = -10.0
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  eta(j,k) = -MIN(0.0,hzero(j,k))
  etan(j,k) = eta(j,k)
END DO
END DO
!XXXXXXXXXXXXXXXXXXX

DO j = 0,ny+1
DO k = 0,nx+1
  h(j,k) = hzero(j,k)+eta(j,k)
! wet = 1 defines "wet" grid cells
! wet = 0 defines "dry" grid cells
  wet(j,k) = 1
  if(h(j,k) < hmin)wet(j,k) = 0
  u(j,k) = 0.
  un(j,k) = 0.
  v(j,k) = 0.
  vn(j,k) = 0.
END DO
END DO

END SUBROUTINE init

!================
SUBROUTINE dyn

! local parameters
REAL :: du(0:ny+1,0:nx+1), dv(0:ny+1,0:nx+1)
REAL :: Rx(0:ny+1,0:nx+1), Ry(0:ny+1,0:nx+1)
REAL :: uu, vv, duu, dvv, speed, uv, vu
REAL :: pgrdx, pgrdy, corx, cory, tx, ty, advx, advy
REAL :: hue, huw, hwp, hwn, hen, hep
REAL :: hvn, hvs, hsp, hsn, hnn, hnp

DO j = 1,ny
DO k = 1,nx
  pgrdx = -dt*g*(eta(j,k+1)-eta(j,k))/dx
  tx = 0.0
  Rx(j,k) = 1.0
  hu = 0.5*(h(j,k)+h(j,k+1))
  uu = u(j,k)
  vu = 0.25*(v(j,k)+v(j,k+1)+v(j-1,k)+v(j-1,k+1))
  speed = SQRT(uu*uu+vu*vu)
  IF(hu > 0.0)THEN
    tx = dt*taux/(rho*hu)
    Rx(j,k) = 1.0+dt*r*speed/hu 
  END IF
  du(j,k) = pgrdx + tx

  pgrdy = -dt*g*(eta(j+1,k)-eta(j,k))/dy
  ty = 0.0
  Ry(j,k) = 1.0
  hv = 0.5*(h(j+1,k)+h(j,k))
  vv = v(j,k)
  uv = 0.25*(u(j,k)+u(j+1,k)+u(j,k-1)+u(j+1,k-1))
  speed = SQRT(uv*uv+vv*vv)
  IF(hv > 0.0)THEN
    ty = dt*tauy/(rho*hv)
    Ry(j,k) = 1.0+dt*r*speed/hv 
  END IF
  dv(j,k) = pgrdy + ty

END DO
END DO

DO j = 1,ny
DO k = 1,nx

! prediction for u
un(j,k) = 0.0
uu = u(j,k)
duu = du(j,k)
IF(wet(j,k)==1) THEN
  IF((wet(j,k+1)==1).or.(duu>0.0)) un(j,k) = (uu+duu)/Rx(j,k)
ELSE
  IF((wet(j,k+1)==1).and.(duu<0.0)) un(j,k) = (uu+duu)/Rx(j,k)
END IF

! prediction for v
vv = v(j,k)
dvv = dv(j,k)
vn(j,k) = 0.0
IF(wet(j,k)==1) THEN
  IF((wet(j+1,k)==1).or.(dvv>0.0)) vn(j,k) = (vv+dvv)/Ry(j,k)
ELSE
  IF((wet(j+1,k)==1).and.(dvv<0.0)) vn(j,k) = (vv+dvv)/Ry(j,k)
END IF

END DO
END DO

! sea level predictor
DO j = 1,ny
DO k = 1,nx
  
 hep = 0.5*(un(j,k)+abs(un(j,k)))*h(j,k)
  hen = 0.5*(un(j,k)-abs(un(j,k)))*h(j,k+1)
  hue = hep+hen
  hwp = 0.5*(un(j,k-1)+abs(un(j,k-1)))*h(j,k-1)
  hwn = 0.5*(un(j,k-1)-abs(un(j,k-1)))*h(j,k)
  huw = hwp+hwn
  
  hnp = 0.5*(vn(j,k)+abs(vn(j,k)))*h(j,k)
  hnn = 0.5*(vn(j,k)-abs(vn(j,k)))*h(j+1,k)
  hvn = hnp+hnn
  hsp = 0.5*(vn(j-1,k)+abs(vn(j-1,k)))*h(j-1,k)
  hsn = 0.5*(vn(j-1,k)-abs(vn(j-1,k)))*h(j,k)
  hvs = hsp+hsn
  etan(j,k) = eta(j,k)-dt*(hue-huw)/dx-dt*(hvn-hvs)/dy

END DO
END DO

END SUBROUTINE dyn

!======================
SUBROUTINE shapiro

!local parameters
REAL :: term1,term2,term3

! 1-order Shapiro filter

DO j = 1,ny
DO k = 1,nx

IF(wet(j,k)==1)THEN
  term1 = (1.0-0.25*eps*(wet(j,k+1)+wet(j,k-1)+ wet(j+1,k)+wet(j-1,k)))*etan(j,k)
  term2 = 0.25*eps*(wet(j,k+1)*etan(j,k+1)+wet(j,k-1)*etan(j,k-1))
  term3 = 0.25*eps*(wet(j+1,k)*etan(j+1,k)+wet(j-1,k)*etan(j-1,k))
  eta(j,k) = term1+term2+term3
ELSE
  eta(j,k) = etan(j,k)
END IF

END DO
END DO

END SUBROUTINE shapiro

END MODULE sub