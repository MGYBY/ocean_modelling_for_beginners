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
g = 9.81
r = 1.e-3
rho = 1028.0

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
RETURN
END SUBROUTINE init

!================
SUBROUTINE dyn

! local parameters
REAL :: du(0:ny+1,0:nx+1), dv(0:ny+1,0:nx+1)
REAL :: Rx(0:ny+1,0:nx+1), Ry(0:ny+1,0:nx+1)
REAL :: uu, vv, duu, dvv, speed, uv, vu
REAL :: pgrdx, pgrdy, tx, ty
REAL :: div, div1, div2
REAL :: advx(0:ny+1,0:nx+1), advy(0:ny+1,0:nx+1)

! calculate the nonlinear terms for u-momentum equation
DO j = 0,ny+1
DO k = 0,nx
  CuP(j,k) = 0.25*(u(j,k)+u(j,k+1)+abs(u(j,k))+abs(u(j,k+1)))*dt/dx
  CuN(j,k) = 0.25*(u(j,k)+u(j,k+1)-abs(u(j,k))-abs(u(j,k+1)))*dt/dx
  CvP(j,k) = 0.25*(v(j,k)+v(j,k+1)+abs(v(j,k))+abs(v(j,k+1)))*dt/dy
  CvN(j,k) = 0.25*(v(j,k)+v(j,k+1)-abs(v(j,k))-abs(v(j,k+1)))*dt/dy
  Cu(j,k) = 0.5*abs(u(j,k)+u(j,k+1))*dt/dx
  Cv(j,k) = 0.5*abs(v(j,k)+v(j,k+1))*dt/dy
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  B(j,k)=u(j,k)
END DO
END DO

mode = 3
CALL advect

DO j = 1,ny
DO k = 1,nx
  div1 = 0.5*(u(j,k+1)-u(j,k-1))/dx
  div2 = 0.5*(v(j,k)+v(j,k+1)-v(j-1,k)-v(j-1,k+1))/dy
  div = dt*B(j,k)*(div1+div2)  
  advx(j,k)= BN(j,k)+div
END DO
END DO

! calculate the nonlinear terms for v-momentum equation
DO j = 0,ny
DO k = 0,nx+1
  CuP(j,k) = 0.25*(u(j,k)+u(j+1,k)+abs(u(j,k))+abs(u(j+1,k)))*dt/dx
  CuN(j,k) = 0.25*(u(j,k)+u(j+1,k)-abs(u(j,k))-abs(u(j+1,k)))*dt/dx
  CvP(j,k) = 0.25*(v(j,k)+v(j+1,k)+abs(v(j,k))+abs(v(j+1,k)))*dt/dy
  CvN(j,k) = 0.25*(v(j,k)+v(j+1,k)-abs(v(j,k))-abs(v(j+1,k)))*dt/dy
  Cu(j,k) = 0.5*abs(u(j,k)+u(j+1,k))*dt/dx
  Cv(j,k) = 0.5*abs(v(j,k)+v(j+1,k))*dt/dy
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  B(j,k)=v(j,k)
END DO
END DO

mode = 3
CALL advect

DO j = 1,ny
DO k = 1,nx
  div1 = 0.5*(u(j,k)+u(j+1,k)-u(j,k-1)-u(j+1,k-1))/dx
  div2 = 0.5*(v(j+1,k)-v(j-1,k))/dy
  div = dt*B(j,k)*(div1+div2)  
  advy(j,k)= BN(j,k) + div
END DO
END DO

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
  du(j,k) = pgrdx + tx + advx(j,k)

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
  dv(j,k) = pgrdy + ty + advy(j,k)

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
! calculate the nonlinear terms for v-momentum equation
DO j = 0,ny+1
DO k = 0,nx+1
  CuP(j,k) = 0.5*(un(j,k)+abs(un(j,k)))*dt/dx
  CuN(j,k) = 0.5*(un(j,k)-abs(un(j,k)))*dt/dx
  CvP(j,k) = 0.5*(vn(j,k)+abs(vn(j,k)))*dt/dy
  CvN(j,k) = 0.5*(vn(j,k)-abs(vn(j,k)))*dt/dy
  Cu(j,k) = abs(un(j,k))*dt/dx
  Cv(j,k) = abs(vn(j,k))*dt/dy
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
B(j,k)=h(j,k)
END DO
END DO

mode = 3
CALL advect

DO j = 1,ny
DO k = 1,nx
  
   etan(j,k) = eta(j,k)+BN(j,k)

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

SUBROUTINE advect
! local parameters
REAL :: RxP(0:ny+1,0:nx+1), RxN(0:ny+1,0:nx+1)
REAL :: RyP(0:ny+1,0:nx+1), RyN(0:ny+1,0:nx+1)
REAL :: dB, term1, term2, term3, term4
REAL :: BwP, BwN, BeP, BeN, BsP, BsN, BnP, BnN 

DO j = 0,ny+1
DO k = 0,nx+1
  RxP(j,k) = 0.0
  RxN(j,k) = 0.0
  RyP(j,k) = 0.0
  RyN(j,k) = 0.0
END DO
END DO

DO j = 1,ny
DO k = 1,nx
  dB =  B(j,k+1)-B(j,k)
  IF(ABS(dB) > 0.0) RxP(j,k) = (B(j,k)-B(j,k-1))/dB
  dB =  B(j+1,k)-B(j,k)
  IF(ABS(dB) > 0.0) RyP(j,k) = (B(j,k)-B(j-1,k))/dB
END DO
END DO

DO j = 1,ny
DO k = 0,nx-1
  dB =  B(j,k+1)-B(j,k)
  IF(ABS(dB) > 0.0) RxN(j,k) = (B(j,k+2)-B(j,k+1))/dB
END DO
END DO

DO j = 0,ny-1
DO k = 1,nx
  dB =  B(j+1,k)-B(j,k)
  IF(ABS(dB) > 0.0) RyN(j,k) = (B(j+2,k)-B(j+1,k))/dB
END DO
END DO   

DO j = 1,ny
DO k = 1,nx

BwP = B(j,k-1)+0.5*PSI(RxP(j,k-1),Cu(j,k-1),mode)*(1.0-CuP(j,k-1))*(B(j,k)-B(j,k-1))
BwN = B(j,k)-0.5*PSI(RxN(j,k-1),Cu(j,k-1),mode)*(1.0+CuN(j,k-1))*(B(j,k)-B(j,k-1))
BeP = B(j,k)+0.5*PSI(RxP(j,k),Cu(j,k),mode)*(1.0-CuP(j,k))*(B(j,k+1)-B(j,k))
BeN = B(j,k+1)-0.5*PSI(RxN(j,k),Cu(j,k),mode)*(1.0+CuN(j,k))*(B(j,k+1)-B(j,k))

BsP = B(j-1,k)+0.5*PSI(RyP(j-1,k),Cv(j-1,k),mode)*(1.0-CvP(j-1,k))*(B(j,k)-B(j-1,k))
BsN = B(j,k)-0.5*PSI(RyN(j-1,k),Cv(j-1,k),mode)*(1.0+CvN(j-1,k))*(B(j,k)-B(j-1,k))
BnP = B(j,k)+0.5*PSI(RyP(j,k),Cv(j,k),mode)*(1.0-CvP(j,k))*(B(j+1,k)-B(j,k))
BnN = B(j+1,k)-0.5*PSI(RyN(j,k),Cv(j,k),mode)*(1.0+CvN(j,k))*(B(j+1,k)-B(j,k))

term1 = CuP(j,k-1)*BwP+CuN(j,k-1)*BwN
term2 = CuP(j,k)*BeP+CuN(j,k)*BeN
term3 = CvP(j-1,k)*BsP+CvN(j-1,k)*BsN
term4 = CvP(j,k)*BnP+CvN(j,k)*BnN

BN(j,k) = term1-term2+term3-term4

END DO
END DO
RETURN
END SUBROUTINE advect

REAL FUNCTION psi(r,cfl,mmode)

! input parameters

REAL, INTENT(IN) :: r, cfl  
INTEGER, INTENT(IN) :: mmode

! local parameters 
REAL :: term1, term2, term3

IF(mmode == 1) psi = 0. 
IF(mmode == 2) psi = 1.
IF(mmode == 3)THEN
  term1 = MIN(2.0*r,1.0)
  term2 = MIN(r,2.0)
  term3 = MAX(term1,term2)
  psi = MAX(term3,0.0)
END IF
IF(mmode == 4)THEN
  psi = 0.0
  IF(r > 0.0)THEN
    IF(r > 1.0) THEN
      psi = MIN(r, 2.0/(1.-cfl))
    ELSE
      IF(cfl > 0.0) psi = MIN(2.0*r/cfl,1)
    END IF
  END IF
END IF

RETURN
END FUNCTION psi 

END MODULE sub