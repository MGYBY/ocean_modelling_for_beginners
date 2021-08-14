MODULE sub

USE param

CONTAINS

!=======================
SUBROUTINE init

hmin = 0.10

! grid parameters
dx = 1000.0
dy = 1000.0
dt = 20.0

! physical parameters
g = 9.81
f = -1.e-4
r = 0.0
rho = 1028.0

! parameter for semi-implicit treatment of Coriolis force
beta = 0.5*f*dt
beta = beta*beta

! choice of TVD advection scheme
!======================================
! mode = 1 => upstream scheme
! mode = 2 => Lax-Wendroff scheme
! mode = 3 => Superbee scheme
! mode = 4 => Super-C scheme
!=======================================

mode = 3

! read bathymetry file
OPEN(10,file ='topo.dat',form='formatted',recl=100000)
DO j = 0,ny+1
  READ(10,'(153F12.6)')(hzero(j,k),k=0,nx+1)
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  eta(j,k) = -MIN(0.0,hzero(j,k))
  etan(j,k) = eta(j,k)
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  h(j,k) = hzero(j,k)+eta(j,k)
! wet = 1 defines "wet" grid cells
! wet = 0 defines "dry" grid cells
  wet(j,k) = 1
  if(h(j,k) < hmin) wet(j,k) = 0
  u(j,k) = 0.
  un(j,k) = 0.
  v(j,k) = 0.
  vn(j,k) = 0.
END DO
END DO

! output of initial eta distribution
OPEN(10,file ='eta0.dat',form='formatted',recl=1000000)
  DO j = 1,ny
    WRITE(10,'(151F12.6)')(eta(j,k),k=1,nx)
  END DO
CLOSE(10)

! output of initial layer thickness distribution
OPEN(10,file ='h0.dat',form='formatted',recl=10000000)
  DO j = 1,ny
    WRITE(10,'(151F12.6)')(hzero(j,k),k=1,nx)
  END DO
CLOSE(10)

END SUBROUTINE init

!================
SUBROUTINE dyn

! local parameters
REAL :: du(0:ny+1,0:nx+1), dv(0:ny+1,0:nx+1)
REAL :: ustar(0:ny+1,0:nx+1), vstar(0:ny+1,0:nx+1)
REAL :: Rx(0:ny+1,0:nx+1), Ry(0:ny+1,0:nx+1), rr
REAL :: uu, vv, duu, dvv, speed, uv, vu
REAL :: pgrdx, pgrdy, corx, cory, tx, ty
REAL :: div, div1, div2
REAL :: advx(0:ny+1,0:nx+1), advy(0:ny+1,0:nx+1)
REAL :: hu, hv


! calculate the nonlinear terms for u-momentum equation
DO j = 0,ny+1
DO k = 0,nx
  CuP(j,k) = 0.25*(u(j,k)+u(j,k+1)+2.0*UGEO+abs(u(j,k)+UGEO)+abs(u(j,k+1)+UGEO))*dt/dx
  CuN(j,k) = 0.25*(u(j,k)+u(j,k+1)+2.0*UGEO-abs(u(j,k)+UGEO)-abs(u(j,k+1)+UGEO))*dt/dx
  CvP(j,k) = 0.25*(v(j,k)+v(j,k+1)+abs(v(j,k))+abs(v(j,k+1)))*dt/dy
  CvN(j,k) = 0.25*(v(j,k)+v(j,k+1)-abs(v(j,k))-abs(v(j,k+1)))*dt/dy
  Cu(j,k) = 0.5*abs(u(j,k)+u(j,k+1)+2.0*UGEO)*dt/dx
  Cv(j,k) = 0.5*abs(v(j,k)+v(j,k+1))*dt/dy 
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  B(j,k) = u(j,k)+UGEO
END DO
END DO

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
  CuP(j,k) = 0.25*(u(j,k)+u(j+1,k)+2.0*UGEO+abs(u(j,k)+UGEO)+abs(u(j+1,k)+UGEO))*dt/dx
  CuN(j,k) = 0.25*(u(j,k)+u(j+1,k)+2.0*UGEO-abs(u(j,k)+UGEO)-abs(u(j+1,k)+UGEO))*dt/dx
  CvP(j,k) = 0.25*(v(j,k)+v(j+1,k)+abs(v(j,k))+abs(v(j+1,k)))*dt/dy
  CvN(j,k) = 0.25*(v(j,k)+v(j+1,k)-abs(v(j,k))-abs(v(j+1,k)))*dt/dy
  Cu(j,k) = 0.5*abs(u(j,k)+u(j+1,k)+2.0*UGEO)*dt/dx
  Cv(j,k) = 0.5*abs(v(j,k)+v(j+1,k))*dt/dy
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  B(j,k) = v(j,k)
END DO
END DO

CALL advect

DO j = 1,ny
DO k = 1,nx
  div1 = 0.5*(u(j,k)+u(j+1,k)-u(j,k-1)-u(j+1,k-1))/dx
  div2 = 0.5*(v(j+1,k)-v(j-1,k))/dy
  div = dt*B(j,k)*(div1+div2)  
  advy(j,k)= BN(j,k) + div
END DO
END DO

! Step 1: prediction of u & v without Coriolis force

DO j = 1,ny
DO k = 1,nx
  rr = r
!
! Strong friction is required (applied to velocity fluctuations only) 
! to eliminate the built-up of wave disturbances near the right boundary
! 
  IF(k > nx-50)rr = 50.0*REAL(k-nx+50)/50.0
!
  pgrdx = -dt*g*(eta(j,k+1)-eta(j,k))/dx
  tx = 0.0
  Rx(j,k) = 1.0
  hu = 0.5*(h(j,k)+h(j,k+1))
  uu = u(j,k)+UGEO
  vu = 0.25*(v(j,k)+v(j,k+1)+v(j-1,k)+v(j-1,k+1))
  speed = SQRT(uu*uu+vu*vu)
  IF(hu > 0.0)THEN
    tx = dt*taux/(rho*hu)
    Rx(j,k) = 1.0+dt*rr*speed/hu 
  END IF

  ustar(j,k) = u(j,k)+ pgrdx + tx + advx(j,k)

  pgrdy = -dt*g*(eta(j+1,k)-eta(j,k))/dy
  ty = 0.0
  Ry(j,k) = 1.0
  hv = 0.5*(h(j+1,k)+h(j,k))
  vv = v(j,k)
  uv = 0.25*(u(j,k)+u(j+1,k)+u(j,k-1)+u(j+1,k-1))
  uv = uv+UGEO
  speed = SQRT(uv*uv+vv*vv)
  IF(hv > 0.0)THEN
    ty = dt*tauy/(rho*hv)
    Ry(j,k) = 1.0+dt*rr*speed/hv 
  END IF

  vstar(j,k) = v(j,k) + pgrdy + ty + advy(j,k)

END DO
END DO

! Step 2: Semi-implicit treatment of Coriolis force

DO j = 1,ny
DO k = 1,nx

vu = 0.25*(v(j,k)+v(j,k+1)+v(j-1,k)+v(j-1,k+1))
corx = dt*f*vu
du(j,k) = (ustar(j,k)-beta*u(j,k)+corx)/(1.0+beta)-u(j,k)

uv = 0.25*(u(j,k)+u(j+1,k)+u(j,k-1)+u(j+1,k-1))
cory = -dt*f*uv
dv(j,k) = (vstar(j,k)-beta*v(j,k)+cory)/(1.0+beta)-v(j,k)

END DO
END DO

! Step 3: final prediction of u & v including the flooding algorithm

DO j = 1,ny
DO k = 1,nx

un(j,k) = 0.0
uu = u(j,k)
duu = du(j,k)
IF(wet(j,k)==1) THEN
  IF((wet(j,k+1)==1).or.(duu>0.0)) un(j,k) = (uu+duu)/Rx(j,k)
ELSE
  IF((wet(j,k+1)==1).and.(duu<0.0)) un(j,k) = (uu+duu)/Rx(j,k)
END IF

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

! flow disturbances vanish at inflow boundary
DO j = 1,ny
 un(j,1) = 0.0
 vn(j,1) = 0.0
END DO

! zero-gradient boundary conditions
DO j = 0,ny+1
  un(j,0) = un(j,1)
  un(j,nx+1) = un(j,nx)
  vn(j,0) = vn(j,1)
  vn(j,nx+1) = vn(j,nx)
END DO

DO k = 0,nx+1
  un(0,k) = un(1,k)
  un(ny+1,k) = un(ny,k)
  vn(0,k) = vn(1,k)
  vn(ny+1,k) = vn(ny,k)
END DO

! Eulerian tracer predictor
DO j = 0,ny+1
DO k = 0,nx+1
  CuP(j,k) = 0.5*(u(j,k)+UGEO+abs(u(j,k)+UGEO))*dt/dx
  CuN(j,k) = 0.5*(u(j,k)+UGEO-abs(u(j,k)+UGEO))*dt/dx
  CvP(j,k) = 0.5*(v(j,k)+abs(v(j,k)))*dt/dy
  CvN(j,k) = 0.5*(v(j,k)-abs(v(j,k)))*dt/dy
  Cu(j,k) = abs(u(j,k)+UGEO)*dt/dx 
  Cv(j,k) = abs(v(j,k))*dt/dy 
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  B(j,k) = T(j,k)
END DO
END DO

CALL advect

DO j = 1,ny
DO k = 1,nx
  div1 = (u(j,k)-u(j,k-1))/dx
  div2 = (v(j,k)-v(j-1,k))/dy
  div = dt*B(j,k)*(div1+div2)  
  TN(j,k)= T(j,k)+BN(j,k)+div
END DO
END DO

! boundary conditions
DO j = 0,ny+1
  TN(j,0) = TN(j,1)
  TN(j,nx+1) = TN(j,nx)
END DO

DO k = 0,nx+1
  TN(0,k) = TN(1,k)
  TN(ny+1,k) = TN(ny,k)
END DO


! sea level predictor

DO j = 0,ny+1
DO k = 0,nx+1
  CuP(j,k) = 0.5*(un(j,k)+UGEO+abs(un(j,k)+UGEO))*dt/dx
  CuN(j,k) = 0.5*(un(j,k)+UGEO-abs(un(j,k)+UGEO))*dt/dx
  CvP(j,k) = 0.5*(vn(j,k)+abs(vn(j,k)))*dt/dy
  CvN(j,k) = 0.5*(vn(j,k)-abs(vn(j,k)))*dt/dy
  Cu(j,k) = abs(un(j,k)+UGEO)*dt/dx 
  Cv(j,k) = abs(vn(j,k))*dt/dy 
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  B(j,k) = h(j,k)
END DO
END DO

CALL advect

DO j = 1,ny
DO k = 1,nx
  etan(j,k) = eta(j,k)+BN(j,k)
END DO
END DO

! damping layer
! elimination of sea-level anomalies near upstream boundary
DO j = 1,ny
DO k = 1,10
  etan(j,k) = etan(j,k)*( 1.0-1.0*REAL(11-k)/10. )
END DO
END DO

DO j = 1,ny
  etan(j,0) = 0.0
  etan(j,nx+1) = etan(j,nx) 
END DO

DO k = 0,nx+1
  etan(0,k) = etan(1,k)
  etan(ny+1,k) = etan(ny,k)
END DO

! updating of sea level
DO j = 0,ny+1
DO k = 0,nx+1
  eta(j,k) = etan(j,k)
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  h(j,k) = hzero(j,k) + eta(j,k)
  wet(j,k) = 1
  IF(h(j,k)<hmin)wet(j,k) = 0
  u(j,k) = un(j,k)
  v(j,k) = vn(j,k)
  T(j,k) = TN(j,k)
END DO
END DO

END SUBROUTINE dyn

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

END SUBROUTINE advect

REAL FUNCTION psi(r2,cfl,mmode)

! input parameters

REAL, INTENT(IN) :: r2, cfl  
INTEGER, INTENT(IN) :: mmode

! local parameters 
REAL :: term1, term2, term3

IF(mmode == 1) psi = 0. 
IF(mmode == 2) psi = 1.
IF(mmode == 3)THEN
  term1 = MIN(2.0*r2,1.0)
  term2 = MIN(r2,2.0)
  term3 = MAX(term1,term2)
  psi = MAX(term3,0.0)
END IF
IF(mmode == 4)THEN
  psi = 0.0
  IF(r2 > 0.0)THEN
    IF(r2 > 1.0) THEN
      psi = MIN(r2, 2.0/(1.-cfl))
    ELSE
      IF(cfl > 0.0) psi = MIN(2.0*r2/cfl,1)
    END IF
  END IF
END IF

RETURN

END FUNCTION psi 

END MODULE sub