MODULE sub
USE param

CONTAINS

!******************
SUBROUTINE init
! local parameters
REAL :: htot, hini(nz)

hmin = 0.1 ! minimum layer thickness

! grid parameters
dx = 10000.0
dy = 10000.0
dt = 10.0

! physical parameters
g = 9.81

! variable Coriolis parameter & wind-stress forcing
DO j = 0,ny+1
  f(j) = 1.0e-4 + 2.0e-5*REAL(j)/REAL(ny+1)
  taux(j) = -0.05*COS(REAL(j-1)/REAL(ny-1)*PI)
END DO
tauy = 0.0 

r = 1.e-3  ! bottom-drag coefficient
ah = 500.0 ! horizontal eddy viscosity

!**** choice of lateral friction scheme
! no-slip condition => slip = 2
! semi-slip condition => slip = 1
! full-slip condition => slip = 0

slip = 2.0

!*** choice of TVD advection scheme
mode = 3 ! Superbee limiter

! set bathymetry
DO j = 0,ny+1
DO k = 0,nx+1
  htotal(j,k) = 1000.0 
END DO
END DO
 
! land boundaries
DO k = 0,nx+1
  htotal(0,k) = -10.0
  htotal(ny+1,k) = -10.0
END DO

DO j = 0,ny+1
  htotal(j,0) = -10.0
  htotal(j,nx+1) = -10.0
END DO

! undisturbed layer thicknesses & interface displacements
hini(1) = 200.0
hini(2) = 800.0

DO j = 0,ny+1
DO k = 0,nx+1
 htot = htotal(j,k)
 DO i = 1,nz
 hzero(i,j,k) = MAX(AMIN1(hini(i),htot),0.0)
 eta(i,j,k) = MAX(0.0,-htot)
 htot = htot - hini(i)
 END DO
END DO
END DO

! layer densities
rho(0) = 0.0 ! air density ignored
rho(1) = 1025.0
rho(2) = 1030.0

! boundary values for dp and eta
DO j = 0,ny+1
DO k = 0,nx+1
  dp(0,j,k) = 0.0 ! air pressure ignored
  eta(nz+1,j,k) = 0.0 ! sea floor is rigid
END DO
END DO

! store initial interface displacements
DO i = 1,nz+1
DO j = 0,ny+1
DO k = 0,nx+1
  eta0(i,j,k) = eta(i,j,k)
END DO
END DO
END DO

! layer thicknesses, wet\dry pointers and velocities
DO i = 1,nz
DO j = 0,ny+1
DO k = 0,nx+1
  h(i,j,k) = hzero(i,j,k)
  wet(i,j,k) = 1
  if(h(i,j,k) < hmin) wet(i,j,k) = 0
  u(i,j,k) = 0.
  un(i,j,k) = 0.
  v(i,j,k) = 0.
  vn(i,j,k) = 0.
END DO
END DO
END DO

! output of bathymetry distribution
OPEN(10,file ='h0.dat',form='formatted')
  DO j = 1,ny
    WRITE(10,'(101F12.6)')(htotal(j,k),k=1,nx)
  END DO
CLOSE(10)

END SUBROUTINE init

!*******************
SUBROUTINE dyn

! local parameters
REAL :: du(nz,0:ny+1,0:nx+1), dv(nz,0:ny+1,0:nx+1)
REAL :: ustar(nz,0:ny+1,0:nx+1), vstar(nz,0:ny+1,0:nx+1)
REAL :: Rx(nz,0:ny+1,0:nx+1), Ry(nz,0:ny+1,0:nx+1)
REAL :: uu, vv, duu, dvv, speed, uv, vu, vh
REAL :: pgrdx, pgrdy, corx, cory, tx, ty, fm
REAL :: div, div1, div2, deta
REAL :: advx(nz,0:ny+1,0:nx+1), advy(nz,0:ny+1,0:nx+1)
REAL :: diffu, diffv, term1, term2, term3, term4, hu, hv
REAL :: s1, s2, h1

! calculate the nonlinear terms for u-momentum equation
!===========
DO i = 1,nz
!===========

DO j = 0,ny+1
DO k = 0,nx
  CuP(j,k) = 0.25*(u(i,j,k)+u(i,j,k+1)+abs(u(i,j,k))+abs(u(i,j,k+1)))*dt/dx
  CuN(j,k) = 0.25*(u(i,j,k)+u(i,j,k+1)-abs(u(i,j,k))-abs(u(i,j,k+1)))*dt/dx
  CvP(j,k) = 0.25*(v(i,j,k)+v(i,j,k+1)+abs(v(i,j,k))+abs(v(i,j,k+1)))*dt/dy
  CvN(j,k) = 0.25*(v(i,j,k)+v(i,j,k+1)-abs(v(i,j,k))-abs(v(i,j,k+1)))*dt/dy
  Cu(j,k) = 0.5*abs(u(i,j,k)+u(i,j,k+1))*dt/dx
  Cv(j,k) = 0.5*abs(v(i,j,k)+v(i,j,k+1))*dt/dy
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  B(j,k) = u(i,j,k)
END DO
END DO

CALL advect

DO j = 1,ny
DO k = 1,nx
  div1 = 0.5*(u(i,j,k+1)-u(i,j,k-1))/dx
  div2 = 0.5*(v(i,j,k)+v(i,j,k+1)-v(i,j-1,k)-v(i,j-1,k+1))/dy
  div = dt*B(j,k)*(div1+div2)  
  advx(i,j,k)= BN(j,k)+div
END DO
END DO

!==========
END DO
!==========

! calculate the nonlinear terms for v-momentum equation
!===========
DO i = 1,nz
!===========

DO j = 0,ny
DO k = 0,nx+1
  CuP(j,k) = 0.25*(u(i,j,k)+u(i,j+1,k)+abs(u(i,j,k))+abs(u(i,j+1,k)))*dt/dx
  CuN(j,k) = 0.25*(u(i,j,k)+u(i,j+1,k)-abs(u(i,j,k))-abs(u(i,j+1,k)))*dt/dx
  CvP(j,k) = 0.25*(v(i,j,k)+v(i,j+1,k)+abs(v(i,j,k))+abs(v(i,j+1,k)))*dt/dy
  CvN(j,k) = 0.25*(v(i,j,k)+v(i,j+1,k)-abs(v(i,j,k))-abs(v(i,j+1,k)))*dt/dy
  Cu(j,k) = 0.5*abs(u(i,j,k)+u(i,j+1,k))*dt/dx
  Cv(j,k) = 0.5*abs(v(i,j,k)+v(i,j+1,k))*dt/dy
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  B(j,k) = v(i,j,k)
END DO
END DO

CALL advect

DO j = 1,ny
DO k = 1,nx
  div1 = 0.5*(u(i,j,k)+u(i,j+1,k)-u(i,j,k-1)-u(i,j+1,k-1))/dx
  div2 = 0.5*(v(i,j+1,k)-v(i,j-1,k))/dy
  div = dt*B(j,k)*(div1+div2)  
  advy(i,j,k)= BN(j,k) + div
END DO
END DO

!==========
END DO
!==========

!===========
DO i = 1,nz
!===========

! calculate dynamic pressure
DO j = 0,ny+1
DO k = 0,nx+1
  dp(i,j,k) = dp(i-1,j,k)+(rho(i)-rho(i-1))*g*eta(i,j,k)  
END DO
END DO

! Step 1: prediction of u & v without Coriolis force

DO j = 1,ny
DO k = 1,nx
  pgrdx = -dt*(dp(i,j,k+1)-dp(i,j,k))/rho(i)/dx
  tx = 0.0
  diffu = 0.0
  hu = 0.5*(h(i,j,k)+h(i,j,k+1))
  Rx(i,j,k) = 1.0
  uu = u(i,j,k)
  vu = 0.25*(v(i,j,k)+v(i,j,k+1)+v(i,j-1,k)+v(i,j-1,k+1))
  speed = SQRT(uu*uu+vu*vu)
  IF(hu > 0.0)THEN
    term1 = h(i,j,k+1)*(u(i,j,k+1)-u(i,j,k))/dx
    term2 = h(i,j,k)*(u(i,j,k)-u(i,j,k-1))/dx
! slip conditions
    s1 = 1.0
    s2 = 1.0
    h1 = h(i,j+1,k) + h(i,j+1,k+1) 
    IF( h1 < hmin ) THEN
      s1 = 0.0
      s2 = slip
      h1 = h(i,j,k) + h(i,j,k+1) 
    END IF
    term3 = 0.25*(h(i,j,k)+h(i,j,k+1)+h1)*(s1*u(i,j+1,k)-s2*u(i,j,k))/dy
    s1 = 1.0
    s2 = 1.0
    h1 = h(i,j-1,k) + h(i,j-1,k+1) 
    IF( h1 < hmin ) THEN
      s1 = 0.0
      s2 = slip
      h1 = h(i,j,k) + h(i,j,k+1) 
    END IF
    term4 = 0.25*(h(i,j,k)+h(i,j,k+1)+h1)*(s2*u(i,j,k)-s1*u(i,j-1,k))/dy
    diffu = dt*ah*((term1-term2)/dx+(term3-term4)/dy)/hu
    IF(i == 1) tx = ad*dt*taux(j)/(rho(1)*hu)
    IF(i == nz)Rx(i,j,k) = 1.0+dt*r*speed/hu 
  END IF

  ustar(i,j,k) = u(i,j,k) + pgrdx + tx + advx(i,j,k) + diffu

  pgrdy = -dt*(dp(i,j+1,k)-dp(i,j,k))/rho(i)/dy
  vv = v(i,j,k)
  uv = 0.25*(u(i,j,k)+u(i,j+1,k)+u(i,j,k-1)+u(i,j+1,k-1))
  ty = 0.0
  Ry(i,j,k) = 1.0
  hv = 0.5*(h(i,j+1,k)+h(i,j,k))
  speed = SQRT(uv*uv+vv*vv)
  diffv = 0.0
  IF(hv > 0.0)THEN
! slip conditions
    s1 = 1.0
    s2 = 1.0
    h1 = h(i,j,k+1)+h(i,j+1,k+1) 
    IF( h1 < hmin ) THEN
      s1 = 0.0
      s2 = slip
      h1 = h(i,j,k)+h(i,j+1,k)
    END IF
    term1 = 0.25*( h(i,j,k)+h(i,j+1,k)+ h1)*(s1*v(i,j,k+1)-s2*v(i,j,k))/dx
    s1 = 1.0
    s2 = 1.0
    h1 = h(i,j,k-1)+h(i,j+1,k-1) 
    IF( h1 < hmin ) THEN
      s1 = 0.0
      s2 = slip
      h1 = h(i,j,k)+h(i,j+1,k)
    END IF  
    term2 = 0.25*( h(i,j,k)+h(i,j+1,k)+h1)*(s2*v(i,j,k)-s1*v(i,j,k-1))/dx
    term3 = h(i,j+1,k)*(v(i,j+1,k)-v(i,j,k))/dy
    term4 = h(i,j,k)*(v(i,j,k)-v(i,j-1,k-1))/dy

    diffv = dt*ah*((term1-term2)/dx+(term3-term4)/dy)/hu  
    IF(i == 1)ty = dt*tauy/(rho(1)*hv)
    IF(i==nz)Ry(i,j,k) = 1.0+dt*r*speed/hv 
   END IF

   vstar(i,j,k) = v(i,j,k) + pgrdy + ty + advy(i,j,k) + diffv 
END DO
END DO

! Step 2: Semi-implicit treatment of Coriolis force

DO j = 1,ny
DO k = 1,nx

fm = f(j)
beta = 0.5*fm*dt
beta = beta*beta

vu = 0.25*(v(i,j,k)+v(i,j,k+1)+v(i,j-1,k)+v(i,j-1,k+1))
corx = dt*fm*vu

du(i,j,k) = (ustar(i,j,k)-beta*u(i,j,k)+corx)/(1.0+beta)-u(i,j,k)

fm = 0.5*(f(j)+f(j+1))
beta = 0.5*fm*dt
beta = beta*beta

uv = 0.25*(u(i,j,k)+u(i,j+1,k)+u(i,j,k-1)+u(i,j+1,k-1))
cory = -dt*fm*uv

dv(i,j,k) = (vstar(i,j,k)-beta*v(i,j,k)+cory)/(1.0+beta)-v(i,j,k)

END DO
END DO

! Step 3: final prediction of u & v including the flooding algorithm

DO j = 1,ny
DO k = 1,nx

! prediction for u
un(i,j,k) = 0.0
uu = u(i,j,k)
duu = du(i,j,k)
IF(wet(i,j,k)==1) THEN
  IF((wet(i,j,k+1)==1).or.(duu>0.0)) un(i,j,k) = (uu+duu)/Rx(i,j,k)
ELSE
  IF((wet(i,j,k+1)==1).and.(duu<0.0)) un(i,j,k) = (uu+duu)/Rx(i,j,k)
END IF

! prediction for v
vv = v(i,j,k)
dvv = dv(i,j,k)
vn(i,j,k) = 0.0
IF(wet(i,j,k)==1) THEN
  IF((wet(i,j+1,k)==1).or.(dvv>0.0)) vn(i,j,k) = (vv+dvv)/Ry(i,j,k)
ELSE
  IF((wet(i,j+1,k)==1).and.(dvv<0.0)) vn(i,j,k) = (vv+dvv)/Ry(i,j,k)
END IF

END DO
END DO

!========
END DO
!========

! boundary conditions
DO i = 1,nz
DO j = 1,ny
  un(i,j,nx) = 0.0
END DO
DO k = 1,nx
  vn(i,ny,k) = 0.0
END DO
END DO 

! layer-thickness change predictor
!============
DO i = 1,nz
!============

DO j = 0,ny+1
DO k = 0,nx+1
  CuP(j,k) = 0.5*(un(i,j,k)+abs(un(i,j,k)))*dt/dx
  CuN(j,k) = 0.5*(un(i,j,k)-abs(un(i,j,k)))*dt/dx
  CvP(j,k) = 0.5*(vn(i,j,k)+abs(vn(i,j,k)))*dt/dy
  CvN(j,k) = 0.5*(vn(i,j,k)-abs(vn(i,j,k)))*dt/dy
  Cu(j,k) = abs(un(i,j,k))*dt/dx
  Cv(j,k) = abs(vn(i,j,k))*dt/dy
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  B(j,k) = h(i,j,k)
END DO
END DO

CALL advect

DO j = 1,ny
DO k = 1,nx
  dhdt(i,j,k) = BN(j,k)
END DO
END DO

!=========
END DO
!=========

! update interface displacements
DO j = 1,ny
DO k = 1,nx
  deta = 0.0
  DO i = nz,1,-1
    deta = deta + dhdt(i,j,k)
    etan(i,j,k)= eta(i,j,k)+dt*deta
  END DO
END DO
END DO

DO i = 1,nz
DO j = 1,ny
DO k = 1,nx
  eta(i,j,k)= etan(i,j,k)
END DO
END DO
END DO

! update layer thicknesses, lateral velocities and wet/dry pointers
DO j = 1,ny
DO k = 1,nx
   DO i = 1,nz
    h(i,j,k) = hzero(i,j,k)+eta(i,j,k)-eta(i+1,j,k)-eta0(i,j,k)+eta0(i+1,j,k)
    u(i,j,k) = un(i,j,k)
    v(i,j,k) = vn(i,j,k)
    wet(i,j,k) = 1
    if(h(i,j,k)<hmin)wet(i,j,k) = 0
  END DO
END DO
END DO

RETURN
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

term1 = (1.0-CuP(j,k-1))*(B(j,k)-B(j,k-1))
BwP = B(j,k-1)+0.5*PSI(RxP(j,k-1),Cu(j,k-1),mode)*term1
term1 = (1.0+CuN(j,k-1))*(B(j,k)-B(j,k-1))
BwN = B(j,k)-0.5*PSI(RxN(j,k-1),Cu(j,k-1),mode)*term1
term1 = (1.0-CuP(j,k))*(B(j,k+1)-B(j,k))
BeP = B(j,k)+0.5*PSI(RxP(j,k),Cu(j,k),mode)*term1
term1 = (1.0+CuN(j,k))*(B(j,k+1)-B(j,k))  
BeN = B(j,k+1)-0.5*PSI(RxN(j,k),Cu(j,k),mode)*term1
term1 = (1.0-CvP(j-1,k))*(B(j,k)-B(j-1,k))
BsP = B(j-1,k)+0.5*PSI(RyP(j-1,k),Cv(j-1,k),mode)*term1
term1 = (1.0+CvN(j-1,k))*(B(j,k)-B(j-1,k))  
BsN = B(j,k)-0.5*PSI(RyN(j-1,k),Cv(j-1,k),mode)*term1
term1 = (1.0-CvP(j,k))*(B(j+1,k)-B(j,k)) 
BnP = B(j,k)+0.5*PSI(RyP(j,k),Cv(j,k),mode)*term1
term1 = (1.0+CvN(j,k))*(B(j+1,k)-B(j,k)) 
BnN = B(j+1,k)-0.5*PSI(RyN(j,k),Cv(j,k),mode)*term1

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