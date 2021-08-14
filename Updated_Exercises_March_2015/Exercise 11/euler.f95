PROGRAM advect

!*****************************************!
!* Eulerian advection using TVD schemes  *!
!*                                       *!
!*=======================================*!
!* mode = 1 => upstream scheme           *!
!* mode = 2 => Lax-Wendroff scheme 	 *!
!* mode = 3 => Superbee scheme           *!
!* mode = 4 => Super-C scheme            *!
!========================================*!                  *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!

! IMPORTANT INFORMATION
! This code requires the files "uF.dat", "vF.dat", and "h0.dat" 
! from Exercise 10 as input data

INTEGER, PARAMETER :: nx = 51	
INTEGER, PARAMETER :: ny = 51	

REAL :: u(0:ny+1,0:nx+1), v(0:ny+1,0:nx+1)
REAL :: CuP(0:ny+1,0:nx+1), CuN(0:ny+1,0:nx+1)
REAL :: CvP(0:ny+1,0:nx+1), CvN(0:ny+1,0:nx+1)
REAL :: Cu(0:ny+1,0:nx+1), Cv(0:ny+1,0:nx+1)
REAL :: RxP(0:ny+1,0:nx+1), RxN(0:ny+1,0:nx+1)
REAL :: RyP(0:ny+1,0:nx+1), RyN(0:ny+1,0:nx+1)
REAL :: B(0:ny+1,0:nx+1), Bn(0:ny+1,0:nx+1)
REAL :: dt, dx, dy
REAL :: dB, term1, term2, term3, term4, div
REAL :: BwP, BwN, BeP, BeN, BsP, BsN, BnP, BnN 
INTEGER :: i,j,k, mode

mode = 3

! grid parameters
dx = 100.0
dy = 100.0
dt = 200.0

!read steady-state velocity field from Exercise 10

OPEN(10,file ='uF.dat',form='formatted')
DO j = 1,ny
  READ(10,'(51F12.6)')(u(j,k),k=1,nx)
END DO
CLOSE(10)

OPEN(10,file ='vF.dat',form='formatted')
DO j = 1,ny
  READ(10,'(51F12.6)')(v(j,k),k=1,nx)
END DO
CLOSE(10)

! set boundary values 
DO k = 0,nx+1
  u(0,k) = 0.0
  u(ny+1,k) = 0.0
  v(0,k) = 0.0
  v(ny+1,k) = 0.0
END DO

DO j = 0,ny+1
  u(j,0) = 0.0
  u(j,nx+1) = 0.0
  v(j,0) = 0.0
  v(j,nx+1) = 0.0
END DO

! set initial concentration field
DO j = 0,ny+1
DO k = 0,nx+1
  B(j,k) = 0.0
END DO
END DO

DO j = 40,45
DO k = 10,15
  B(j,k) = 1.0
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  Bn(j,k) = B(j,k)
END DO
END DO

! runtime parameters
ntot = INT(24.*3600./dt)

! output parameter
nout = INT(0.5*3600./dt)

! open files for output
OPEN(10,file ='B.dat',form='formatted')

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = REAL(n)*dt

DO j = 0,ny+1
DO k = 0,nx+1
  CuP(j,k) = 0.5*(u(j,k)+abs(u(j,k)))*dt/dx
  CuN(j,k) = 0.5*(u(j,k)-abs(u(j,k)))*dt/dx
  CvP(j,k) = 0.5*(v(j,k)+abs(v(j,k)))*dt/dy
  CvN(j,k) = 0.5*(v(j,k)-abs(v(j,k)))*dt/dy
  Cu(j,k) = abs(u(j,k))*dt/dx
  Cv(j,k) = abs(v(j,k))*dt/dx
END DO
END DO

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

div = dt*B(j,k)*( (u(j,k)-u(j,k-1) )/dx + ( v(j,k)-v(j-1,k) )/dy)  

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

Bn(j,k) = B(j,k)+term1-term2+term3-term4+div

END DO
END DO

! updating

DO j = 1,ny
DO k = 1,nx
  B(j,k) = Bn(j,k)
END DO
END DO

! data output
IF(MOD(n,nout)==0)THEN
  DO j = 1,ny
    WRITE(10,'(51F12.6)')(B(j,k),k=1,nx)
  END DO
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)

ENDIF

END DO ! end of iteration loop

END PROGRAM advect

REAL FUNCTION psi(r,cfl,mmode)

! input parameters

REAL, INTENT(IN) :: r, cfl  
INTEGER, INTENT(IN) :: mmode 

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

END FUNCTION psi 