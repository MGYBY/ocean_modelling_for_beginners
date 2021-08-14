PROGRAM traject

!*****************************************!
!* Lagrangian float prediction           *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!

!**************************************************
! Compile with g95 lagrange.f95 random.f95
! The files "uF.dat" and "vF.dat" from Exercise 10
! are required for this code.
!**************************************************

USE random

INTEGER, PARAMETER :: nx = 51	
INTEGER, PARAMETER :: ny = 51	
INTEGER, PARAMETER :: ntra = 3000 ! number of floats

REAL :: u(0:ny+1,0:nx+1), v(0:ny+1,0:nx+1), uu, vv
REAL :: tra(ntra,2) ! tracer coordinates
REAL :: dt, dx, dy, randm, xlen, ylen, xpos, ypos
INTEGER :: i, j, k, ist, jpos, kpos

! grid parameters
dx = 100.0
dy = 100.0
dt = 10.0

! read steady-state velocity field from Exercise 10

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

! initial random locations of tracers

!initialisation of random function
ist = -1
randm = ran3(ist)
xlen = REAL(nx)*dx - 10.*dx
ylen = REAL(ny)*dy - 10.*dy

DO i = 1,ntra
  xpos = 5.*dx+ran3(ist)*xlen
  ypos = 5.*dy+ran3(ist)*ylen
  tra(i,1) = ypos
  tra(i,2) = xpos
END DO

! runtime parameters
ntot = INT(24.*3600./dt)

! output parameter
nout = INT(5.*60./dt)

! open files for output
OPEN(10,file ='TRx.dat',form='formatted',recl = 10000000)
OPEN(20,file ='TRy.dat',form='formatted',recl = 10000000)

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = REAL(n)*dt

DO i = 1,ntra

! locate grid cell of tracer
jpos = INT( (tra(i,1)+0.5*dy)/dy) + 1
kpos = INT( (tra(i,2)+0.5*dx)/dx) + 1

! velocity interpolated to nearest h gridpoint
uu = 0.5*(u(jpos,kpos-1)+u(jpos,kpos))
vv = 0.5*(v(jpos-1,kpos)+v(jpos,kpos))

! change of location
tra(i,1) = tra(i,1)+dt*vv
tra(i,2) = tra(i,2)+dt*uu

END DO

! data output
IF(MOD(n,nout)==0)THEN
  WRITE(10,'(5000F12.6)')(tra(i,2)/1000.,i=1,ntra)
  WRITE(20,'(5000F12.6)')(tra(i,1)/1000.,i=1,ntra)
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)
ENDIF

END DO ! end of iteration loop

END PROGRAM traject

