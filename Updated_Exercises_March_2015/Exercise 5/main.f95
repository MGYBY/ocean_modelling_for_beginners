PROGRAM wave1D

!*****************************************!
!* 1d shallow-Water model                *!
!*                                       *!
!* including:                            *!
!* - horizontal pressure-gradient force  *!
!* - Shapiro filter                      *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!

USE param
USE sub

! local parameters
REAL :: hmax,time,dtmax
REAL :: Apadle,Tpadle,pi,c,lambda
REAL :: dist, dwdz
INTEGER :: n,ntot,nout,mode

!**********
CALL INIT  ! initialisation
!**********

! set local parameters

! set epsilon for Shapiro filter
eps = 0.0

! parameters for wave padle
Apadle = 1.0   ! amplitude 
Tpadle = 20.0 ! period
pi = 4.*atan(1.)

! mode = 1.  ! use mode = 1 for wave paddle
mode = 2.  ! use mode = 2 for dam-break scenario

IF(mode==2)THEN
DO k = 51-5,51+5
  eta(k) = 1.  
END DO
END IF
!XXXXXXXXXXXXXXXXXXXX

! runtime parameters
ntot = 200.0/dt

! output parameter
nout = 2.0/dt

! determine maximum water depth
hmax = 0.
DO k = 0,nx+1
  hmax = MAX(hmax,h(k))
END DO

! maximum phase speed
c = SQRT(g*hmax)

! determine stability parameter
lambda = dt*SQRT(g*hmax)/dx

IF(lambda > 1)THEN
  WRITE(6,*) "This will not work. Do you know why?"   
  STOP
END IF

! open files for output
OPEN(10,file ='eta.dat',form='formatted')
OPEN(20,file ='u.dat',form='formatted')
OPEN(30,file ='w.dat',form='formatted')
!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = REAL(n)*dt

IF(mode==1)THEN
!wave padle at k = 51
  eta(51) = Apadle*sin(2.0*pi*time/Tpadle) 
END IF

! call predictor
CALL dyn

! updating including Shapiro filter
CALL shapiro

DO k = 0,nx+1
  h(k) = hzero(k) + eta(k)
  u(k) = un(k)
END DO

! calculate w at 5 m below sea surface
DO k = 1,nx
  dist = h(k)-5.0;
  w(k) = -(u(k)-u(k-1))/dx*dist
END DO

w(0) = 0.
w(nx+1) = 0.

! calculate w at u-gridpoints

DO k = 1,nx
  wu(k) = 0.5*(W(k)+w(k+1))
END DO

! data output
IF(MOD(n,nout)==0)THEN
  WRITE(10,'(101F12.6)')(eta(k),k=1,nx)
  WRITE(20,'(101F12.6)')(u(k),k=1,nx)
  WRITE(30,'(101F12.6)')(wu(k),k=1,nx)
  WRITE(6,*)"Data output at time = ",time
ENDIF

END DO ! end of iteration loop

END PROGRAM wave1D
