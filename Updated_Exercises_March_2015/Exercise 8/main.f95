PROGRAM wave2D

!*****************************************!
!* 2d shallow-Water model                *!
!*                                       *!
!* including:                            *!
!* - horizontal pressure-gradient force  *!
!* - Shapiro filter                      *!
!* - flooding algorithm                  *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!

USE param
USE sub

! local parameters
REAL :: hmax,time,dtmax
REAL :: c,lambda
INTEGER :: n,ntot,nout

!**********
CALL INIT  ! initialisation
!**********

! output of initial eta distribution
OPEN(10,file ='eta0.dat',form='formatted')
  DO j = 1,ny
    WRITE(10,'(101F12.6)')(eta(j,k),k=1,nx)
  END DO
CLOSE(10)

! output of initial layer thickness distribution
OPEN(10,file ='h0.dat',form='formatted')
  DO j = 1,ny
    WRITE(10,'(101F12.6)')(hzero(j,k),k=1,nx)
  END DO
CLOSE(10)

! set local parameters

! set epsilon for Shapiro filter
eps = 0.05

! runtime parameters
ntot = 1000

! output parameter
nout = 5

! determine maximum water depth
hmax = 0.

DO j = 1,ny
DO k = 1,nx
  hmax = MAX(hmax,h(j,k))
END DO
END DO

! maximum phase speed
c = SQRT(2*g*hmax)

! determine stability parameter
lambda = dt*SQRT(g*hmax)/MIN(dx,dy)

IF(lambda > 1)THEN
  WRITE(6,*) "This will not work. Do you know why?"   
  STOP
END IF

! open files for output
OPEN(10,file ='eta.dat',form='formatted')
OPEN(20,file ='h.dat',form='formatted')
OPEN(30,file ='u.dat',form='formatted')
OPEN(40,file ='v.dat',form='formatted')

DO j = 26,26
DO k = 26,26
eta(j,k) = 1.0  
END DO
END DO

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = REAL(n)*dt

! call predictor
CALL dyn

! updating including Shapiro filter

CALL shapiro

DO j = 0,ny+1
DO k = 0,nx+1
  h(j,k) = hzero(j,k) + eta(j,k)
  wet(j,k) = 1
  IF(h(j,k)<hmin)wet(j,k) = 0
  u(j,k) = un(j,k)
  v(j,k) = vn(j,k)
END DO
END DO

! data output
IF(MOD(n,nout)==0)THEN
  DO j = 1,ny
    WRITE(10,'(101F12.6)')(eta(j,k),k=1,nx)
    WRITE(20,'(101F12.6)')(h(j,k),k=1,nx)
    WRITE(30,'(101F12.6)')(u(j,k),k=1,nx)
    WRITE(40,'(101F12.6)')(v(j,k),k=1,nx)
  END DO
  WRITE(6,*)"Data output at time = ",time
ENDIF

END DO ! end of iteration loop

END PROGRAM wave2D
