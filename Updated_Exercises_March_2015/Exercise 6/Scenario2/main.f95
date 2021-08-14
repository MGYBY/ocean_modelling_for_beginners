PROGRAM wave1D

!*****************************************!
!* 1d shallow-Water model                *!
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
REAL :: time
INTEGER :: n,ntot,nout

!**********
CALL INIT  ! initialisation
!**********
! output of initial eta distribution
OPEN(10,file ='eta0.dat',form='formatted')
  WRITE(10,'(101F12.6)')(eta(k),k=1,nx)
CLOSE(10)

! output of initial layer thickness distribution
OPEN(10,file ='h0.dat',form='formatted')
  WRITE(10,'(101F12.6)')(hzero(k),k=1,nx)
CLOSE(10)

! set local parameters

! set epsilon for Shapiro filter
eps = 0.05

! runtime parameters
ntot = 200.0/dt

! output parameter
nout = 2.0/dt

! open files for output
OPEN(10,file ='eta.dat',form='formatted')
OPEN(20,file ='u.dat',form='formatted')

DO k = 2,20
  eta(k) = eta(k)+1.
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

DO k = 0,nx+1
  h(k) = hzero(k) + eta(k)
  wet(k) = 1
  if(h(k)<hmin)wet(k) = 0
  u(k) = un(k)
END DO

! data output
IF(MOD(n,nout)==0)THEN
  WRITE(10,'(101F12.6)')(eta(k),k=1,nx)
  WRITE(20,'(101F12.6)')(u(k),k=1,nx)
  WRITE(6,*)"Data output at time = ",time
ENDIF

END DO ! end of iteration loop

END PROGRAM wave1D
