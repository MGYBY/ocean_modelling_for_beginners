PROGRAM wind

!*****************************************!
!* 2d shallow-Water model                *!
!*                                       *!
!* including:                            *!
!* - wind-stress forcing 		 *!
!* - semi-implicit bottom friction       *!
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
ntot = INT(10.*24.*3600./dt)

! output parameter
nout = INT(2.*3600./dt)

! open files for output
OPEN(10,file ='eta.dat',form='formatted')
OPEN(20,file ='h.dat',form='formatted')
OPEN(30,file ='u.dat',form='formatted')
OPEN(40,file ='v.dat',form='formatted')

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = REAL(n)*dt

! write(6,*)"time (days)", time/(24.*3600.)

taux = 0.0
tauy = 0.2*MIN(time/(1.*24.*3600.),1.0)

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
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)

ENDIF

END DO ! end of iteration loop

! open files for output
OPEN(11,file ='etaF.dat',form='formatted')
OPEN(21,file ='hF.dat',form='formatted')
OPEN(31,file ='uF.dat',form='formatted')
OPEN(41,file ='vF.dat',form='formatted')

 DO j = 1,ny
    WRITE(11,'(101F12.6)')(eta(j,k),k=1,nx)
    WRITE(21,'(101F12.6)')(h(j,k),k=1,nx)
    WRITE(31,'(101F12.6)')(u(j,k),k=1,nx)
    WRITE(41,'(101F12.6)')(v(j,k),k=1,nx)
END DO

WRITE(6,*)"Final data output at time = ",time/(24.*3600.)


END PROGRAM wind
