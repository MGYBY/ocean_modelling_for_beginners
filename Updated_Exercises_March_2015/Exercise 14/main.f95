PROGRAM wind

!*****************************************!
!* 2d shallow-Water model                *!
!*                                       *!
!* including:                            *!
!* - wind-stress forcing 		 *!
!* - semi-implicit bottom friction       *!
!* - nonlinear terms                     *!
!* - horizontal pressure-gradient force  *!
!* - lateral momentum diffusion/friction *!
!* - Shapiro filter                      *!
!* - flooding algorithm                  *!
!* - TDV advection scheme 		 *!
!* - Eulerian tracer prediction          *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!

!***********************************************************
! bathymetry must be supplied in terms of a file "topo.dat"
!***********************************************************

USE param
USE sub

! local parameters
REAL :: time
INTEGER :: n,ntot,nout

!================================================
! mode = 3 => Superbee scheme is used throughout
!================================================

mode = 3

!**********
CALL INIT  ! initialisation
!**********

! initialisation of Eulerian tracer

DO j = 0,ny+1
DO k = 0,nx+1
  T(j,k) = 0.
  TN(j,k) = 0.
END DO
END DO

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
ntot = INT(2.*24.*3600./dt)

! output parameter
nout = INT(900./dt)

! open files for output
OPEN(10,file ='eta.dat',form='formatted')
OPEN(20,file ='h.dat',form='formatted')
OPEN(30,file ='u.dat',form='formatted')
OPEN(40,file ='v.dat',form='formatted')
OPEN(50,file ='T.dat',form='formatted')

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = REAL(n)*dt

write(6,*)"time (days)", time/(24.*3600.)

taux = 0.2*MIN(time/(1.*24.*3600.),1.0)
tauy = 0.0

DO j = 1,ny
  T(j,0) = 0.0
  T(J,1) = 0.0
END DO

DO j = 25,27
  T(j,0) = 1.0
  T(J,1) = 1.0
END DO

! call predictor
CALL dyn

! updating including Shapiro filter

CALL shapiro

! cyclic boundary conditions
DO j = 1,ny
  eta(j,0) = eta(j,nx)
  eta(j,nx+1) = eta(j,1)
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

! data output
IF(MOD(n,nout)==0)THEN
  DO j = 1,ny
    WRITE(10,'(101F12.6)')(eta(j,k),k=1,nx)
    WRITE(20,'(101F12.6)')(h(j,k),k=1,nx)
    WRITE(30,'(101F12.6)')(u(j,k),k=1,nx)
    WRITE(40,'(101F12.6)')(v(j,k),k=1,nx)
    WRITE(50,'(101F12.6)')(T(j,k),k=1,nx)
  END DO

  WRITE(6,*)"Data output at time = ",time/(24.*3600.)

ENDIF

END DO ! end of iteration loop

END PROGRAM wind

