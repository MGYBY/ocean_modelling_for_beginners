PROGRAM kelvin

!*****************************************!
!* 2d shallow-Water model                *!
!*                                       *!
!* including:                            *!
!* - wind-stress forcing		 *!
!* - semi-implicit bottom friction       *!
!* - nonlinear terms                     *!
!* - horizontal pressure-gradient force  *!
!* - Coriolis force (semi-implicit)      *!
!* - Shapiro filter                      *!
!* - flooding algorithm                  *!
!* - TDV advection scheme 		 *!
!* - Eulerian tracer prediction          *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!
!*
!* Lateral momentum diffusion/friction
!* is not included.
!*
!*****************************************

USE param
USE sub

! local parameters
REAL :: time
REAL :: ad
INTEGER :: n,ntot,nout

!======================================
! mode = 1 => upstream scheme
! mode = 2 => Lax-Wendroff scheme
! mode = 3 => Superbee scheme
! mode = 4 => Super-C scheme
!=======================================

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
OPEN(10,file ='eta0.dat',form='formatted',recl=1000000)
  DO j = 1,ny
    WRITE(10,'(201F12.6)')(eta(j,k),k=1,nx)
  END DO
CLOSE(10)

! output of initial layer thickness distribution
OPEN(10,file ='h0.dat',form='formatted',recl=10000000)
  DO j = 1,ny
    WRITE(10,'(201F12.6)')(hzero(j,k),k=1,nx)
  END DO
CLOSE(10)

! set local parameters

! set epsilon for Shapiro filter
eps = 0.05

! runtime parameters
ntot = INT(12.*3600./dt)

! output parameter
nout = INT(10.*60./dt)

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

taux = 0.0
tauy = 0.0

! sea-level forcing
ad = 1.0*SIN(time/(2.*3600.)*2.*PI)
eta(6,2) = ad

! call predictor
CALL dyn

! data output
IF(MOD(n,nout)==0)THEN
  DO j = 1,ny
    WRITE(10,'(201F12.6)')(eta(j,k),k=1,nx)
    WRITE(20,'(201F12.6)')(h(j,k),k=1,nx)
    WRITE(30,'(201F12.6)')(u(j,k),k=1,nx)
    WRITE(40,'(201F12.6)')(v(j,k),k=1,nx)
    WRITE(50,'(201F12.6)')(T(j,k),k=1,nx)
  END DO

  WRITE(6,*)"Data output at time = ",time/(24.*3600.)

ENDIF

END DO ! end of iteration loop

END PROGRAM kelvin

