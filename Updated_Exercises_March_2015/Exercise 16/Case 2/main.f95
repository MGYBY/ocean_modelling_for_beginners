PROGRAM steer

!*****************************************!
!* 2d shallow-Water model                *!
!*                                       *!
!* including:                            *!
!* - wind-stress forcing		 *!
!* - semi-implicit bottom friction       *!
!* - nonlinear terms                     *!
!* - horizontal pressure-gradient force  *!
!* - Coriolis force (semi-implicit)      *!
!* - flooding algorithm                  *!
!* - TDV advection scheme 		 *!
!* - Eulerian tracer prediction          *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!
!*
!* Lateral momentum diffusion/friction
!*     is not included.
!*
!* Bathymetry file "topo.dat" must be
!*    supplied in local folder.
!******************************************
USE param
USE sub

! local parameters
REAL :: time
REAL :: ad
INTEGER :: n,ntot,nout

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

! runtime parameters
ntot = INT(20.*24.*3600./dt)

! output parameter
nout = INT(6.0*3600./dt)

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

ad = MIN(time/(5.*24.*3600.),1.0)

UGEO = ad*0.1

DO j = 0,ny+1
  T(j,1) = 0.0
  T(j,0) = 0.0
END DO

DO j = 39,41
  T(j,1) = 1.0
  T(j,0) = 1.0
END DO

CALL dyn

! data output
IF(MOD(n,nout)==0)THEN
  DO j = 1,ny
    WRITE(10,'(151F12.6)')(eta(j,k),k=1,nx)
    WRITE(20,'(151F12.6)')(h(j,k),k=1,nx)
    WRITE(30,'(151F12.6)')(u(j,k)+UGEO,k=1,nx)
    WRITE(40,'(151F12.6)')(v(j,k),k=1,nx)
    WRITE(50,'(151F12.6)')(T(j,k),k=1,nx)
  END DO

  WRITE(6,*)"Data output at time = ",time/(24.*3600.)

ENDIF

END DO ! end of iteration loop

END PROGRAM steer

