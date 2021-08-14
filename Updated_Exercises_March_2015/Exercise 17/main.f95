PROGRAM instab

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
!* The file "random.f95" needs to be compiled as well.
!* This contains the random-number creator and
!*     is used in the file "sub.f95".
!******************************************************

USE param
USE sub

! local parameters
REAL :: time
INTEGER :: n,ntot,nout

!**********
CALL INIT  ! initialisation
!**********

! initialisation of Eulerian tracer

DO k = 0,nx+1
DO j = 0,25
  T(j,k) = 0.
  TN(j,k) = 0.
END DO
DO j = 26,ny+1
  T(j,k) = 1.
  TN(j,k) = 1.
END DO
END DO

! initialise ambient geostrophic flow
! forcing
DO j = 0,22
  UGEO(j) = -0.2
END DO
DO j = 23,28
  UGEO(j) = -0.2 + 0.4*REAL(j-22)/7.
END DO
DO j = 29,ny+1
  UGEO(j) = 0.2
END DO

! runtime parameters
ntot = INT(2*24.*3600./dt)

! output parameter
nout = INT(1800./dt)

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

time = real(n)*dt

write(6,*)"time (days)", time/(24.*3600.)

taux = 0.0
tauy = 0.0

!call predictor
CALL dyn

! data output
IF(MOD(n,nout)==0)THEN
  DO j = 1,ny
    WRITE(10,'(201F12.6)')(eta(j,k),k=1,nx)
    WRITE(20,'(201F12.6)')(h(j,k),k=1,nx)
    WRITE(30,'(201F12.6)')(u(j,k)+UGEO(j),k=1,nx)
    WRITE(40,'(201F12.6)')(v(j,k),k=1,nx)
    WRITE(50,'(201F12.6)')(T(j,k),k=1,nx)
  END DO
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)
ENDIF

END DO ! end of iteration loop

END PROGRAM instab