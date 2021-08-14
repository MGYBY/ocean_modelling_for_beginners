PROGRAM multi

!*****************************************!
!* 2d multi-layer shallow-Water model    *!
!*                                       *!
!* including:                            *!
!* - wind-stress forcing		 *!
!* - semi-implicit bottom friction       *!
!* - nonlinear terms                     *!
!* - horizontal pressure-gradient force  *!
!* - Coriolis force (beta plane)	 *!
!* - Lateral momentum diffusion/friction *!
!* - flooding algorithm                  *!
!* - TDV advection scheme 		 *!
!* - Eulerian tracer prediction          *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!

USE param
USE sub
USE random

! local parameters
REAL :: etN, etS
INTEGER :: n, ntot, nout, noutra

!**********
CALL INIT  ! initialisation
!**********

! runtime parameters
ntot = INT(20.*24.*3600./dt)

! output parameter
nout = INT(6.*3600./dt)

! initial tracer distribution
DO i = 1,nz
DO j = 0,ny+1
DO k = 0,nx+1
  T(i,j,k) = 0.
  TN(i,j,k) = 0.
END DO
END DO
END DO

DO i = 1,nz
DO j = 26,ny
DO k = 0,nx+1
  T(i,j,k) = 1.
  TN(i,j,k) = 1.
END DO
END DO
END DO

! initial surface pressure field
DO k = 0,nx+1
  DO j = 26-5,26+5
    eta(1,j,k) = -0.05*(SIN(REAL(j-26)*0.5*PI/5.))
  END DO
  DO j = 0,20
    eta(1,j,k) = +0.05
  END DO
  DO j = 32,ny+1
    eta(1,j,k) = -0.05
  END DO
END DO

! initial structure of density interface
DO k = 0,nx+1
DO j = 0,ny+1
  eta(2,j,k) = -eta(1,j,k)*rho(2)/(rho(2)-rho(1))
END DO
END DO

! initial layer configuration
DO j = 0,ny+1
DO k = 0,nx+1
DO i = 1,nz
  h(i,j,k) = hzero(i,j,k)+eta(i,j,k)-eta(i+1,j,k)-eta0(i,j,k)+eta0(i+1,j,k)
  wet(i,j,k) = 1
  if(h(i,j,k)<hmin)wet(i,j,k) = 0
END DO
END DO
END DO

! initial geostrophic flow field in upper layer
DO j = 1,ny
DO k = 1,nx
  etN = 0.5*(eta(1,j+1,k)+eta(1,j+1,k+1))
  etS = 0.5*(eta(1,j-1,k)+eta(1,j-1,k+1))
  u(1,j,k) = -0.5*g*(etN-etS)/dy/f(1)
END DO
END DO

DO j = 1,ny
  u(1,j,0) = u(1,j,1) 
  u(1,j,nx+1) = u(1,j,nx)
END DO 

! add small random perturbations to surface pressure field
ist = -1
randm = ran3(ist)

DO j = 0,ny+1
DO k = 0,nx+1
  eta(1,j,k) = eta(1,j,k)+0.005*ran3(ist)
END DO
END DO

! open files for output 
OPEN(10,file ='eta1.dat',form='formatted')
OPEN(20,file ='h1.dat',form='formatted')
OPEN(30,file ='u1.dat',form='formatted')
OPEN(40,file ='v1.dat',form='formatted')
OPEN(50,file ='T1.dat',form='formatted')

OPEN(11,file ='eta2.dat',form='formatted')
OPEN(21,file ='h2.dat',form='formatted')
OPEN(31,file ='u2.dat',form='formatted')
OPEN(41,file ='v2.dat',form='formatted')
OPEN(51,file ='T2.dat',form='formatted')

! output of initial configuration
DO j = 1,ny
  WRITE(10,'(101F12.6)')(eta(1,j,k),k=1,nx)
  WRITE(20,'(101F12.6)')(h(1,j,k),k=1,nx)
  WRITE(30,'(101F12.6)')(u(1,j,k),k=1,nx)
  WRITE(40,'(101F12.6)')(v(1,j,k),k=1,nx)
  WRITE(50,'(101F12.6)')(T(1,j,k),k=1,nx)
  WRITE(11,'(101F12.6)')(eta(2,j,k),k=1,nx)
  WRITE(21,'(101F12.6)')(h(2,j,k),k=1,nx)
  WRITE(31,'(101F12.6)')(u(2,j,k),k=1,nx)
  WRITE(41,'(101F12.6)')(v(2,j,k),k=1,nx)
  WRITE(51,'(101F12.6)')(T(2,j,k),k=1,nx)
END DO

! output parameter
noutra = INT(1.*3600./dt)

!---------------------------
! simulation loop
!---------------------------

DO n = 1,ntot

time = REAL(n)*dt
write(6,*)"time (hours)", time/(3600.)
ad = 0.0 ! adjustment not used in this exercise

! call prognostic equations
CALL dyn

! data output
IF(MOD(n,nout)==0)THEN
  DO j = 1,ny
    WRITE(10,'(101F12.6)')(eta(1,j,k),k=1,nx)
    WRITE(20,'(101F12.6)')(h(1,j,k),k=1,nx)
    WRITE(30,'(101F12.6)')(u(1,j,k),k=1,nx)
    WRITE(40,'(101F12.6)')(v(1,j,k),k=1,nx)
    WRITE(50,'(101F12.6)')(T(1,j,k),k=1,nx)
    WRITE(11,'(101F12.6)')(eta(2,j,k),k=1,nx)
    WRITE(21,'(101F12.6)')(h(2,j,k),k=1,nx)
    WRITE(31,'(101F12.6)')(u(2,j,k),k=1,nx)
    WRITE(41,'(101F12.6)')(v(2,j,k),k=1,nx)
    WRITE(51,'(101F12.6)')(T(2,j,k),k=1,nx)

  END DO
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)

ENDIF

END DO ! end of iteration loop

END PROGRAM multi
