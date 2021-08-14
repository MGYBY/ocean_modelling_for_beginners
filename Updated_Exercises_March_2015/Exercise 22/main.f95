PROGRAM plume

!*****************************************!
!* 2d reduced-gravity plume model        *!
!*                                       *!
!* including:                            *!
!* - semi-implicit bottom friction       *!
!* - nonlinear terms                     *!
!* - horizontal pressure-gradient force  *!
!* - Coriolis force (f plane)		 *!
!* - flooding algorithm                  *!
!* - TDV advection scheme 		 *!
!* - Eulerian tracer prediction 	 *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!

USE param
USE sub

! local parameters
INTEGER :: n, ntot, nout

!**********
CALL INIT  ! initialisation
!**********

! runtime parameters
ntot = INT(5.*24.*3600./dt)

! output parameter
nout = INT(1.*3600./dt)

! initial tracer distributions (disabled here)
DO i = 1,nz
DO j = 0,ny+1
DO k = 0,nx+1
  T(i,j,k) = 0.
  TN(i,j,k) = 0.
END DO
END DO
END DO

! open files for output
OPEN(10,file ='eta20.dat',form='formatted')
DO j = 1,ny
  WRITE(10,'(101F12.6)')(eta0(2,j,k),k=1,nx)
END DO
CLOSE(10)

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
    WRITE(10,'(101F12.6)')(eta(1,j,k)-eta0(1,j,k),k=1,nx)
    WRITE(20,'(101F12.6)')(h(1,j,k),k=1,nx)
    WRITE(30,'(101F12.6)')(u(1,j,k),k=1,nx)
    WRITE(40,'(101F12.6)')(v(1,j,k),k=1,nx)
    WRITE(50,'(101F12.6)')(T(1,j,k),k=1,nx)
    WRITE(11,'(101F12.6)')(eta(2,j,k)-eta0(2,j,k),k=1,nx)
    WRITE(21,'(101F12.6)')(h(2,j,k),k=1,nx)
    WRITE(31,'(101F12.6)')(u(2,j,k),k=1,nx)
    WRITE(41,'(101F12.6)')(v(2,j,k),k=1,nx)
    WRITE(51,'(101F12.6)')(T(2,j,k),k=1,nx)
  END DO
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)

ENDIF

END DO ! end of iteration loop

END PROGRAM plume
