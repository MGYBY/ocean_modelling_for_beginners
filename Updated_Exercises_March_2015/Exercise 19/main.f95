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
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!
!
! This model is applied with 2 layers only.
!******************************************

USE param
USE sub

! local parameters
REAL :: time
INTEGER :: n, ntot, nout

!**********
CALL INIT  ! initialisation
!**********

! runtime parameters
ntot = INT(100.*24.*3600./dt)

! output parameter
nout = INT(60.*3600./dt)

! open files for output
OPEN(10,file ='eta1.dat',form='formatted')
OPEN(20,file ='h1.dat',form='formatted')
OPEN(30,file ='u1.dat',form='formatted')
OPEN(40,file ='v1.dat',form='formatted')

OPEN(11,file ='eta2.dat',form='formatted')
OPEN(21,file ='h2.dat',form='formatted')
OPEN(31,file ='u2.dat',form='formatted')
OPEN(41,file ='v2.dat',form='formatted')

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = REAL(n)*dt

ad = MIN(time/(50.*24.*3600.),1.0)

write(6,*)"time (days)", time/(24.*3600.)

! call prognostic equations
CALL dyn

! data output
IF(MOD(n,nout)==0)THEN
  DO j = 1,ny
    WRITE(10,'(101F12.6)')(eta(1,j,k),k=1,nx)
    WRITE(20,'(101F12.6)')(h(1,j,k),k=1,nx)
    WRITE(30,'(101F12.6)')(u(1,j,k),k=1,nx)
    WRITE(40,'(101F12.6)')(v(1,j,k),k=1,nx)
    WRITE(11,'(101F12.6)')(eta(2,j,k),k=1,nx)
    WRITE(21,'(101F12.6)')(h(2,j,k),k=1,nx)
    WRITE(31,'(101F12.6)')(u(2,j,k),k=1,nx)
    WRITE(41,'(101F12.6)')(v(2,j,k),k=1,nx)
  END DO
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)

ENDIF

END DO ! end of iteration loop

END PROGRAM multi
