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
!* - Lagrangian float prediction         *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!
!*****************************************!
!
! This model is applied with 2 layers only
!******************************************

USE param
USE sub
USE random

! local parameters
REAL :: x, y, radius, angle, dro, term1
INTEGER :: n, ntot, nout, noutra

!**********
CALL INIT  ! initialisation
!**********

! runtime parameters
ntot = INT(5.*24.*3600./dt)

! output parameter
nout = INT(2.*3600./dt)

! initial low-density surface patch
DO j = 1,ny
DO k = 1,nx
  y = REAL(26-j)
  x = REAL(26-k) 
  radius = SQRT(x*x+y*y)
  IF(radius < 10.0)eta(2,j,k) = eta(2,j,k)-100.0
END DO
END DO

DO j = 1,ny
DO k = 1,nx
DO i = 1,nz
  term1 = eta(i,j,k)-eta(i+1,j,k)-eta0(i,j,k)+eta0(i+1,j,k)
  h(i,j,k) = hzero(i,j,k)+ term1 
  wet(i,j,k) = 1 
  if(h(i,j,k) < hmin)wet(i,j,k) = 0
END DO
END DO
END DO

!initialisation of random function
ist = -1
randm = ran3(ist)

! add Lagrangian floats to both layers
! in concentric circle segments 

DO i = 1,nz

DO itra = 1,200
  radius = 10.0-5.0*ran3(ist)
  angle = 2.0*PI*ran3(ist)
  xpos = (26.0+radius*COS(angle))*dx
  ypos = (26.0+radius*SIN(angle))*dy
  tra(i,itra,1) = xpos
  tra(i,itra,2) = ypos
END DO

DO itra = 201,400
  radius = 15.0-5.0*ran3(ist)
  angle = 2.0*PI*ran3(ist)
  xpos = (26.0+radius*COS(angle))*dx
  ypos = (26.0+radius*SIN(angle))*dy
  tra(i,itra,1) = xpos
  tra(i,itra,2) = ypos
END DO

DO itra = 401,600
  radius = 20.0-5.0*ran3(ist)
  angle = 2.0*PI*ran3(ist)
  xpos = (26.0+radius*COS(angle))*dx
  ypos = (26.0+radius*SIN(angle))*dy
  tra(i,itra,1) = xpos
  tra(i,itra,2) = ypos
END DO

END DO

! open files for output
! top layer
OPEN(10,file ='eta1.dat',form='formatted')
OPEN(20,file ='h1.dat',form='formatted')
OPEN(30,file ='u1.dat',form='formatted')
OPEN(40,file ='v1.dat',form='formatted')
! bottom layer
OPEN(11,file ='eta2.dat',form='formatted')
OPEN(21,file ='h2.dat',form='formatted')
OPEN(31,file ='u2.dat',form='formatted')
OPEN(41,file ='v2.dat',form='formatted')

! output of initial distributions
DO j = 1,nx
  WRITE(10,'(51F12.6)')(eta(1,j,k),k=1,nx)
  WRITE(20,'(51F12.6)')(h(1,j,k),k=1,nx)
  WRITE(30,'(51F12.6)')(u(1,j,k),k=1,nx)
  WRITE(40,'(51F12.6)')(v(1,j,k),k=1,nx)
  WRITE(11,'(51F12.6)')(eta(2,j,k),k=1,nx)
  WRITE(21,'(51F12.6)')(h(2,j,k),k=1,nx)
  WRITE(31,'(51F12.6)')(u(2,j,k),k=1,nx)
  WRITE(41,'(51F12.6)')(v(2,j,k),k=1,nx)
END DO

! output parameter for float locations
noutra = INT(1.*3600./dt)

! open files for output
OPEN(60,file ='TRx1.dat',form='formatted',recl = 10000000)
OPEN(70,file ='TRy1.dat',form='formatted',recl = 10000000)
OPEN(61,file ='TRx2.dat',form='formatted',recl = 10000000)
OPEN(71,file ='TRy2.dat',form='formatted',recl = 10000000)

!---------------------------
! simulation loop
!---------------------------
DO n = 1,ntot

time = real(n)*dt
write(6,*)"time (hours)", time/(3600.)

! gradual adjustment of density contrast
ad = MIN(time/(48.*3600.),1.0)
dro = 1.0*ad
rho(2) = 1028.0
rho(1) = rho(2)-dro

! call prognostic equations
CALL dyn

! start Lagrangian float prediction

IF(ad == 1.0)THEN

DO i = 1,nz

DO itra = 1,ntra

! locate grid cell of tracer
jpos = INT((tra(i,itra,1)+0.5*dy)/dy) + 1
kpos = INT((tra(i,itra,2)+0.5*dx)/dx) + 1
utra = 0.5*(u(i,jpos,kpos)+u(i,jpos,kpos-1))
vtra = 0.5*(v(i,jpos,kpos)+v(i,jpos-1,kpos))
! change of location
tra(i,itra,1) = tra(i,itra,1)+dt*vtra
tra(i,itra,2) = tra(i,itra,2)+dt*utra

END DO

END DO

! output of tracer locations
IF(MOD(n,noutra)==0)THEN
    WRITE(60,'(900F12.6)')(tra(1,itra,2)/1000.,itra=1,ntra)
    WRITE(70,'(900F12.6)')(tra(1,itra,1)/1000.,itra=1,ntra)
    WRITE(61,'(900F12.6)')(tra(2,itra,2)/1000.,itra=1,ntra)
    WRITE(71,'(900F12.6)')(tra(2,itra,1)/1000.,itra=1,ntra)
ENDIF
END IF

! output or other data
IF(MOD(n,nout)==0)THEN
  DO j = 1,ny
    WRITE(10,'(51F12.6)')(eta(1,j,k),k=1,nx)
    WRITE(20,'(51F12.6)')(h(1,j,k),k=1,nx)
    WRITE(30,'(51F12.6)')(u(1,j,k),k=1,nx)
    WRITE(40,'(51F12.6)')(v(1,j,k),k=1,nx)
    WRITE(11,'(51F12.6)')(eta(2,j,k),k=1,nx)
    WRITE(21,'(51F12.6)')(h(2,j,k),k=1,nx)
    WRITE(31,'(51F12.6)')(u(2,j,k),k=1,nx)
    WRITE(41,'(51F12.6)')(v(2,j,k),k=1,nx)
  END DO
  WRITE(6,*)"Data output at time = ",time/(24.*3600.)

ENDIF

END DO ! end of iteration loop

END PROGRAM multi