PROGRAM multi

!*****************************************!
!* 1d multi-layer shallow-Water model    *!
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
REAL :: Apadle, Tpadle, pi
REAL :: hmax
INTEGER :: n, ntot, nout

!---- initialisation ----
CALL INIT 
!------------------------
! output of initial layer thickness distributions
OPEN(10,file ='h0.dat',form='formatted')
WRITE(10,'(101F12.6)')(htotal(k),k=1,nx)
CLOSE(10)

! automatic setting of time step (10% below CFL threshold) 
hmax = 100.0 ! total water depth
dt = 0.9*dx/sqrt(g*hmax)
write(6,*)"Time step = ", dt," seconds"
pause

! set epsilon for Shapiro filter
eps = 0.05

! parameters for wave padle
Apadle = 1.   ! amplitude in metres
! Tpadle = 10.0 ! period in seconds CASE 1
 Tpadle = 2.0*3600.0 ! period in seconds CASE 2

pi = 4.*atan(1.)

! runtime parameters
ntot = INT(10.*Tpadle/dt)

! output parameter
nout = INT(Tpadle/10./dt)

! open files for output
OPEN(10,file ='h.dat',form='formatted')
OPEN(20,file ='u.dat',form='formatted')

! initial configuration
DO i = 1,nz
  WRITE(10,'(101F12.6)')(h(i,k),k=1,nx)
  WRITE(20,'(101F12.6)')(u(i,k),k=1,nx)
END DO

! header file
OPEN(33,file ='header.txt',form='formatted')
  WRITE(33,*)dt*nout, dx, nx, nz, hmax
CLOSE(33)

!---- simulation loop ----
DO n = 1,ntot
!-------------------------
time = REAL(n)*dt

DO i = 1,nz
  eta(i,1) = Apadle*sin(2.0*pi*time/Tpadle)
END DO

!---- prognostic equations ----
CALL dyn
!------------------------------
! data output
IF(MOD(n,nout)==0)THEN
  DO i = 1,nz
    WRITE(10,'(101F12.6)')(h(i,k),k=1,nx)
    WRITE(20,'(101F12.6)')(u(i,k),k=1,nx)
  END DO
  WRITE(6,*)"Data output at time = ",time
ENDIF

!---- end of iteration loop ----
END DO
!-------------------------------

END PROGRAM multi