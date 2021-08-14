PROGRAM buoy

!*************************************************
!* Simulation of oscillations of a buoyant object 
!*
!* Author: J. Kaempf, 2008
!*************************************************

INTEGER :: n    ! time level
INTEGER :: ntot ! total number of iterations

REAL :: z       ! location at time level n
REAL :: zn      ! location at time level n+1
REAL :: w       ! vertical speed at time level n
REAL :: wn      ! vertical speed at time level n+1
REAL :: rhosea  ! seawater density
REAL :: rho     ! density of your object
REAL :: N2      ! stability freqency squared
REAL :: g       ! gravity
REAL :: bf      ! buoyancy force
REAL :: dt      ! time step
REAL :: time    ! time counter
REAL :: zz      ! another z for output purposes
REAL :: r ! friction parameter

! Initialization of variables and parameters

z = -80.0       ! initial location is 50 m below sea surface
w = 0.0         ! no vertical speed at time zero
dt = 1.0        ! value of time step; dt = 1 second here
rho = 1025.5    ! density of your object
g = 9.81        ! gravity
N2 = 1.0e-4     ! stability frequency squared of ambient ocean
r = 0.00         ! friction parameter

! Total simulation time is 1 hour. How many iteration steps relate to this?

ntot = 3600./dt

! Open file for output
open(10, file = 'output.txt', form = 'formatted', status = 'unknown')

! write initial time, location, speed and relative density to this file
write(10,'(4F12.4)')time, Z, W, RHO-1000.

!****** start of iteration *******
DO n = 1,ntot   
!*********************************

rhosea = density(Z,N2)       ! determine ambient density at current location
bf = -g*(rho-rhosea)/rho  ! calculate buoyancy force
wn = (w+dt*bf)/(1.+r*dt)   ! predict new vertical speed
zn = z+dt*wn              ! predict new location
zn = MIN(zn,0.0)         ! location constrained by sea surface
zn = MAX(zn,-100.0)      ! location constrained by seafloor

! update Z and W for next time step

w = wn
z = zn 

time = REAL(n)*dt ! counting time as it progresses

IF(MOD(n,10) == 0) THEN
   write(10,'(4F12.4)')time, Z, W, RHO-1000.
END IF

!****** end of iteration *******
END DO
!*******************************

WRITE(6,*)" *** Simulation completed *** "

END PROGRAM buoy

!-----------------------------

REAL FUNCTION density(zin,N2in)

! input parameters

REAL, INTENT(IN) :: zin  ! vertical location
REAL, INTENT(IN) :: N2in ! stability frequency squared

! local parameters

REAL :: rhos    ! surface density
REAL :: g       ! gravity

rhos = 1025.0
g = 9.81

! return value

density = rhos*(1.-N2in/G*zin) 

END FUNCTION density
