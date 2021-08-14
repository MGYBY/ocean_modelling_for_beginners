PROGRAM corio

!*******************************************
! Predictions of the path of a water
! parcel subject to inertial oscillations.
!
! Author: J. Kaempf, 2008
!*******************************************

REAL :: u,v,un,vn,x,y,xn,yn
REAL :: uzero, vzero, du, dv, ustar, vstar
REAL :: dt,fre,f,pi,alpha,beta
INTEGER :: n,ntot,mode

! ambient flow
uzero = 0.05
vzero = 0.05

! initial relative speed and location
u = 0.1
v = 0.0
x = 0.
y = 0.

pi = 4.*atan(1.)
freq = -2.*pi/(24.*3600.)
f = 2*freq
dt = 6.*24.*3600./120.
alpha = f*dt
beta = 0.25*alpha*alpha

ntot = 120;

mode = 2   ! choose between mode 1 and mode 2

IF(mode == 1)THEN
  OPEN(10,FILE='output1.txt',FORM='formatted')
ELSE
  OPEN(10,FILE='output2.txt',FORM='formatted')
END IF

WRITE(10,*)freq,dt,ntot

!**** start of iteration loop ****
DO n = 1,ntot
!*********************************

time = REAL(n)*dt
du = 0.0
dv = 0.0

IF(n == 40) THEN
 du = 0.0
 dv = -0.3
END IF

IF(n == 80) THEN
 du = 0.0
 dv = 0.1
END IF

ustar = u + du
vstar = v + dv

IF (mode == 1) THEN
  un = (ustar*(1-beta)+alpha*vstar)/(1+beta)
  vn = (vstar*(1-beta)-alpha*ustar)/(1+beta)
ELSE
  un = cos(alpha)*ustar+sin(alpha)*vstar
  vn = cos(alpha)*vstar-sin(alpha)*ustar
END IF

xn = x + dt*(un+uzero)/1000.0
yn = y + dt*(vn+vzero)/1000.0

u = un
v = vn
x = xn
y = yn

WRITE(10,*)x,y,time

!**** end of iteration loop ****
END DO
!*******************************

END PROGRAM corio




