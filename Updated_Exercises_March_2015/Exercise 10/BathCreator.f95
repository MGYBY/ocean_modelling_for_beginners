!*************************************
! This code is an example of how 
! to create variable bathymetry used
! for in the shallow-water model.
!
! Author: J. Kaempf, 2008
!*************************************

PROGRAM topo

INTEGER, PARAMETER :: nx = 51	
INTEGER, PARAMETER :: ny = 51	
REAL :: h2(0:ny+1,0:nx+1), h1(0:ny+1,0:nx+1)
REAL :: diff, rad
INTEGER :: n,ntot, jis, kis

! smoothing parameter
diff = 0.02

! set coarse batymetry 
DO j = 0,ny+1
DO k = 0,nx+1
  h1(j,k) = 1.0 
END DO
END DO
 
DO j = 3,ny-2
DO k = 3,nx-2
  h1(j,k) = 20.0 
END DO
END DO

! land boundaries
DO k = 0,nx+1
  h1(0,k) = -10.0
  h1(1,k) = -10.0
  h1(ny+1,k) = -10.0
  h1(ny,k) = -10.0
END DO

DO j = 0,ny+1
  h1(j,0) = -10.0
  h1(j,1) = -10.0
  h1(j,nx+1) = -10.0
  h1(j,nx) = -10.0
END DO

! add circular island
jis = 26 ! centre of island
kis = 36

DO j = 1,ny
DO k = 1,nx
  rad = SQRT( (REAL(j-jis))**2 + (REAL(k-kis))**2 )
  IF(rad < 5.0) h1(j,k) = 1.0 ! shallow water around island 
  IF(rad < 4.0) h1(j,k) =  -0.2 ! island with 20 cm elevation
END DO
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  h2(j,k) = h1(j,k)
END DO
END DO

! open output file
OPEN(10,file ='topo.dat',form='formatted')

! runtime parameters
ntot = 1000

!---------------------------
! smoothing loop
!---------------------------
DO n = 1,ntot

DO j = 1,ny
DO k = 1,nx

  hh1 = h1(j,k)

  IF(h1(j,k)<= 1.0)THEN
     h2(j,k) = hh1 ! excluded from smoothing
  ELSE
     dhe = h1(j,k+1)-hh1
     IF(h1(j,k+1)< 0.0)dhe = 0.0
     dhw = hh1-h1(j,k-1)
     IF(h1(j,k-1)< 0.0)dhw = 0.0
     dhn = h1(j+1,k)-hh1
     IF(h1(j+1,k)< 0.0)dhn = 0.0
     dhs = hh1-h1(j-1,k)
     IF(h1(j-1,k)< 0.0)dhs = 0.0
     h2(j,k) = hh1 + diff*(dhe-dhw+dhn-dhs)
  ENDIF

END DO
END DO

! update for next iteration step

DO j = 1,ny
DO k = 1,nx
  h1(j,k) = h2(j,k)
END DO
END DO

END DO ! end of iteration loop

DO j = 1,ny
  WRITE(10,'(51F12.6)')(h2(j,k),k=1,nx)
END DO

END PROGRAM topo
