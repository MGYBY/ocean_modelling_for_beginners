PROGRAM topo

INTEGER(4), PARAMETER :: nx = 201	
INTEGER(4), PARAMETER :: ny = 51	
REAL :: h2(0:ny+1,0:nx+1), h1(0:ny+1,0:nx+1)
REAL :: diff, rad, xx
INTEGER :: n,ntot, jis, kis, kc(ny)

! smoothing parameter
diff = 0.2

! define coastline

DO j = 0,ny+1
xx = REAL(j)/REAL(ny+1)
kc(j) = 170 + 30*xx
END DO


! set coarse batymetry 
DO j = 0,ny+1
DO k = 0,nx+1
  h1(j,k) = 30.0
  if(k>kc(j))h1(j,k) = -2.0 
END DO
  h1(j,kc(j)) = 2.0
END DO

!
DO k = 0,nx+1
  h1(0,k) = h1(1,k)
  h1(ny+1,k) = h1(ny,k)
END DO

DO j = 0,ny+1
   h1(j,0) = h1(j,1)
   h1(j,nx+1) = -10.0
END DO

DO j = 0,ny+1
DO k = 0,nx+1
  h2(j,k) = h1(j,k)
END DO
END DO

! open output file
OPEN(10,file ='topo.dat',form='formatted')

! runtime parameters
ntot = 100000

!---------------------------
! smoothing loop
!---------------------------
DO n = 1,ntot

DO j = 1,ny
DO k = 1,kc(j)-1

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

!
DO k = 0,nx+1
  h1(0,k) = h1(1,k)
  h1(ny+1,k) = h1(ny,k)
END DO

!DO j = 0,ny+1
!   h1(j,0) = h1(j,1)
!END DO

END DO ! end of iteration loop

DO j = 0,ny+1
  WRITE(10,'(203F12.6)')(h1(j,k),k=0,nx+1)
END DO

END PROGRAM topo
