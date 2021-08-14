MODULE RANDOM

CONTAINS

!=======================================================================

REAL FUNCTION ran3(iseed)
!***********************************************************************

!    *RAN3*       RANDOM GENERATOR RETURNING A RANDOM NUMBER BETWEEN
!                 0 AND 1
!        REFERENCE - Press W.H., Flannery B.P., Teukolsky S.A. and Vetterling
!                    W.T., 1989. Numerical Recipes. The art of scientific
!                    computing. Cambridge University Press, Cambridge, 702 pp.
!***********************************************************************

!*    ARGUMENTS

INTEGER, INTENT(IN OUT)                     :: iseed


!------------------------------------------------------------------------
! (RE)INITIALISES THE RANDOM GENERATOR IF ISEED < 0 (IN)
!------------------------------------------------------------------------

INTEGER :: i, iff, ii, inext, inextp, k
INTEGER ::  mj, mk
INTEGER :: ma(55)

INTEGER, PARAMETER :: mbig=1000000000
INTEGER, PARAMETER :: mseed=161803398
INTEGER, PARAMETER :: mz=0
REAL, PARAMETER :: fac=1.0/mbig

SAVE inext, inextp, ma
DATA iff /0/

IF (iseed < 0.OR.iff == 0) THEN
  iff = 1
  mj = mseed - IABS(iseed)
  mj = MOD(mj,mbig)
  ma(55) = mj
  mk = 1
  DO  i=1,54
    ii = MOD(21*i,55)
    ma(ii) = mk
    mk = mj - mk
    IF (mk < mz) mk = mk + mbig
    mj = ma(ii)
  END DO
  DO  k=1,4
    DO  i=1,55
      ma(i) = ma(i) - ma(1+MOD(i+30,55))
      IF (ma(i) < mz) ma(i) = ma(i) + mbig
    END DO
  END DO
  inext = 0
  inextp = 31
  iseed = 1
END IF
inext = inext + 1
IF (inext == 56) inext = 1
inextp = inextp + 1
IF (inextp == 56) inextp = 1
mj = ma(inext)-ma(inextp)
IF (mj < mz) mj = mj+mbig
ma(inext) = mj
ran3 = mj*fac


RETURN

END FUNCTION ran3

END MODULE RANDOM