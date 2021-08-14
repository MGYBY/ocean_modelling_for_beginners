PROGRAM decay

!**********************************************
!* Fortran code dealing with the decay problem. 
!*
!* Author: J. Kaempf, 2008
!**********************************************

! Text after pronounciation marks are treated as comments.
! Every FORTRAN program starts with the "program" statement.
! The program name, "decay" here, should be different from
! the file name. Lower case or upper case doesn't matter in
! FORTRAN.
! 
! Now comes the declaration section.
! Every symbol used later on needs to be declared.

INTEGER :: n    ! time level
INTEGER :: ntot ! total number of iterations
INTEGER :: nout ! output every nout's iteration step
INTEGER :: mode ! use either explicit or implicit scheme 

REAL :: C       ! concentration at time level n
REAL :: CN      ! "N" for new, concentration at next time level n+1
REAL :: CTRUE   ! exact analytical solution for comparison
REAL :: CZERO   ! Initial value of C
REAL :: dt      ! time step
REAL :: kappa   ! decay constant
REAL :: time    ! time counter
REAL :: fac     ! factor used in the scheme

! Now comes the initialization of variables and parameters.

CZERO = 100.0   ! initial concentration is 100
C = CZERO       ! initial value is used for exact solution
kappa = 0.0001  ! value of decay constant; you can also write this as kappa = 1.0e-4
dt = 3600.      ! value of time step; dt = 60 seconds here

! Here is an example of an IF-statement that writes a message to the command screen,
! if the stability condition is violated.

mode = 2        ! mode = 1 (explicit scheme); mode = 2 (implicit scheme)

fac = 1.0-dt*kappa

IF(mode == 2) fac = 1.0/(1.+dt*kappa) ! notice the two equal signs in the if-statement

IF(mode == 1)THEN
   if(fac<=0.0)write(6,*)'STABITY CRITERION ALERT: REDUCE TIME STEP'
END IF

! Total simulation time is 24 hours. How many interation steps is this?

ntot = 24.0*3600./dt  ! Remember that an hour has 60 minutes of 60 seconds each

! Data output at every hour of the simulation.

nout = 1.0*3600./dt

! Open file for output.
! The first (unit) number is a referene for output (see below).
! The "file" statement specifies the desired file name. 
! The statement "formatted" means ascii output.
! The status "unknown" implies new creation of a file if this does not exist,
! otherwise an existing file will be overwritten.
! Be careful not to overwrite files that you might need in the future.   

IF(mode == 1)THEN
  open(10, file = 'output1.txt', form = 'formatted', status = 'unknown')
END IF

IF(mode == 2)THEN
  open(10, file = 'output2.txt', form = 'formatted', status = 'unknown')
END IF

! Write initial time and concentration to this file.

WRITE(10,*)0,100.0,100.0 

! Now comes the start of the iteration loop.
! It means that the statements between the DO-loop start and 
! the corresponding "END DO" statement are repeated for n = 1 to
! n = ntot at steps of 1. In this case, statements are repeated
! a total number of 864000 times. Try to do this on a piece of paper... 

!****** Start of iteration *******
DO n = 1,ntot   
!*********************************

CN = C*fac ! prediction for next time step

time = REAL(n)*dt ! time counter

CTRUE = CZERO*exp(-kappa*time)  ! exact analytical solution

C = CN  ! updating for upcoming time step

! Data output if the counter i is a constant integer multiple of nout.
! For instance, nout = 30 produces outputs at i = 30, 60. 90, ....
! Note that two equal signs (==) have to be used in the IF statement.

IF(mod(i,nout) == 0)THEN

   ! Output to unit 10. The associated file name is given above. The "*" creates an automated format.
   ! Here we produce three output columns: the first is the time in hours, the second our prediction, and the
   ! third the exact solution. IF-statements with more than 1 line have to use a "THEN" and an "END IF".

   WRITE(10,*)time/3600.0, C, CTRUE

   ! We can also write a note a message to the screen here.

   write(6,*)"Data output at time = ",time/3600.0," hours"

END IF

!****** End of iteration *******
END DO
!*******************************

STOP
END PROGRAM decay