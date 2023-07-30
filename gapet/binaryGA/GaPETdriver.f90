!! GA is inherently a maximization algorithm
!! Need to change the flag based on the specific problem
PROGRAM GaPETdriver
IMPLICIT NONE

INCLUDE 'mpif.h'
!columns are contiguous so pass columns instead of rows
!https://fortran-lang.org/en/learn/best_practices/multidim_arrays/

! Simulation variables
INTEGER :: flag, iter

! GA variables
INTEGER :: nvar, npop, ngen, nbit, tbit
REAL(KIND=8) :: best_val_cur
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: popul, bounds
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: values

! selection variables
INTEGER :: idx
REAL(KIND=8) :: best_val, avgfit
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: best_point
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: fit, fit_prob, prob_per_string, c_prob

! crossover variables
REAL(KIND=8) :: crossrate

! mutation variables
REAL(KIND=8) :: mutrate

! MPI
INTEGER :: errcode, rank, np

!! -- initialize MPI 
CALL MPI_INIT(errcode)
CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,errcode)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD,np,errcode)

!! experimenting with np = npop, other option os npop+1

! OPEN(UNIT=10,FILE='parameters.txt', &
! 	FORM='FORMATTED',STATUS='OLD',ACTION='READ')
! READ(10,*) params(1,1:)
! CLOSE(10)

! can read from a text file instead
flag = 0                ! 0 for minimization; 1 for maximazation
nvar = 2                ! number of variables in the optimization problem
! should be even for current implementation                
npop = 40               ! population size (number of points in one generation)        
ngen = 150              ! max number of generations (iterations)
nbit = 40               ! number of bits to represent one variable
mutrate = 0.05          ! mutation rate
crossrate = 0.8         ! crossover rate

! bounds for each variable
ALLOCATE(bounds(nvar,2))
! DO ivar = 1:nvar
!     bounds(:,ivar) = 
! ENDDO
bounds(:,1) = 0
bounds(:,2) = 6

tbit = nbit*nvar        ! total bits to represent each point

ALLOCATE(popul(npop,tbit))

! intitialize population randomly
CALL random_number(popul)
popul(:,:) = NINT(popul(:,:))       ! round off the randomized population to 0's or 1's

! Recording best values
IF (flag == 0) THEN
    best_val = 10**15
ELSEIF (flag == 1) THEN
    best_val = -10**15
ENDIF

iter = 1
DO WHILE (iter <= ngen)
    ! converting binary representation to real values within bounds
    CALL binary2cont(popul, bounds)

    ! Evaluate the objective function at the current population
    values

    ! Updating the best values and best point
    IF (flag == 0) THEN
        best_val_cur = MINVAL(values)
        idx = MINLOC(values)
        IF (best_val_cur <= best_val) THEN
            best_val = best_val_cur
            best_point = popul(idx,:)
        ENDIF
    ELSEIF (flag == 1) THEN
        best_val_cur = MAXVAL(values)
        idx = MAXLOC(values)
        IF (best_val_cur >= best_val) THEN
            best_val = best_val_cur
            best_point = popul(idx,:)
        ENDIF
    ENDIF

    ! fitness values
    IF (flag==0) THEN
        fit = 1/(1+values)
    ELSEIF (flag == 1) THEN
        fit = values
    ENDIF

    ! calculating cumulative probability for roulette wheel selection
    avgfit = SUM(fit)/npop
    fit_prob = fit/avgfit
    prob_per_string = fit_prob/npop
    c_prob = cumsum(prob_per_string)

    ! roulette wheel selection
    CALL binary_sel(popul, values, flag)

    iter = iter+1
ENDDO

CALL MPI_BARRIER(MPI_COMM_WORLD,errcode)
CALL MPI_FINALIZE(errcode)

END PROGRAM GaPETdriver