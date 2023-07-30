SUBROUTINE binary2cont(popul, bounds, pool)
IMPLICIT NONE

! INCLUDE mpif.h

real(kind=8), intent(in) :: popul(:,:), bounds(:,:)
real(kind=8), intent(out), allocatable, dimension(:,:) :: pool

real(kind=8), allocatable, dimension(:,:) :: rePool
real(kind=8), allocatable, dimension(:) :: quant, tempPool, par

integer :: npop, tbit, nbit, nvar, tvar

integer :: i, j !, count

npop = size(popul, dim=1)
tbit = size(popul, dim=2)
nvar = size(bounds, dim=1)

nbit = tbit/nvar
tvar = nvar*npop

allocate(rePool(tvar,nbit))
allocate(quant(nbit), tempPool(tvar), par(tvar))
allocate(pool(npop,nvar))

rePool = reshape(transpose(popul), (/nbit, tvar/))
! count = 0
! do i = 1, tvar
!     tempPool(i,1:nbit) = popul(count+1:count+nbit)
!     count = count+nbit
! end do

do j = nbit-1, 0, -1
    quant(j) = 2**j/(2**nbit-1)
end do

tempPool = matmul(rePool,quant)

do i = 1, nvar
    par(i:2:npop*nvar) = ((tempPool(i:2:npop*nvar))*(bounds(i,2)-bounds(i,1)) + bounds(i,1))
end do

pool = reshape(par,(/nvar, npop/))

END SUBROUTINE  binary2cont