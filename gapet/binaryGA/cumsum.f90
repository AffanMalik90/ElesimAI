FUNCTION cumsum(avgprob)
integer :: i
real(kind=8), intent(in) :: avgprob(:)
real(kind=8) :: cumsum(size(avgprob))

cumsum(:) = [(sum(avgprob(1:i)), i = 1, size(avgprob))]

END FUNCTION cumsum