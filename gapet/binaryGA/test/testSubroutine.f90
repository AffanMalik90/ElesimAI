PROGRAM testSubroutine
IMPLICIT NONE

!! binary2cont variables
real(kind=8), dimension(3,4) :: popul
real(kind=8), dimension(2,2) :: bounds
real(kind=8), dimension(2,3) :: pool

!! binary2cont
popul = reshape((/ 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0 /), shape(popul))
bounds = reshape((/ 0, 6, 0, 6/), shape(bounds))

write(*,*) popul

! CALL binary2cont(popul,bounds,pool)

END PROGRAM testSubroutine