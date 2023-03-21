!
! Module to convert rads to degrees and degrees to rads. pi is also declared and public
!
module math_module
implicit none

integer, parameter :: realkind = 8
real(realkind), parameter :: pi = 4.*atan(1.)

private
public pi, deg_to_rad, rad_to_deg, cross_product, vector_lenght


contains

! Converting degrees to radian
real(realkind) function deg_to_rad(deg) result(rad)
real(realkind) :: deg

rad = deg * (pi/180)

end function deg_to_rad

! Converting radian to degrees
real(realkind) function rad_to_deg(rad) result(deg)
real(realkind) :: rad

deg = rad * (180/pi)

end function rad_to_deg

! Calculating the cross product of two 3D vectors
function cross_product(vec1, vec2) result(result)
real(realkind), intent(in) :: vec1(3), vec2(3)
real(realkind) :: result(3)

result(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
result(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
result(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)

end function 

! Calculating the length of a vector
real(realkind) function vector_lenght(vector)
real(realkind), intent(in) :: vector(3)

vector_lenght = sqrt(dot_product(vector, vector))

end function

end module