!
! Module to convert rads to degrees and degrees to rads. pi is also declared and public
!
module deg_rad
implicit none

integer, parameter :: realkind = 8
real(realkind), parameter :: pi = 4.*atan(1.)

private
public pi, deg_to_rad, rad_to_deg



contains

real(realkind) function deg_to_rad(deg) result(rad)
real(realkind) :: deg

rad = deg * (pi/180)

end function deg_to_rad



real(realkind) function rad_to_deg(rad) result(deg)
real(realkind) :: rad

deg = rad * (180/pi)

end function rad_to_deg

end module