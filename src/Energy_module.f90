module energy_module
use atom_module
implicit none

private
public stretch_energy

integer, parameter :: realkind = 8
real(realkind), parameter :: k_CC = 310.
real(realkind), parameter :: k_CH = 340.
real(realkind), parameter :: R_CC = 1.526
real(realkind), parameter :: R_CH = 1.090

contains


real(realkind) function stretch_energy(bonds) result(E)
type (bond) :: bonds(:)
real(realkind) :: E_bond, E_sum
integer :: i
E_sum = 0.

do i = 1,size(bonds) - count(bonds%type == 'EE')
    if (bonds(i)%type == 'CC') then
        E_bond = k_CC * (bonds(i)%length - R_CC)**2
    else if (bonds(i)%type == 'CH') then
        E_bond = k_CH * (bonds(i)%length - R_CH)**2
    end if
    E_sum = E_sum + E_bond
enddo

E = E_sum

end function stretch_energy


real(realkind) function bending_energy(bonds, atoms) result(E)
type(bond) :: bonds(:)
type (atom) :: atoms(:)



end function bending_energy



end module