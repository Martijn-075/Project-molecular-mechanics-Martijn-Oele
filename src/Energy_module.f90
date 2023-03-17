module energy_module
use atom_module
implicit none

private
public stretch_energy, bending_energy

integer, parameter :: realkind = 8
real(realkind), parameter :: k_CC = 310.
real(realkind), parameter :: k_CH = 340.
real(realkind), parameter :: R_CC = 1.526
real(realkind), parameter :: R_CH = 1.090

real(realkind), parameter :: k_CC_CC = 40.0
real(realkind), parameter :: k_CC_CH = 50.0
real(realkind), parameter :: k_CH_CH = 35.0
real(realkind), parameter :: theta_0_SP3 = 109.50


contains


real(realkind) function stretch_energy(mol) result(E)
type (molecule) :: mol
real(realkind) :: E_bond, E_sum
integer :: i
E_sum = 0.

do i = 1,size(mol%bonds) - count(mol%bonds%type == 'EE')
    if (mol%bonds(i)%type == 'CC') then
        E_bond = k_CC * (mol%bonds(i)%length - R_CC)**2
    else if (mol%bonds(i)%type == 'CH') then
        E_bond = k_CH * (mol%bonds(i)%length - R_CH)**2
    end if
    E_sum = E_sum + E_bond
enddo

E = E_sum

end function stretch_energy



real(realkind) function bending_energy(mol) result(E)
type(molecule) :: mol
real(realkind) :: E_bend, E_sum
integer :: i

E_sum = 0.

do i = 1, size(mol%angles)
    if (mol%angles(i)%bonds(1) == 'CH' .and. mol%angles(i)%bonds(2) == 'CH') then
        E_bend = k_CH_CH * (mol%angles(i)%angle - theta_0_SP3)**2
    else if (mol%angles(i)%bonds(1) == 'CC' .and. mol%angles(i)%bonds(2) == 'CH') then
        E_bend = k_CC_CH * (mol%angles(i)%angle - theta_0_SP3)**2
    else if (mol%angles(i)%bonds(1) == 'CC' .and. mol%angles(i)%bonds(2) == 'CC') then
        E_bend = k_CC_CC * (mol%angles(i)%angle - theta_0_SP3)**2
    end if

    E_sum = E_sum + E_bend
enddo

E = E_sum

end function bending_energy



end module