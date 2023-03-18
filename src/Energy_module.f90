module energy_module
use atom_module
implicit none

private
public stretch_energy, bending_energy, van_der_waals_energy, electrostatic_energy

integer, parameter :: realkind = 8
! Stretching parameters
real(realkind), parameter :: k_CC = 310.
real(realkind), parameter :: k_CH = 340.
real(realkind), parameter :: R_CC = 1.526
real(realkind), parameter :: R_CH = 1.090
! Bending parameters
real(realkind), parameter :: k_CC_CC = 40.0
real(realkind), parameter :: k_CC_CH = 50.0
real(realkind), parameter :: k_CH_CH = 35.0
real(realkind), parameter :: theta_0_SP3 = 109.50
! Electrostatic parameters
real(realkind), parameter :: q_H = 0.078
real(realkind), parameter :: q_C = -0.344
real(realkind), parameter :: coulombe_constant = 331.1908 !! Needs to be checked 2.146639897e16, eV*A/e2 converted to kcal/mol
! Van der waals parameters
real(realkind), parameter :: A_vdw_HH = 0.0157 * (2 * 1.4870)**12
real(realkind), parameter :: B_vdw_HH = 0.0157 * (2 * 1.4870)**6
real(realkind), parameter :: A_vdw_CH = sqrt(0.0157 * 0.1094) * (1.4870 + 1.9080)**12
real(realkind), parameter :: B_vdw_CH = sqrt(0.0157 * 0.1094) * (1.4870 + 1.9080)**6
real(realkind), parameter :: A_vdw_CC = 0.1094 * (2 * 1.9080)**12
real(realkind), parameter :: B_vdw_CC = 0.1094 * (2 * 1.9080)**6

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



real(realkind) function van_der_waals_energy(mol) result(E)
type (molecule), intent(inout) :: mol
real(realkind) :: E_non_bonding, E_sum
integer :: i, j

E_sum = 0.

do i = 1, size(mol%distance, 1) - 1
    do j = i + 1, size(mol%distance, 1)
        if (mol%bonding(j,i)) cycle

        if (mol%atoms(j)%element == 'H' .and. mol%atoms(i)%element == 'H') then
            E_non_bonding = (A_vdw_HH / mol%distance(j,i)**12) - (B_vdw_HH / mol%distance(j,i)**6)
        else if (mol%atoms(j)%element == 'C' .and. mol%atoms(i)%element == 'C') then
            E_non_bonding = (A_vdw_CC / mol%distance(j,i)**12) - (B_vdw_CC / mol%distance(j,i)**6)
        else 
            E_non_bonding = (A_vdw_CH / mol%distance(j,i)**12) - (B_vdw_CH / mol%distance(j,i)**6)
        end if
        E_sum = E_sum + E_non_bonding
    enddo
enddo

E = E_sum

end function van_der_waals_energy


real(realkind) function electrostatic_energy(mol) result(E)
type (molecule), intent(inout) :: mol
real(realkind) :: E_electrostatic, E_sum
integer :: i, j

E_sum = 0.

do i = 1, size(mol%distance, 1) - 1
    do j = i + 1, size(mol%distance, 1)
        if (mol%bonding(j,i)) cycle

        if (mol%atoms(j)%element == 'H' .and. mol%atoms(i)%element == 'H') then
            E_electrostatic = ((q_H * q_H) / mol%distance(j,i)) * coulombe_constant
        else if (mol%atoms(j)%element == 'C' .and. mol%atoms(i)%element == 'C') then
            E_electrostatic = ((q_C * q_C) / mol%distance(j,i)) * coulombe_constant
        else 
            E_electrostatic = ((q_H * q_C) / mol%distance(j,i)) * coulombe_constant
        end if
        E_sum = E_sum + E_electrostatic
    enddo
enddo

E = E_sum

end function electrostatic_energy

end module