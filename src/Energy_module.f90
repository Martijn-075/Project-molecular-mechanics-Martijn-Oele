!
! Energy module
! Module contianing the energy functions (bonding and nonbonding) to calculate the molecular forcefield energy
! Author: Martijn Oele (GitHub: Martijn-075)
!
module energy_module
use molecule_module
use math_module
use constant_module
implicit none

private
public stretch_energy, bending_energy, van_der_waals_energy, electrostatic_energy, torsion_energy, forcefield_energy


contains

! The total forcield energy function calling all sub energy functions
real(realkind) function forcefield_energy(mol) result(E)
type (molecule), intent(inout) :: mol

E = stretch_energy(mol) + bending_energy(mol) + torsion_energy(mol) + electrostatic_energy(mol) + van_der_waals_energy(mol)

end function

! The strech energy function. The bonds_angle subroutine must be called before this function
real(realkind) function stretch_energy(mol) result(E)
type (molecule) :: mol
real(realkind) :: E_bond, E_sum
integer :: i

E_sum = 0.

! Calculating the strech energy for every bond based on the bond type
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

! Calculating the bending energy for all conecting bonds
real(realkind) function bending_energy(mol) result(E)
type(molecule) :: mol
real(realkind) :: E_bend, E_sum
integer :: i

E_sum = 0.

! Calculating the bending energy based on the two bond types
do i = 1, size(mol%angles)
    if (mol%angles(i)%bonds(1)%type == 'CH' .and. mol%angles(i)%bonds(2)%type == 'CH') then
        E_bend = k_CH_CH * (mol%angles(i)%angle - theta_0_SP3)**2
    else if (mol%angles(i)%bonds(1)%type == 'CC' .and. mol%angles(i)%bonds(2)%type == 'CC') then
        E_bend = k_CC_CC * (mol%angles(i)%angle - theta_0_SP3)**2
    else
        E_bend = k_CC_CH * (mol%angles(i)%angle - theta_0_SP3)**2
    end if

    E_sum = E_sum + E_bend
enddo

E = E_sum

end function bending_energy

! Calcualting the torsion energy. There must be an CC bond else the torsion energy = 0.0
real(realkind) function torsion_energy(mol) result(E)
type (molecule), intent(inout) :: mol 
real(realkind) :: E_torsion, E_sum
integer :: i

E_sum = 0.

! Checking if there is a CC bond that can function as a centeral bond
! If CC bond present the torsion energy is calculated
if (count(mol%bonds%type == 'CC') > 0) then
    do i = 1,size(mol%torsion_angles)
        E_torsion = V2_CC * (1 + cos(deg_to_rad(n_CC * mol%torsion_angles(i)%angle - gamma_CC)))
        E_sum = E_sum + E_torsion
    enddo
end if

E = E_sum

end function torsion_energy

! Calculating the non bonding van der waals energy
real(realkind) function van_der_waals_energy(mol) result(E)
type (molecule), intent(inout) :: mol
real(realkind) :: E_non_bonding, E_sum
integer :: i, j

E_sum = 0.

! Calcualting the nonbonding van der waals energy based on the two interacting atoms
! Only non bonding atoms at least 4.5 A (3 CC bonds) and not further than 9 A are taken into account
do i = 1, size(mol%distance, 1) - 1
    do j = i + 1, size(mol%distance, 1)
        if (mol%bonding(j,i) .or. mol%distance(j,i) < 4.5 .or. mol%distance(j,i) > 9) cycle

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

! Calcualting the nonbonding electrastatic energy
real(realkind) function electrostatic_energy(mol) result(E)
type (molecule), intent(inout) :: mol
real(realkind) :: E_electrostatic, E_sum
integer :: i, j

E_sum = 0.

! Calcualting the electrastatic energy based on the two interacting atoms
! Only non bonding atoms at least 4.5 A (3 CC bonds) and not further than 9 A are taken into account
do i = 1, size(mol%distance, 1) - 1
    do j = i + 1, size(mol%distance, 1)
        if (mol%bonding(j,i) .or. mol%distance(j,i) < 4.5 .or. mol%distance(j,i) > 9) cycle

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