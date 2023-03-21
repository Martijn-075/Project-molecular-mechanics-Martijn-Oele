!
! Module containing the metropolis algorithm to minimize the energy of a molecular conformation  
! Author: Martijn Oele (GitHub: Martijn-075)
!
module metropolis_module
use molecule_module
use energy_module
use constant_module
implicit none

private
public random_atom_metropolis, metropolis

real(realkind), parameter :: kb = 1.987204259e-3
real(realkind), parameter :: T = 293

contains

! Moving all the atoms at random whitin a certain radius 
subroutine random_atom_metropolis(mol, r)
type (molecule), intent(inout) :: mol
real(realkind), intent(in), optional :: r
real(realkind) :: r_local
real(realkind), allocatable :: q(:,:)
integer :: i

! Checking if r is specified by the user
if (present(r)) then
    r_local = r
else 
    r_local = 0.001
end if


allocate(q(size(mol%atoms),3))

! Creating a random array with real numbers between -1 and 1
call random_number(q)
q = (q - 0.5) * 2

q = q * r_local

! Adding the random real numbers to the atoms cordinates
do i = 1,size(mol%atoms)
    mol%atoms(i)%cords = mol%atoms(i)%cords + q(i,:)
enddo

end subroutine random_atom_metropolis

! The actuale metropolis algorithm to minimize the enrgy by randomly moving the atoms whitin a certain range
subroutine metropolis(mol)
type (molecule), intent(inout) :: mol
type (molecule) :: old_mol
real(realkind) :: energy, old_energy, delta_energy, starting_energy, P, Pa
integer :: i

! Getting the energy of the unchnaged molecule 
call create_molecule(mol)
old_energy = forcefield_energy(mol)
starting_energy = old_energy

print *, 'Stretching energy', stretch_energy(mol)
print *, 'Bending energy', bending_energy(mol)
print *, 'Torsion energy', torsion_energy(mol)
print *, 'Electrstatic energy', electrostatic_energy(mol)
print *, 'Van der waals energy', van_der_waals_energy(mol)
print *, ''

! The metropolis algorithm. If the new cordinates are rejected 1000 times in a row it is assumed that an minimal energy is found
i = 1
do while (i < 1000)
    ! Creating a new moelcule with the atoms randomly moved
    old_mol = mol
    call delete_molecule(mol)
    call random_atom_metropolis(mol)
    call create_molecule(mol)
    energy = forcefield_energy(mol)
    delta_energy = energy - old_energy

    ! Checking if the energy of the new (randomly moved) molecule is lower than the old molecule
    if (delta_energy < 0.) then
        old_energy = energy
        i = 1
    else 
        call random_number(P)
        Pa = min(1.,exp((-1 * delta_energy) / kb * T))

        ! If the energy of the new molecule is higher than that of the old molecule there is still a change that the new molecule is accepted based on the boltzzman distribution
        if (p < Pa) then
            old_energy = energy
            i = 1
        else
            mol = old_mol
            energy = old_energy
            i = i + 1
        end if
    end if
enddo


print *, 'Starting energy', starting_energy
print *, 'minimized energy', energy
print *, 'Energy reduced by', energy - starting_energy
print *, 'Energy reduced %', (energy - starting_energy) / starting_energy * 100


end subroutine metropolis

end module