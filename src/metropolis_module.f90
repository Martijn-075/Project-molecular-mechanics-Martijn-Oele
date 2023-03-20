module metropolis_module
use atom_module
use energy_module
implicit none

integer, parameter :: realkind = 8

private
public random_atom_metropolis, metropolis

real(realkind), parameter :: kb = 1.987204259e-3
real(realkind), parameter :: T = 293

contains

subroutine random_atom_metropolis(mol, r)
type (molecule), intent(inout) :: mol
real(realkind), intent(in), optional :: r
real(realkind) :: r_local
real(realkind), allocatable :: q(:,:)
integer :: i

if (present(r)) then
    r_local = r
else 
    r_local = 0.001
end if


allocate(q(size(mol%atoms),3))

call random_number(q)
q = (q - 0.5) * 2

q = q * r_local

do i = 1,size(mol%atoms)
    mol%atoms(i)%cords = mol%atoms(i)%cords + q(i,:)
enddo

end subroutine random_atom_metropolis



subroutine metropolis(mol)
type (molecule), intent(inout) :: mol
type (molecule) :: old_mol
real(realkind) :: energy, old_energy, delta_energy, starting_energy, P, Pa
integer :: i

call create_molecule(mol)
old_energy = forcefield_energy(mol)
starting_energy = old_energy

i = 1

do while (i < 1000)
    old_mol = mol
    call delete_molecule(mol)
    call random_atom_metropolis(mol)
    call create_molecule(mol)
    energy = forcefield_energy(mol)
    delta_energy = energy - old_energy


    if (delta_energy < 0.) then
        old_energy = energy
        i = 1
    else 
        call random_number(P)
        Pa = min(1.,exp((-1 * delta_energy) / kb * T))

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