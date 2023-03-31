!
! Metropolis module 
! Module containing the metropolis algorithm to minimize the energy of a molecular conformation  
! Author: Martijn Oele (GitHub: Martijn-075)
!
module metropolis_module
use molecule_module
use energy_module
use constant_module
implicit none

private
public random_atom_metropolis, metropolis, minimize_energy

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
    r_local = 0.0001
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
subroutine metropolis(mol, T, r)
type (molecule), intent(inout) :: mol
type (molecule) :: old_mol
real(realkind) :: energy, old_energy, delta_energy, starting_energy, P, Pa, T_local
real(realkind), optional :: T, r
integer :: i, j, k

! Checking if T is given by the user
if (present(T)) then
    T_local = T
else 
    T_local = 293
end if 

! Getting the energy of the unchnaged molecule 
call create_molecule(mol)
old_energy = forcefield_energy(mol)
starting_energy = old_energy

! The metropolis algorithm. If the new cordinates are rejected 1000 times in a row it is assumed that an minimal energy is found
i = 0
j = 0
k = 0
do while (i < 1000)
    ! Creating a new molecule with the atoms randomly moved
    old_mol = mol
    call delete_molecule(mol)
    call random_atom_metropolis(mol, r)
    call create_molecule(mol)
    energy = forcefield_energy(mol)
    delta_energy = energy - old_energy

    ! Checking if the energy of the new (randomly moved) molecule is lower than the old molecule
    if (delta_energy < 0.) then
        old_energy = energy
        i = 1
    else 
        call random_number(P)
        Pa = min(1.,exp((-1 * delta_energy) / kb * T_local))

        ! If the energy of the new molecule is higher than that of the old molecule there is still a change that the new molecule is accepted based on the boltzzman distribution
        if (p < Pa) then
            old_energy = energy
            i = 1
            k = k + 1
        else
            mol = old_mol
            energy = old_energy
            i = i + 1
        end if
    end if

    ! If not converged whitin 1.000.000
    if (j > 1000000) then
        print *, 'Not converged within 1.000.000 iterations. Try making r bigger and/or T lower'
    end if
    j = j + 1

enddo

mol%minimized_energy = energy

print '(a,x,es12.6,a)', 'Starting energy:', starting_energy, ' kcal/mol'
print '(a,x,es12.6,a)', 'Minimized energy:', energy, ' kcal/mol'
print '(a,x,es13.6,a)', 'Energy reduced by:', energy - starting_energy, ' kcal/mol'
print '(a,x,f5.1,a)', 'Energy reduced %:', (energy - starting_energy) / starting_energy * 100, '%'
print '(a,x,i5)', 'Number of iterations:', j
print '(a,x,i1)', 'Number of acceptance with higer energy:', k

end subroutine metropolis


subroutine minimize_energy
type (molecule) :: mol
character(32) :: filename
character(1) :: ans
real(realkind) :: r, T


print *, 'File name of the molecule (.xyz)?'
read (*,*) filename
    
call read_atom(mol, filename)

print *, 'Do you want to give custom parameters for r and T? (standard r = 0.0001 and T = 293 K)'
read (*,*) ans
if (ans == 'y') then
    print *, 'T?'
    read (*,*) T
    print *, 'r?'
    read (*,*) r
    call metropolis(mol, T, r)
else 
    call metropolis(mol)
end if


print *, 'Do you want to save the minimized energy configuration (y/n)'
read (*,*) ans

if (ans == 'y') then
    print *, 'File name(.xyz)?'
    read (*,*) filename
    call write_atom(mol, filename)
end if

end subroutine minimize_energy

end module