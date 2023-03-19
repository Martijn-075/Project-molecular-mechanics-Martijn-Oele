program main
use atom_module
use energy_module
use deg_rad
implicit none

type (molecule) :: mol
integer i, j



call read_atom(mol, 'c4h10.xyz')

call print_atom(mol%atoms)

call bonds_atom(mol)

call angle_bonds(mol)


do i = 1,size(mol%bonds)
    print *, mol%bonds(i)
enddo

print *, 'Stretch energy: ', stretch_energy(mol)

print *, 'Bending energy: ', bending_energy(mol)

print *, 'Van der waals: ', van_der_waals_energy(mol)

print *, 'Electrostatic: ', electrostatic_energy(mol)

do i = 1,size(mol%angles)
    print *, mol%angles(i)%bonds, mol%angles(i)%atoms_indicies, mol%angles(i)%angle
enddo


end program