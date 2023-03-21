!
! Author: Martijn Oele (GitHub: Martijn-075)
!
program main
use molecule_module
use energy_module
use metropolis_module
use constant_module
implicit none

type (molecule) :: mol

call read_atom(mol, 'c6h12.xyz')


call metropolis(mol)

print *, ""
print *, 'Stretching energy', stretch_energy(mol)
print *, 'Bending energy', bending_energy(mol)
print *, 'Torsion energy', torsion_energy(mol)
print *, 'Electrstatic energy', electrostatic_energy(mol)
print *, 'Van der waals energy', van_der_waals_energy(mol)

end program