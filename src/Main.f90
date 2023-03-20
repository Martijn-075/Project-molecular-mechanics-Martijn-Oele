program main
use atom_module
use energy_module
use metropolis_module
use deg_rad
implicit none

type (molecule) :: mol

call read_atom(mol, 'c6h12.xyz')

call metropolis(mol)

end program