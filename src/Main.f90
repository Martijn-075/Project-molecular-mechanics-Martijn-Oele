program main
use atom_module
implicit none
type(atom), allocatable :: atoms(:)
type (bond), allocatable :: bonds(:)
integer i


call read_atom(atoms, 'c4h10.xyz')

call bonds_atom(atoms, bonds)


call print_atom(atoms)

do i = 1,size(bonds)
    print *, bonds(i)%link(1), bonds(i)%link(2), bonds(i)%length, bonds(i)%type
    
enddo

end program