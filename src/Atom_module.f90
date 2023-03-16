module atom_module
implicit none

integer, parameter :: realkind = 8

private
public atom, bond, bond_angle, molecule, read_atom, write_atom, print_atom, bonds_atom

type atom
! Setting the element in this project only C (carbon) or H (hydrogen)
character(1) :: element = 'E'
! The xyz coordinates of the atom in the molecule
real(realkind) :: cords(3) = 0.
end type

type bond
! Setting the indices of the atoms that are bonding
integer :: link(2) = 0
! Setting the type of bond in this project only CC or CH
character(2) :: type = 'EE'
! Setting the calculated distance of the bond in A
real(realkind) :: length = 0.
! Setting the vector for the bond
real(realkind) :: vector(3) = 0
end type

type bond_angle
! Setting the types of bonds over which the angle is calulated
character(2) :: bonds(2) = (/'EE', 'EE'/)
! Setting the indicies of the atoms over wich the angle is calcualted
integer :: atom_indicies(3)
! The actual angle of the bonds
real(realkind) :: angle = 0.
end type

type molecule
type (atom), allocatable :: atoms(:)
type (bond), allocatable :: bonds(:)
type (bond_angle), allocatable :: angles(:)
end type


contains


subroutine read_atom(mol, filename)
type (molecule) :: mol
character(len=*), intent(in) :: filename
integer :: i, iu, ios, n

open(newunit=iu, file=filename, status='old', iostat=ios)
read(iu,*) n
read(iu,*)
allocate(mol%atoms(n))

if (ios == 0) then
    do i = 1,n
        read(iu,*) mol%atoms(i)%element, mol%atoms(i)%cords
    enddo
else
    print *, 'error', ios
end if
close(iu)

end subroutine read_atom



subroutine write_atom(mol, filename)
type (molecule) :: mol
character(len=*), intent(in) :: filename
integer i, iu

open(newunit=iu, file=filename, status='replace')
write(iu,*) size(mol%atoms)
write(iu,*)

do i = 1, size(mol%atoms)
    write(iu,*) mol%atoms(i)
enddo

close(iu)
end subroutine



subroutine print_atom(atoms)
type (atom) :: atoms(:)
integer :: i

do i = 1,size(atoms)
    print *, i, atoms(i)
enddo

end subroutine print_atom



subroutine bonds_atom(mol)
type (molecule) :: mol
integer :: i, j, k
real(realkind) :: CC, CH, distance

allocate(mol%bonds(size(mol%atoms) - 1))

CC = 1.535
CH = 1.094

k = 1
do i = 1, size(mol%atoms) - 1
    do j = i + 1, size(mol%atoms)
        distance = sqrt((mol%atoms(i)%cords(1) - mol%atoms(j)%cords(1))**2 + (mol%atoms(i)%cords(2) - mol%atoms(j)%cords(2))**2 &
        + (mol%atoms(i)%cords(3) - mol%atoms(j)%cords(3))**2)
        if (mol%atoms(i)%element == 'C' .and. mol%atoms(j)%element == 'C' .and. abs(CC - distance) < 0.1) then
            mol%bonds(k)%link(1) = i
            mol%bonds(k)%link(2) = j
            mol%bonds(k)%length = distance
            mol%bonds(k)%vector = mol%atoms(i)%cords - mol%atoms(j)%cords
            mol%bonds(k)%type = 'CC'
            k = k + 1
        else if (mol%atoms(i)%element == 'C' .and. mol%atoms(j)%element == 'H' .and. abs(CH - distance) < 0.1) then
            mol%bonds(k)%link(1) = i
            mol%bonds(k)%link(2) = j
            mol%bonds(k)%length = distance
            mol%bonds(k)%vector = mol%atoms(i)%cords - mol%atoms(j)%cords
            mol%bonds(k)%type = 'CH'
            k = k + 1
        end if
    enddo
enddo

end subroutine bonds_atom



subroutine angle_bonds(mol)
type (molecule) :: mol
integer :: n, i, k

n = count(mol%atoms%element == 'C') * 6

allocate(mol%angles(n))

k = 1




end subroutine angle_bonds


end module