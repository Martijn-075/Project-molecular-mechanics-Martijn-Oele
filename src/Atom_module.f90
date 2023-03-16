module atom_module
implicit none

integer, parameter :: realkind = 8

private
public atom, bond, read_atom, write_atom, print_atom, bonds_atom

type atom
! Setting the element in this project only C (carbon) or H (hydrogen)
character(1) :: element = 'E'
! The xyz coordinates of the atom in the molecule
real(realkind) :: x, y, z = 0.
end type

type bond
! Setting the indices of the atoms that are bonding
integer :: link(2) = 0
! Setting the type of bond in this project only CC or CH
character(2) :: type = 'EE'
! Setting the calculated distance of the bond in A
real(realkind) :: length = 0.
end type

contains


subroutine read_atom(atoms, filename)
type (atom), allocatable, intent(inout) :: atoms(:)
character(len=*), intent(in) :: filename
integer :: i, iu, ios, n

open(newunit=iu, file=filename, status='old', iostat=ios)
read(iu,*) n
read(iu,*)
allocate(atoms(n))

if (ios == 0) then
    do i = 1,n
        read(iu,*) atoms(i)%element, atoms(i)%x, atoms(i)%y, atoms(i)%z
    enddo
else
    print *, 'error', ios
end if
close(iu)

end subroutine read_atom



subroutine write_atom(atoms, filename)
type (atom), intent(in) :: atoms(:)
character(len=*), intent(in) :: filename
integer i, iu

open(newunit=iu, file=filename, status='replace')
write(iu,*) size(atoms)
write(iu,*)

do i = 1, size(atoms)
    write(iu,*) atoms(i)
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



subroutine bonds_atom(atoms, bonds)
type (atom) :: atoms(:)
type (bond), allocatable, intent(inout) :: bonds(:)
integer :: i, j, k
real(realkind) :: CC, CH, distance

allocate(bonds(size(atoms) - 1))

CC = 1.535
CH = 1.094

k = 1
do i = 1, size(atoms) - 1
    do j = i + 1, size(atoms)
        distance = sqrt((atoms(i)%x - atoms(j)%x)**2 + (atoms(i)%y - atoms(j)%y)**2 + (atoms(i)%z - atoms(j)%z)**2)
        if (atoms(i)%element == 'C' .and. atoms(j)%element == 'C' .and. abs(CC - distance) < 0.1) then
            bonds(k)%link(1) = i
            bonds(k)%link(2) = j
            bonds(k)%length = distance
            bonds(k)%type = 'CC'
            k = k + 1
        else if (atoms(i)%element == 'C' .and. atoms(j)%element == 'H' .and. abs(CH - distance) < 0.1) then
            bonds(k)%link(1) = i
            bonds(k)%link(2) = j
            bonds(k)%length = distance
            bonds(k)%type = 'CH'
            k = k + 1
        end if
    enddo
enddo

end subroutine bonds_atom


end module