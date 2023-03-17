module atom_module
use deg_rad
implicit none

integer, parameter :: realkind = 8

private
public atom, bond, bond_angle, molecule, read_atom, write_atom, print_atom, bonds_atom, angle_bonds

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
end type

type bond_angle
! Setting the types of bonds over which the angle is calulated
character(2) :: bonds(2) = (/'EE', 'EE'/)
! Setting the indicies of the atoms over wich the angle is calcualted
integer :: atom_indicies(3) = 0
! The actual angle of the bonds
real(realkind) :: angle = 0.
end type

type molecule
type (atom), allocatable :: atoms(:)
type (bond), allocatable :: bonds(:)
type (bond_angle), allocatable :: angles(:)
real(realkind), allocatable :: distance(:,:)
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
real(realkind) :: CC, CH

allocate(mol%bonds(size(mol%atoms) - 1))

allocate(mol%distance(size(mol%atoms), size(mol%atoms)))

do i = 1, size(mol%atoms)
    do j = 1,size(mol%atoms)
        mol%distance(i,j) = sqrt((mol%atoms(i)%cords(1) - mol%atoms(j)%cords(1))**2 + &
        (mol%atoms(i)%cords(2) - mol%atoms(j)%cords(2))**2 + (mol%atoms(i)%cords(3) - mol%atoms(j)%cords(3))**2)
    enddo
enddo

CC = 1.535
CH = 1.094

k = 1
do i = 1, size(mol%atoms) - 1
    do j = i + 1, size(mol%atoms)
        if (mol%atoms(i)%element == 'C' .and. mol%atoms(j)%element == 'C' .and. abs(CC - mol%distance(i,j)) < 0.1) then
            mol%bonds(k)%link(1) = i
            mol%bonds(k)%link(2) = j
            mol%bonds(k)%length = mol%distance(i,j)
            mol%bonds(k)%type = 'CC'
            k = k + 1
        else if (mol%atoms(i)%element == 'C' .and. mol%atoms(j)%element == 'H' .and. abs(CH - mol%distance(i,j)) < 0.1) then
            mol%bonds(k)%link(1) = i
            mol%bonds(k)%link(2) = j
            mol%bonds(k)%length = mol%distance(i,j)
            mol%bonds(k)%type = 'CH'
            k = k + 1
        end if
    enddo
enddo

end subroutine bonds_atom



subroutine angle_bonds(mol)
type (molecule) :: mol
type (bond) :: bonds_holder(4,count(mol%atoms%element == 'C'))
real(realkind) :: A, B, C
integer :: carbon_indicies(count(mol%atoms%element == 'C'))
integer :: n, i, j, k, l

n = count(mol%atoms%element == 'C') * 6

allocate(mol%angles(n))

k = 1
do j =1,size(mol%atoms)
    if (mol%atoms(j)%element == 'C') then
        carbon_indicies(k) = j
        k = k + 1
    end if
enddo


do i = 1, size(carbon_indicies)
    k = 1
    do j = 1,size(mol%bonds)
        if (mol%bonds(j)%link(1) == carbon_indicies(i) .or. mol%bonds(j)%link(2) == carbon_indicies(i)) then
            bonds_holder(k,i) = mol%bonds(j)
            k = k + 1
        end if
    enddo
end do


k = 1
do l = 1,size(bonds_holder, 2)
    do i = 1,3
        do j = i + 1,4
            mol%angles(k)%bonds(1) = bonds_holder(i,l)%type
            mol%angles(k)%bonds(2) = bonds_holder(j,l)%type

            if (bonds_holder(i,l)%link(1) == bonds_holder(j,l)%link(1)) then
                mol%angles(k)%atom_indicies(1) = bonds_holder(i,l)%link(2)
                mol%angles(k)%atom_indicies(2) = bonds_holder(i,l)%link(1)
                mol%angles(k)%atom_indicies(3) = bonds_holder(j,l)%link(2)
            else if (bonds_holder(i,l)%link(1) == bonds_holder(j,l)%link(2)) then
                mol%angles(k)%atom_indicies(1) = bonds_holder(i,l)%link(2)
                mol%angles(k)%atom_indicies(2) = bonds_holder(i,l)%link(1)
                mol%angles(k)%atom_indicies(3) = bonds_holder(j,l)%link(1)
            else if (bonds_holder(i,l)%link(2) == bonds_holder(j,l)%link(1)) then
                mol%angles(k)%atom_indicies(1) = bonds_holder(i,l)%link(1)
                mol%angles(k)%atom_indicies(2) = bonds_holder(i,l)%link(2)
                mol%angles(k)%atom_indicies(3) = bonds_holder(j,l)%link(2)
            else
                mol%angles(k)%atom_indicies(1) = bonds_holder(i,l)%link(1)
                mol%angles(k)%atom_indicies(2) = bonds_holder(i,l)%link(2)
                mol%angles(k)%atom_indicies(3) = bonds_holder(j,l)%link(1)
           end if

            A = mol%distance(mol%angles(k)%atom_indicies(1), mol%angles(k)%atom_indicies(2))
            B = mol%distance(mol%angles(k)%atom_indicies(2), mol%angles(k)%atom_indicies(3))
            C = mol%distance(mol%angles(k)%atom_indicies(1), mol%angles(k)%atom_indicies(3))
            mol%angles(k)%angle =rad_to_deg(acos((A**2 + B**2 - C**2) / (2 * A * B)))
            k = k + 1
        enddo
    enddo
enddo

end subroutine angle_bonds



end module