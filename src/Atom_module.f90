module atom_module
use deg_rad
implicit none

integer, parameter :: realkind = 8

private
public atom, bond, bond_angle, molecule, read_atom, write_atom, print_atom, create_molecule, delete_molecule

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
real(realkind) :: vector(3) = 0.
end type

type bond_angle
! Setting the types of bonds over which the angle is calulated
character(2) :: bonds_type(2) = (/'EE', 'EE'/)
integer :: atoms_indicies(3)
type (bond) :: bonds(2)
! The actual angle of the bonds
real(realkind) :: angle = 0.
end type

type torsion
type (bond) :: end_bonds(2)
type (bond) :: center_bond
real(realkind) :: angle = 0.
end type



type molecule
type (atom), allocatable :: atoms(:)
type (bond), allocatable :: bonds(:)
type (bond_angle), allocatable :: angles(:)
type (torsion), allocatable :: torsion_angles(:)
real(realkind), allocatable :: distance(:,:)
logical, allocatable :: bonding(:,:)
end type


contains


subroutine create_molecule(mol)
type (molecule), intent(inout) :: mol

call bonds_atom(mol)
call angle_bonds(mol)
call angle_torsion(mol)

end subroutine create_molecule


subroutine delete_molecule(mol)
type (molecule), intent(inout) :: mol

if (count(mol%bonds%type == 'CC') > 0) deallocate(mol%torsion_angles)

deallocate(mol%bonds)
deallocate(mol%bonding)
deallocate(mol%distance)
deallocate(mol%angles)

end subroutine


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

!! Originaly used size - 1 but to keep competibility with cyclo compounds removed no known bugs
allocate(mol%bonds(size(mol%atoms)))


allocate(mol%distance(size(mol%atoms), size(mol%atoms)))
mol%distance = 0.
allocate(mol%bonding(size(mol%atoms), size(mol%atoms)))
mol%bonding = .false.

!! Can be optimized
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
            mol%bonds(k)%length = mol%distance(j,i)
            mol%bonds(k)%vector = mol%atoms(i)%cords - mol%atoms(j)%cords
            mol%bonding(j,i) = .true.
            mol%bonds(k)%type = 'CC'
            k = k + 1
        else if (mol%atoms(i)%element == 'C' .and. mol%atoms(j)%element == 'H' .and. abs(CH - mol%distance(i,j)) < 0.1) then
            mol%bonds(k)%link(1) = i
            mol%bonds(k)%link(2) = j
            mol%bonds(k)%length = mol%distance(j,i)
            mol%bonds(k)%vector = mol%atoms(i)%cords - mol%atoms(j)%cords
            mol%bonding(j,i) = .true.
            mol%bonds(k)%type = 'CH'
            k = k + 1
        end if
    enddo
enddo

end subroutine bonds_atom



subroutine angle_bonds(mol)
type (molecule) :: mol
type (bond) :: bonds_holder(4,count(mol%atoms%element == 'C'))
integer :: carbon_indicies(count(mol%atoms%element == 'C'))
integer :: n, i, j, k, l

allocate(mol%angles(count(mol%atoms%element == 'C') * 6))

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
do l = 1,size(carbon_indicies)
    do i = 1,3
        do j = i + 1,4
            mol%angles(k)%bonds_type(1) = bonds_holder(i,l)%type
            mol%angles(k)%bonds_type(2) = bonds_holder(j,l)%type

            mol%angles(k)%atoms_indicies(1) = bonds_holder(i,l)%link(2)
            mol%angles(k)%atoms_indicies(2) = bonds_holder(i,l)%link(1)
            mol%angles(k)%atoms_indicies(3) = bonds_holder(j,l)%link(2)

            mol%angles(k)%bonds(1) = bonds_holder(i,l)
            mol%angles(k)%bonds(2) = bonds_holder(j,l)

            mol%angles(k)%angle = rad_to_deg(acos(dot_product(mol%angles(k)%bonds(1)%vector, mol%angles(k)%bonds(2)%vector) &
            / (mol%angles(k)%bonds(1)%length * mol%angles(k)%bonds(2)%length)))

            if (mol%angles(k)%angle <= 90.) mol%angles(k)%angle = 180.0 - mol%angles(k)%angle

            k = k + 1
        enddo
    enddo
enddo

end subroutine angle_bonds



subroutine angle_torsion(mol)
type (molecule), intent(inout) :: mol
type (bond), allocatable :: end_bonds_holder(:,:,:), CC_bonds_holder(:)
real(realkind) :: n1(3), n2(3)
integer :: i, j, k, l
if (count(mol%bonds%type == 'CC') > 0) then

allocate(mol%torsion_angles(count(mol%bonds%type == 'CC') * 9))
allocate(end_bonds_holder(3, 2, count(mol%bonds%type == 'CC')))
allocate(CC_bonds_holder(count(mol%bonds%type == 'CC')))

k = 1
do i = 1, size(mol%bonds)
        if (mol%bonds(i)%type == 'CC') then
            CC_bonds_holder(k) = mol%bonds(i)
            k = k + 1
        end if
enddo

do i = 1, size(CC_bonds_holder)
    k = 1
    l = 1
    do j = 1,size(mol%bonds)
        if ((mol%bonds(j)%link(1) == CC_bonds_holder(i)%link(1) .or. mol%bonds(j)%link(2) == CC_bonds_holder(i)%link(1)) .and. &
        (mol%bonds(j)%vector(1) /= CC_bonds_holder(i)%vector(1) .or. mol%bonds(j)%vector(2) /= CC_bonds_holder(i)%vector(2))) then
            end_bonds_holder(k,1,i) = mol%bonds(j)
            k = k + 1
        end if

        if ((mol%bonds(j)%link(1) == CC_bonds_holder(i)%link(2) .or. mol%bonds(j)%link(2) == CC_bonds_holder(i)%link(2)) .and. &
        (mol%bonds(j)%vector(1) /= CC_bonds_holder(i)%vector(1) .or. mol%bonds(j)%vector(2) /= CC_bonds_holder(i)%vector(2))) then
            end_bonds_holder(l,2,i) = mol%bonds(j)
            l = l + 1
        end if
    enddo
enddo

l = 1
do i = 1, size(CC_bonds_holder)
    do j = 1,3
        do k = 1,3
            mol%torsion_angles(l)%center_bond = CC_bonds_holder(i)

            mol%torsion_angles(l)%end_bonds(1) = end_bonds_holder(j,1,i)
            mol%torsion_angles(l)%end_bonds(2) = end_bonds_holder(k,2,i)

            l = l + 1
        enddo 
    enddo
enddo


do i = 1, size(mol%torsion_angles)
    n1 = cross_product(mol%torsion_angles(i)%end_bonds(1)%vector, mol%torsion_angles(i)%center_bond%vector)
    n2 = cross_product(mol%torsion_angles(i)%end_bonds(2)%vector, mol%torsion_angles(i)%center_bond%vector)

    mol%torsion_angles(i)%angle = rad_to_deg(acos(dot_product(n1, n2) / (vector_distance(n1) * vector_distance(n2))))
enddo

end if
end subroutine angle_torsion



end module