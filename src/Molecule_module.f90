!
! Molecule module
! Module containing the data types and subroutines for creating an molecule from atoms inclusing the bonds, bond angles and torsion angles
! Author: Martijn Oele (GitHub: Martijn-075)
!
 module molecule_module
use math_module
use constant_module
implicit none

private
public atom, bond, molecule, read_atom, write_atom, create_molecule, delete_molecule

type atom
    ! The element of the atom (in this project only C (carbon) or H (hydrogen))
    character(1) :: element = 'E'
    ! The xyz coordinates of the atom in the molecule
    real(realkind) :: cords(3) = 0.
end type

type bond
    ! The indices of the atoms that are bonding
    integer :: link(2) = 0
    ! The type of bond in this project only CC or CH
    character(2) :: type = 'EE'
    ! The calculated distance of the bond in (A) angstrom 
    real(realkind) :: length = 0.
    ! the vecotr of the bond (3D)
    real(realkind) :: vector(3) = 0.
end type

type bond_angle
    ! The bonds that have an angle
    type (bond) :: bonds(2)
    ! The angle of the bonds
    real(realkind) :: angle = 0.
end type

type torsion
    ! The bonds at the ends (1 and 3) that make the angle
    type (bond) :: end_bonds(2)
    ! The central bond that conect the end bonds. In this project always a CC bond
    type (bond) :: center_bond
    ! The angle the end bonds make acourding to the newman projection
    real(realkind) :: angle = 0.
end type



type molecule
    ! The atoms in the molecule
    type (atom), allocatable :: atoms(:)
    ! The bonds in the molecule
    type (bond), allocatable :: bonds(:)
    ! The angles between the bonds in the molecule
    type (bond_angle), allocatable :: angles(:)
    ! The torsion angles in the molecule
    type (torsion), allocatable :: torsion_angles(:)
    ! The distances between all atoms in the molecule
    real(realkind), allocatable :: distance(:,:)
    ! Suporting the distance array and holds if two molecules are form a bond
    logical, allocatable :: bonding(:,:)
    ! The minimized energy of the molecule
    real(realkind) :: minimized_energy = 0.
end type


contains

! Create a molecule after it is read in. This is needed before the forcefield energy function can be called
subroutine create_molecule(mol)
type (molecule), intent(inout) :: mol

call bonds_atom(mol)
call angle_bonds(mol)
call angle_torsion(mol)

end subroutine create_molecule

! Deallocate everything from the molecule except the atoms. Neededd to recreate the molecule after the atoms are changed
subroutine delete_molecule(mol)
type (molecule), intent(inout) :: mol

if (count(mol%bonds%type == 'CC') > 0) deallocate(mol%torsion_angles)

deallocate(mol%bonds)
deallocate(mol%bonding)
deallocate(mol%distance)
deallocate(mol%angles)

end subroutine

! Read in the atom information from an xyz file
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

! Writes the atoms information to a xyz file
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

write(iu,*)
write(iu,*) 'Minimized energy', mol%minimized_energy, 'kcal/mol'
close(iu)
end subroutine

! Establishes the bonds between atoms based on distance
subroutine bonds_atom(mol)
type (molecule) :: mol
integer :: i, j, k

allocate(mol%bonds(size(mol%atoms)))


allocate(mol%distance(size(mol%atoms), size(mol%atoms)))
mol%distance = 0.
allocate(mol%bonding(size(mol%atoms), size(mol%atoms)))
mol%bonding = .false.

! Creates the distnace array which is used for cheking if atoms are bonded but also for nonbonding interactions 
do i = 1, size(mol%atoms)
    do j = 1, size(mol%atoms)
        mol%distance(j,i) = sqrt((mol%atoms(i)%cords(1) - mol%atoms(j)%cords(1))**2 + &
        (mol%atoms(i)%cords(2) - mol%atoms(j)%cords(2))**2 + (mol%atoms(i)%cords(3) - mol%atoms(j)%cords(3))**2)
    enddo
enddo


! Checking if atoms are bonding
k = 1
do i = 1, size(mol%atoms) - 1
    do j = i + 1, size(mol%atoms)
        if (mol%atoms(i)%element == 'C' .and. mol%atoms(j)%element == 'C' .and. &
        abs(ideal_CC_length - mol%distance(i,j)) < 0.1) then
            mol%bonds(k)%link(1) = i
            mol%bonds(k)%link(2) = j
            mol%bonds(k)%length = mol%distance(j,i)
            mol%bonds(k)%vector = mol%atoms(i)%cords - mol%atoms(j)%cords
            mol%bonding(j,i) = .true.
            mol%bonding(i,j) = .true.
            mol%bonds(k)%type = 'CC'
            k = k + 1
        else if (mol%atoms(i)%element == 'C' .and. mol%atoms(j)%element == 'H' .and. &
        abs(ideal_CH_length - mol%distance(i,j)) < 0.1) then
            mol%bonds(k)%link(1) = i
            mol%bonds(k)%link(2) = j
            mol%bonds(k)%length = mol%distance(j,i)
            mol%bonds(k)%vector = mol%atoms(i)%cords - mol%atoms(j)%cords
            mol%bonding(j,i) = .true.
            mol%bonding(i,j) = .true.
            mol%bonds(k)%type = 'CH'
            k = k + 1
        end if
    enddo
enddo

end subroutine bonds_atom

! Calculating the angels between neighboring bonds (in this project must be an carbon atom) 
subroutine angle_bonds(mol)
type (molecule) :: mol
type (bond) :: bonds_holder(4,count(mol%atoms%element == 'C'))
integer :: carbon_indicies(count(mol%atoms%element == 'C'))
integer :: i, j, k, l

allocate(mol%angles(count(mol%atoms%element == 'C') * 6))

! Getting the atom indicies of the carbon atoms (is the center atom to which the two bonds tha make an angle are conected)
k = 1
do j =1,size(mol%atoms)
    if (mol%atoms(j)%element == 'C') then
        carbon_indicies(k) = j
        k = k + 1
    end if
enddo

! Puts for every carbon atom the 4 bonds that are conected to that atom in an bonds holder array
do i = 1, size(carbon_indicies)
    k = 1
    do j = 1,size(mol%bonds)
        if (mol%bonds(j)%link(1) == carbon_indicies(i) .or. mol%bonds(j)%link(2) == carbon_indicies(i)) then
            bonds_holder(k,i) = mol%bonds(j)
            k = k + 1
        end if
    enddo
end do

! Creates (6) unique paires of 2 bonds for each carbon atom (4 bonds) and calculates the angle between them using vectors
k = 1
do l = 1,size(carbon_indicies)
    do i = 1,3
        do j = i + 1,4
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

! Calculating the torsion angles. In this project the center bond must be an CC
subroutine angle_torsion(mol)
type (molecule), intent(inout) :: mol
type (bond), allocatable :: end_bonds_holder(:,:,:), CC_bonds_holder(:)
real(realkind) :: n1(3), n2(3)
integer :: i, j, k, l

! If there are 0 bonds over which an torsion angle can be calcualted the subroutine is skipped
if (count(mol%bonds%type == 'CC') > 0) then

allocate(mol%torsion_angles(count(mol%bonds%type == 'CC') * 9))
!! look intyo allocating at beggining 
allocate(end_bonds_holder(3, 2, count(mol%bonds%type == 'CC')))
allocate(CC_bonds_holder(count(mol%bonds%type == 'CC')))

! Extracting all the CC (center) bonds
k = 1
do i = 1, size(mol%bonds)
        if (mol%bonds(i)%type == 'CC') then
            CC_bonds_holder(k) = mol%bonds(i)
            k = k + 1
        end if
enddo

! Getting all the the (2 x 3) (end) bonds of the CC bond per C atom
do i = 1, size(CC_bonds_holder)
    k = 1
    l = 1
    do j = 1,size(mol%bonds)
        !! can be optimised
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

! Makking (9) unique end bonds pairs
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

! Calcualting the angle according to the newman projection 
do i = 1, size(mol%torsion_angles)
    n1 = cross_product(mol%torsion_angles(i)%end_bonds(1)%vector, mol%torsion_angles(i)%center_bond%vector)
    n2 = cross_product(mol%torsion_angles(i)%end_bonds(2)%vector, mol%torsion_angles(i)%center_bond%vector)

    mol%torsion_angles(i)%angle = rad_to_deg(acos(dot_product(n1, n2) / (vector_lenght(n1) * vector_lenght(n2))))
enddo

end if
end subroutine angle_torsion

end module