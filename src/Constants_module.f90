!
! Module containing all constants (parameters)
! Author: Martijn Oele (GitHub: Martijn-075)
!
module constant_module
implicit none

public

! The bits for the real numbers (8 = double precision)
integer, parameter :: realkind = 8

! Molecule_module
real(realkind), parameter :: ideal_CC_length = 1.535
real(realkind), parameter :: ideal_CH_length = 1.094

! Energy_module
! Stretching parameters
real(realkind), parameter :: k_CC = 310.
real(realkind), parameter :: k_CH = 340.
real(realkind), parameter :: R_CC = 1.526
real(realkind), parameter :: R_CH = 1.090
! Bending parameters
real(realkind), parameter :: k_CC_CC = 40.0
real(realkind), parameter :: k_CC_CH = 50.0
real(realkind), parameter :: k_CH_CH = 35.0
real(realkind), parameter :: theta_0_SP3 = 109.50
! Electrostatic parameters
real(realkind), parameter :: q_H = 0.078
real(realkind), parameter :: q_C = -0.344
real(realkind), parameter :: coulombe_constant = 331.1908 !! Needs to be checked 2.146639897e16, eV*A/e2 converted to kcal/mol
! Van der waals parameters
real(realkind), parameter :: A_vdw_HH = 0.0157 * (2 * 1.4870)**12
real(realkind), parameter :: B_vdw_HH = 0.0157 * (2 * 1.4870)**6
real(realkind), parameter :: A_vdw_CH = sqrt(0.0157 * 0.1094) * (1.4870 + 1.9080)**12
real(realkind), parameter :: B_vdw_CH = sqrt(0.0157 * 0.1094) * (1.4870 + 1.9080)**6
real(realkind), parameter :: A_vdw_CC = 0.1094 * (2 * 1.9080)**12
real(realkind), parameter :: B_vdw_CC = 0.1094 * (2 * 1.9080)**6
! Torsion angle
real(realkind), parameter :: n_CC = 3.0
real(realkind), parameter :: V2_CC = 1.40
real(realkind), parameter :: gamma_CC = 0.

end module