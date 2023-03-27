!
! constants module
! Module containing all constants (parameters)
! Author: Martijn Oele (GitHub: Martijn-075)
!
module constant_module
implicit none

public

! The bytes for the real numbers (8 = double precision)
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
real(realkind), parameter :: coulombe_constant = 0.6246348496 !! converted from n m**2 / C**2 to A Kcal / mol C**2
! Van der waals parameters
real(realkind), parameter :: A_vdw_HH = 0.0157 * (2 * 1.4870)**12
real(realkind), parameter :: B_vdw_HH = 2. * 0.0157 * (2 * 1.4870)**6
real(realkind), parameter :: A_vdw_CH = sqrt(0.0157 * 0.1094) * (1.4870 + 1.9080)**12
real(realkind), parameter :: B_vdw_CH = 2. * sqrt(0.0157 * 0.1094) * (1.4870 + 1.9080)**6
real(realkind), parameter :: A_vdw_CC = 0.1094 * (2 * 1.9080)**12
real(realkind), parameter :: B_vdw_CC = 2. * 0.1094 * (2 * 1.9080)**6
! Torsion angle
real(realkind), parameter :: n_CC = 3.0
real(realkind), parameter :: V2_CC = 1.40
real(realkind), parameter :: gamma_CC = 0.

! Metropolis module
real(realkind), parameter :: kb = 1.987204259e-3

! Math module
real(realkind), parameter :: pi = 4.*atan(1.)

end module