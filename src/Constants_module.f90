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
real(realkind), parameter :: ideal_CC_length = 1.535 ! A
real(realkind), parameter :: ideal_CH_length = 1.094 ! A

! Energy_module
! Stretching parameters
real(realkind), parameter :: k_CC = 310. ! kcal / mol A**2
real(realkind), parameter :: k_CH = 340. ! kcal / mol A**2
real(realkind), parameter :: R_CC = 1.526 ! A
real(realkind), parameter :: R_CH = 1.090 ! A
! Bending parameters
real(realkind), parameter :: k_CC_CC = 40.0 ! kcal / mol rad**2
real(realkind), parameter :: k_CC_CH = 50.0 ! kcal / mol rad**2
real(realkind), parameter :: k_CH_CH = 35.0 ! kcal / mol rad**2
real(realkind), parameter :: theta_0_SP3 = 109.50 ! deg
! Electrostatic parameters
real(realkind), parameter :: q_H = 0.078 ! C
real(realkind), parameter :: q_C = -0.344 ! C
real(realkind), parameter :: coulombe_constant = 0.6246348496 ! A Kcal / mol C**2
! Van der waals parameters
real(realkind), parameter :: A_vdw_HH = 0.0157 * (2 * 1.4870)**12
real(realkind), parameter :: B_vdw_HH = 2. * 0.0157 * (2 * 1.4870)**6
real(realkind), parameter :: A_vdw_CH = sqrt(0.0157 * 0.1094) * (1.4870 + 1.9080)**12
real(realkind), parameter :: B_vdw_CH = 2. * sqrt(0.0157 * 0.1094) * (1.4870 + 1.9080)**6
real(realkind), parameter :: A_vdw_CC = 0.1094 * (2 * 1.9080)**12
real(realkind), parameter :: B_vdw_CC = 2. * 0.1094 * (2 * 1.9080)**6
! Torsion angle
real(realkind), parameter :: n_CC = 3.0
real(realkind), parameter :: V2_CC = 1.40 ! kcal / mol
real(realkind), parameter :: gamma_CC = 0. ! deg

! Metropolis module
real(realkind), parameter :: kb = 1.987204259e-3 ! kcal / mol K

! Math module
real(realkind), parameter :: pi = 4.*atan(1.)

end module