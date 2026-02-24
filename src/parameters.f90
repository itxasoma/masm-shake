MODULE GLOBAL
  IMPLICIT NONE
  ! System size
  integer:: N = 0, M = 0
  integer:: nmol = 0
  integer:: atoms_per_mol = 0

  ! Parameters (reduced units)
  double precision:: rho   = 0.d0
  double precision:: cutoff = 2.5d0
  double precision:: L     = 0.d0
  double precision:: a     = 0.d0

  ! Some cutoff radii (from MoMo)
  double precision:: rc(4) = [1.5d0, 2.d0, 2.5d0, 3.d0]

END MODULE GLOBAL
