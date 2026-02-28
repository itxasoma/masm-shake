MODULE force
  use GLOBAL
  IMPLICIT NONE
contains

  SUBROUTINE pbc(rij, L)
    IMPLICIT NONE
    double precision, intent(in)    :: L
    double precision, intent(inout) :: rij(3)
    integer:: i
    do i = 1, 3
      if (rij(i) .ge. (L/2.d0)) then
        rij(i) = rij(i) - L
      elseif (rij(i) .lt. (-L/2.d0)) then
        rij(i) = rij(i) + L
      end if
    end do
  END SUBROUTINE

  SUBROUTINE find_force_lj(cutoff, pos, Upot, F)
    double precision, intent(in):: cutoff
    double precision, intent(in):: pos(:,:)
    double precision, allocatable, intent(out):: F(:,:)
    double precision, intent(out):: Upot
    double precision:: cutoff2, cf4, cf6, cf12, shift
    double precision:: rij(3), d2, d4, d6, d8, d12, d14
    double precision:: pref
    integer:: i, j, mol_i, mol_j

    cutoff2 = cutoff*cutoff
    cf4 = cutoff2*cutoff2
    cf6 = cutoff2*cf4
    cf12 = cf6*cf6
    shift = 4.0d0*(1.0d0/cf12 - 1.0d0/cf6)
    Upot = 0.0d0

    if (allocated(F)) deallocate(F)
    allocate(F(N,3))
    F(:,:) = 0.0d0

    do i = 1, N
      do j = i+1, N

        ! Skip intramolecular interactions
        if (atoms_per_mol > 0) then
          mol_i = (i-1)/atoms_per_mol
          mol_j = (j-1)/atoms_per_mol
          if (mol_i == mol_j) cycle
        endif

        rij = pos(i,:) - pos(j,:)
        call pbc(rij, L)
        d2 = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)

        if (d2 .lt. cutoff2) then
          d4  = d2*d2
          d6  = d4*d2
          d8  = d4*d4
          d12 = d6*d6
          d14 = d12*d2

          Upot = Upot + 4.0d0*(1.0d0/d12 - 1.0d0/d6) - shift

          pref = (48.d0/d14 - 24.d0/d8)
          F(i,:) = F(i,:) + pref * rij(:)
          F(j,:) = F(j,:) - pref * rij(:)
        endif
      enddo
    enddo
  END SUBROUTINE

END MODULE
