MODULE potentials
  use GLOBAL
  IMPLICIT NONE
contains
! Periodic boundary conditions for a cutoff < L/2
    SUBROUTINE pbc1(x, L)
        double precision, intent(in):: L
        double precision:: x
        if (x.ge.(L/2.d0)) then
            x = x - L
        elseif (x.lt.(-L/2.d0)) then
            x = x + L
        endif
    END SUBROUTINE

! eps = sigma = 1 (reduced units)
    SUBROUTINE lj_energy(cutoff, pos, Upot)
      double precision, intent(in):: cutoff
      double precision, allocatable, intent(in):: pos(:,:)
      double precision, intent(out):: Upot
      double precision:: cutoff2, cf4, cf6, cf12, shift
      double precision:: dx, dy, dz, d2, d4, d6, d12
      integer:: i, j
      cutoff2 = cutoff*cutoff
      cf4 = cutoff2*cutoff2
      cf6 = cutoff2*cf4
      cf12 = cf6*cf6
      shift = 4.0d0*(1.0d0/cf12 - 1.0d0/cf6)
      Upot = 0.0d0
      do i = 1, N
        do j = i+1, N
          dx = pos(i,1)-pos(j,1); dy = pos(i,2)-pos(j,2); dz = pos(i,3)-pos(j,3)
          d2 = dx*dx + dy*dy + dz*dz
          if ((d2).lt.cutoff2) then
            d4=d2*d2; d6=d4*d2; d12=d6*d6
            Upot = Upot + 4.0d0*(1.0d0/d12 - 1.0d0/d6) - shift
          endif
        enddo
      enddo
    END SUBROUTINE

! with pbc:
    SUBROUTINE lj_energy_pbc(cutoff, pos, Upot)
      double precision, intent(in):: cutoff
      double precision, allocatable, intent(in):: pos(:,:)
      double precision, intent(out):: Upot
      double precision:: cutoff2, cf4, cf6, cf12, shift
      double precision:: dx, dy, dz, d2, d4, d6, d12
      integer:: i, j
      cutoff2 = cutoff*cutoff
      cf4 = cutoff2*cutoff2
      cf6 = cutoff2*cf4
      cf12 = cf6*cf6
      shift = 4.0d0*(1.0d0/cf12 - 1.0d0/cf6)
      Upot = 0.0d0
      do i = 1, N
        do j = i+1, N
          dx = pos(i,1)-pos(j,1); dy = pos(i,2)-pos(j,2); dz = pos(i,3)-pos(j,3)
          call pbc1(dx, L); call pbc1(dy, L); call pbc1(dz, L)
          d2 = dx*dx + dy*dy + dz*dz
          if ((d2).lt.cutoff2) then
            d4=d2*d2; d6=d4*d2; d12=d6*d6
            Upot = Upot + 4.0d0*(1.0d0/d12 - 1.0d0/d6) - shift
          endif
        enddo
      enddo
    END SUBROUTINE

END MODULE
