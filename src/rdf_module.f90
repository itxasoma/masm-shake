MODULE rdf_module
  use GLOBAL
  use force, only: pbc
  IMPLICIT NONE

  integer :: nhis = 0
  integer :: ngr  = 0
  double precision :: delg = 0.d0
  double precision, allocatable :: hist(:)

contains

  SUBROUTINE rdf_init(nhis_in, box)
    IMPLICIT NONE
    integer, intent(in) :: nhis_in
    double precision, intent(in) :: box

    nhis = nhis_in
    delg = (box/2.d0) / dble(nhis)

    if (allocated(hist)) deallocate(hist)
    allocate(hist(nhis))
    call rdf_reset()
  END SUBROUTINE rdf_init


  SUBROUTINE rdf_reset()
    IMPLICIT NONE
    if (allocated(hist)) hist(:) = 0.d0
    ngr = 0
  END SUBROUTINE rdf_reset


  SUBROUTINE rdf_sample(pos)
    IMPLICIT NONE
    double precision, intent(in) :: pos(:,:)

    integer :: ic, is, jc, js
    integer :: i, j, ibin
    double precision :: rij(3), rr2, dist
    double precision :: rmax

    ngr = ngr + 1
    rmax = L/2.d0

    do ic = 1, nmol-1
      do is = 1, atoms_per_mol
        i = (ic-1)*atoms_per_mol + is

        do jc = ic+1, nmol
          do js = 1, atoms_per_mol
            j = (jc-1)*atoms_per_mol + js

            rij(:) = pos(j,:) - pos(i,:)
            call pbc(rij, L)

            rr2 = rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3)
            dist = sqrt(rr2)

            if (dist .lt. rmax) then
              ibin = int(dist/delg) + 1
              if (ibin <= nhis) hist(ibin) = hist(ibin) + 2.d0
            end if

          end do
        end do
      end do
    end do

  END SUBROUTINE rdf_sample


  SUBROUTINE rdf_write(filename)
    IMPLICIT NONE
    character(len=*), intent(in) :: filename

    integer :: u, i
    double precision :: pi, vol, rho_inter, total_atoms
    double precision :: r_bin, del_r, vol_shell, ideal_count, g_val

    pi = 4.d0*atan(1.d0)

    if (ngr <= 0) then
      write(*,*) 'RDF: no samples collected.'
      return
    end if

    del_r = delg
    vol = L**3

    rho_inter   = dble((nmol-1)*atoms_per_mol) / vol
    total_atoms = dble(nmol*atoms_per_mol)

    open(newunit=u, file=filename, status='replace', action='write')
      write(u,'(A)') '# r   g(r)'
      do i = 1, nhis
        r_bin = (dble(i)-0.5d0)*del_r

        vol_shell = (4.d0*pi/3.d0) * ( (r_bin+0.5d0*del_r)**3 - (r_bin-0.5d0*del_r)**3 )
        ideal_count = total_atoms * rho_inter * vol_shell

        if (ideal_count > 1.d-14) then
          g_val = hist(i) / (dble(ngr) * ideal_count)
        else
          g_val = 0.d0
        end if

        write(u,'(2(1X,F20.10))') r_bin, g_val
      end do
    close(u)
  END SUBROUTINE rdf_write

END MODULE rdf_module
