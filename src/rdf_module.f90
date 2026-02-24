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
    delg = box / (2.d0*dble(nhis)) ! bin width, so that the last bin ends at L/2
    ngr  = 0

    if (allocated(hist)) deallocate(hist)
    allocate(hist(nhis))
    hist(:) = 0.d0
  END SUBROUTINE

  SUBROUTINE rdf_sample(pos)
    IMPLICIT NONE
    double precision, intent(in) :: pos(:,:)
    integer :: i, j, ig
    double precision :: rij(3), r

    ngr = ngr + 1

    do i = 1, N-1
      do j = i+1, N
        rij(:) = pos(i,:) - pos(j,:)
        call pbc(rij, L)
        r = sqrt(rij(1)*rij(1) + rij(2)*rij(2) + rij(3)*rij(3))

        if (r < L/2.d0) then
          ig = int(r/delg) + 1
          if (ig >= 1 .and. ig <= nhis) hist(ig) = hist(ig) + 2.d0
        end if
      end do
    end do
  END SUBROUTINE

  SUBROUTINE rdf_write(filename, rho_in)
    IMPLICIT NONE
    character(len=*), intent(in) :: filename
    double precision, intent(in) :: rho_in
    integer :: u, i
    double precision :: pi, r_in, r_out, r_center, shell_vol, nid, g

    pi = 4.d0*atan(1.d0)

    if (ngr <= 0) then
      write(*,*) 'RDF: no samples collected.'
      return
    end if

    open(newunit=u, file=filename, status='replace', action='write')
      write(u,'(A)') '# r   g(r)'
      do i = 1, nhis
        r_in     = dble(i-1)*delg
        r_out    = dble(i)*delg
        r_center = (r_in + r_out)/2.d0

        shell_vol = (4.d0*pi/3.d0) * (r_out**3 - r_in**3)
        nid = shell_vol * rho_in

        if (nid > 1.d-14) then
          g = hist(i) / (dble(ngr) * dble(N) * nid)
        else
          g = 0.d0
        end if

        write(u,'(2(1X,F20.10))') r_center, g
      end do
    close(u)
  END SUBROUTINE

END MODULE rdf_module

