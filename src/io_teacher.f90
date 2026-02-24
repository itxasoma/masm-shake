MODULE io_teacher
  use GLOBAL
  IMPLICIT NONE
contains

  SUBROUTINE read_input_dades(filename, nsteps, dt_ps, tauT_ps, nmolecules, Tref_K, sigma_A, epsil_K, mass_gmol, r0_A)
    IMPLICIT NONE
    character(len=*), intent(in) :: filename
    integer, intent(out) :: nsteps, nmolecules
    double precision, intent(out) :: dt_ps, tauT_ps, Tref_K, sigma_A, epsil_K, mass_gmol, r0_A
    integer :: u

    open(newunit=u, file=filename, status='old', action='read')
      read(u,*) nsteps
      read(u,*) dt_ps, tauT_ps
      read(u,*) nmolecules, Tref_K
      read(u,*) sigma_A, epsil_K
      read(u,*) mass_gmol
      read(u,*) r0_A
    close(u)
  END SUBROUTINE

  SUBROUTINE read_conf_data(filename, nmolecules, posA, velAps, boxA)
    IMPLICIT NONE
    character(len=*), intent(in) :: filename
    integer, intent(in) :: nmolecules
    double precision, allocatable, intent(out) :: posA(:,:), velAps(:,:)
    double precision, intent(out) :: boxA
    integer :: u, ic, is, i
    integer, parameter :: natoms = 3

    if (allocated(posA)) deallocate(posA)
    if (allocated(velAps)) deallocate(velAps)

    allocate(posA(3*nmolecules,3))
    allocate(velAps(3*nmolecules,3))

    open(newunit=u, file=filename, status='old', action='read')
      do ic = 1, nmolecules
        do is = 1, natoms
          i = (ic-1)*natoms + is
          read(u,*) posA(i,1), posA(i,2), posA(i,3)
          read(u,*) velAps(i,1), velAps(i,2), velAps(i,3)
        end do
      end do
      read(u,*) boxA
    close(u)
  END SUBROUTINE

  SUBROUTINE to_reduced_units(nmolecules, sigma_A, epsil_K, mass_gmol, dt_ps, tauT_ps, Tref_K, posA, velAps, boxA, &
                             dt, tauT, Tref, pos, vel, box_red, uvel_Aps, utime_ps)
    IMPLICIT NONE
    integer, intent(in) :: nmolecules
    double precision, intent(in) :: sigma_A, epsil_K, mass_gmol
    double precision, intent(in) :: dt_ps, tauT_ps, Tref_K
    double precision, intent(in) :: posA(:,:), velAps(:,:), boxA
    double precision, intent(out) :: dt, tauT, Tref, box_red, uvel_Aps, utime_ps
    double precision, allocatable, intent(out) :: pos(:,:), vel(:,:)

    double precision :: rgas
    integer :: i

    rgas = 8.314472673d0  ! J/(mol*K)
    utime_ps = sigma_A * sqrt(mass_gmol/epsil_K) * sqrt(10.d0/rgas)
    uvel_Aps = sigma_A / utime_ps

    box_red = boxA / sigma_A
    dt   = dt_ps   / utime_ps
    tauT = tauT_ps / utime_ps
    Tref = Tref_K  / epsil_K

    if (allocated(pos)) deallocate(pos)
    if (allocated(vel)) deallocate(vel)
    allocate(pos(size(posA,1),3))
    allocate(vel(size(velAps,1),3))

    do i = 1, size(posA,1)
      pos(i,:) = posA(i,:) / sigma_A
      vel(i,:) = velAps(i,:) / uvel_Aps
    end do
  END SUBROUTINE

END MODULE io_teacher

