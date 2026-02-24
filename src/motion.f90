MODULE motion
  use GLOBAL
  use force
  IMPLICIT NONE
contains

  SUBROUTINE kinetic_energy(vel, kin)
    IMPLICIT NONE
    double precision, intent(in):: vel(:,:)
    double precision, intent(out):: kin
    integer:: i
    kin = 0.d0
    do i = 1, N
      kin = kin + 0.5d0*(vel(i,1)**2 + vel(i,2)**2 + vel(i,3)**2)
    enddo
  END SUBROUTINE

  SUBROUTINE temperature_from_kin(kin, nf, temp)
    IMPLICIT NONE
    double precision, intent(in):: kin
    integer, intent(in):: nf
    double precision, intent(out):: temp
    temp = 2.d0*kin / dble(nf)
  END SUBROUTINE

  SUBROUTINE apply_berendsen(dt, tauT, Tref, nf, vel, lambda)
    IMPLICIT NONE
    double precision, intent(in):: dt, tauT, Tref
    integer, intent(in):: nf
    double precision, intent(inout):: vel(:,:)
    double precision, intent(out):: lambda
    double precision:: kin, temp

    call kinetic_energy(vel, kin)
    call temperature_from_kin(kin, nf, temp)

    ! Avoid division by zero if something exploded
    if (temp <= 1.d-14) then
      lambda = 1.d0
      return
    endif

    lambda = sqrt(1.d0 + (dt/tauT) * (Tref/temp - 1.d0))
    vel(:,:) = vel(:,:) * lambda
  END SUBROUTINE

  SUBROUTINE apply_pbc_all(pos)
    IMPLICIT NONE
    double precision, intent(inout):: pos(:,:)
    integer :: i
    do i = 1, N
      pos(i,:) = modulo(pos(i,:), L)
    enddo
  END SUBROUTINE

  SUBROUTINE time_step_VelocityVerlet_NVT(dt, cutoff, tauT, Tref, nf, pos, vel, Upot, kin, temp, lambda)
    IMPLICIT NONE
    double precision, intent(in)    :: dt, cutoff, tauT, Tref
    integer, intent(in)             :: nf
    double precision, intent(inout) :: pos(:,:), vel(:,:)
    double precision, intent(out)   :: Upot, kin, temp, lambda
    double precision, allocatable   :: F(:,:), Fnew(:,:)

    call find_force_lj(cutoff, pos, Upot, F)

    vel(:,:) = vel(:,:) + 0.5d0*F(:,:)*dt
    pos(:,:) = pos(:,:) + vel(:,:)*dt
    call apply_pbc_all(pos)

    call find_force_lj(cutoff, pos, Upot, Fnew)
    vel(:,:) = vel(:,:) + 0.5d0*Fnew(:,:)*dt

    ! Thermostat (Berendsen)
    call apply_berendsen(dt, tauT, Tref, nf, vel, lambda)

    call kinetic_energy(vel, kin)
    call temperature_from_kin(kin, nf, temp)

    if (allocated(F))    deallocate(F)
    if (allocated(Fnew)) deallocate(Fnew)
  END SUBROUTINE

  SUBROUTINE write_xyz(unit, pos, names)
    IMPLICIT NONE
    integer, intent(in)         :: unit
    double precision, intent(in):: pos(:,:)
    character(len=*), intent(in):: names(:)
    integer:: i
    write(unit,*) N
    write(unit,*)  ! blank comment line
    do i = 1, N
      write(unit,'(A,3(1X,F20.10))') trim(names(i)), pos(i,1), pos(i,2), pos(i,3)
    enddo
  END SUBROUTINE

END MODULE
