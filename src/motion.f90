MODULE motion
  use GLOBAL
  use force
  IMPLICIT NONE
contains

! Computes the kinetic energy for all three velocity
! components of N particles. vel (vector) --> kin (scalar)
    SUBROUTINE kinetic_energy(vel, kin)
        IMPLICIT NONE
        double precision, intent(in):: vel(:,:)
        double precision, intent(out):: kin
        integer :: i
        kin = 0.d0
        do i = 1, N
          kin = kin + 0.5d0*(vel(i,1)**2 + vel(i,2)**2 + vel(i,3)**2)
        enddo
    END SUBROUTINE

! Applies periodic boundary conditions.
    SUBROUTINE apply_pbc_all(pos)
        IMPLICIT NONE
        double precision, intent(inout):: pos(:,:)
        integer :: i
        do i = 1, N
          pos(i,:) = mod(pos(i,:), L)
        enddo
    END SUBROUTINE

! Euler algorithm
    SUBROUTINE time_step_Euler_pbc(dt, cutoff, pos, vel, Upot, kin)
        IMPLICIT NONE
        double precision, intent(in):: dt, cutoff
        double precision, intent(inout):: pos(:,:), vel(:,:)
        double precision, intent(out):: Upot, kin
        double precision, allocatable:: F(:,:)
        ! Calls the LJ force
        call find_force_lj(cutoff, pos, Upot, F)
        pos(:,:) = pos(:,:) + vel(:,:)*dt + 0.5d0*F(:,:)*dt*dt
        vel(:,:) = vel(:,:) + F(:,:)*dt
        call apply_pbc_all(pos)
        call kinetic_energy(vel, kin)
        if (allocated(F)) deallocate(F)
    END SUBROUTINE

! Verlet algorithm
    SUBROUTINE time_step_Verlet(dt, cutoff, pos, vel, Upot, kin)
        IMPLICIT NONE
        double precision, intent(in):: dt, cutoff
        double precision, intent(inout):: pos(:,:), vel(:,:)
        double precision, intent(out):: Upot, kin
        double precision, allocatable:: F(:,:), Fnew(:,:)
        ! Find the position and apply pbc
        call find_force_lj(cutoff, pos, Upot, F)
        pos(:,:) = pos(:,:) + vel(:,:)*dt + 0.5d0*dt*dt*F(:,:)
        call apply_pbc_all(pos)
        ! Find kin with velocities
        call find_force_lj(cutoff, pos, Upot, Fnew)
        vel(:,:) = vel(:,:) + 0.5d0*(F(:,:) + Fnew(:,:))*dt
        call kinetic_energy(vel, kin)

        if (allocated(F))    deallocate(F)
        if (allocated(Fnew)) deallocate(Fnew)
    END SUBROUTINE

! velocity Verlet algorithm
    SUBROUTINE time_step_VelocityVerlet(dt, cutoff, pos, vel, Upot, kin)
        IMPLICIT NONE
        double precision, intent(in):: dt, cutoff
        double precision, intent(inout):: pos(:,:), vel(:,:)
        double precision, intent(out):: Upot, kin
        double precision, allocatable:: F(:,:), Fnew(:,:)

        call find_force_lj(cutoff, pos, Upot, F)
        vel(:,:) = vel(:,:) + 0.5d0*F(:,:)*dt
        pos(:,:) = pos(:,:) + vel(:,:)*dt
        call apply_pbc_all(pos)

        call find_force_lj(cutoff, pos, Upot, Fnew)
        vel(:,:) = vel(:,:) + 0.5d0*Fnew(:,:)*dt
        call kinetic_energy(vel, kin)

        if (allocated(F))    deallocate(F)
        if (allocated(Fnew)) deallocate(Fnew)
    END SUBROUTINE

! Write down the trajectory coordinates
    SUBROUTINE write_xyz(unit, pos, names)
        IMPLICIT NONE
        integer, intent(in):: unit
        double precision, intent(in):: pos(:,:)
        character(len=*), intent(in):: names(:)
        integer:: i
        write(unit,*) N
        write(unit,*)
        do i = 1, N
          write(unit,'(A,3(1X,F20.10))') trim(names(i)), pos(i,1), pos(i,2), pos(i,3)
        enddo
    END SUBROUTINE

! Distance between two particles when applying pbc
    FUNCTION distance12(pos) result(d)
        IMPLICIT NONE
        double precision, intent(in):: pos(:,:)
        double precision:: d, rij(3)
        rij(:) = pos(2,:) - pos(1,:)
        call pbc(rij, L)
        d = sqrt(rij(1)**2 + rij(2)**2 + rij(3)**2)
    END FUNCTION

END MODULE