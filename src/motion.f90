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
    double precision, intent(in):: dt, cutoff, tauT, Tref
    integer, intent(in):: nf
    double precision, intent(inout):: pos(:,:), vel(:,:)
    double precision, intent(out):: Upot, kin, temp, lambda
    double precision, allocatable:: F(:,:), Fnew(:,:)

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

! HERE, WE ADD SHAKE: CONSTRAINTS ON BONDS
! ---------------------------------------------------------------
! SHAKE: corregeix posicions provisionals per mantenir el triangle
! equilàter de cada molècula (r_ij = d per als 3 parells).
!
! pos_old : posicions a temps t (no es modifiquen)
! pos_pro : posicions provisionals t+dt (es corregeixen in-place)
! dt : pas de temps (u.r.)
! d : costat del triangle (u.r.)
!
! Referència: Ryckaert, Ciccotti & Berendsen, J.Comput.Phys 23, 327 (1977)
! ---------------------------------------------------------------
SUBROUTINE shake(pos_old, pos_pro, dt, d)
  IMPLICIT NONE
  double precision, intent(in):: pos_old(:,:)
  double precision, intent(inout):: pos_pro(:,:)
  double precision, intent(in):: dt, d
  integer, parameter:: max_iter = 1000
  double precision, parameter:: tol = 1.d-8
  integer:: ic, k, iter
  integer:: a, b
  integer:: i1, i2, i3
  integer:: pairs(3,2)
  double precision:: r_ab(3), rp_ab(3)
  double precision:: rp2, dot_rp_r, mu, dsq
  logical:: converged
  
  dsq = d * d

    do ic = 1, nmol
      ! Índexos dels 3 àtoms de la molècula
      i1 = (ic-1)*atoms_per_mol + 1
      i2 = (ic-1)*atoms_per_mol + 2
      i3 = (ic-1)*atoms_per_mol + 3

      ! Els 3 parells del triangle
      pairs(1,1) = i1 ; pairs(1,2) = i2 ! parell 1-2
      pairs(2,1) = i1 ; pairs(2,2) = i3 ! parell 1-3
      pairs(3,1) = i2 ; pairs(3,2) = i3 ! parell 2-3
    
      do iter = 1, max_iter
        converged = .true.
        do k = 1, 3
          a = pairs(k,1)
          b = pairs(k,2)
          ! Vector provisional r_pij(t+dt) = r_pb - r_pa 
          rp_ab(:) = pos_pro(b,:) - pos_pro(a,:)
          call pbc(rp_ab, L)
          rp2 = rp_ab(1)**2 + rp_ab(2)**2 + rp_ab(3)**2

          if (abs(rp2 - dsq) > tol) then
            converged = .false.

            ! Vector antic r_ij(t) = r_b(t) - r_a(t) (mínim imatge, fix durant iter)
            r_ab(:) = pos_old(b,:) - pos_old(a,:)
            call pbc(r_ab, L)
            dot_rp_r = rp_ab(1)*r_ab(1) + rp_ab(2)*r_ab(2) + rp_ab(3)*r_ab(3)

            ! Multiplicador de Lagrange (m=1 en u.r., 1/ma+1/mb = 2)
            ! mu = (rp_ij^2 - d^2)/(4*(1/ma+1/mb)*dt^2 * rp_ij·r_ij(t))
            mu = (rp2 - dsq)/(4.d0 * 2.d0 * dt*dt * dot_rp_r)

            ! Correcció de posicions:
            ! r_a = r_pa + 2*mu*dt^2/ma * r_ab(t)
            ! r_b = r_pb - 2*mu*dt^2/mb * r_ab(t)
            pos_pro(a,:) = pos_pro(a,:) + 2.d0 * mu * dt*dt * r_ab(:)
            pos_pro(b,:) = pos_pro(b,:) - 2.d0 * mu * dt*dt * r_ab(:)
          endif

        enddo ! k (parells)
        if (converged) exit
      enddo ! iter
      if (.not. converged) &
      write(*,'(A,I6,A)') 'WARNING: SHAKE did not converge for molecule ', ic, '!'
    enddo ! ic (molècules)
END SUBROUTINE shake


! ---------------------------------------------------------------
! Velocity-Verlet NVT with Berendsen + SHAKE
! If d <= 0.d0, no SHAKE (d=0 means no constraint))
! ---------------------------------------------------------------
SUBROUTINE time_step_VelocityVerlet_NVT_shake(dt, cutoff, tauT, Tref, nf, d, pos, vel, Upot, kin, temp, lambda)
IMPLICIT NONE
double precision, intent(in):: dt, cutoff, tauT, Tref, d
integer, intent(in):: nf
double precision, intent(inout):: pos(:,:), vel(:,:)
double precision, intent(out):: Upot, kin, temp, lambda
double precision, allocatable:: F(:,:), Fnew(:,:)
double precision, allocatable:: pos_old(:,:)

  ! Save positions at time t (SHAKE)
  allocate(pos_old(N,3))
  pos_old(:,:) = pos(:,:)
  
  call find_force_lj(cutoff, pos, Upot, F)
  ! Half time step: v(t+dt/2) = v(t) + F(t)*dt/2
  vel(:,:) = vel(:,:) + 0.5d0*F(:,:)*dt

  ! Provisional positions: r_pro = r(t) + v(t+dt/2)*dt
  pos(:,:) = pos(:,:) + vel(:,:)*dt
  call apply_pbc_all(pos)

  ! SHKE: correct positions
  if (d > 0.d0) call shake(pos_old, pos, dt, d)

  ! New forces with corrected positions
  call find_force_lj(cutoff, pos, Upot, Fnew)

  ! Second half step: v(t+dt) = v(t+dt/2) + F(t+dt)*dt/2
  vel(:,:) = vel(:,:) + 0.5d0*Fnew(:,:)*dt

  ! Berendsen 
  call apply_berendsen(dt, tauT, Tref, nf, vel, lambda)
  call kinetic_energy(vel, kin)
  call temperature_from_kin(kin, nf, temp)
  
  if (allocated(F)) deallocate(F)
  if (allocated(Fnew)) deallocate(Fnew)
  if (allocated(pos_old)) deallocate(pos_old)
END SUBROUTINE time_step_VelocityVerlet_NVT_shake


END MODULE
