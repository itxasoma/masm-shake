PROGRAM main
  use GLOBAL
  use crystal
  use motion
  use force
  IMPLICIT NONE

  integer :: step, nsteps
  double precision :: dt, t
  double precision, allocatable :: pos(:,:), vel(:,:)
  double precision :: Upot, kin
  character(len=8), allocatable :: names(:)

  ! 4. System parameters
  L = 10.d0
  cutoff = 2.5d0  ! < L/2

  ! Initialize two particles
  call two_particles(1.01d0, pos)
  allocate(vel(N,3)); vel(:,:) = 0.d0
  allocate(names(N)); names = 'A'

! -------------------------------
! Single run: Velocity-Verlet
! -------------------------------
  dt = 1.0d-3
  nsteps = 10000
  t = 0.d0

  open(unit=10, file='trajectory_velVerlet.xyz',         status='replace', action='write')
  open(unit=20, file='energies_velVerlet.dat',           status='replace', action='write')
  write(*,*) '--- Running Velocity-Verlet (single run) ---'
  do step = 0, nsteps
     if (mod(step,10) == 0) call write_xyz(10, pos, names)
     write(20,'(F12.6,3(1X,F20.10))') t, Upot, kin, Upot+kin
     call time_step_VelocityVerlet(dt, cutoff, pos, vel, Upot, kin)
     t = t + dt
  enddo
  close(20)
  close(10)

! -------------------------------
! Single run: Euler
! -------------------------------
  call two_particles(1.01d0, pos)
  vel(:,:) = 0.d0
  t = 0.d0

  open(unit=30, file='trajectory_Euler.xyz',             status='replace', action='write')
  open(unit=40, file='energies_euler.dat',               status='replace', action='write')
  write(*,*) '--- Running Euler (single run) ---'
  do step = 0, nsteps
     if (mod(step,10) == 0) call write_xyz(30, pos, names)
     write(40,'(F12.6,3(1X,F20.10))') t, Upot, kin, Upot+kin
     call time_step_Euler_pbc(dt, cutoff, pos, vel, Upot, kin)
     t = t + dt
  enddo
  close(30)
  close(40)

! -------------------------------
! Hand-written dt sweep: Velocity-Verlet
! -------------------------------
  write(*,*) '--- Sweep dt: Velocity-Verlet ---'

  ! dt = 1E-4
  call two_particles(1.01d0, pos); vel(:,:) = 0.d0; t = 0.d0; dt = 1.0d-4
  open(unit=41, file='energies_vV_dt_-1E-04.dat',        status='replace', action='write')
  open(unit=45, file='trajectory_vV_dt_-1E-04.xyz',      status='replace', action='write')
  do step = 0, nsteps
     if (mod(step,10) == 0) call write_xyz(45, pos, names)
     write(41,'(F12.6,3(1X,F20.10))') t, Upot, kin, Upot+kin
     call time_step_VelocityVerlet(dt, cutoff, pos, vel, Upot, kin)
     t = t + dt
  enddo
  close(41); close(45)

  ! dt = 1E-3
  call two_particles(1.01d0, pos); vel(:,:) = 0.d0; t = 0.d0; dt = 1.0d-3
  open(unit=42, file='energies_vV_dt_-1E-03.dat',        status='replace', action='write')
  open(unit=46, file='trajectory_vV_dt_-1E-03.xyz',      status='replace', action='write')
  do step = 0, nsteps
     if (mod(step,10) == 0) call write_xyz(46, pos, names)
     write(42,'(F12.6,3(1X,F20.10))') t, Upot, kin, Upot+kin
     call time_step_VelocityVerlet(dt, cutoff, pos, vel, Upot, kin)
     t = t + dt
  enddo
  close(42); close(46)

  ! dt = 1E-2
  call two_particles(1.01d0, pos); vel(:,:) = 0.d0; t = 0.d0; dt = 1.0d-2
  open(unit=43, file='energies_vV_dt_-1E-02.dat',        status='replace', action='write')
  open(unit=47, file='trajectory_vV_dt_-1E-02.xyz',      status='replace', action='write')
  do step = 0, nsteps
     if (mod(step,10) == 0) call write_xyz(47, pos, names)
     write(43,'(F12.6,3(1X,F20.10))') t, Upot, kin, Upot+kin
     call time_step_VelocityVerlet(dt, cutoff, pos, vel, Upot, kin)
     t = t + dt
  enddo
  close(43); close(47)

  ! dt = 1E-1
  call two_particles(1.01d0, pos); vel(:,:) = 0.d0; t = 0.d0; dt = 1.0d-1
  open(unit=44, file='energies_vV_dt_-1E-01.dat',        status='replace', action='write')
  open(unit=48, file='trajectory_vV_dt_-1E-01.xyz',      status='replace', action='write')
  do step = 0, nsteps
     if (mod(step,10) == 0) call write_xyz(48, pos, names)
     write(44,'(F12.6,3(1X,F20.10))') t, Upot, kin, Upot+kin
     call time_step_VelocityVerlet(dt, cutoff, pos, vel, Upot, kin)
     t = t + dt
  enddo
  close(44); close(48)

! -------------------------------
! Hand-written dt sweep: Euler
! -------------------------------
  write(*,*) '--- Sweep dt: Euler ---'

  ! dt = 1E-4
  call two_particles(1.01d0, pos); vel(:,:) = 0.d0; t = 0.d0; dt = 1.0d-4
  open(unit=51, file='energies_E_dt_-1E-04.dat',         status='replace', action='write')
  open(unit=55, file='trajectory_E_dt_-1E-04.xyz',       status='replace', action='write')
  do step = 0, nsteps
     if (mod(step,10) == 0) call write_xyz(55, pos, names)
     write(51,'(F12.6,3(1X,F20.10))') t, Upot, kin, Upot+kin
     call time_step_Euler_pbc(dt, cutoff, pos, vel, Upot, kin)
     t = t + dt
  enddo
  close(51); close(55)

  ! dt = 1E-3
  call two_particles(1.01d0, pos); vel(:,:) = 0.d0; t = 0.d0; dt = 1.0d-3
  open(unit=52, file='energies_E_dt_-1E-03.dat',         status='replace', action='write')
  open(unit=56, file='trajectory_E_dt_-1E-03.xyz',       status='replace', action='write')
  do step = 0, nsteps
     if (mod(step,10) == 0) call write_xyz(56, pos, names)
     write(52,'(F12.6,3(1X,F20.10))') t, Upot, kin, Upot+kin
     call time_step_Euler_pbc(dt, cutoff, pos, vel, Upot, kin)
     t = t + dt
  enddo
  close(52); close(56)

  ! dt = 1E-2
  call two_particles(1.01d0, pos); vel(:,:) = 0.d0; t = 0.d0; dt = 1.0d-2
  open(unit=53, file='energies_E_dt_-1E-02.dat',         status='replace', action='write')
  open(unit=57, file='trajectory_E_dt_-1E-02.xyz',       status='replace', action='write')
  do step = 0, nsteps
     if (mod(step,10) == 0) call write_xyz(57, pos, names)
     write(53,'(F12.6,3(1X,F20.10))') t, Upot, kin, Upot+kin
     call time_step_Euler_pbc(dt, cutoff, pos, vel, Upot, kin)
     t = t + dt
  enddo
  close(53); close(57)

  ! dt = 1E-1
  call two_particles(1.01d0, pos); vel(:,:) = 0.d0; t = 0.d0; dt = 1.0d-1
  open(unit=54, file='energies_E_dt_-1E-01.dat',         status='replace', action='write')
  open(unit=58, file='trajectory_E_dt_-1E-01.xyz',       status='replace', action='write')
  do step = 0, nsteps
     if (mod(step,10) == 0) call write_xyz(58, pos, names)
     write(54,'(F12.6,3(1X,F20.10))') t, Upot, kin, Upot+kin
     call time_step_Euler_pbc(dt, cutoff, pos, vel, Upot, kin)
     t = t + dt
  enddo
  close(54); close(58)


   call system("gnuplot *.gnu")
END PROGRAM
