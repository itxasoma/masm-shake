PROGRAM main
  use GLOBAL
  use force
  use motion
  use io_teacher
  use rdf_module
  IMPLICIT NONE

  integer :: step, nsteps, nmolecules
  integer :: nf
  integer :: eq_steps, traj_stride, sample_stride
  double precision :: dt_ps, tauT_ps, Tref_K, sigma_A, epsil_K, mass_gmol, r0_A
  double precision :: dt, tauT, Tref
  double precision :: box_A, uvel_Aps, utime_ps
  double precision :: t
  double precision, allocatable :: posA(:,:), velA(:,:)
  double precision, allocatable :: pos(:,:),  vel(:,:)
  double precision :: Upot, kin, Etot, temp, lambda
  character(len=2), allocatable :: names(:)

  ! -------------------------
  ! Read inputs/files
  ! -------------------------
  call read_input_dades('exercicishake.dades', nsteps, dt_ps, tauT_ps, nmolecules, Tref_K, sigma_A, epsil_K, mass_gmol, r0_A)
  call read_conf_data('conf.data', nmolecules, posA, velA, box_A)

  call to_reduced_units(nmolecules, sigma_A, epsil_K, mass_gmol, dt_ps, tauT_ps, Tref_K, posA, velA, box_A, &
                        dt, tauT, Tref, pos, vel, L, uvel_Aps, utime_ps)

  N = 3*nmolecules
  nmol = nmolecules
  atoms_per_mol = 3

  cutoff = 2.5d0
  rho = dble(N) / (L**3)

  nf = 6*nmolecules - 3   

  allocate(names(N))
  names(:) = 'Ar'

  ! Controls
  traj_stride = 10
  eq_steps = max(0, min(200, nsteps/5))
  sample_stride = 1

  t = 0.d0

  ! -------------------------
  ! Outputs
  ! -------------------------
  open(unit=10, file='trajectory.xyz', status='replace', action='write')
  open(unit=20, file='energies_T.dat', status='replace', action='write')
  write(20,'(A)') '# t  Upot  K  Etot  T  lambda'

  ! RDF
  call rdf_init(200, L)

  write(*,*) '--- Running NVT (Berendsen) Velocity-Verlet ---'
  write(*,*) 'N=', N, ' L=', L, ' rho=', rho
  write(*,*) 'dt=', dt, ' tauT=', tauT, ' Tref=', Tref, ' nsteps=', nsteps
  write(*,*) 'Equil steps (no RDF sampling)=', eq_steps

  do step = 1, nsteps

    call time_step_VelocityVerlet_NVT(dt, cutoff, tauT, Tref, nf, pos, vel, Upot, kin, temp, lambda)
    t = t + dt
    Etot = Upot + kin

    if (mod(step, traj_stride) == 0) call write_xyz(10, pos, names)
    write(20,'(F12.6,5(1X,F20.10))') t, Upot, kin, Etot, temp, lambda

    if (step > eq_steps) then
      if (mod(step, sample_stride) == 0) call rdf_sample(pos)
    endif
    
  enddo

  close(10)
  close(20)

  call rdf_write('gr_ArAr.dat', rho)

  write(*,*) 'Done.'
  write(*,*) 'Wrote: trajectory.xyz, energies_T.dat, gr_ArAr.dat'

END PROGRAM
