MODULE crystal
  use GLOBAL
  IMPLICIT NONE
contains

! Only two concatenated particles, not a crystal, at a distance
! along x axis named dis.
    SUBROUTINE two_particles(dis, pos)
        use GLOBAL
        IMPLICIT NONE
        double precision, intent(in) :: dis
        double precision, allocatable, intent(out) :: pos(:,:)
        N = 2
        if (allocated(pos)) deallocate(pos)
        allocate(pos(N,3))
        pos(:,:) = 0.d0
        pos(2,1) = dis
    END SUBROUTINE



! SIMPLE CUBIC
    SUBROUTINE sc_lattice(M, rho, filename, pos)
        IMPLICIT NONE
        integer, intent(in):: M
        double precision, intent(in):: rho
        character(len=17), intent(in):: filename
        double precision, allocatable, intent(out):: pos(:,:)
        integer:: nx, ny, nz, i
        N = M**3 ! Number of particles
        L = (dble(N)/rho)**(1.d0/3.d0) ! Lenght of the box (size of system)
        a = L/M ! Lattice spacing
        if (allocated(pos)) deallocate(pos)
        allocate(pos(N,3))
        i = 1
        do nx = 0, M-1
            do ny = 0, M-1
                do nz = 0, M-1
                pos(i,1) = nx
                pos(i,2) = ny
                pos(i,3) = nz
                i = i + 1
                enddo
            enddo
        enddo
        pos(:,:) = pos(:,:)*a

    ! write conf_sc for atomic positions
    ! N
    ! (blank)
    ! A x1 x2 x3
    ! A x2 x2 x3
    ! ...
    open(100, file=filename)
    write(100,*) N
    write(100,*) ""
    do i = 1, N
        write(100,*) "SC", pos(i,1), pos(i,2), pos(i,3)
    enddo
    close(100)
    END SUBROUTINE

! FACE-CENTERED CUBIC
    SUBROUTINE fcc_lattice(M, rho, filename, pos)
        IMPLICIT NONE
        integer, intent(in):: M
        double precision, intent(in):: rho
        character(len=18), intent(in):: filename
        double precision, allocatable, intent(out):: pos(:,:)
        integer:: nx, ny, nz, i
        N = 4*M**3
        L = (dble(N)/rho)**(1.d0/3.d0)
        a = L/dble(M)
        if (allocated(pos)) deallocate(pos)
        allocate(pos(N,3))
        i = 1
        do nx = 0, M-1
            do ny = 0, M-1
                do nz = 0, M-1
                ! (0,0,0)
                pos(i,1) = dble(nx)
                pos(i,2) = dble(ny)
                pos(i,3) = dble(nz)
                i = i + 1
                ! (0,1/2,1/2)
                pos(i,1) = dble(nx)
                pos(i,2) = dble(ny) + 0.5d0
                pos(i,3) = dble(nz) + 0.5d0
                i = i + 1
                ! (1/2,0,1/2)
                pos(i,1) = dble(nx) + 0.5d0
                pos(i,2) = dble(ny)
                pos(i,3) = dble(nz) + 0.5d0
                i = i + 1
                ! (1/2,1/2,0)
                pos(i,1) = dble(nx) + 0.5d0
                pos(i,2) = dble(ny) + 0.5d0
                pos(i,3) = dble(nz)
                i = i + 1
                enddo
            enddo
        enddo
        pos(:,:) = pos(:,:)*a
    ! Write file
    open(200, file=filename)
    write(200,*) N
    write(200,*) ""
    do i = 1, N
        write(200,*) "FCC", pos(i,1), pos(i,2), pos(i,3)
    enddo
    close(200)
    END SUBROUTINE


! BODY-CENTERED CUBIC
    SUBROUTINE bcc_lattice(M, rho, filename, pos)
        IMPLICIT NONE
        integer, intent(in):: M
        double precision, intent(in):: rho
        character(len=18), intent(in):: filename
        double precision, allocatable, intent(out):: pos(:,:)
        integer:: nx, ny, nz, i
        N = 2*M**3
        L = (dble(N)/rho)**(1.d0/3.d0)
        a = L/dble(M)
        if (allocated(pos)) deallocate(pos)
        allocate(pos(N,3))
        i = 1
        do nx = 0, M-1
            do ny = 0, M-1
                do nz = 0, M-1
                ! (0,0,0)
                pos(i,1) = dble(nx)
                pos(i,2) = dble(ny)
                pos(i,3) = dble(nz)
                i = i + 1
                ! (1/2,1/2,1/2)
                pos(i,1) = dble(nx) + 0.5d0
                pos(i,2) = dble(ny) + 0.5d0
                pos(i,3) = dble(nz) + 0.5d0
                i = i + 1
                enddo
            enddo
        enddo
        pos(:,:) = pos(:,:)*a
    ! Write file
    open(300, file=filename)
    write(300,*) N
    write(300,*) ""
    do i = 1, N
        write(300,*) "BCC", pos(i,1), pos(i,2), pos(i,3)
    enddo
    close(300)
    END SUBROUTINE


END MODULE
