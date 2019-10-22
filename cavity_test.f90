program cavity_flow
    use precision_mod, only: wp
    use d2q9_mod, only: LatticeD2Q9
    use vtk_mod
    use accelerate_mod, only: ramp_with_time

    real(wp) :: reynolds
    integer  :: nx
    integer  :: ny

    real(wp) :: uwall, ulid
    real(wp) :: magic
    real(wp) :: ni
    real(wp) :: tau
    real(wp) :: omega, omega_2, dx

    integer :: tmax
    integer :: tplot
    integer :: tconv
    integer :: t, stream_nconv, tf

    real(wp) :: tol, stream_tol, stream_omega
    type(LatticeD2Q9) :: fluid
    real(wp), allocatable :: density(:,:), velocity(:,:,:), old_vel(:,:,:)
    real(wp) :: tstart, tfin
    integer :: itstart, itend, count_rate
    logical :: finish

    character(len=100) :: ctrl_file, output_folder

    ! read parameters
    if (command_argument_count() == 1) then
        call get_command_argument(1,ctrl_file)
        call read_control_file(ctrl_file)
    else
        ! default parameters
        reynolds = 100._wp
        nx = 120
        ny = 120
        uwall = 0.1_wp
        magic = 0.25_wp
        tmax = 1000000
        tplot = 10000
        tconv = 2000
        tol = 1.e-6_wp
        stream_nconv = 5000
        stream_tol = 1.e-6_wp
        stream_omega = 1.4_wp
        output_folder = "out"
    end if

    ! create folder for results
    call system("mkdir -p "//trim(adjustl(output_folder)))

    ! allocate result arrays
    allocate(density(nx,ny))
    allocate(velocity(nx,ny,2))
    allocate(old_vel(nx,ny,2))

    ! step size in real units
    dx = 1.0_wp/real(nx,wp)

    ! kinematic viscosity
    ni = uwall*real(nx,wp)/reynolds
    
    ! relaxation time and frequency
    tau = (3.0_wp*ni + 0.5_wp)
    omega = 1.0_wp/tau

    ! TRT magic parameter
    omega_2 = (magic/(tau-0.5_wp) + 0.5_wp)**(-1)

    ! initialize as empty
    old_vel = 0._wp

    ! initialize fluid
    fluid = LatticeD2Q9(nx,ny,omega)


    ! time until which we accelerate
    tf = int(50.*real(nx,wp)/uwall)

    
    call cpu_time(tstart)
    call system_clock(itstart,count_rate=count_rate)

    ! set fluid to set_equilibrium
    call fluid%set_equilibrium(1.0_wp)

    ! TIME STEPPING
    do t = 1, tmax

        ! COLLISION STEP
        ! call fluid%collide_and_stream_trt(omega_2)
        call fluid%collide_and_stream_risc()
        ! call fluid%collide_and_stream()

        ! set boundary conditions
        ulid = uwall*ramp_with_time(t,tf,2.0_wp)
        call fluid%cavity_bc(ulid)

        ! swap arrays
        call fluid%prepare_next_step()

        ! export results
        if (mod(t,tplot) == 0) then 
            density = fluid%get_density_all()
            velocity = fluid%get_velocity_all()
            ! call write_density("density_",density,t)
            ! call write_velocity("velocity_",velocity,t)
            call write_density_matrix("densmat",density,t)
            call write_velocity_matrix("velmat",velocity,t)
            print *, "Time = ", t, ulid
        end if

        !  check convergence
        if (mod(t,tconv) == 0) then
            velocity = fluid%get_velocity_all()
            call check_convergence(velocity,old_vel,tol,finish)
            old_vel = velocity
            if (finish) then
                print *, "Convergence reached at time = ", t
                exit
            end if
        end if
    end do

    call cpu_time(tfin)
    call system_clock(itend)
    print *, "MLUPS cpu_time = ", real(nx,wp)*real(ny,wp)*real(t,wp)/(tfin-tstart)/1.e6_wp
    print *, "MLUPS system_clock = ", real(nx,wp)*real(ny,wp)*real(t,wp)/(real(itend - itstart,wp)/real(count_rate,wp))/1.e6_wp
    
    ! get final results
    density = fluid%get_density_all()
    velocity = fluid%get_velocity_all()

    ! save as matrices
    call write_density_matrix("dens",density)
    call write_velocity_matrix("vel",velocity)
    !call write_fluid_vtk(nx,ny,density,velocity(:,:,1),velocity(:,:,2))
    !call write_fluid_vtk_binary(nx,ny,density,velocity(:,:,1),velocity(:,:,2))

    ! peform post processing
    call cavity_post_process(velocity,uwall,dx)

    ! END

contains
    
    subroutine read_control_file(filename)
        character(len=*), intent(in) :: filename

        integer, parameter :: fu = 56
        character(len=100) :: buffer, label
        integer :: pos
        integer :: ios = 0
        integer :: line = 0
        
        open(fu,file=filename,status='old')
        do while (ios == 0)
            read(fu,'(A)',iostat=ios) buffer
            if (ios == 0) then
                line = line + 1

                ! Find the first instance of whitespace. Split label and data
                pos = scan(buffer,' ')
                label = buffer(1:pos)
                buffer = buffer(pos+1:)

                select case (label)
                case("reynolds")
                    read(buffer,*,iostat=ios) reynolds
                    print *, "Read reynolds = ", reynolds
                case("size")
                    read(buffer,*,iostat=ios) nx
                    ny = nx
                    print *, "Read nx = ", nx
                    print *, "Read ny = ", ny
                case("uwall")
                    read(buffer,*,iostat=ios) uwall
                    print *, "Read uwall = ", uwall
                    if (uwall > 0.3_wp) stop "U is too large in lattice units"
                case("magic")
                    read(buffer,*,iostat=ios) magic
                    print *, "Read magic = ", magic
                case("tmax")
                    read(buffer,*,iostat=ios) tmax
                    print *, "Read tmax = ", tmax
                case("tplot")
                    read(buffer,*,iostat=ios) tplot
                    print *, "Read tplot = ", tplot
                case("tconv")
                    read(buffer,*,iostat=ios) tconv
                    print *, "Read tconv = ", tconv
                case("conv_tol")
                    read(buffer,*,iostat=ios) tol
                    print *, "Read convergence tolerance = ", tol
                case("stream_nconv")
                    read(buffer,*,iostat=ios) stream_nconv
                    print *, "Read stream SOR convergence interval = ",stream_nconv
                case("stream_tol")
                    read(buffer,*,iostat=ios) stream_tol
                    print *, "Read stream SOR function tolerance = ", stream_tol
                case("stream_omega")
                    read(buffer,*,iostat=ios) stream_omega
                    print *, "Read stream SOR omega = ", stream_omega
                case ("folder")
                    read(buffer,*,iostat=ios) output_folder
                    print *, "Output folder = ", trim(output_folder)
                case default
                    print *, 'Skipping invalid label at line', line
                end select
            end if
        end do
        close(fu)
    end subroutine read_control_file
                

    subroutine check_convergence(vel_new,vel_old,tol,is_exit)
        real(wp), intent(in) :: vel_new(:,:,:)
        real(wp), intent(in) :: vel_old(:,:,:)
        real(wp), intent(in) :: tol
        logical, intent(out) :: is_exit
        
        real(wp) :: norm_val

        norm_val = sum(norm2(vel_new(:,:,:)-vel_old(:,:,:),dim=3))/sum(norm2(vel_new(:,:,:),dim=3))

        print *, "convergence criterion = ", norm_val

        if (norm_val < tol) then
            is_exit = .true.
        else
            is_exit = .false.
        end if
    end subroutine check_convergence


    subroutine write_density(filename,density,t)
        character(len=*), intent(in) :: filename
        real(wp), intent(in) :: density(:,:)
        integer, intent(in) :: t

        integer :: i, j
        character(len=100) :: tnum

        write(tnum,*) t
        tnum = adjustl(tnum)

        open(22,file=trim(output_folder)//"/"//trim(filename)//trim(tnum)//'.out')

        do i = 1, size(density,dim=1)
            do j = 1, size(density,dim=2)
                write(22,*) i, j, density(i,j)
            end do
            write(22,*)
        end do

        close(22)
    end subroutine write_density


    subroutine write_velocity(filename,velocity,t)
        character(len=*), intent(in) :: filename
        real(wp), intent(in) :: velocity(:,:,:)
        integer, intent(in) :: t

        integer :: i, j
        character(len=100) :: tnum

        write(tnum,*) t
        tnum = adjustl(tnum)

        open(22,file=trim(output_folder)//"/"//trim(filename)//trim(tnum)//'.out')

        do i = 1, size(velocity,dim=1)
            do j = 1, size(velocity,dim=2)
                write(22,*) i, j, velocity(i,j,1), velocity(i,j,2)
            end do
            write(22,*)
        end do

        close(22)
    end subroutine write_velocity


    subroutine write_density_matrix(filename,density,t)
        character(len=*), intent(in) :: filename
        real(wp), intent(in) :: density(:,:)
        integer, intent(in), optional :: t

        integer :: i, j
        character(len=100) :: tnum

        if (present(t)) then
            write(tnum,*) t
            tnum = adjustl(tnum)
            open(22,file=trim(output_folder)//"/"//trim(filename)//trim(tnum)//".out")
        else
            open(22,file=trim(output_folder)//"/"//trim(filename)//".out")
        end if

        do j = 1, size(density,dim=2)
            write(22,*) (density(i,j),i=1,size(density,dim=1))
        end do

        close(22)

    end subroutine write_density_matrix

    subroutine write_velocity_matrix(filename,velocity,t)
        character(len=*), intent(in) :: filename
        real(wp), intent(in) :: velocity(:,:,:)
        integer, intent(in), optional :: t

        integer :: i, j
        character(len=100) :: tnum

        if (present(t)) then
            write(tnum,*) t
            tnum = adjustl(tnum)
            open(22,file=trim(output_folder)//"/"//trim(filename)//"x"//trim(tnum)//".out")
            open(24,file=trim(output_folder)//"/"//trim(filename)//"y"//trim(tnum)//".out")
        else
            open(22,file=trim(output_folder)//"/"//trim(filename)//"x"//".out")
            open(24,file=trim(output_folder)//"/"//trim(filename)//"y"//".out")
        end if

        do j = 1, size(density,dim=2)
            write(22,*) (velocity(i,j,1),i=1,size(velocity,dim=1))
            write(24,*) (velocity(i,j,2),i=1,size(velocity,dim=1))
        end do

        close(22)
        close(24)

    end subroutine write_velocity_matrix


    subroutine cavity_post_process(velocity,uwall,dx,nconv,tolconv,omegaconv)
        real(wp), intent(in) :: velocity(nx,ny,2)
        real(wp), intent(in) :: uwall, dx
        integer, intent(in), optional :: nconv
        real(wp), intent(in), optional :: tolconv, omegaconv
        real(wp) :: stream(nx+1,ny+1)
        real(wp) :: vorticity(nx+1,ny+1)
        real(wp) :: U(nx+1,ny+1), V(nx+1,ny+1), sf(nx+1,ny+1)
        integer :: i, j, niter
        integer :: my_nconv = 10000
        real(wp) :: my_omegaconv = 1.2_wp
        real(wp) :: my_tolconv = 1.0e-6_wp
        real(wp) :: conv, stemp

        if (present(nconv)) my_nconv = nconv
        if (present(tolconv)) my_tolconv = tolconv
        if (present(omegaconv)) my_omegaconv = omegaconv


        print *, "interpolating velocity field"
        ! linear interpolate velocity onto new grid
        V = 0.0_wp
        U = 0.0_wp
        U(:,ny+1) = uwall
        do i = 2, nx
            do j = 2, ny
                U(i,j) = 0.25_wp*(velocity(i-1,j-1,1)+velocity(i,j-1,1)+velocity(i,j,1)+velocity(i-1,j,1))
                V(i,j) = 0.25_wp*(velocity(i-1,j-1,2)+velocity(i,j-1,2)+velocity(i,j,2)+velocity(i-1,j,2))
            end do
        end do

        ! normalize velocity field to region [0, 1]
        U = U/uwall
        V = V/uwall

        print *, "calculating vorticity field"
        vorticity = 0.0_wp
        ! calculate vorticity at internal points
        do i = 2, nx
            do j = 2, ny
                vorticity(i,j) = 0.5_wp*(V(i+1,j)-V(i-1,j)) - 0.5_wp*(U(i,j+1)-U(i,j-1))
            end do
        end do
        ! left wall
        do j = 2, ny
            vorticity(1,j) = -1.5_wp*V(1,j) + 2.0_wp*V(2,j) - 0.5_wp*V(3,j)
        end do
        ! right wall
        do j = 2, ny
            vorticity(nx+1,j) = 1.5_wp*V(nx+1,j) - 2.0_wp*V(nx,j) + 0.5_wp*V(nx-1,j)
        end do
        ! top wall
        do i = 2, nx
            vorticity(i,ny+1) = -(1.5_wp*U(i,ny+1) - 2.0_wp*U(i,ny) + 0.5_wp*U(i,ny-1))
        end do
        ! bottom wall
        do i = 2, nx
            vorticity(i,1) = -(-1.5_wp*U(i,1) + 2.0_wp*U(i,2) - 0.5_wp*U(i,3))
        end do
        vorticity = vorticity/dx

        stream = 0.0_wp
        sf = 0.0_wp
        niter = 0
        ! calculate stream function with SOR
        outer : do
            do i = 2, nx
                do j = 2, ny
                    stemp = 0.25_wp*(stream(i+1,j)+stream(i-1,j)+stream(i,j+1)+stream(i,j-1) + dx**2*vorticity(i,j))
                    stream(i,j) = my_omegaconv*stemp + (1.0_wp-my_omegaconv)*stream(i,j)
                end do
            end do
            niter = niter + 1
            if (mod(niter,my_nconv) == 0) then
                conv = sum(abs(stream-sf))/sum(abs(stream))
                print *, "niter = ",niter
                print *, "Convergence criterion = ", conv
                sf = stream
                if (conv < my_tolconv) then
                    print *, "Stream function has converged"
                    exit outer
                end if
            end if
        end do outer

        ! compare stream function to velocity field

        ! velocity
        open(22,file=trim(output_folder)//"/ppvelx.out")
        open(23,file=trim(output_folder)//"/ppvely.out")
        open(24,file=trim(output_folder)//"/stream.out")
        open(25,file=trim(output_folder)//"/vortic.out")

        do j = 1, ny+1
            write(22,*) (U(i,j),i=1,nx+1)
            write(23,*) (V(i,j),i=1,nx+1)
            write(24,*) (stream(i,j),i=1,nx+1)
            write(25,*) (vorticity(i,j),i=1,nx+1)
        end do

        close(22)
        close(23)
        close(24)
        close(25)

        print *, "Finished post-processing"

    end subroutine cavity_post_process

end program