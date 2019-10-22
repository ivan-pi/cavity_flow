module d2q9_mod
    use precision_mod, only: wp
    use omp_lib
    implicit none

    private
    public LatticeD2Q9

        ! lattice weights
    real(wp), parameter :: w(0:8) = [4.0_wp/9.0_wp,&
                                     1.0_wp/9.0_wp, 1.0_wp/9.0_wp, 1.0_wp/9.0_wp, 1.0_wp/9.0_wp,&
                                     1.0_wp/36.0_wp,1.0_wp/36.0_wp,1.0_wp/36.0_wp,1.0_wp/36.0_wp]

        ! lattice speed of sound    
    real(wp), parameter :: cs = 1.0_wp/sqrt(3.0_wp)
    real(wp), parameter :: inv_cs2 = 1.0_wp/cs**2
    real(wp), parameter :: inv_cs4 = 1.0_wp/cs**4

    real(wp), allocatable, target :: f1(:,:,:)
    real(wp), allocatable, target :: f2(:,:,:)

    type :: LatticeD2Q9
        integer :: nx
        integer :: ny
        real(wp) :: omega
        real(wp) :: one_min_omega
        real(wp) :: flag
        real(wp), pointer, contiguous  :: fpre(:,:,:) => null()
        real(wp), pointer, contiguous :: fpost(:,:,:) => null()
    contains
        procedure :: set_equilibrium
        procedure :: collide_and_stream
        procedure :: collide_and_stream_risc
        procedure :: collide_and_stream_trt
        procedure :: prepare_next_step
        procedure :: get_density_at
        procedure :: get_density_all
        procedure :: get_velocity_at
        procedure :: get_velocity_all
        procedure :: cavity_bc
    end type LatticeD2Q9

    interface LatticeD2Q9
        module procedure :: new_lattice
    end interface LatticeD2Q9

contains

    function new_lattice(nx,ny,omega) result (lattice)
        integer, intent(in) :: nx
        integer, intent(in) :: ny
        real(wp), intent(in) :: omega
        type(LatticeD2Q9) :: lattice

            ! save parameters
        lattice%nx = nx
        lattice%ny = ny
        lattice%omega = omega
        lattice%one_min_omega = 1.0_wp - omega

            ! allocate array space
        allocate(f1(0:nx+1,0:ny+1,0:8))
        allocate(f2(0:nx+1,0:ny+1,0:8))

            ! connect pre and post collision pointers
            ! to storage arrays
        lattice%fpre => f1
        lattice%fpost => f2

            ! set flag for pointer switch
        lattice%flag = 0
    end function new_lattice
    
    !
    ! SET EQUILIBRIUM
    !
    subroutine set_equilibrium(self,rho)
        class(LatticeD2Q9), intent(inout) :: self
        real(wp), intent(in), optional :: rho
        real(wp) :: my_rho
        integer :: i, j

        my_rho = 1.0_wp
        if (present(rho)) my_rho = rho

        do i = 1, self%nx
            do j = 1, self%ny
                self%fpre(i,j,0) = w(0)*my_rho
                self%fpre(i,j,1) = w(1)*my_rho
                self%fpre(i,j,2) = w(2)*my_rho
                self%fpre(i,j,3) = w(3)*my_rho
                self%fpre(i,j,4) = w(4)*my_rho
                self%fpre(i,j,5) = w(5)*my_rho
                self%fpre(i,j,6) = w(6)*my_rho
                self%fpre(i,j,7) = w(7)*my_rho
                self%fpre(i,j,8) = w(8)*my_rho
            end do
        end do

    end subroutine set_equilibrium

    !
    ! COLLIDE AND STREAM
    !
    subroutine collide_and_stream(self)
        class(LatticeD2Q9), intent(inout) :: self
        real(wp) :: rho, rhoinv, ux, uy, usqr, ne(0:8)
        integer :: i, j

        associate(fpre => self%fpre,fpost => self%fpost, omega => self%omega, one_min_omega => self%one_min_omega)

        do i = 1, self%nx
            do j = 1, self%ny
                    ! macroscopic variables
                rho = sum(fpre(i,j,:))
                rhoinv = 1.0_wp/rho
                ux = rhoinv*(fpre(i,j,1)+fpre(i,j,5)+fpre(i,j,8)-fpre(i,j,3)-fpre(i,j,6)-fpre(i,j,7))
                uy = rhoinv*(fpre(i,j,2)+fpre(i,j,5)+fpre(i,j,6)-fpre(i,j,4)-fpre(i,j,7)-fpre(i,j,8))
                usqr = ux*ux + uy*uy

                    ! non-equlibrium part of the distribution function
                ne(0)  = 0.5_wp*w(0)*rho*(2.0_wp - inv_cs2*usqr)
                ne(1)  = 0.5_wp*w(1)*rho*(2.0_wp + 2.0_wp*inv_cs2*ux + inv_cs4*ux*ux - inv_cs2*usqr)
                ne(2)  = 0.5_wp*w(2)*rho*(2.0_wp + 2.0_wp*inv_cs2*uy + inv_cs4*uy*uy - inv_cs2*usqr)
                ne(3)  = 0.5_wp*w(3)*rho*(2.0_wp - 2.0_wp*inv_cs2*ux + inv_cs4*ux*ux - inv_cs2*usqr)
                ne(4)  = 0.5_wp*w(4)*rho*(2.0_wp - 2.0_wp*inv_cs2*uy + inv_cs4*uy*uy - inv_cs2*usqr)
                ne(5)  =        w(5)*rho*(1.0_wp + inv_cs2*(ux + uy) + inv_cs4*ux*uy + inv_cs2*usqr)
                ne(6)  =        w(6)*rho*(1.0_wp - inv_cs2*(ux - uy) - inv_cs4*ux*uy + inv_cs2*usqr)
                ne(7)  =        w(7)*rho*(1.0_wp - inv_cs2*(ux + uy) + inv_cs4*ux*uy + inv_cs2*usqr)
                ne(8)  =        w(8)*rho*(1.0_wp + inv_cs2*(ux - uy) - inv_cs4*ux*uy + inv_cs2*usqr)

                    ! perform collide and stream
                fpost(i,j,0)     = one_min_omega*fpre(i,j,0) + omega*ne(0)
                fpost(i+1,j,1)   = one_min_omega*fpre(i,j,1) + omega*ne(1)
                fpost(i,j+1,2)   = one_min_omega*fpre(i,j,2) + omega*ne(2)
                fpost(i-1,j,3)   = one_min_omega*fpre(i,j,3) + omega*ne(3)
                fpost(i,j-1,4)   = one_min_omega*fpre(i,j,4) + omega*ne(4)
                fpost(i+1,j+1,5) = one_min_omega*fpre(i,j,5) + omega*ne(5)
                fpost(i-1,j+1,6) = one_min_omega*fpre(i,j,6) + omega*ne(6)
                fpost(i-1,j-1,7) = one_min_omega*fpre(i,j,7) + omega*ne(7)
                fpost(i+1,j-1,8) = one_min_omega*fpre(i,j,8) + omega*ne(8)
            end do
        end do

        end associate
    end subroutine collide_and_stream


    subroutine collide_and_stream_risc(self)
        class(LatticeD2Q9), intent(inout) :: self
        real(wp) :: rho(self%nx), rhoinv, ux(self%nx), uy(self%nx), usqr(self%nx), ne(self%nx,0:8)
        integer :: i, j

        !$omp parallel shared(self) private(j,i,rho,rhoinv,ux,uy,usqr,ne)
        associate(fpre => self%fpre,fpost => self%fpost, omega => self%omega, one_min_omega => self%one_min_omega)

        !$omp do schedule(static)
        do j = 1, self%ny
            do i = 1, self%nx
                    ! macroscopic variables
                rho(i) = sum(fpre(i,j,:))
                rhoinv = 1.0_wp/rho(i)
                ux(i) = rhoinv*(fpre(i,j,1)+fpre(i,j,5)+fpre(i,j,8)-fpre(i,j,3)-fpre(i,j,6)-fpre(i,j,7))
                uy(i) = rhoinv*(fpre(i,j,2)+fpre(i,j,5)+fpre(i,j,6)-fpre(i,j,4)-fpre(i,j,7)-fpre(i,j,8))
                usqr(i) = ux(i)*ux(i) + uy(i)*uy(i)
            end do

            do i = 1, self%nx
                    ! non-equlibrium part of the distribution function
                ne(i,0)  = 0.5_wp*w(0)*rho(i)*(2.0_wp - inv_cs2*usqr(i))
                ne(i,1)  = 0.5_wp*w(1)*rho(i)*(2.0_wp + 2.0_wp*inv_cs2*ux(i) + inv_cs4*ux(i)*ux(i) - inv_cs2*usqr(i))
                ne(i,2)  = 0.5_wp*w(2)*rho(i)*(2.0_wp + 2.0_wp*inv_cs2*uy(i) + inv_cs4*uy(i)*uy(i) - inv_cs2*usqr(i))
                ne(i,3)  = 0.5_wp*w(3)*rho(i)*(2.0_wp - 2.0_wp*inv_cs2*ux(i) + inv_cs4*ux(i)*ux(i) - inv_cs2*usqr(i))
                ne(i,4)  = 0.5_wp*w(4)*rho(i)*(2.0_wp - 2.0_wp*inv_cs2*uy(i) + inv_cs4*uy(i)*uy(i) - inv_cs2*usqr(i))
                ne(i,5)  =        w(5)*rho(i)*(1.0_wp + inv_cs2*(ux(i) + uy(i)) + inv_cs4*ux(i)*uy(i) + inv_cs2*usqr(i))
                ne(i,6)  =        w(6)*rho(i)*(1.0_wp - inv_cs2*(ux(i) - uy(i)) - inv_cs4*ux(i)*uy(i) + inv_cs2*usqr(i))
                ne(i,7)  =        w(7)*rho(i)*(1.0_wp - inv_cs2*(ux(i) + uy(i)) + inv_cs4*ux(i)*uy(i) + inv_cs2*usqr(i))
                ne(i,8)  =        w(8)*rho(i)*(1.0_wp + inv_cs2*(ux(i) - uy(i)) - inv_cs4*ux(i)*uy(i) + inv_cs2*usqr(i))
            end do

                    ! perform collide and stream
            do i = 1, self%nx
                fpost(i,j,0)     = one_min_omega*fpre(i,j,0) + omega*ne(i,0)
            end do
            do i = 1, self%nx
                fpost(i+1,j,1)   = one_min_omega*fpre(i,j,1) + omega*ne(i,1)
            end do
            do i = 1, self%nx
                fpost(i,j+1,2)   = one_min_omega*fpre(i,j,2) + omega*ne(i,2)
            end do
            do i = 1, self%nx    
                fpost(i-1,j,3)   = one_min_omega*fpre(i,j,3) + omega*ne(i,3)
            end do
            do i = 1, self%nx
                fpost(i,j-1,4)   = one_min_omega*fpre(i,j,4) + omega*ne(i,4)
            end do
            do i = 1, self%nx
                fpost(i+1,j+1,5) = one_min_omega*fpre(i,j,5) + omega*ne(i,5)
            end do
            do i = 1, self%nx
                fpost(i-1,j+1,6) = one_min_omega*fpre(i,j,6) + omega*ne(i,6)
            end do
            do i = 1, self%nx
                fpost(i-1,j-1,7) = one_min_omega*fpre(i,j,7) + omega*ne(i,7)
            end do
            do i = 1, self%nx
                fpost(i+1,j-1,8) = one_min_omega*fpre(i,j,8) + omega*ne(i,8)
            end do
        end do
        !$omp end do

        end associate
        !$omp end parallel
    end subroutine collide_and_stream_risc


    subroutine collide_and_stream_trt(self,omega_2)
        class(LatticeD2Q9), intent(inout) :: self
        real(wp), intent(in) :: omega_2
        real(wp) :: rho(self%nx), rhoinv, ux(self%nx), uy(self%nx), usqr(self%nx), ne(self%nx,0:8)
        real(wp) :: fp, fm, fpeq, fmeq
        integer :: i, j

        !$omp parallel shared(self) private(j,i,rho,rhoinv,ux,uy,usqr,ne)
        associate(fpre => self%fpre,fpost => self%fpost, omega => self%omega, one_min_omega => self%one_min_omega)

        !$omp do schedule(static)
        do j = 1, self%ny
            do i = 1, self%nx
                    ! macroscopic variables
                rho(i) = sum(fpre(i,j,:))
                rhoinv = 1.0_wp/rho(i)
                ux(i) = rhoinv*(fpre(i,j,1)+fpre(i,j,5)+fpre(i,j,8)-fpre(i,j,3)-fpre(i,j,6)-fpre(i,j,7))
                uy(i) = rhoinv*(fpre(i,j,2)+fpre(i,j,5)+fpre(i,j,6)-fpre(i,j,4)-fpre(i,j,7)-fpre(i,j,8))
                usqr(i) = ux(i)*ux(i) + uy(i)*uy(i)
            end do

            do i = 1, self%nx
                    ! non-equlibrium part of the distribution function
                ne(i,0)  = 0.5_wp*w(0)*rho(i)*(2.0_wp - inv_cs2*usqr(i))
                ne(i,1)  = 0.5_wp*w(1)*rho(i)*(2.0_wp + 2.0_wp*inv_cs2*ux(i) + inv_cs4*ux(i)*ux(i) - inv_cs2*usqr(i))
                ne(i,2)  = 0.5_wp*w(2)*rho(i)*(2.0_wp + 2.0_wp*inv_cs2*uy(i) + inv_cs4*uy(i)*uy(i) - inv_cs2*usqr(i))
                ne(i,3)  = 0.5_wp*w(3)*rho(i)*(2.0_wp - 2.0_wp*inv_cs2*ux(i) + inv_cs4*ux(i)*ux(i) - inv_cs2*usqr(i))
                ne(i,4)  = 0.5_wp*w(4)*rho(i)*(2.0_wp - 2.0_wp*inv_cs2*uy(i) + inv_cs4*uy(i)*uy(i) - inv_cs2*usqr(i))
                ne(i,5)  =        w(5)*rho(i)*(1.0_wp + inv_cs2*(ux(i) + uy(i)) + inv_cs4*ux(i)*uy(i) + inv_cs2*usqr(i))
                ne(i,6)  =        w(6)*rho(i)*(1.0_wp - inv_cs2*(ux(i) - uy(i)) - inv_cs4*ux(i)*uy(i) + inv_cs2*usqr(i))
                ne(i,7)  =        w(7)*rho(i)*(1.0_wp - inv_cs2*(ux(i) + uy(i)) + inv_cs4*ux(i)*uy(i) + inv_cs2*usqr(i))
                ne(i,8)  =        w(8)*rho(i)*(1.0_wp + inv_cs2*(ux(i) - uy(i)) - inv_cs4*ux(i)*uy(i) + inv_cs2*usqr(i))
            end do

                    ! perform collide and stream
            do i = 1, self%nx
                fpost(i,j,0)     = one_min_omega*fpre(i,j,0) + omega*ne(i,0)
            end do
            do i = 1, self%nx
                fp = 0.5_wp*(fpre(i,j,1)+fpre(i,j,3))
                fm = 0.5_wp*(fpre(i,j,1)-fpre(i,j,3))
                fpeq = 0.5_wp*(ne(i,1)+ne(i,3))
                fmeq = 0.5_wp*(ne(i,1)-ne(i,3))
                fpost(i+1,j,1)   = fpre(i,j,1) - omega*(fp-fpeq) - omega_2*(fm-fmeq)
            end do
            do i = 1, self%nx
                fp = 0.5_wp*(fpre(i,j,2)+fpre(i,j,4))
                fm = 0.5_wp*(fpre(i,j,2)-fpre(i,j,4))
                fpeq = 0.5_wp*(ne(i,2)+ne(i,4))
                fmeq = 0.5_wp*(ne(i,2)-ne(i,4))
                fpost(i,j+1,2)   = fpre(i,j,2) - omega*(fp-fpeq) - omega_2*(fm-fmeq)
            end do
            do i = 1, self%nx
                fp = 0.5_wp*(fpre(i,j,3)+fpre(i,j,1))
                fm = 0.5_wp*(fpre(i,j,3)-fpre(i,j,1))
                fpeq = 0.5_wp*(ne(i,3)+ne(i,1))
                fmeq = 0.5_wp*(ne(i,3)-ne(i,1))    
                fpost(i-1,j,3)   = fpre(i,j,3) - omega*(fp-fpeq) - omega_2*(fm-fmeq)
            end do
            do i = 1, self%nx
                fp = 0.5_wp*(fpre(i,j,4)+fpre(i,j,2))
                fm = 0.5_wp*(fpre(i,j,4)-fpre(i,j,2))
                fpeq = 0.5_wp*(ne(i,4)+ne(i,2))
                fmeq = 0.5_wp*(ne(i,4)-ne(i,2))
                fpost(i,j-1,4)   = fpre(i,j,4) - omega*(fp-fpeq) - omega_2*(fm-fmeq)
            end do
            do i = 1, self%nx
                fp = 0.5_wp*(fpre(i,j,5)+fpre(i,j,7))
                fm = 0.5_wp*(fpre(i,j,5)-fpre(i,j,7))
                fpeq = 0.5_wp*(ne(i,5)+ne(i,7))
                fmeq = 0.5_wp*(ne(i,5)-ne(i,7))
                fpost(i+1,j+1,5) = fpre(i,j,5) - omega*(fp-fpeq) - omega_2*(fm-fmeq)
            end do
            do i = 1, self%nx
                fp = 0.5_wp*(fpre(i,j,6)+fpre(i,j,8))
                fm = 0.5_wp*(fpre(i,j,6)-fpre(i,j,8))
                fpeq = 0.5_wp*(ne(i,6)+ne(i,8))
                fmeq = 0.5_wp*(ne(i,6)-ne(i,8))
                fpost(i-1,j+1,6) = fpre(i,j,6) - omega*(fp-fpeq) - omega_2*(fm-fmeq)
            end do
            do i = 1, self%nx
                fp = 0.5_wp*(fpre(i,j,7)+fpre(i,j,5))
                fm = 0.5_wp*(fpre(i,j,7)-fpre(i,j,5))
                fpeq = 0.5_wp*(ne(i,7)+ne(i,5))
                fmeq = 0.5_wp*(ne(i,7)-ne(i,5))
                fpost(i-1,j-1,7) = fpre(i,j,7) - omega*(fp-fpeq) - omega_2*(fm-fmeq)
            end do
            do i = 1, self%nx
                fp = 0.5_wp*(fpre(i,j,8)+fpre(i,j,6))
                fm = 0.5_wp*(fpre(i,j,8)-fpre(i,j,6))
                fpeq = 0.5_wp*(ne(i,8)+ne(i,6))
                fmeq = 0.5_wp*(ne(i,8)-ne(i,6))
                fpost(i+1,j-1,8) = fpre(i,j,8) - omega*(fp-fpeq) - omega_2*(fm-fmeq)
            end do
        end do
        !$omp end do

        end associate
        !$omp end parallel
    end subroutine collide_and_stream_trt


    subroutine prepare_next_step(self)
        class(LatticeD2Q9), intent(inout) :: self

        if (self%flag == 0) then
            self%fpre => f2
            self%fpost => f1
            self%flag = 1
        else
            self%fpre => f1
            self%fpost => f2
            self%flag = 0
        end if

    end subroutine prepare_next_step

    
    function get_density_at(self,x,y,to_pressure) result(rho)
        class(LatticeD2Q9), intent(in) :: self
        integer, intent(in) :: x, y
        logical, intent(in), optional :: to_pressure
        logical :: my_to_pressure = .false.
        real(wp) :: rho

            ! save pressure flag
        if (present(to_pressure)) my_to_pressure = to_pressure

            ! calculate density at position x, y 
        rho = sum(self%fpre(x,y,:))

            ! convert to pressure
        if (my_to_pressure) then
            rho = cs**2*rho
        end if

    end function get_density_at


    function get_density_all(self,to_pressure) result(rho)
        class(LatticeD2Q9), intent(in) :: self
        logical, intent(in), optional :: to_pressure
        real(wp) :: rho(self%nx,self%ny)
        logical :: my_to_pressure = .false.
        integer :: i, j

            ! save pressure flag
        if (present(to_pressure)) my_to_pressure = to_pressure

            ! calculate densities 
        do i = 1, self%nx
            do j = 1, self%ny
                rho(i,j) = sum(self%fpre(i,j,:))
            end do
        end do

            ! convert to pressure
        if (my_to_pressure) then
            rho = cs**2*rho
        end if

    end function get_density_all


    function get_velocity_at(self,x,y) result(velocity)
        class(LatticeD2Q9), intent(in) :: self
        integer, intent(in) :: x, y
        real(wp) :: velocity(2)
        real(wp) :: rho

        associate(fpre => self%fpre)
            
            ! get density
        rho = self%get_density_at(x,y)
            ! calculate velocity
        velocity(1) = (fpre(x,y,1)+fpre(x,y,5)+fpre(x,y,8)-fpre(x,y,3)-fpre(x,y,6)-fpre(x,y,7))/rho
        velocity(2) = (fpre(x,y,2)+fpre(x,y,5)+fpre(x,y,6)-fpre(x,y,4)-fpre(x,y,7)-fpre(x,y,8))/rho

        end associate
    end function get_velocity_at


    function get_velocity_all(self) result(velocity)
        class(LatticeD2Q9), intent(in) :: self
        real(wp) :: velocity(self%nx,self%ny,2)
        integer :: i, j
        real(wp) :: rho
        associate(fpre => self%fpre)

        do i = 1, self%nx
            do j = 1, self%ny            
                    ! get density
                rho = sum(self%fpre(i,j,:))
                    ! calculate velocity
                velocity(i,j,1) = (fpre(i,j,1)+fpre(i,j,5)+fpre(i,j,8)-fpre(i,j,3)-fpre(i,j,6)-fpre(i,j,7))/rho
                velocity(i,j,2) = (fpre(i,j,2)+fpre(i,j,5)+fpre(i,j,6)-fpre(i,j,4)-fpre(i,j,7)-fpre(i,j,8))/rho
            end do
        end do

        end associate
    end function get_velocity_all

    subroutine cavity_bc(self,uwall)
        class(LatticeD2Q9), intent(inout) :: self
        real(wp), intent(in) :: uwall
        integer :: x, y
        real(wp) :: rho

        associate(fpost => self%fpost, fpre => self%fpre, nx => self%nx, ny => self%ny)

            ! top wall
        do x = 1, nx
            rho = sum(fpre(x,ny,:))
            fpost(x,ny,4) = fpost(x,ny+1,2)
            fpost(x,ny,7) = fpost(x+1,ny+1,5) - 2.0_wp*inv_cs2*w(7)*rho*uwall
            fpost(x,ny,8) = fpost(x-1,ny+1,6) + 2.0_wp*inv_cs2*w(8)*rho*uwall
        end do

            ! bottom wall
        do x = 1, nx
            fpost(x,1,2) = fpost(x,  0,4)
            fpost(x,1,6) = fpost(x+1,0,8)
            fpost(x,1,5) = fpost(x-1,0,7)
        end do
        
            ! left wall
        do y = 1, ny
            fpost(1,y,1) = fpost(0,y  ,3)
            fpost(1,y,5) = fpost(0,y-1,7)
            fpost(1,y,8) = fpost(0,y+1,6)
        end do

            ! right wall
        do y = 1, ny
            fpost(nx,y,3) = fpost(nx+1,y  ,1)
            fpost(nx,y,6) = fpost(nx+1,y-1,8)
            fpost(nx,y,7) = fpost(nx+1,y+1,5)
        end do

        ! ! bottom left corner
        ! fpost(1,1,1) = fpost(0,1,3)
        ! fpost(1,1,2) = fpost(1,0,4)
        ! fpost(1,1,5) = fpost(0,0,7) 

        ! ! top left corner
        ! fpost(1,ny,1) =  fpost(0,ny,3)
        ! fpost(1,ny,4) =  fpost(1,ny+1,2)
        ! fpost(1,ny,8) =  fpost(0,ny+1,6)

        ! ! bottom right corner
        ! fpost(nx,1,2) = fpost(nx,0,4)
        ! fpost(nx,1,3) = fpost(nx+1,1,1)
        ! fpost(nx,1,6) = fpost(nx+1,0,8)

        ! ! top right corner
        ! fpost(nx,ny,3) = fpost(nx+1,ny,1)
        ! fpost(nx,ny,4) = fpost(nx,ny+1,2)
        ! fpost(nx,ny,7) = fpost(nx+1,ny+1,5)

        end associate
    end subroutine cavity_bc
end module