module VTK_mod
    use precision_mod, only: wp
    implicit none
    
    private
    public :: write_fluid_vtk, write_vtk_2dscalar, write_fluid_vtk_binary

contains


    subroutine write_fluid_vtk(nx,ny,rho,ux,uy,t)
        integer, intent(in) :: nx, ny
        real(wp), intent(in) :: rho(nx,ny)
        real(wp), intent(in) :: ux(nx,ny)
        real(wp), intent(in) :: uy(nx,ny)
        integer, intent(in), optional :: t
        
        character(len=24) :: tnum
        integer :: i, j 
        write(tnum,*) t
        tnum = adjustl(tnum)

        call system("mkdir -p vtk")

        if (present(t)) then
            open(10,file="vtk/fluid_t"//trim(tnum)//".vtk")
        else
            open(10,file="vtk/fluid.vtk")
        end if

        write(10,"(A)") "# vtk DataFile Version 3.0"
        write(10,"(A)") "fluid"
        write(10,"(A)") "ASCII"
        write(10,"(A)") "DATASET RECTILINEAR_GRID"
        write(10,"(A,I5,I5,I5)") "DIMENSIONS ", nx, ny, 1
        write(10,"(A,I5,A)") "X_COORDINATES ", nx, " float"
        do i = 1, nx
            write(10,"(ES12.5)") real(i-1,wp)
        end do
        write(10,"(A,I5,A)") "Y_COORDINATES ", ny, " float"
        do j = 1, ny
            write(10,"(ES12.5)") real(j-1,wp)
        end do
        write(10,"(A,I5,A)") "Z_COORDINATES ", 1, " float"
        write(10,"(ES12.5)") 0.0_wp
        write(10,"(A,I5)") "POINT_DATA ", nx*ny
        write(10,"(A)") "SCALARS density float 1"
        write(10,"(A)") "LOOKUP_TABLE default"
        do j = 1, ny
            do i = 1, nx
                write(10,"(ES12.5)") rho(i,j)
            end do
        end do
        write(10,*) "VECTORS velocity float"
        do j = 1, ny
            do i = 1, nx
                write(10,"(3ES12.5)") ux(i,j), uy(i,j), 0.0_wp
            end do
        end do
        
        close(10)   

    end subroutine write_fluid_vtk

    function int_to_string(my_int) result(my_string)
        integer, intent(in) :: my_int
        character(len=8) :: my_string
        
        write(my_string,'(i8)') my_int
        
    end function int_to_string


    subroutine write_fluid_vtk_binary(nx,ny,rho,ux,uy,t)
        integer, intent(in) :: nx
        integer, intent(in) :: ny
        real(wp), intent(in) :: rho(nx,ny), ux(nx,ny), uy(nx,ny)
        integer, intent(in), optional :: t
        character(len=1) :: lf = achar(10) ! line feed character
        character(len=8) :: cx, cy, cz
        real(wp) :: vel(3,nx*ny)
        integer :: vtk = 10
        
        character(len=24) :: tnum
        integer :: i, j, k 
        
        do j = 1, ny
            do i = 1, nx
                k = (j-1)*ny+i
                vel(1,k) = ux(i,j)
                vel(2,k) = uy(i,j)
                vel(3,k) = 0.0_wp
            end do
        end do

        write(tnum,*) t
        tnum = adjustl(tnum)

        call system("mkdir -p vtk")

        if (present(t)) then
            open(unit=vtk,file="vtk/fluid_t"//trim(tnum)//".vtk",access="stream",convert="BIG_ENDIAN")
        else
            open(unit=vtk,file="vtk/fluid.vtk",access="stream",convert="BIG_ENDIAN")
        end if

        cx = int_to_string(nx)
        cy = int_to_string(ny)
        cz = int_to_string(1)

        write(vtk) '# vtk DataFile Version 3.0'//lf
        write(vtk) 'fluid'//lf
        write(vtk) 'BINARY'//lf
        write(vtk) 'DATASET RECTILINEAR_GRID'//lf
        write(vtk) 'DIMENSIONS '//cx//cy//cz//lf
        write(vtk) 'X_COORDINATES '//cx//' DOUBLE'//lf
        write(vtk) (real(i-1,wp),i=1,nx), lf
        write(vtk) 'Y_COORDINATES '//cy//' DOUBLE'//lf
        write(vtk) (real(j-1,wp),j=1,ny)
        write(vtk) 'Z_COORDINATES '//cz//' DOUBLE'//lf
        write(vtk) 0.0_wp, lf
        write(vtk) "POINT_DATA", int_to_string(nx*ny), lf
        write(vtk) "SCALARS density DOUBLE 1"//lf
        write(vtk) "LOOKUP_TABLE default",lf
        write(vtk) rho,lf
        write(vtk) "VECTORS velocity DOUBLE 3"//lf
        write(vtk) ((vel(i,j),i=1,3),j=1,nx*ny),lf

        close(vtk)        
        
    end subroutine write_fluid_vtk_binary


    subroutine write_vtk_2dscalar(filename,fieldname,x,y,scalar)
        character(len=*), intent(in) :: filename, fieldname
        real(wp), intent(in) :: x(:), y(:)
        real(wp), intent(in) :: scalar(:,:)
        character(len=1) :: lf = achar(10) ! line feed character
        character(len=8) :: nx, ny, nz
        integer :: vtk = 10
        
        nx = int_to_string(size(x))
        ny = int_to_string(size(y))
        nz = int_to_string(1)

        open(unit=vtk,file=filename,access='stream',convert='BIG_ENDIAN')
        
        write(vtk) '# vtk DataFile Version 3.0'//lf
        write(vtk) fieldname//lf
        write(vtk) 'BINARY'//lf
        write(vtk) 'DATASET RECTILINEAR_GRID'//lf
        write(vtk) 'DIMENSIONS '//nx//ny//nz//lf
        write(vtk) 'X_COORDINATES '//nx//' DOUBLE'//lf
        write(vtk) x, lf
        write(vtk) 'Y_COORDINATES '//ny//' DOUBLE'//lf
        write(vtk) y, lf
        write(vtk) 'Z_COORDINATES '//nz//' DOUBLE'//lf
        write(vtk) 0.0d0, lf
        write(vtk) "POINT_DATA", int_to_string(size(x)*size(y)), lf
        write(vtk) "SCALARS density DOUBLE 1"//lf
        write(vtk) "LOOKUP_TABLE default",lf
        write(vtk) scalar,lf
        
        close(vtk)        
        
    end subroutine write_vtk_2dscalar


end module VTK_mod