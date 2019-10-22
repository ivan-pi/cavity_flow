module accelerate_mod

  use precision_mod, only: wp
  use iso_c_binding, only: c_double
  implicit none
  private

  public :: ramp_with_time

  interface 
    real(c_double) function my_ramp(x,p) bind(C,name="my_ramp")
      import c_double
      real(c_double), intent(in), value :: x
      real(c_double), intent(in), value :: p
    end function
  end interface

contains

  real(wp) function ramp_with_time(t,tf,p)
    integer, intent(in) :: t, tf
    real(wp), intent(in) :: p
    real(wp) :: x

    if (t < tf) then
      x = 2.0_wp*real(t,wp)/real(tf,wp) - 1.0_wp
      ramp_with_time = my_ramp(x,p)
    else
      ramp_with_time = 1.0_wp
    end if
  end function

end module