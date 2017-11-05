module donnees
  implicit none

contains
  
  function alpha(x)
    real*8, intent(in) :: x
    real*8 :: alpha

    alpha = 1.

  end function alpha

  function f(t,x)
    real*8,intent(in) :: t,x
    real*8 :: f

    f = 0.

  end function f

  function u0(x)
    real*8, intent(in) :: x
    real*8 :: u0

    u0 = 3.*sin(x)

  end function u0
    
end module donnees
