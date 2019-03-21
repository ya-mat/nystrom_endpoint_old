! -*- coding: utf-8 -*-
module mod_force_raise
  implicit none
contains
!--------------------------------------------
  subroutine force_raise()
    implicit none

    integer :: zero

    zero = 0
    write(*,*) '# force raise !!'
    write(*,*) '# error', 0/zero
    stop

  end subroutine force_raise
!--------------------------------------------
end module mod_force_raise
