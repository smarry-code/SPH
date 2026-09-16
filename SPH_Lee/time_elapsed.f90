subroutine time_elapsed(s)
!===============================================================
! The standard Fortran 90 routine RTC is used to calculate the 
! elapsed CPU time since the start of the program.
! 
! Note: 
!   This implementation is compiler-specific (Intel Fortran).
!   For portability, consider using CPU_TIME instead.
!===============================================================

    use dfport  ! Модуль специфичный для компилятора Intel
    implicit none
    real(8), intent(out) :: s  ! Выходное значение: время в секундах
    
    ! Вызов функции RTC() из модуля dfport
    s = rtc()
    
end subroutine time_elapsed