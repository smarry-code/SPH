subroutine time_print
!===============================================================
! TIME_PRINT    Print out the current date and time.
!
! Notes:
!   The standard Fortran 90 routine DATE_AND_TIME is used to get
!   the current date and time strings.
!===============================================================

    implicit none
    integer, parameter :: output = 6  ! Стандартный вывод (экран)
    
    ! Локальные переменные
    character(len = 8)  :: datstr     ! Строка даты (ГГГГММДД)
    character(len = 10) :: timstr     ! Строка времени (ЧЧММСС.МС)
    
    ! Получение текущей даты и времени
    call date_and_time(datstr, timstr)
    
    ! Форматированный вывод
    write(output, "(/a)") "*******************************************************"
    write(output, "(a)") "  SPH Simulation - Start Time"
    write(output, "(a)") "*******************************************************"
    write(output, "(a)") " Date: " // datstr(7:8) // "/" // datstr(5:6) // "/" // datstr(1:4)
    write(output, "(a)") " Time: " // timstr(1:2) // ":" // timstr(3:4) // ":" // timstr(5:10)
    write(output, "(a)") "*******************************************************"
    write(output, *)

end subroutine time_print