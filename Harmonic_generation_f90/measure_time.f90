! modules to measure execution Wall-time.

module timer_mod

    implicit none
    integer :: tick_rate, tick_start
contains
    subroutine tick()
        call system_clock(count_rate = tick_rate)
        call system_clock(count = tick_start)
    end subroutine tick

    subroutine tock(elapsed)
        real(kind=8), intent(out) :: elapsed
        integer :: tick_end
        call system_clock(count = tick_end)
        elapsed = real(tick_end - tick_start, 8) / real(tick_rate, 8)
    end subroutine tock

end module timer_mod

