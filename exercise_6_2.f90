! exercise_6_2.f90
! Ushijima p156
module subprogs
        implicit none
        contains
                subroutine gauss_jordan(a0, x, b, n)
                        integer, intent(in) :: n
                        real(8), intent(in) :: a0(n, n), b(n)
                        real(8), intent(out) :: x(n)
                        integer i, k
                        real(8) ar, a(n, n)

                        a(:, :) = a0(:, :)
                        x(:) = b(:)
                        do k =1, n
                                if (a(k, k) == 0.0d0 ) stop 'pivot = 0'
                                ar = 1.0d0 / a(k, k)
                                a(k, k) = 1.0d0
                                a(k, k + 1: n) = ar * a(k, k + 1 : n)
                                x(k) = ar * x(k)
                                do i = 1, n
                                        if (i /= k) then
                                                a(i, k + 1: n) = a(i, k + 1: n) - a(i, k) * a(k, k + 1: n)
                                                x(i) = x(i) - a(i, k) * x(k)
                                                a(i, k) = 0.0d0
                                        end if
                                end do
                        end do
                end subroutine gauss_jordan

                subroutine set_random_ab(a, b, x, n)
                        integer :: i, j
                        real :: uni_rand_num 
                        real(8), intent(inout), allocatable :: a(:, :), b(:), x(:)
                        integer, intent(inout) :: n

                        call random_seed
                        write(*, *) "Input n:"
                        read(*, *)n
                        allocate(a(1:n, 1:n))
                        allocate(b(1:n))
                        allocate(x(1:n))
                        do i = 1, n
                                do j = 1, n
                                        call random_number(uni_rand_num)
                                        a(i, j) = uni_rand_num
                                end do
                                call random_number(uni_rand_num)
                                b(i) = uni_rand_num
                        end do
                end subroutine set_random_ab
end module subprogs

program main
        use subprogs
        implicit none
        real(8), allocatable :: a(:, :), b(:), x(:), r(:)
        integer n
        real(8) t1, t2

        call cpu_time(t1)   
        call set_random_ab(a, b, x, n)
        call gauss_jordan(a, x, b, n)
        allocate(r(n))
        r(:) = b(:) - matmul(a, x)
        write(*, *) 'Gauss-Jordan error = ', dot_product(r, r)
        call cpu_time(t2)
        write(*, *) "cpu time =", t2 - t1
        deallocate(a, b, x)
end program main
