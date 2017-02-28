program main
      integer, parameter :: n = 4, nsteps = 10000
      real, parameter :: a = 0.0, b =100.0
      real ::  x(n)
      x = (/0.0, 0.0, 4.0, 8.0/)
      h = (b - a)/nsteps
      print*, h
      call rk4sys(n,h,x,nsteps)
end program main

subroutine xpsys(n,x,f)
      real, dimension (n) ::  x, f
      integer n
      real :: C
      C = 5.2
      f(1) = x(3)
      f(2) = x(4)
      f(3) = -x(1) - 4 * x(1) ** 3 - 4.0 * C * x(1) * x(2) ** 2
      f(4) = -x(2) - 4 * x(2) ** 3 - 4.0 * C * x(1) ** 2 * x(2)
end subroutine xpsys

subroutine rk4sys(n,h,x,nsteps)
      real ::  x(n)
      real, allocatable :: y(:), f(:,:)
      integer :: i, k, n, poin
      real :: h, xold
      print *,0,x
      allocate (y(n), f(n,4))
      open(10,file='phase.csv')
      open(20,file='poincare.csv')
      open(30,file='phase2.csv')
      open(40,file='phase3.csv')
      open(50,file='phase4.csv')
      open(60,file='3d.csv')

      poin = 0
out:  do k = 1, nsteps
        call xpsys(n,x,f(1,1))
in1:    do i = 1,n
          y(i) = x(i) + 0.5*h*f(i,1)
        end do in1
        call xpsys(n,y,f(1,2))
in2:    do i = 1,n
          y(i) = x(i) + 0.5*h*f(i,2)
        end do in2
        call xpsys(n,y,f(1,3))
in3:    do i = 1,n
          y(i) = x(i) + h*f(i,3)
        end do in3
        call xpsys(n,y,f(1,4))
in4:    do i = 1,n
          x(i) = x(i) + (h/6.0)* (f(i,1) + 2.0*(f(i,2) + f(i,3)) + f(i,4))
        end do in4

        print *, k, x
        !if(k>1000) write(10,*) x(1), x(2)
        !if(k>1000) write(30,*) x(1), x(3)
        !if(k>1000) write(40,*) x(2), x(4)
        !if(k>1000) write(50,*) x(3), x(4)
        if(k>1000) write(60,*) x(1), x(2), x(3)
        if((xold < 0.d0 .and. x(1) > 0.d0 .and. k>1000) .or. &
           (xold > 0.d0 .and. x(1) < 0.d0 .and. k>1000)) then
          write(20,*) x(2), x(3), x(4)
          poin = poin + 1
          !print*, poin, k
        end if
        xold = x(1)
        if(poin >= 100000) exit
      end do out
end subroutine rk4sys
