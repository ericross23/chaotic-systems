program test
  use fft_mod
  implicit none
  complex(kind=dp), dimension(8192) :: datax1, datax2
  double precision, dimension(8192) :: wx1, wx2
  integer :: i
  double precision :: t_f, t, dt, dw, Fs
  integer, parameter :: n = 4, nsteps = 8192
  real, parameter :: a = 0.0, b = 100.0
  real ::  x(n), ham, C, h
  x = (/0.0, 0.0, 4.0, 4.0/)
  h = (b - a)/nsteps

  open(10,file='fftinx1.csv')
  open(20,file='fftinx2.csv')
  open(30,file='fftoutx1.csv')
  open(40,file='fftoutx2.csv')



  dt = h
  Fs = 1.0 / dt
  !print*, Fs

  wx1(1) = 0
  wx2(1) = 0

  call rk4sys(n,h,x,nsteps,datax1,datax2)

  t = 0.d0
  do i = 1, size(datax1)
    write(10,*) t, realpart(datax1(i))
    t = t + dt
    if(i>1) wx1(i) = (i-1) * Fs / size(datax1)
  end do

  t = 0.d0
  do i = 1, size(datax2)
    write(20,*) t, realpart(datax2(i))
    t = t + dt
    if(i>1) wx2(i) = (i-1) * Fs / size(datax2)
  end do

  call fft(datax1)
  call fft(datax2)

  do i=1, size(datax1)
    write(30,*) wx1(i), sqrt(realpart(datax1(i))**2 + imagpart(datax1(i))**2)
    write(40,*) wx2(i), sqrt(realpart(datax2(i))**2 + imagpart(datax2(i))**2)
    write(*,*) wx1(i), sqrt(realpart(datax1(i))**2 + imagpart(datax1(i))**2)
  end do

end program test

subroutine xpsys(n,x,f)
      real, dimension (n) ::  x, f
      integer n
      real :: C
      C = -0.5
      f(1) = x(3)
      f(2) = x(4)
      f(3) = -x(1) - 4 * x(1) ** 3 - 4.0 * C * x(1) * x(2) ** 2
      f(4) = -x(2) - 4 * x(2) ** 3 - 4.0 * C * x(1) ** 2 * x(2)
end subroutine xpsys

subroutine rk4sys(n,h,x,nsteps,datax1,datax2)
    use fft_mod
      real ::  x(n)
      real, allocatable :: y(:), f(:,:)
      integer :: i, k, n, poin, nsteps
      real :: h
      complex(kind=dp), dimension(8192) :: datax1, datax2

      print *,0,x
      allocate (y(n), f(n,4))

      poin = 0
out:  do k = 1, nsteps+1000
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
        if(k>1000)datax1(k-1000) = x(1)
        if(k>1000)datax2(k-1000) = x(2)
      end do out
end subroutine rk4sys
