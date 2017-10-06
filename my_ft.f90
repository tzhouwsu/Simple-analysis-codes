! here I want to calculate fourier transformation of a function f(x) = (f_r(x), f_i(x))
! the output is g(w) = int_(0 to inf) {f(x)*exp(iwt)dt}, use summation as integral
! by Tiecheng Zhou

module iFT_mod
  implicit none

  type func_data   ! this is the data of initial function f(t)
    integer :: n  ! this is the number of data points
    real :: dt   ! this is the time separation 
    complex, dimension(:), allocatable :: fx
    complex, dimension(:), allocatable :: x

   contains
    procedure :: get_obj  ! read data from a file
    procedure :: del_obj  ! delete it
    procedure :: prt_obj  ! print it
  end type func_data

  interface func_data
    procedure :: init_obj  ! constructror
  end interface func_data

 contains
  type(func_data) function init_obj()
    integer :: ni = 100000  ! default array size for data point in f(t)
    real :: dti = 0.01d0  ! default time separation of data point in f(t)
    allocate(init_obj%fx(ni))
    allocate(init_obj%x(ni))
    init_obj%n = ni
    init_obj%dt = dti
  end function init_obj

  subroutine get_obj(fdata,fname)   ! read the original function data from a file
    class(func_data), intent(inout) :: fdata
    character(LEN=100), intent(in):: fname
    real :: dtj, tempdtj, temptj, temprj, tempij
    integer :: num, iread, dsize
    logical :: iexist

    dsize = size(fdata%x)
    num = 0;
    inquire(file=fname, exist=iexist)
    if(.not. iexist) then
       print *, 'Error: cannot find the input-file'
    else
      open(unit=11, file=fname)
      num = 1
      iread = 0
      do
        read(11,*,IOSTAT=iread) temptj, temprj, tempij  ! read tj (pure real), and f(tj) (real, and imaginary part)
        if(iread /= 0) then ! at the end of file , or there is some error
          exit
        end if
        fdata%x(num) = cmplx(temptj,0.0d0)    ! tj is always pure real
        fdata%fx(num) = cmplx(temprj,tempij) ! f(tj) is a complex number
        if(num .eq. 1) then
          dtj = real(fdata%x(num))
        else if(num .eq. 2) then
          dtj = real(fdata%x(num))-real(fdata%x(num-1))
        else
          tempdtj = real(fdata%x(num)) - real(fdata%x(num-1))
          if( abs(dtj-tempdtj) .ge. 0.01d-2 ) then
             print *, 'Error: t is not evenly sampled', num, dtj, tempdtj
             exit
          end if
        end if
        num = num + 1
        if(num > dsize) then
          print *, 'data exceeds the default arraysize'
        end if
      end do

      if(iread > 0) then
        print *, 'Error: something is wrong in the input-file'
      else if (iread < 0) then   ! reach the end of the file
        fdata%n = num - 1  
        fdata%dt = dtj
      end if

    end if
  end subroutine get_obj

  subroutine del_obj(ddata)
    class(func_data) :: ddata
    deallocate(ddata%x)
    deallocate(ddata%fx)
  end subroutine del_obj

  subroutine prt_obj(pdata)
    class(func_data), intent(in) :: pdata
    integer :: pi
    pi = 1
    do while(pi .le. pdata%n)
      print *, real(pdata%x(pi)),aimag(pdata%x(pi)),real(pdata%fx(pi)),aimag(pdata%fx(pi))
      pi = pi + 1
    end do
  end subroutine prt_obj

  


end module iFT_mod


program myFTcode
use iFT_mod
  implicit none
  type(func_data) :: ft, fw
  character(LEN=100) :: finput = 'input'   ! name of the input-file, that has f(t)
  character(LEN=100) :: foutput = 'output'  ! the output-file
  real, parameter :: CONVERT = 2.99793d-2   ! conversion factor, 1 cm-1 = 2.99793*10+10 Hz (s-1) = 2.99793*10-2 (ps-1)
  real, parameter :: PI = 3.1415926535 
  integer :: ss=1

  ft=func_data()
  fw=func_data()
  call get_obj(ft,finput)
!  call prt_obj(ft)
!  print *, CONVERT
!  do while (ss .le. 10)
!     print *, ss
!     ss = ss + 1
!  end do

  call get_FT_func(fw,ft)     ! exam 0 to 4000 cm-1
!  call prt_obj(fw)

  open(unit = 13, file=foutput)
    ss = 1
    do while (ss .le. fw%n)
      write(13,*) real(fw%x(ss)),aimag(fw%x(ss)),real(fw%fx(ss)),aimag(fw%fx(ss)), &
         & real(fw%fx(ss))*real(fw%fx(ss))+aimag(fw%fx(ss))*aimag(fw%fx(ss))
      ss = ss + 1
    end do
  close(13)


 contains
   subroutine get_FT_func(fw,ft,dw,nw)
      class(func_data), intent(in) :: ft
      class(func_data), intent(inout) :: fw
      real, intent(in), optional :: dw
      integer, intent(in), optional :: nw
      integer :: ni, nj
      real :: wi, tj
      real :: rpart, ipart

 !     print *, CONVERT
      if(present(dw)) then
        fw%dt = dw * CONVERT ! if you enter dw in cm-1, convert to ps-1
      else
        fw%dt = 1.0d0/(ft%n*ft%dt)   ! resolution in frequency domain, in unit of ps-1
      end if
      if(present(nw)) then
        fw%n = nw
      else
        fw%n = ft%n
      end if

      ni = 1
      do while(ni .le. fw%n)  ! loop over each w data point in fw, the fourier transformed function
        wi = ni*fw%dt*2.0d0*PI ! sampled w, it should be pure real, and in unit of rad*ps-1
!        fw%x(ni) = cmplx(ni*fw%dt/CONVERT, 0.0d0)  ! output g(w), w in unit of cm-1
        fw%x(ni) = cmplx(ni*fw%dt, 0.0d0)  ! output g(w), w in unit of ps-1
        rpart = 0.0d0
        ipart = 0.0d0   ! initialize it
        nj = 1
        do while (nj .le. ft%n)
          tj = nj*ft%dt  ! in unit of ps, so that wi*tj is in unit of rad
!          rpart = rpart + (cos(wi*tj)*real(ft%fx(nj)) - sin(wi*tj)*aimag(ft%fx(nj)))*ft%dt
!          ipart = ipart + (cos(wi*tj)*aimag(ft%fx(nj)) + sin(wi*tj)*real(ft%fx(nj)))*ft%dt
          rpart = rpart + (cos(wi*tj)*real(ft%fx(nj)) - sin(wi*tj)*aimag(ft%fx(nj)))*ft%dt*sqrt(1-(nj+0.0d0)/(ft%n+0.0d0)) ! multiply a damping func
          ipart = ipart + (cos(wi*tj)*aimag(ft%fx(nj)) + sin(wi*tj)*real(ft%fx(nj)))*ft%dt*sqrt(1-(nj+0.0d0)/(ft%n+0.0d0))
          nj = nj + 1
        end do
        fw%fx(ni) = cmplx(rpart,ipart)
        ni = ni + 1
      end do

   end subroutine get_FT_func
   

end program myFTcode


