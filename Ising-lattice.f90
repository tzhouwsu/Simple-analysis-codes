program Ising2d

   integer, dimension(100,100) :: lattice, avglattice
   real :: irand
   real, dimension(100,100) :: jrand
   integer :: k,m,n
   real, dimension(10000) :: trand,orand
   real :: tk,tm
   integer, dimension(10000) :: torder
   integer, dimension(8) :: current_time
   integer, dimension(10000) :: iseed,jseed
   integer, dimension(100) :: left,right
   real :: H_field, J_coup, total_E, delta_E, sum_neighbs
   integer :: step

   H_field = 0.0
   J_coup = 1d-2   ! in unit of kT


!  initialize the lattice spin
   call DATE_AND_TIME(values=current_time)
   do k =1,10000
     jseed(k) = step * (current_time(8)+current_time(7)*1000)
   end do
   call RANDOM_SEED(put=jseed)
   call RANDOM_NUMBER(trand)
   do k = 1,100
     do m = 1,100
       if(trand(k*100+m) .ge. 0.5) then
         lattice(k,m) = -1
       else
         lattice(k,m) = 1
       endif
       avglattice(k,m) = 0
     end do
!     call sleep(1)
   end do

! define the left, right neighbors
   do k = 1,100
     if(k .eq. 1) then
       left(k) = 100
       right(k) = k+1
     elseif (k .eq. 100) then
       left(k) = k-1
       right(k) = 1
     else
       left(k) = k-1
       right(k) = k+1
     end if
   end do

! calculate the initial total energy
   total_E = 0.0
   do k = 1,100
     do m = 1,100
       sum_neighbs = 0.0

       sum_neighbs = lattice(left(k),m) + lattice(right(k),m) &
              & + lattice(k,left(m)) + lattice(k,right(m))
       total_E = total_E + lattice(k,m)*(-1.0*J_coup*sum_neighbs-1.0*H_field)
     end do
   end do


   do step = 1,1000

     call DATE_AND_TIME(values=current_time)
     do k =1,10000
       iseed(k) = step * (current_time(8) + current_time(7)*1000)
     end do
     call RANDOM_SEED(put=iseed)
     call RANDOM_NUMBER(trand)  ! this is to determine the order of updating the matrix
     torder(1) = 1
     orand(1) = trand(1)
     do k = 2,10000      ! sort the array to get the order-list in torder
       tk = trand(k)
       do m = 1,k-1
         tm = orand(m)
         if(tk .gt. tm) then
           exit
         end if
       end do
       if(m .lt. k) then 
         do n=k-1,m
           orand(n+1) = orand(n)
           torder(n+1) = torder(n)
         enddo
         orand(m) = trand(k)
         torder(m) = k
       else
         orand(k) = trand(k)
         torder(k) = k
       end if
     end do

!    create an matrix of random number, where it may be used for updating spins
     call DATE_AND_TIME(values=current_time)
     do k =1,10000
       jseed(k) = step * (current_time(8)+current_time(7)*1000)
     end do
     call RANDOM_SEED(put=jseed)
     call RANDOM_NUMBER(trand)

     do k = 1,10000   ! loop over the order-list  to update the lattice
       m = int(torder(k)/100) + 1
       n = mod(torder(k),100)
       if(n .eq. 0) then
         n = 100
         m = m - 1
       end if

!  update the lattice at site (m,n)
       sum_neighbs = 0.0
       sum_neighbs = 0.0 + lattice(m,right(n)) + lattice(m,left(n)) &
                & + lattice(right(m),n) + lattice(left(m),n)
       delta_E = (-2.0*lattice(m,n)) * (-1.0*J_coup*sum_neighbs-1.0*H_field)
       if(delta_E .le. 0.0) then
         lattice(m,n) = -1 * lattice(m,n)
       else
         if(trand(k) .lt. exp(-1.0*delta_E)) then
           lattice(m,n) = -1 * lattice(m,n)
         end if
       end if

     end do ! this is the end of k in order-list, i.e. the lattice has been fully updated

! calculate the total energy
   total_E = 0.0
   do k = 1,100
     do m = 1,100
       sum_neighbs = 0.0

       sum_neighbs = lattice(left(k),m) + lattice(right(k),m) &
              & + lattice(k,left(m)) + lattice(k,right(m))
       total_E = total_E + lattice(k,m)*(-1.0*J_coup*sum_neighbs-1.0*H_field)
       avglattice(k,m) = avglattice(k,m) + lattice(k,m)
     end do
   end do

     print *, step,torder(1),total_E

   end do


    call my_print_lattice(lattice,100)
    call my_print_lattice(avglattice,100)


contains

real function get_rand_num(i,num)
   integer, intent(in) :: i
   integer, intent(in) :: num
   integer :: iseed(num)
   real :: tempx(num)
   integer, dimension(8) :: current_time

   call DATE_AND_TIME(values=current_time)
!   print *,(current_time(k), k=1,8)
   iseed(1) = current_time(8) + 1000*current_time(7)
   call RANDOM_SEED(PUT=iseed)
!   call RANDOM_SEED()
   call RANDOM_NUMBER(HARVEST = tempx)
   get_rand_num = tempx(i)
end function get_rand_num


subroutine my_print_lattice(array,num)
   integer, intent(in) :: num
   integer, dimension(num,num), intent(in) :: array
   integer :: tempi, tempj
   do tempi = 1,num
      print *, (array(tempi,tempj),tempj=1,num)
   end do

end subroutine my_print_lattice

end program Ising2d



