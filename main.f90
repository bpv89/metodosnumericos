
program progre
    implicit none
real :: dx,l,c
integer :: i,j,k
real :: nx
real, allocatable :: m(:,:)
real, allocatable :: r(:)
real, allocatable :: u(:)
real, allocatable :: m1(:,:)
real, allocatable :: r1(:)



write (*,*) "Tamanho de Delta x"
read *,dx
l = 1
nx = l / dx
allocate(m(nx,nx))
allocate(r(nx))
allocate(m1(nx,nx+1))
allocate(r1(nx))
allocate(u(nx))


do i=1 , nx

        if (i-1 > 0 .and. i+1 <= nx) then
        m(i,i-1) = 1
        m(i,i+1) =1
        m(i,i) = dx * dx -2
        r(i)= -dx*dx*dx*i

            else if (i-1 <= 0 .and. i+1 <= nx) then
            m(i,i) = dx * dx -2
            m(i,i+1) = 1
            r(i)= -dx*dx*dx*i

            else
            m(i,i-1) = 1
            m(i,i) = dx * dx -1
            r(i)= (-dx*dx*dx*i)-dx

        endif

end do

!Substituição de gaus
do i=1, nx
    do j=1, nx
        m1(i,j) = m(i,j)
    end do
    m1(i,nx+1) = r(i)
end do
 do i = 1, nx-1
    do k = i+1,nx
        c = m1(k,i)/m1(i,i)
        do j =i , nx+1
            m1(k,j)= m1(k,j)- c* m1(i,j)
        end do
    end do
 end do

u(nx) = m1(nx,nx+1)/m1(nx,nx)
do i = nx-1,1,-1
    c=0
    do j = i+1,nx
        c=c+ m1(i,j)* u(j)
    end do
    u(i) = (m1(i,nx+1)-c)/m1(i,i)
end do
! Fim Substituição de Gauss

open (10,file ="ResultadosProgressiva.csv" )

do i = 1, ubound(m, 1)
    write (10,*) m(i, :), r(i) , u(i)
end do

write (*,*) "Resultados impressos no arquivo"

close(10)


end program progre

