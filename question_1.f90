 !============================================================
    !AM5850 Tutorial_1 Question_1
    ! Solutions to a system of linear equations A*x=b
    ! Method: the basic elimination (simple Gauss elimination)
    !===========================================================
program linear_system
    implicit none
    integer,parameter::n=3
    integer::i
    double precision a(n,n), b(n), x(n)
    data a/1,2,3,1,-3,1,-2,1,4/
    data b/1,-8,7/
    CALL gauss_1(a,b,x,n)

    do i=1,n
        print 10,"x",i,"=",x(i)
    enddo
    10 format(2x,a,i1,a,f10.5)
end program linear_system

subroutine gauss_1(a,b,x,n)
    implicit none
    integer n
    double precision a(n,n), b(n), x(n)
    double precision c
    integer i, j, k
    !step 1: forward elimination
    do k=1, n-1
        ! using the first row as a standard with a constant multiplier 'c' to substract with below rows
        do i=k+1,n   
            !going through all the below rows from 2 to n
            c=a(i,k)/a(k,k) ! defining a constant c
            a(i,k) = 0.0
            b(i)=b(i)- c*b(k)
            do j=k+1,n
                 ! along the row  substitution except first coloumn
                a(i,j) = a(i,j)-c*a(k,j)
             end do
        end do 
    end do
    !step 2: back substitution
    x(n) = b(n)/a(n,n)
    do i=n-1,1,-1 !from last but one row to first
        c=0.0
        do j=i+1,n
            c= c + a(i,j)*x(j)
        end do
        x(i) = (b(i)- c)/a(i,i)
    end do
end subroutine gauss_1