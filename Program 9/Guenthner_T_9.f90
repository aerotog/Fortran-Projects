PROGRAM cont_mech

!Declare variables and arrays
implicit none
real::I,E,L,p0
integer::k,h,z,p,r,j
real,allocatable::A(:,:),B(:),x(:),y(:),S(:),y2(:),C(:)

!Open and start data txt
OPEN(3,File='Guenthner_T_9.txt')
WRITE(3,*)''
WRITE(3,*)'Written by Timothy Guenthner'

!Format codes for later in program
10 FORMAT (1x,A9, 1x, i3.1)
11 FORMAT(7x,a1,11x,a12,18x,a7)
12 FORMAT (3x,f5.1,11x,f12.8,15x,f12.8)

!Define known variables
I=10.4
E=30000000
L=100
p0=10

!Start do loop for main part of code
DO k=1,4
        !Get stepsize input from user
        WRITE(*,*) 'Stepsize?'
        READ(*,*)z

        !Find array size paramters        
        h=L/z
        p=h-1
        r=h+1

        !Allocate arrays
        ALLOCATE(A(p,p))
        ALLOCATE(B(p))
        ALLOCATE(x(r))
        ALLOCATE(y(r))
        ALLOCATE(S(r))
        ALLOCATE(y2(r))
        ALLOCATE(C(r))

        !Fill in matrix A
        A(1,1)=-2
        DO j=2,p
                A(j,j)=-2
                A(j,j-1)=1
                A(j-1,j)=1
        END DO

        !Fill in array x
        x(1)=0
        x(r)=0
        DO j=2,r
                x(j)=x(j-1)+z
        END DO

        !Use known variables along with array x to calculate B
        DO j=1,p
                B(j)=((z**2)/(E*I))*(p0/2)*(x(j+1))*(x(j+1)-L)
        END DO

        !Find exact solution to given function for displacement
        DO j =1,r
                y2(j)=(p0/(24*E*I))*(x(j)**4-2*L*(x(j))**3+x(j)*L**3)
        END DO

        !Call Gaussian Elimination and solution subroutines to solve system
        Call GE (A,B,p)
        Call Sol (A,B,p,S)

        !Put solutions into final answer array 
        C(1)=0
        C(r)=0 
        DO j=1,p
                C(j+1)=S(j)
        END DO

        !Write data to txt
        WRITE(3,*)''
        WRITE(3,10)'Stepsize:', z
        WRITE(3,11)'x', 'Approximate', 'Exact'
        DO j=1,r
                 WRITE(3,12)x(j),C(j),y2(j)
        END DO

        !Deallocate arrays for use in next loop
        DEALLOCATE(A)
        DEALLOCATE(B)
        DEALLOCATE(x)
        DEALLOCATE(y)
        DEALLOCATE(S)
        DEALLOCATE(C)
        DEALLOCATE(y2)

END DO

END PROGRAM cont_mech

!===================================================


SUBROUTINE GE (A,B,p)
implicit none
integer, intent(inout)::p
real, Intent(inout)::B(p),A(p,p)
Integer::i, j, k

DO i = 1, p
        DO j=(i+1),p
                A(i,j)= A(i,j)/A(i,i)
        END DO
        DO j = (i+1),p
                DO k = (i+1),p
                        A(k,j) = A(k,j)-A(i,j)*A(k,i)
                END DO
                B(j) = B(j) - A(i,j)*B(i)
        END DO
END DO

DO i =1,p-1
        A(i,(i+1):p) = 0
END DO

END SUBROUTINE GE

!===================================================

SUBROUTINE Sol (A,B,p,S)
implicit none
integer,intent(inout)::p
real, intent(in)::A(p,p), B(p)
real, intent(out)::S(p)
integer::i, j

S(p)= B(p)/A(p,p)

DO i = p,1,-1
        S(i) = B(i)
        DO j = p,(i+1),-1
                S(i) = S(i) - S(j)*A(j,i)
        END DO
        S(i) = S(i)/A(i,i)
END DO

END SUBROUTINE sol


