PROGRAM RK4

!Declare variables and open data txt
implicit none
integer::i
real::x(15), y(15),t(15),h,K1x,K2x,K3x,K4x,K1y,K2y,K3y,K4y

Open(3,File='Guenthner_T_8.txt')

!Define atepsize and initial values
h = 0.1
x(1) = 1
y(1) = 1

!Formats for code later on
10 FORMAT(1x,a4,4x,a7,7x,a7)
11 FORMAT(2x,i2.0,3x,f12.5,3x,f12.5)
12 FORMAT (8x,a3,15x,a3,12x,a3,12x,a3)
13 FORMAT (1x,4(f13.5,3x))
14 FORMAT (2x,a2,3x,f12.5,2x,f12.5)

!Start data sheet
WRITE(3,*)''
WRITE(3,*)'Written by Timothy Guenthner'
WRITE(3,*)''
WRITE(3,10)'step','x','y'
WRITE(3,14)'1',x(1),y(1)
WRITE(3,*)''
WRITE(3,*)''

!Initiate do loop for main part of code
DO i = 1,14
        !Define inital t
        t(1)=0
        !Call subroutine to find k values
        CALL Kval (h,x(i),y(i),t(i),K1x,K2x,K3x,K4x,K1y,K2y,K3y,K4y)

        !Use calculated k values to perform Runge Kutta
        x(i+1) = x(i) + (K1x + 2*K2x + 2*K3x + K4x)/6
        y(i+1) = y(i) + (K1y + 2*K2y + 2*K3y + K4y)/6
        t(i+1) = t(i)+h

        !Write to data; include h, x, y and k values
        WRITE(3,*)'================================================================'
        WRITE(3,*)''
        WRITE(3,10)'stepsize','x','y'        
        WRITE(3,11)i+1,x(i+1),y(i+1)
        WRITE(3,*)''
        WRITE(3,12)'K1x','K2x','K3x','K4x'
        WRITE(3,13)K1x,K2x,K3x,K4x
        WRITE(3,*)''
        WRITE(3,12)'K1y','K2y','K3y','K4y'
        WRITE(3,13)K1y,K2y,K3y,K4y
        WRITE(3,*)''        
        WRITE(3,*)''


END DO

END PROGRAM RK4

!=====================================

SUBROUTINE Kval (h,x,y,t,K1x,K2x,K3x,K4x,K1y,K2y,K3y,K4y)

!Declare varibles and determine input/ouput
implicit none
real,intent(in)::h,x,y,t
real,intent(out)::K1x,K2x,K3x,K4x,K1y,K2y,K3y,K4y

!Use given equations to calculate each k(1 to 4 for x and y)
K1x = (x*y+t)*h
K1y = (t*y+x)*h
K2x = ((x+K1x/2)*(y+K1y/2)+(t+h/2))*h
K2y = ((t+h/2)*(y+K1y/2)+(x+K1x/2))*h
K3x = ((x+K2x/2)*(y+K2y/2)+(t+h/2))*h
K3y = ((t+h/2)*(y+K2y/2)+(x+K2x/2))*h
K4x = ((x+K3x)*(y+K3y)+(t+h))*h
K4y = ((t+h)*(y+K3y)+(x+K3x))*h

END SUBROUTINE Kval

