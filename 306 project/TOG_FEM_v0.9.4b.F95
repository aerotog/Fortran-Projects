PROGRAM FEM
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!PROGRAM WRITTEN BY TIMOTHY GUENTHNER
!ALL VARIABLES MUST BE DECLARED
IMPLICIT NONE
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!HOUSEKEEPING
!DECLARE VARIABLES
REAL::LENGTH(100),THETA(100),E(100),A(100),IYY(100),F(101,3),PX(100),PZ(100),kb(6,6)
REAL::STIFF(303,303),FORCE(303),GLOBK(6,6),TRANS_TRANS(6,6),TRANS(6,6),FDIST(6)
REAL::GLOBFDIST(6),STIFF_REACT(303,303),Q(303),GLOBREACT(303),QB(6),FORCEINT(6)
REAL::QSUB(6),FORCE_INT(6),DENS(100)
INTEGER::NUMEL,NUMNODES,ELCONNECT(100,2),NUMBC,BC(303,2)
INTEGER::i,g

!INITIALIZE TXT DOCUMENTS
OPEN(1,FILE='FEM_DATA.txt')
OPEN(2,FILE='FEM_OUTPUT.txt')
OPEN(3,FILE='DEBUG.txt')
!/HOUSEKEEPING
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!MAIN BODY OF PROGRAM
!GATHER INPUT INFO
WRITE(*,*)'USE GRAVITY? (1 - yes, 0 - no)'
READ(*,*)g
CALL INPUT(LENGTH,THETA,E,A,IYY,F,NUMEL,NUMNODES,ELCONNECT,NUMBC,BC,PX,PZ,DENS,g)

!INITILIZE MATRICES
CALL INITIALIZE(kb,STIFF,NUMEL,NUMNODES,FORCE)

!LOOP FOR EACH ELEMENT, BUILDING GLOBAL K AND Q MATRICES
DO i=1,NUMEL
        CALL ELESTIFF(i,E,A,IYY,LENGTH,kb)
        CALL TRANSFORMK(i,kb,GLOBK,THETA,TRANS_TRANS,TRANS)
        CALL ELEFORCE(i,LENGTH,PX,PZ,FDIST)
        CALL TRANSFORMQ(i,TRANS_TRANS,FDIST,GLOBFDIST)
        CALL ASSEMBLEKQ(i,ELCONNECT,STIFF,STIFF_REACT,GLOBK,FORCE,F,GLOBFDIST)
END DO

!APPLY THE BOUNDARY CONDITIONS
CALL APPLYBC(NUMBC,BC,STIFF,FORCE,NUMNODES)

!SOLVE THE SYSTEM OF EQS
CALL SOLVE(STIFF,FORCE,NUMEL,Q,NUMNODES)

!SOLUTION SECTION OF DATA OUT TXT
WRITE(2,*)'~~~~~~~DATA OUTPUT~~~~~~~'
!CALCULATE AND OUTPUT DATA ON A PER ELEMENT BASIS
DO i=1,NUMEL
        CALL ELESTIFF(i,E,A,IYY,LENGTH,kb)
        CALL TRANSFORMK(i,kb,GLOBK,THETA,TRANS_TRANS,TRANS)
        CALL ELEFORCE(i,LENGTH,PX,PZ,FDIST)
        CALL REACTIONS(i,NUMNODES,STIFF_REACT,Q,GLOBREACT,TRANS,kb,QB&
                &,FDIST,FORCE_INT,ELCONNECT,QSUB,NUMEL)
        CALL ELEOUTPUT(i,QSUB,FORCE_INT)
END DO

!OUTPUT DATA PER GLOBAL BASIS
CALL GLOBOUTPUT(Q,FORCE,GLOBREACT)


!/MAIN BODY OF PROG
END PROGRAM FEM

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!SUBROUTINE: INPUT

SUBROUTINE INPUT(LENGTH,THETA,E,A,IYY,F,NUMEL,NUMNODES,ELCONNECT,NUMBC,BC,PX,PZ,DENS,g)
IMPLICIT NONE
REAL::X(100),Z(100),DELNX,DELNZ,LENGTH(100),THETA(100)
REAL::E(100),A(100),IYY(100),PX(100),PZ(100),F(101,3),DENS(100)
INTEGER::NUMEL,NUMNODES,i,j,ND1,ND2,ELCONNECT(100,2),NUMBC,BC(303,2),g

!NEED # OF ELEMENTS AND NODES
!# OF ELEMENTS
READ(1,*)NUMEL
!# OF NODES
READ(1,*)NUMNODES

!BUILD TRUSS SYSTEM VIA COORDINATES
DO i=1,NUMNODES
        !INPUT X AND Z COORDINATES OF ELEMENT i
        READ(1,*)X(i),Z(i)
END DO

!DEFINE CONNECTIVITY OF ELEMENTS
DO i=1,NUMEL
        !DEFINE ELEMENT CONNETIVITY STARTING WITH ELEMENT i
        READ(1,*)(ELCONNECT(i,j),j=1,2)
END DO

!REASSIGN NODES TO PROG VARIABLES AND CALC. ANGLE AND LENGTH
DO i=1,NUMEL
        ND1=ELCONNECT(i,1)
        ND2=ELCONNECT(i,2)
        DELNX=X(ND2)-X(ND1)
        DELNZ=Z(ND2)-Z(ND1)
        LENGTH(i)=SQRT(DELNX**2+DELNZ**2)
        THETA(i)=ATAN2(DELNZ,DELNX)
END DO

!READ AND DEFINE INPUT MATERIAL PROPERTIES
DO i=1,NUMEL
        READ(1,*)E(i),A(i),IYY(i),DENS(i)
END DO

!DEFINE ANY DISTRIBUTED LOADS
DO i=1,NUMEL
        READ(1,*)PX(i),PZ(i)
        IF(g.EQ.1) THEN
                PX(i)=PX(i)-A(i)*DENS(i)*SIN(THETA(i))
                PZ(i)=PZ(i)-A(i)*DENS(i)*COS(THETA(i))
        END IF
END DO

!DEFINE NODAL FORCES
DO i=1,NUMNODES
        READ(1,*)(F(i,j),j=1,3)
END DO

!DEFINE BOUNDARY CONDITIONS
!FIRST READ # OF BCs
READ(1,*)NUMBC
!DEFINE BC SUCH THAT BC(*,1) IS THE NODE AND BC(*,2) IS THE RESTRAINED DoF
DO i=1,NUMBC
        READ(1,*)(BC(i,j),j=1,2)
END DO

!ECHO INPUT
3 FORMAT(I3)
4 FORMAT(I3, 3X, F20.10, 3X, F20.10)
5 FORMAT(I3, 3X, F20.10, 3X, F20.10, 3X, F20.10)
6 FORMAT(I3, 15X, I3)
WRITE(2,*)'~~~~~~~ECHO INPUT~~~~~~~'
WRITE(2,*)
IF(g.EQ.1) THEN
        WRITE(2,*)'Gravity is ON!'
ELSE
        WRITE(2,*)'Gravity is OFF!'
END IF
WRITE(2,*)
WRITE(2,*)'NUMBER OF ELEMENTS:'
WRITE(2,3)NUMEL
WRITE(2,*)
WRITE(2,*)
WRITE(2,*)'NUMBER OF NODES:'
WRITE(2,3)NUMNODES
WRITE(2,*)
WRITE(2,*)
WRITE(2,*)"Nodes and location in global X-Z plane"
WRITE(2,*)'NODE NO.     X-Coordinate           Z-Coordinate'
DO i=1,NUMNODES
        WRITE(2,4)i,X(i),Z(i)
END DO
WRITE(2,*)
WRITE(2,*)
WRITE(2,*)'Element Connectivity Table:'
WRITE(2,*)
WRITE(2,*)'              ','LOCAL NODE','   ','GLOBAL NODE'
DO i=1,NUMEL
        WRITE(2,*)'Element',i,':'
        DO J=1,2
                WRITE(2,*)'Connectivity',j,ELCONNECT(i,j)
        END DO
        WRITE(2,*)
END DO
WRITE(2,*)
WRITE(2,*)
WRITE(2,*)'Calculated Length and Angle of Elements:'
DO i=1,NUMEL
        WRITE(2,*)'Length of element',i,':',LENGTH(i)
        WRITE(2,*)'Angle of element',i,':',THETA(i)
        WRITE(2,*)
END DO
WRITE(2,*)
WRITE(2,*)'Material Properties:'
WRITE(2,*)"Element      Modulus                    Area            Moment of inertia"
DO i=1,NUMEL
        WRITE(2,5)i,E(i),A(i),IYY(i)
END DO
WRITE(2,*)
WRITE(2,*)
WRITE(2,*)'Distributed Loads:'
WRITE(2,*)'Element      1-Direction           2-Direction'
DO i=1,NUMEL
  WRITE(2,4)i,PX(i),PZ(i)
END DO
WRITE(2,*)
WRITE(2,*)
WRITE(2,*)'Applied forces at nodes:'
WRITE(2,*)'Node No.     1-Direction            2-Direction            3-Direction'
DO I=1,NUMNODES
  WRITE(2,5)I,(F(I,J),J=1,3)
END DO
WRITE(2,*)
WRITE(2,*)
WRITE(2,*)'Number of boundary conditions:'
WRITE(2,3)NUMBC
WRITE(2,*)
WRITE(2,*)
WRITE(2,*)'Boundary conditions:'
WRITE(2,*)'Node No.         Direction'
DO I=1,NUMBC
  WRITE(2,6)(BC(I,J),J=1,2)
END DO  
WRITE(2,*)
WRITE(2,*)
WRITE(2,*)'<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>'
WRITE(2,*)'<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>'
WRITE(2,*)
WRITE(2,*)
END

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!SUBROUTINE: INITIALIZE

SUBROUTINE INITIALIZE(kb,STIFF,NUMEL,NUMNODES,FORCE)
IMPLICIT NONE
REAL::kb(6,6),STIFF(303,303),FORCE(303)
INTEGER::i,j,NUMEL,NUMNODES

!INITIALIZE LOCAL K
DO i=1,6
        DO j=1,6
                kb(i,j)=0
        END DO
END DO

!INITIALIZE STIFFNESS MATRIX
DO i=1,303       
        DO j=1,303
                STIFF(i,j)=0
        END DO
END DO

!INITIALIZE FORCE MATRIX
DO I=1,303
  FORCE(i)=0
END DO

!WRITE(3,*)kb
!WRITE(3,*)STIFF
!WRITE(3,*)FORCE
END

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!SUBROUTINE: ELEMENT STIFFNESS

SUBROUTINE ELESTIFF(i,E,A,IYY,LENGTH,kb)
IMPLICIT NONE
REAL::kb(6,6),E(100),A(100),LENGTH(100),IYY(100)
INTEGER::i,j,k

!BUILD kb MATRIX
DO j=1,6
        DO k=1,6
                kb(j,k)=0
        END DO
END DO
kb(1,1)=E(i)*A(i)/LENGTH(i)
kb(1,4)=-E(i)*A(i)/LENGTH(i)
kb(2,2)=12*E(i)*IYY(i)/LENGTH(i)**3
kb(2,3)=6*E(i)*IYY(i)/LENGTH(i)**2
kb(2,5)=-12*E(i)*IYY(i)/LENGTH(i)**3
kb(2,6)=6*E(i)*IYY(i)/LENGTH(i)**2
kb(3,2)=6*E(i)*IYY(i)/LENGTH(i)**2
kb(3,3)=4*E(i)*IYY(i)/LENGTH(i)
kb(3,5)=-6*E(i)*IYY(i)/LENGTH(i)**2
kb(3,6)=2*E(i)*IYY(i)/LENGTH(i)
kb(4,1)=-E(i)*A(i)/LENGTH(i)
kb(4,4)=E(i)*A(i)/LENGTH(i)
kb(5,2)=-12*E(i)*IYY(i)/LENGTH(i)**3
kb(5,3)=-6*E(i)*IYY(i)/LENGTH(i)**2
kb(5,5)=12*E(i)*IYY(i)/LENGTH(i)**3
kb(5,6)=-6*E(i)*IYY(i)/LENGTH(i)**2
kb(6,2)=6*E(i)*IYY(i)/LENGTH(i)**2
kb(6,3)=2*E(i)*IYY(i)/LENGTH(i)
kb(6,5)=-6*E(i)*IYY(i)/LENGTH(i)**2
kb(6,6)=4*E(i)*IYY(i)/LENGTH(i)


!WRITE(3,*)kb
END


!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!SUBROUTINE:TRANSFORM K

SUBROUTINE TRANSFORMK(i,kb,GLOBK,THETA,TRANS_TRANS,TRANS)
IMPLICIT NONE
REAL::kb(6,6),GLOBK(6,6),THETA(100),TRANS(6,6),TRANS_TRANS(6,6),INTER(6,6)
INTEGER::i,j,k

!BUILD TRANSFORMATION MATRIX
DO j=1,6
        DO k=1,6
                TRANS(j,k)=0
        END DO
END DO
TRANS(1,1)=COS(THETA(i))
TRANS(1,2)=SIN(THETA(i))
TRANS(2,1)=-SIN(THETA(i))
TRANS(2,2)=COS(THETA(i))
TRANS(3,3)=1
TRANS(4,4)=COS(THETA(i))
TRANS(4,5)=SIN(THETA(i))
TRANS(5,4)=-SIN(THETA(i))
TRANS(5,5)=COS(THETA(i))
TRANS(6,6)=1

!BUILD TRANSPOSE OF TRANSFORMATION MATRIX
TRANS_TRANS=TRANSPOSE(TRANS)

!TRANSFORM kb TO GLOBAL COORDs
INTER=MATMUL(TRANS_TRANS,kb)
GLOBK=MATMUL(INTER,TRANS)

!WRITE(3,*)GLOBK

END

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!SUBROUTINE: ELEMENT FORCE, F

SUBROUTINE ELEFORCE(i,LENGTH,PX,PZ,FDIST)
IMPLICIT NONE
REAL::LENGTH(100),PX(100),PZ(100),FDIST(6)
INTEGER::i,j

FDIST(1)=PX(i)*LENGTH(i)/2
FDIST(2)=PZ(i)*LENGTH(i)/2
FDIST(3)=PZ(i)*LENGTH(i)**2/12
FDIST(4)=PX(i)*LENGTH(i)/2
FDIST(5)=PZ(i)*LENGTH(i)/2
FDIST(6)=-PZ(i)*LENGTH(i)**2/12  

!WRITE(3,*)FDIST

END

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!SUBROUTINE: TRANSFORM Q

SUBROUTINE TRANSFORMQ(i,TRANS_TRANS,FDIST,GLOBFDIST)
IMPLICIT NONE
REAL::TRANS_TRANS(6,6),FDIST(6),GLOBFDIST(6)
INTEGER::i

GLOBFDIST=MATMUL(TRANS_TRANS,FDIST)

!WRITE(3,*)GLOBFDIST

END

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!SUBROUTINE: ASSEMBLE FINAL GLOBAL K AND Q MATRICES

SUBROUTINE ASSEMBLEKQ(i,ELCONNECT,STIFF,STIFF_REACT,GLOBK,FORCE,F,GLOBFDIST)
IMPLICIT NONE
REAL::STIFF(303,303),STIFF_REACT(303,303),GLOBK(6,6),FORCE(303),F(101,3)
REAL::GLOBFDIST(6)
INTEGER::i,j,k,ELCONNECT(100,2),ND1,ND2

ND1=ELCONNECT(i,1)
ND2=ELCONNECT(i,2)

!MOVE LOCAL K MATRIX FIELDS TO GLOBAL K MATRIX, 4 NESTED DO LOOPS PER ELEMENT
DO j=1,3
        DO k=1,3
                STIFF(j+(ND1-1)*3,k+(ND1-1)*3)=STIFF(j+(ND1-1)*3,k+(ND1-1)*3)+GLOBK(j,k)
                IF(STIFF(j+(ND1-1)*3,k+(ND1-1)*3) < 0.0001 .and. &
                        &STIFF(j+(ND1-1)*3,k+(ND1-1)*3) > -0.0001) &
                        &STIFF(j+(ND1-1)*3,k+(ND1-1)*3) = 0
        END DO

END DO

DO j=1,3
        DO k=1,3
                STIFF(j+(ND1-1)*3,k+(ND2-1)*3)=STIFF(j+(ND1-1)*3,k+(ND2-1)*3)+GLOBK(j,k+3)
                IF(STIFF(j+(ND1-1)*3,k+(ND2-1)*3) < 0.0001 .and. &
                        &STIFF(j+(ND1-1)*3,k+(ND2-1)*3) > -0.0001) &
                        &STIFF(j+(ND1-1)*3,k+(ND2-1)*3) = 0
        END DO
END DO

DO j=1,3
        DO k=1,3
                STIFF(j+(ND2-1)*3,k+(ND1-1)*3)=STIFF(j+(ND2-1)*3,k+(ND1-1)*3)+GLOBK(j+3,k)
                IF(STIFF(j+(ND2-1)*3,k+(ND1-1)*3) < 0.0001 .and. &
                        &STIFF(j+(ND2-1)*3,k+(ND1-1)*3) > -0.0001) &
                        &STIFF(j+(ND2-1)*3,k+(ND1-1)*3) = 0
        END DO
END DO

DO j=1,3
        DO k=1,3
                STIFF(j+(ND2-1)*3,k+(ND2-1)*3)=STIFF(j+(ND2-1)*3,&
                        &k+(ND2-1)*3)+GLOBK(j+3,k+3)
                IF(STIFF(j+(ND2-1)*3,k+(ND2-1)*3) < 0.0001 .and. &
                        &STIFF(j+(ND2-1)*3,k+(ND2-1)*3) > -0.0001) &
                        &STIFF(j+(ND2-1)*3,k+(ND2-1)*3) = 0
        END DO
END DO
!END ELEMENT TO GLOBAL MOVE

!BUILD NEW STIFF MATRIX FOR LATER USE
DO j = 1, 303
        DO k = 1, 303
                STIFF_REACT(j,k) = STIFF(j,k)
        END DO
END DO

!INCORPORATE DISTRIBUTED LOAD TO TOTAL GLOBAL FORCE
!POSSIBLE FORCE CORRECTION
FORCE(ELCONNECT(i,1)*3-2)=GLOBFDiST(1)+FORCE(ELCONNECT(i,1)*3-2)
FORCE(ELCONNECT(i,1)*3-1)=GLOBFDiST(2)+FORCE(ELCONNECT(i,1)*3-1)
FORCE(ELCONNECT(i,1)*3)=GLOBFDiST(3)+FORCE(ELCONNECT(i,1)*3)
FORCE(ELCONNECT(i,2)*3-2)=GLOBFDiST(4)+FORCE(ELCONNECT(i,2)*3-2)
FORCE(ELCONNECT(i,2)*3-1)=GLOBFDiST(5)+FORCE(ELCONNECT(i,2)*3-1)
FORCE(ELCONNECT(i,2)*3)=GLOBFDiST(6)+FORCE(ELCONNECT(i,2)*3)
FORCE(ELCONNECT(i,1)*3-2)=F(ND1,1)+FORCE(ELCONNECT(i,1)*3-2)
FORCE(ELCONNECT(i,1)*3-1)=F(ND1,2)+FORCE(ELCONNECT(i,1)*3-1)
FORCE(ELCONNECT(i,1)*3-0)=F(ND1,3)+FORCE(ELCONNECT(i,1)*3-0)
FORCE(ELCONNECT(i,2)*3-2)=F(ND2,1)+FORCE(ELCONNECT(i,2)*3-2)
FORCE(ELCONNECT(i,2)*3-1)=F(ND2,2)+FORCE(ELCONNECT(i,2)*3-1)
FORCE(ELCONNECT(i,2)*3-0)=F(ND2,3)+FORCE(ELCONNECT(i,2)*3-0) 
 

!WRITE(3,*)STIFF
!WRITE(3,*)FORCE

END

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!SUBROUTINE: APPLY THE BOUNDARY CONDITIONS TO K AND Q MATRICES

SUBROUTINE APPLYBC(NUMBC,BC,STIFF,FORCE,NUMNODES)
IMPLICIT NONE
REAL::STIFF(303,303),FORCE(303)
INTEGER::i,j,NUMBC,BC(303,2),NUMNODES,IDOF,NUMQ,BOUNDDOF

NUMQ=NUMNODES*3

!DOUBLE CHECK THIS
DO i=1,NUMBC
        IDOF=(BC(i,1)-1)*3+BC(i,2)
        DO j=1,NUMQ
                STIFF(IDOF,j)=0.0
                STIFF(j,IDOF)=0.0
        END DO
        STIFF(IDOF,IDOF)=1.0
        FORCE(IDOF)=0
END DO

!WRITE(3,*)STIFF
!WRITE(3,*)FORCE

END

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!SUBROUTINE: SOLVE SYSTEM VIA GAUSSIAN ELIMINATION, FIND DISPLACEMENTS, Q

SUBROUTINE SOLVE(STIFF,FORCE,NUMEL,Q,NUMNODES)
IMPLICIT NONE
REAL::STIFF(303,303),FORCE(303),Q(303),RATIO,SUM
INTEGER::NUMEL,NUMNODES,j,i,n,nm1,jp1,k,ks

!USING THE CODE FOR GAUSSIAN ELIMINATION PROVIDED BY DR. HAISLER

!FORWARD ELIMINATION
n=NUMNODES*3
nm1=n-1

DO j=1,nm1
        jp1=j+1
        DO i=jp1,n
                RATIO=STIFF(i,j)/STIFF(j,j)
                IF (RATIO.EQ.0.0) go to 10
                DO k=j,n
                        STIFF(i,k)=STIFF(i,k)-RATIO*STIFF(j,k)
                END DO
                FORCE(i)=FORCE(i)-RATIO*FORCE(j)
                10 continue
        END DO
END DO

!BACK SUBSTITUION
Q(n)=FORCE(n)/STIFF(n,n)

DO j=1,nm1        
        i=n-j        
        SUM=FORCE(i)
        ks=i+1
        DO k=ks,n                
                SUM=SUM-STIFF(i,k)*Q(k)                
        END DO
        Q(i)=SUM/STIFF(i,i)
        END DO
!WRITE(3,*)Q

RETURN

END

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!SUBROUTINE: CALCULATE ELEMENT DISPLACEMENTS AND REACTIONS

SUBROUTINE REACTIONS(i,NUMNODES,STIFF_REACT,Q,GLOBREACT,TRANS,kb,QB&
                &,FDIST,FORCE_INT,ELCONNECT,QSUB,NUMEL)
IMPLICIT NONE
REAL::STIFF_REACT(303,303),Q(303),GLOBREACT(303),TRANS(6,6),kb(6,6),QB(6)
REAL::FDIST(6),FORCE_INT(6),QSUB(6),kt(6,6)
INTEGER::i,j,NUMNODES,ELCONNECT(100,2),NUMEL,ND1,ND2

ND1=ELCONNECT(i,1)
ND2=ELCONNECT(i,2)

GLOBREACT=MATMUL(STIFF_REACT,Q)

!FIND DISPLACEMENT OF INDIVIDUAL ELEMENTS FROM GLOBAL DISCPLACEMENTS
DO j=1,3
        QSUB(j)=Q((ND1-1)*3+j)
        QSUB(j+3)=Q((ND2-1)*3+j)
END DO

!ZERO OUT SMALL VALUES FOR EASIER ANALYSIS, NOT TRUE ZEROES DUE TO ERROR
DO j=1,6
        IF(FDIST(j)<0.0001.AND.FDIST(j)>-0.0001) FDIST(j)=0
END DO

!CALCULATE INTERNAL FORCES, ELEMENT REACTIONS, S
kt=MATMUL(kb,TRANS)
QB=MATMUL(kt,QSUB)

DO j=1,6
        FORCE_INT(j)=QB(j)-FDIST(j)
END DO

!WRITE(3,*)QSUB
!WRITE(3,*)FORCE_INT
END

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!SUBROUTINE: WRITE FINAL OUTPUT TO OUTPUT FILE

SUBROUTINE ELEOUTPUT(i,QSUB,FORCE_INT)
IMPLICIT NONE
REAL::QSUB(6),FORCE_INT(6)
INTEGER::i,j

!DECLARE FORMATTING
8 FORMAT(I1,20X,F20.10)

!WRITE FINAL OUTPUT
IF (i.EQ.1) THEN
        WRITE(2,*)'Displacement and internal force 1-2-3 and 4-5-6 correspond'
        WRITE(2,*)'to x-y-z of node 1 and x-y-z of node 2 respectively'
END IF
WRITE(2,*)
WRITE(2,*)'Displacements on Element:',i
WRITE(2,*)'Displacement No.            Displacement'

DO j=1,6

        WRITE(2,8)j,QSUB(j)

END DO

WRITE(2,*)
WRITE(2,*)'Internal Forces on Element:',i
WRITE(2,*)'Force No.                  Force'

DO j=1,6

        WRITE(2,8)j,FORCE_INT(j)

END DO

WRITE(2,*)

END

!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
!<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

!SUBROUTINE: OUTPUT TOTAL DISPLACEMENTS AND REACTIONS

SUBROUTINE GLOBOUTPUT(Q,FORCE,GLOBREACT)
REAL::Q(303),FORCE(303),FORCE_REACT(303),GLOBREACT(303),FORCE_NEW(303)
INTEGER::i

WRITE(2,*)'Global Displacements'
WRITE(2,*)
DO i=1,303
        IF(Q(i)==0)CYCLE
        WRITE(2,102)'Displacement',i,'is',Q(i)
END DO
WRITE(2,*)
WRITE(2,*)'All other displacements are zero'
WRITE(2,*)
WRITE(2,*)
WRITE(2,*)'Reactions at Boundary Conditions'
WRITE(2,*)
DO i=1,303
        IF(FORCE(i).EQ.0.0) THEN
                FORCE_NEW(i)=1
        ELSE
                FORCE_NEW(i)=0
        END IF
        FORCE_REACT(i)=FORCE_NEW(i)*GLOBREACT(i)
        IF(FORCE_REACT(i)<.001.AND.FORCE_REACT(i)>-.001)CYCLE
        WRITE(2,101)'Reaction',i,'is',FORCE_REACT(i)
END DO

100 FORMAT(6(1x,f10.2))
101 FORMAT(1x,a8,1x,i3,1x,a2,2x,f10.3)
102 FORMAT(1x,a12,1x,i3,1x,a11,2x,f10.7)
103 FORMAT(12(1x,f2.0))

END
