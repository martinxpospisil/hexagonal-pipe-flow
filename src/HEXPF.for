	PROGRAM HEXPF

	IMPLICIT NONE
C
C	Variable declaration
C
	INTEGER HexRings, MaxHexRings, MinHexRings, MaxDim, Dim, NCells 
	INTEGER Iter, IterMax, i, j, k, ii
C
C     Calculation mesh definiton
C     - HexRings = no. of hexafonal layers around the centre = k/2
C     - Dim      = true size of the mesh (no. of nodes)
C
	PARAMETER (MaxDim = 101)
	REAL(8) VelZ(MaxDim,MaxDim), VelZ2(MaxDim,MaxDim)
	REAL(8) Omega, Sum, Sum2, StepSize
	REAL(8) Eps, XPos, YPos, PosShift, Scaling
	LOGICAL IsSite(MaxDim,MaxDim), IsBoundary(MaxDim,MaxDim)
C
C     Prepare external files
C
    OPEN(2, FILE = 'vysledky_relaxace.dat', STATUS = 'REPLACE')
	OPEN(1, FILE = 'velocity_profile.dat', STATUS = 'REPLACE')
C
C     Parameters input
C
	WRITE(*,*) 'Zadejte jemonost deleni k/2 (cele cislo <50):'
	READ(*,*) HexRings
	Dim = 2*HexRings+1
	WRITE(*,*) 'Zadejte hodnotu relaxacniho parametru omega:'
	READ(*,*) Omega
C
C	Another parameters
C   - IterMax  = max. number of iterations (halting criterion 1)
C   - Eps      = max. error (halting criterion 2)
C   - StepSize = size of the calculation mesh
C
	IterMax  = 10000
	Eps      = 1.0D-8
	StepSize = (2.0D0/(3.0D0*3.0D0**(0.5D0)))**(0.5D0)/HexRings
C
C
C	Calculation mesh intialization
C
    DO i = 1, MaxDim
	  DO j = 1, MaxDim
	      IsSite(i,j) = .FALSE.
	  END DO
	END DO
C
C     Mark nodes within calculation mesh
C     - T = within, F = without
C
	DO i = 1, Dim
	  DO j = 1, Dim
	    IF (((i+j) .GT. (2*Dim-HexRings)) .OR.
     *        ((i+j) .LT. (2+HexRings))) THEN
	      IsSite(i,j) = .FALSE.
	    ELSE
	      IsSite(i,j) = .TRUE.
	    END IF
	  END DO
	END DO
C
C     Vypocet poctu bunek
	NCells = 0
      DO i = 1, MaxDim
	  DO j = 1, MaxDim
	    IF (IsSite(i,j)) NCells = NCells+1
	  END DO
	END DO
C
C     Mark nodes with boundary condition
C
      DO i = 1, Dim
	  DO j = 1, Dim
          IsBoundary(i,j) = .FALSE.
	  END DO
	END DO

	DO i = 1, Dim
	  DO j = 1, Dim
          IF (IsSite(i,j) .AND.
     *       ((i .EQ. 1) .OR. (i .EQ. Dim) .OR.
     *          (j .EQ. 1) .OR. (j .EQ. Dim) .OR.
     *          ((i+j) .EQ. (2*Dim-HexRings)) .OR.
     *          ((i+j) .EQ. (2+HexRings)))) THEN
              IsBoundary(i,j) = .TRUE.	          
		END IF
	  END DO
	END DO

C
C	Initialization of the velocity profile
C
      DO i = 1, MaxDim
	  DO j = 1, MaxDim
	    VelZ (i,j) = 0.0D0
	    VelZ2(i,j) = 0.0D0
        END DO
	END DO
	Sum  = 0.0D0
	Sum2 = 0.0D0
C
C     Calculation
C
	DO Iter = 1, IterMax
	  Sum2 = 0.0D0
	  DO i = 1, Dim
	    DO j = 1, Dim
	    VelZ2(i,j) = 0.0D0
	      IF (IsSite(i,j) .AND. (.NOT. IsBoundary(i,j))) THEN
              VelZ2(i,j) = VelZ(i-1,j)   + VelZ(i+1,j)
     *                   + VelZ(i,j-1)   + VelZ(i,j+1)
     *                   + VelZ(i-1,j+1) + VelZ(i+1,j-1)
	        VelZ2(i,j) = VelZ2(i,j)/6.0D0 + StepSize**2/4.0D0
	        VelZ2(i,j) = Omega*VelZ2(i,j) + (1.0D0-Omega)*VelZ(i,j)
		  Sum2 = Sum2 + VelZ2(i,j)
	      END IF
	    END DO
        END DO
	  Sum2 = 3.0D0**0.5D0/2.0D0*StepSize**2*Sum2

	  VelZ = VelZ2
	  IF (ABS(Sum-Sum2) .LT. Eps) EXIT
	  Sum = Sum2
	END DO
C
C	Print-out to screen
C 
      WRITE (*,*) 
     *' k  ',
     *'     I       ',      
     *'Pocet bunek  ',
     *'     h       ',
     *'pocet iteraci'
	WRITE (*,10000) 2*HexRings,Sum2,NCells,StepSize,Iter-1
10000	FORMAT(I3,2X,E12.6E2,2X,I5,6X,E12.6E2,2X,I5)
C
C     Print-out to external files
C
C     Output for plotting the density profile in gnuplot
C
	Scaling = 1.0D0
	PosShift = 0.0D0
	DO j = 1, Dim
	  DO i = 1, Dim
	    XPos = i*Scaling + PosShift
	    YPos = j*Scaling*3.0D0**0.5D0/2.0D0
	    IF (IsSite(i,j)) WRITE(1,10003) YPos, XPos, VelZ(i,j)
	  END DO
	  WRITE(1,*)
	  PosShift = PosShift + Scaling/2.0D0
	END DO
	CLOSE(1)
10003 FORMAT(3(E20.5E2))

	END PROGRAM