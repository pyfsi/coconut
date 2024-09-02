#define RAMP |ramp|

C==============================================================================
C======================== GLOBAL DATA =========================================
C==============================================================================
      BLOCK DATA

      INTEGER D,N,S
      PARAMETER (D = |dimension|)
      PARAMETER (N = |arraySize|)
      PARAMETER (S = |surfaces|)

      CHARACTER (LEN=80), DIMENSION(S) :: SURFACEIDS
      COMMON /SURF/ SURFACEIDS
      SAVE /SURF/
      DATA SURFACEIDS /|surfaceIDs|/

      DOUBLE PRECISION LOADNEW(D+1,N,S)
#if RAMP
      DOUBLE PRECISION LOADOLD(D+1,N,S)
#endif
      INTEGER FILLED,
     &        M(S)
      COMMON /TABLE/ LOADNEW,
#if RAMP
     &               LOADOLD,
#endif
     &               FILLED,
     &               M
      SAVE /TABLE/

      INTEGER OP(2,S),ELS(2,N,S)
      COMMON /OPELS/ OP,ELS
      SAVE /OPELS/

      DATA FILLED/-1/

      END

C==============================================================================
C======================== LOOKUP SUBROUTINE ===================================
C==============================================================================
      SUBROUTINE LOOKUP(NOEL,R,K)

      IMPLICIT NONE

      INTEGER N,S
      PARAMETER (N = |arraySize|)
      PARAMETER (S = |surfaces|)

      INTEGER OP(2,S),ELS(2,N,S)
      COMMON /OPELS/ OP,ELS
      SAVE /OPELS/

      INTEGER NOEL,R,K,LEFT,MID,RIGHT

      K = 0
      LEFT  = 1
      RIGHT = OP(1,R)
      DO
         IF (LEFT > RIGHT) RETURN
         MID = (LEFT+RIGHT)/2
         IF (NOEL == ELS(1,MID,R)) THEN
            K = MID
            RETURN
         ELSE IF (NOEL < ELS(1,MID,R)) THEN
            RIGHT = MID-1
         ELSE
            LEFT = MID+1
         END IF
      END DO

      RETURN
      END
     
C==============================================================================
C======================== READDATA SUBROUTINE =================================
C==============================================================================
      SUBROUTINE READDATA(KSTEP)
      
      IMPLICIT NONE
     
      INTEGER D,N,S
      PARAMETER (D = |dimension|)
      PARAMETER (N = |arraySize|)
      PARAMETER (S = |surfaces|)

      DOUBLE PRECISION LOADNEW(D+1,N,S)
#if RAMP
      DOUBLE PRECISION LOADOLD(D+1,N,S)
#endif
      INTEGER FILLED,
     &        M(S)
      COMMON /TABLE/ LOADNEW,
#if RAMP
     &               LOADOLD,
#endif
     &               FILLED,
     &               M
      SAVE /TABLE/

      INTEGER ID
      INTEGER OP(2,S),ELS(2,N,S)
      COMMON /OPELS/ OP,ELS
      SAVE /OPELS/

      CHARACTER(LEN=40) :: FMT_ELEM
      INTEGER UNIT_ELEM

      CHARACTER(LEN=40) :: FMT_LOAD
      CHARACTER(LEN=200) :: FILENAME
      INTEGER I,R,UNIT_LOAD,IOS,KSTEP

!$OMP CRITICAL
      IF (FILLED < 0) THEN

      FMT_LOAD = '(ES27.17E2,|dimension|ES27.17E2)'
      UNIT_LOAD = 100

      ID = 0

      FMT_ELEM = '(I10,BN,I11)'
      UNIT_ELEM = 101

      DO R = 1,S
         WRITE(FILENAME,'(A,A,A,I0,A)') 
     &      '|PWD|',
     &      '/|CSM_dir|/CSM_Time0',
     &      'Surface',(R-1),'Elements.dat'

         OPEN(UNIT=UNIT_ELEM,FILE=FILENAME,STATUS='OLD')
      
         READ(UNIT_ELEM,FMT_ELEM,IOSTAT=IOS) OP(:,R)
         IF (IOS < 0) THEN
            CALL STDB_ABQERR(-3,'USR-error: problem while opening')
         END IF
         IF (OP(1,R) > N) THEN
            CALL STDB_ABQERR(-3,'USR-error: problem with array length')
         END IF
         DO I = 1,OP(1,R)
            READ(UNIT_ELEM,FMT_ELEM,IOSTAT=IOS) ELS(:,I,R)
            IF (IOS < 0) THEN
               CALL STDB_ABQERR(-3,'USR-error: problem while reading')
            END IF
         END DO
         CLOSE(UNIT_ELEM)
      END DO
      
      DO R = 1,S
         WRITE(FILENAME,'(A,A,I0,A,I0,A,I0,A)') 
     &      '|PWD|',
     &      '/|CSM_dir|/CSM_Time',
     &      KSTEP,'Surface',(R-1),'Cpu',ID,'Input.dat'

         OPEN(UNIT=UNIT_LOAD,FILE=FILENAME,STATUS='OLD')
      
         READ(UNIT_LOAD,'(I)',IOSTAT=IOS) M(R)
         IF (IOS < 0) THEN
            CALL STDB_ABQERR(-3,'USR-error: problem while opening')
         END IF
         IF (M(R) > N) THEN
            CALL STDB_ABQERR(-3,'USR-error: problem with array length')
         END IF
         DO I = 1,M(R)
            READ(UNIT_LOAD,FMT_LOAD,IOSTAT=IOS)
     &         LOADNEW(:,I,R)
            IF (IOS < 0) THEN
               CALL STDB_ABQERR(-3,'USR-error: problem while reading')
            END IF
         END DO
         CLOSE(UNIT_LOAD)

#if RAMP
         WRITE(FILENAME,'(A,A,I0,A,I0,A,I0,A)') 
     &      '|PWD|',
     &      '/|CSM_dir|/CSM_Time',
     &      (KSTEP-1),'Surface',(R-1),'Cpu',ID,'Input.dat'

         OPEN(UNIT=UNIT_LOAD,FILE=FILENAME,STATUS='OLD')
      
         READ(UNIT_LOAD,'(I)',IOSTAT=IOS) M(R)
         IF (IOS < 0) THEN
            CALL STDB_ABQERR(-3,'USR-error: problem while opening')
         END IF
         IF (M(R) > N) THEN
            CALL STDB_ABQERR(-3,'USR-error: problem with array length')
         END IF
         DO I = 1,M(R)
            READ(UNIT_LOAD,FMT_LOAD,IOSTAT=IOS)
     &         LOADOLD(:,I,R)
            IF (IOS < 0) THEN
               CALL STDB_ABQERR(-3,'USR-error: problem while reading')
            END IF
         END DO
         CLOSE(UNIT_LOAD)
#endif
      END DO

      FILLED = 1

      END IF
!$OMP END CRITICAL

      RETURN
      END

C==============================================================================
C======================== DLOAD SUBROUTINE ====================================
C==============================================================================

      SUBROUTINE DLOAD(F,KSTEP,KINC,TIME,NOEL,NPT,LAYER,KSPT,
     &   COORDS,JLTYP,SNAME)

      IMPLICIT NONE

      INTEGER D,N,S
      PARAMETER (D = |dimension|)
      PARAMETER (N = |arraySize|)
      PARAMETER (S = |surfaces|)
      CHARACTER (LEN=80), DIMENSION(S) :: SURFACEIDS
      CHARACTER (LEN=89) :: PREPENDED
      COMMON /SURF/ SURFACEIDS
      SAVE /SURF/

      DOUBLE PRECISION LOADNEW(D+1,N,S)
#if RAMP
      DOUBLE PRECISION LOADOLD(D+1,N,S),DT
#endif
      INTEGER FILLED,
     &        M(S)
      COMMON /TABLE/ LOADNEW,
#if RAMP
     &               LOADOLD,
#endif
     &               FILLED,
     &               M
      SAVE /TABLE/

      DOUBLE PRECISION F,TIME(2),COORDS(D)
      CHARACTER(LEN=80) :: SNAME
      INTEGER KSTEP,KINC,NOEL,NPT,LAYER,KSPT,JLTYP,R
      LOGICAL :: FOUND

      INTEGER K, AXIS
      INTEGER OP(2,S),ELS(2,N,S)
      COMMON /OPELS/ OP,ELS
      SAVE /OPELS/
      IF (FILLED < 0) THEN
         CALL READDATA(KSTEP)
      END IF

#if RAMP
      IF (|deltaT| > 0.0) THEN
         DT = |deltaT|
      ELSE
         DT = 1.0
      END IF
#endif

      FOUND  = .FALSE.
      IF (S > 1) THEN
         DO R = 1,S
            PREPENDED = 'ASSEMBLY_' // SURFACEIDS(R)
            IF (SNAME == PREPENDED) THEN
               FOUND = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT. FOUND) THEN
            PRINT *, 'USR-abort: no matching Modelpart found for surfac
     &e.', SNAME
            CALL FLUSH(6)
            CALL STDB_ABQERR(-3,'USR-abort: no matching ModelPart found
     &for surface.', SNAME)
         END IF
      ELSE
         R = 1
      END IF
      CALL LOOKUP(NOEL,R,K)
      K = ELS(2,K,R)+NPT

      IF (JLTYP == 0) THEN
         AXIS = 1
      ELSE
         AXIS = JLTYP - 39
      END IF
#if RAMP
      F = LOADNEW(AXIS,K,R)*(TIME(1)/DT)
     &   +LOADOLD(AXIS,K,R)*(1.0-TIME(1)/DT)
#else
      F = LOADNEW(AXIS,K,R)
#endif

      RETURN
      END
      
C==============================================================================
C======================== UTRACLOAD SUBROUTINE ================================
C==============================================================================

      SUBROUTINE UTRACLOAD(ALPHA,T_USER,KSTEP,KINC,TIME,
     &   NOEL,NPT,COORDS,DIRCOS,JLTYP,SNAME)

      IMPLICIT NONE

      INTEGER D,N,S
      PARAMETER (D = |dimension|)
      PARAMETER (N = |arraySize|)
      PARAMETER (S = |surfaces|)
      CHARACTER (LEN=80), DIMENSION(S) :: SURFACEIDS
      CHARACTER (LEN=89) :: PREPENDED
      COMMON /SURF/ SURFACEIDS
      SAVE /SURF/

      DOUBLE PRECISION LOADNEW(D+1,N,S)
#if RAMP
      DOUBLE PRECISION LOADOLD(D+1,N,S)
#endif
      INTEGER FILLED,
     &        M(S)
      COMMON /TABLE/ LOADNEW,
#if RAMP
     &               LOADOLD,
#endif
     &               FILLED,
     &               M
      SAVE /TABLE/

      DOUBLE PRECISION ALPHA,T_USER(D),TIME(2),COORDS(D),DIRCOS(3,3)
      CHARACTER(LEN=80) :: SNAME
      INTEGER KSTEP,KINC,NOEL,NPT,JLTYP,R
      LOGICAL :: FOUND

      INTEGER L
      INTEGER OP(2,S),ELS(2,N,S)
      COMMON /OPELS/ OP,ELS
      SAVE /OPELS/

      IF (FILLED < 0) THEN
         CALL READDATA(KSTEP)
      END IF

      FOUND  = .FALSE.
      IF (S > 1) THEN
         DO R = 1,S
            PREPENDED = 'ASSEMBLY_' // SURFACEIDS(R)
            IF (SNAME == PREPENDED) THEN
               FOUND = .TRUE.
               EXIT
            END IF
         END DO
         IF (.NOT. FOUND) THEN
            PRINT *, 'USR-abort: no matching Modelpart found for surfac
     &e.', SNAME
            CALL FLUSH(6)
            CALL STDB_ABQERR(-3,'USR-abort: no matching ModelPart found
     &for surface.', SNAME)
         END IF
      ELSE
         R = 1
      END IF

      CALL LOOKUP(NOEL,R,L)
      L = ELS(2,L,R)+NPT

      T_USER = LOADNEW(2:D+1,L,R)

      ALPHA = SQRT(SUM(T_USER**2))
      IF (ALPHA /= 0.0) THEN
         T_USER = T_USER/ALPHA
      ELSE
         T_USER(1) = 1.0
      END IF
      
      RETURN
      END
