C==============================================================================
C======================== GLOBAL DATA =========================================
C==============================================================================
      BLOCK DATA

      INTEGER S
      PARAMETER (S = |surfaces|)
      INTEGER NOEL_PREV(S)
      COMMON /PREV/ NOEL_PREV
      SAVE /PREV/
      DATA NOEL_PREV /S*0/

      CHARACTER (LEN=80), DIMENSION(S) :: SURFACEIDS
      COMMON /SURF/ SURFACEIDS
      SAVE /SURF/
      DATA SURFACEIDS /|surfaceIDs|/

#ifdef MPI
      INTEGER ID,IDENTIFIED
      COMMON /IDENT/ ID,IDENTIFIED
      SAVE /IDENT/

      DATA IDENTIFIED/-1/
#endif

      END

C==============================================================================
C======================== IDENTIFY SUBROUTINE =================================
C==============================================================================
#ifdef MPI
      SUBROUTINE IDENTIFY
      
      IMPLICIT NONE
      
      INCLUDE 'mpif.h'

      INTEGER ID,IDENTIFIED,IERROR
      COMMON /IDENT/ ID,IDENTIFIED
      SAVE /IDENT/
      
      IF (|cpus| > 1) THEN
         CALL MPI_COMM_RANK(MPI_COMM_WORLD,ID,IERROR)
         IF (IERROR < 0) THEN
            CALL STDB_ABQERR(-3,'USR-error: problem while identifying')
         END IF
      ELSE
         ID = 0
      END IF
       
      IDENTIFIED = 1
      
      RETURN
      END
#endif

C==============================================================================
C======================== DLOAD SUBROUTINE ====================================
C==============================================================================

      SUBROUTINE DLOAD(F,KSTEP,KINC,TIME,NOEL,NPT,LAYER,KSPT,
     &   COORDS,JLTYP,SNAME)

      IMPLICIT NONE

      INTEGER D,S
      PARAMETER (D = |dimension|)
      PARAMETER (S = |surfaces|)
      INTEGER NOEL_PREV(S)
      COMMON /PREV/ NOEL_PREV
      SAVE /PREV/
      CHARACTER (LEN=80), DIMENSION(S) :: SURFACEIDS, PREPENDED
      COMMON /SURF/ SURFACEIDS
      SAVE /SURF/

#ifdef MPI
      INTEGER ID,IDENTIFIED
      COMMON /IDENT/ ID,IDENTIFIED
      SAVE /IDENT/
#else 
      INTEGER ID
#endif

      DOUBLE PRECISION F,TIME(2),COORDS(D)
      CHARACTER(LEN=20) :: FMT_FACES
      CHARACTER(LEN=200) :: FILENAME
      CHARACTER(LEN=80) :: SNAME
      INTEGER KSTEP,KINC,NOEL,NPT,LAYER,KSPT,JLTYP,R,UNIT_FACES(S)
      LOGICAL :: FOUND

      FMT_FACES = '(2I6,|dimension|ES27.17E2)'
      UNIT_FACES = (/ (100+R,R=1,S) /)

#ifdef MPI
      IF (IDENTIFIED < 0) THEN
         CALL IDENTIFY
      END IF
#else
      ID = 0
#endif

      FOUND  = .FALSE.
      IF (S > 1) THEN
         DO R = 1,S
            PREPENDED = 'ASSEMBLY_' // SURFACEIDS(R)
            IF (ALL(SNAME == PREPENDED)) THEN
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

      IF (NOEL >= NOEL_PREV(R)) THEN
         WRITE(FILENAME,'(A,A,I0,A,I0,A,I0,A)')
     &      '|PWD|',
     &      '/|CSM_dir|/CSM_Time',
     &      (KSTEP-1),'Surface',(R-1),'Cpu',ID,'Faces.dat'
         OPEN(UNIT=UNIT_FACES(R),FILE=FILENAME,POSITION='APPEND')
         WRITE(UNIT_FACES(R),FMT_FACES) NOEL,NPT,COORDS
         CLOSE(UNIT_FACES(R))
         NOEL_PREV(R) = NOEL
      ELSE IF (NOEL < NOEL_PREV(R)) THEN
         PRINT *, 'USR-abort: end of faces file reached. Normal term
     &ination'
         CALL FLUSH(6)
         CALL STDB_ABQERR(-3,'USR-abort: end of faces file.')
      END IF

      F = 0
      
      RETURN
      END

C==============================================================================
C======================== UTRACLOAD SUBROUTINE ================================
C==============================================================================

      SUBROUTINE UTRACLOAD(ALPHA,T_USER,KSTEP,KINC,TIME,
     &   NOEL,NPT,COORDS,DIRCOS,JLTYP,SNAME)

      IMPLICIT NONE
      
      INTEGER D
      PARAMETER (D = |dimension|)

      DOUBLE PRECISION ALPHA,T_USER(D),TIME(2),COORDS(D),DIRCOS(D,D)
      CHARACTER(LEN=80) :: SNAME
      INTEGER KSTEP,KINC,NOEL,NPT,JLTYP
          
      T_USER(1) = 1
      T_USER(2:) = 0
      ALPHA = 0
      
      RETURN
      END
