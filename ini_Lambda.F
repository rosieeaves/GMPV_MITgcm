#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: INI_Lambda
C     !INTERFACE:
      SUBROUTINE INI_LAMBDA( myThid )

C     !DESCRIPTION: \bv
C     *=================================================================
C     | SUBROUTINE INI_LAMBDA
C     | o Initialize parameterized Lambda
C     *=================================================================
C     \ev

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     == Routine arguments ==
C     myThid ::  Number of this instance of INI_LAMBDA
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j,k  :: Loop counters
C     bi,bj  :: tile indices
      INTEGER bi,bj,i,j,k
CEOP

C--   Initialise parameterized Lambda
C

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1-Oly,sNy+Oly
         DO i=1-Olx,sNx+Olx
           IF (R_low(i,j,bi,bj).NE.0) THEN
            Lambda(i,j,1,bi,bj)=1.E-20*BT_mask(i,j,bi,bj)
           ELSE 
            Lambda(i,j,1,bi,bj)=0.
           END IF 
           r_Lambda(i,j,bi,bj)=5.E-8
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      WRITE(standardMessageUnit,*)
     &      'ini_LAMBDA = ', Lambda

      RETURN
      END
