#include "CPP_OPTIONS.h"

CBOP
C     !ROUTINE: INI_EKE
C     !INTERFACE:
      SUBROUTINE INI_EKE( myThid )

C     !DESCRIPTION: \bv
C     *=================================================================
C     | SUBROUTINE INI_EKE
C     | o Initialize parameterized EKE
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
C     myThid ::  Number of this instance of INI_EKE
      INTEGER myThid

C     !LOCAL VARIABLES:
C     == Local variables ==
C     i,j,k  :: Loop counters
C     bi,bj  :: tile indices
      INTEGER bi,bj,i,j,k
CEOP

C--   Initialise parameterized EKE
C

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1-Oly,sNy+Oly
         DO i=1-Olx,sNx+Olx
          IF (R_low(i,j,bi,bj).NE.0) THEN
            EKE(i,j,1,bi,bj)=1.E-4
          ELSE 
            EKE(i,j,1,bi,bj)=0.
          END IF
          gamma_q(i,j,bi,bj)=0.1
          r_EKE(i,j,bi,bj)=0.
          mu(i,j,bi,bj)=5000.
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      WRITE(standardMessageUnit,*)
     &      'ini_EKE = ', EKE

      RETURN
      END
