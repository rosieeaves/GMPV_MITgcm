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
      _RL gamma_q,r_EKE,r_Lambda,mu
CEOP

C--   Initialise parameterized EKE
C

      DO bj = myByLo(myThid), myByHi(myThid)
       DO bi = myBxLo(myThid), myBxHi(myThid)
        DO j=1-Oly,sNy+Oly
         DO i=1-Olx,sNx+Olx
          EKE(i,j,1,bi,bj)=1.E-4
         ENDDO
        ENDDO
       ENDDO
      ENDDO

      gamma_q = 0.1
      r_Lambda = 5.E-8
      r_K = 0.
      mu = 5000.

      WRITE(standardMessageUnit,*)
     &      'ini_EKE = ', EKE

C      DO bj = myByLo(myThid), myByHi(myThid)
C       DO bi = myBxLo(myThid), myBxHi(myThid)
C        DO j=1-Oly,sNy+Oly
C         DO i=1-Olx,sNx+Olx
C          EKE(i,j,1,bi,bj)=EKE(i,j,1,bi,bj)*maskC(i,j,1,bi,bj)
C         ENDDO
C        ENDDO
C       ENDDO
C      ENDDO

      WRITE(standardMessageUnit,*)
     &      'maskC = ', maskC

      WRITE(standardMessageUnit,*)
     &      'maskInC = ', maskInC

      WRITE(standardMessageUnit,*)
     &      'maskW = ', maskW

      WRITE(standardMessageUnit,*)
     &      'maskInW = ', maskInW

      WRITE(standardMessageUnit,*)
     &      'maskS = ', maskS

      WRITE(standardMessageUnit,*)
     &      'maskInS = ', maskInS

      RETURN
      END
