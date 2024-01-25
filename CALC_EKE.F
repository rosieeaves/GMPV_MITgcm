#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

C--  File CALC_EKE.F:
C--   Contents
C--   o CALC_EKE

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: CALC_EKE
C     !INTERFACE:
      SUBROUTINE CALC_EKE(
     I                     iMin,iMax,jMin,jMax, k, bi, bj,
     I                     myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R CALC_EKE
C     | o Calculates EKE tendency (gEKE) and steps forward K in time. 
C     *==========================================================*
C     \ev

C     !USES:
      IMPLICIT NONE
C     == Global data ==
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"
#include "GRID.h"
#include "DYNVARS.h"
#include "FFIELDS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     iMin,iMax :: Working range of x-index for applying forcing.
C     jMin,jMax :: Working range of y-index for applying forcing.
C     bi,bj     :: Current tile indices
C     myTime    :: Current time in simulation
C     myIter    :: Current iteration number
C     myThid    :: my Thread Id number
      INTEGER iMin, iMax, jMin, jMax
      INTEGER k, bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid

C     !LOCAL VARIABLES:
C     i,j       :: Loop counters
C     gEKEAdv    :: EKE tendency from advection
C     gEKEConv    :: EKE tendence from energy conversion
C     gEKEr    :: EKE tendency from linear drag
C     gEKEDiff    :: EKE tendency from diffusion
      INTEGER i, j
      _RL     gEKEAdv(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     gEKEConv(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     gEKEr(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL     gEKEDiff(1-OLx:sNx+OLx,1-OLy:sNy+OLy)

C    Initialize gEKE arrays with zero values
      DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
          gEKEAdv(i,j) = 0. _d 0
          gEKEConv(i,j) = 0. _d 0
          gEKEr(i,j) = 0. _d 0
          gEKEDiff(i,j) = 0. _d 0
        ENDDO
      ENDDO 

C     CALCULATE ADVECTION TERM
      DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
          gEKEAdv(i,j) = ((EKE(i,j,bi,bj) + EKE(i,j+1,bi,bj))
     &      *(psiVel(i+1,j+1,bi,bj) - psiVel(i,j+1,bi,bj)) 
     &      - (EKE(i,j,bi,bj) + EKE(i+1,j,bi,bj))
     &      *(psiVel(i+1,j+1,bi,bj) - psiVel(i+1,j,bi,bj))
     &      - (EKE(i,j,bi,bj) + EKE(i,j-1,bi,bj))
     &      *(psiVel(i+1,j,bi,bj) - psiVel(i,j,bi,bj)) 
     &      + (EKE(i,j,bi,bj) + EKE(i-1,j,bi,bj))
     &      *(psiVel(i,j+1,bi,bj) - psiVel(i,j,bi,bj)))
        ENDDO
      ENDDO
C     CALCULATE CONVERSION TERM 
C     CALCULATE DRAG TERM 
C     CALUCLATE DIFFUSION TERM     

      DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
          gEKE(i,j,bi,bj) = -gEKEAdv(i,j) 
     &      + gEKEConv(i,j) 
     &      - gEKEr(i,j) 
     &      + gEKEDiff(i,j)
        ENDDO
      ENDDO

      RETURN
      END

