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
     I                     iMin,iMax,jMin,jMax,bi,bj,
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
#include "SURFACE.h"
#include "RESTART.h"
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
      INTEGER bi, bj
      _RL     myTime
      INTEGER myIter
      INTEGER myThid
      INTEGER nanInteger
      LOGICAL nanBoolean

C     !LOCAL VARIABLES:
C     i,j       :: Loop counters
C     gEKEAdv    :: EKE tendency from advection
C     gEKEConv    :: EKE tendence from energy conversion
C     gEKEr    :: EKE tendency from linear drag
C     gEKEDiff    :: EKE tendency from diffusion
      INTEGER i, j
C     YCXC
C      _RL  gEKEAdv(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) 
C      _RL  gEKEr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  gEKEDiff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  gEKE_AB(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  gLambdaAdv(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  gLambdar(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  gLambdaDiff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  gLambda_AB(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  q_YCXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     YGXC
      _RL  EKE_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Lambda_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  dqdx_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  dqdy_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  dpsidx_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  modGradq_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  kappa_q_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  energyConv_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  enstrophyGen_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     YCXG
      _RL  EKE_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Lambda_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  dqdx_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  dqdy_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  dpsidy_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  modGradq_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  kappa_q_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  energyConv_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  enstrophyGen_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)


C    Initialize arrays with zero values
      DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
          gEKEAdv(i,j,bi,bj) = 0. _d 0
          gEKEConv(i,j,bi,bj) = 0. _d 0
          gEKEr(i,j,bi,bj) = 0. _d 0
          gEKEDiff(i,j,bi,bj) = 0. _d 0
          gEKE_AB(i,j,bi,bj) = 0. _d 0

          gLambdaAdv(i,j,bi,bj) = 0. _d 0
          gLambdaGen(i,j,bi,bj) = 0. _d 0
          gLambdar(i,j,bi,bj) = 0. _d 0
          gLambdaDiff(i,j,bi,bj) = 0. _d 0
          gLambda_AB(i,j,bi,bj) = 0. _d 0
        ENDDO
      ENDDO 
C--------------------------------------------------------------
C     CALCULATE EKE & Lambda ADVECTION TERM AND LINEAR DRAG
C--------------------------------------------------------------

      DO j=1,sNy
        DO i=1,sNx
          IF (R_low(i,j,bi,bj).NE.0) THEN

C   H = -1*R_low => 1/H = -1/R_low

            gEKEAdv(i,j,bi,bj) = (-1/(2*R_low(i,j,bi,bj)
     &        *rA(i,j,bi,bj)))
     &        *((EKE(i,j,1,bi,bj) + EKE(i,j+1,1,bi,bj))
     &        *(psiVelInt(i+1,j+1,bi,bj) - psiVelInt(i,j+1,bi,bj)) 
     &        - (EKE(i,j,1,bi,bj) + EKE(i+1,j,1,bi,bj))
     &        *(psiVelInt(i+1,j+1,bi,bj) - psiVelInt(i+1,j,bi,bj))
     &        - (EKE(i,j,1,bi,bj) + EKE(i,j-1,1,bi,bj))
     &        *(psiVelInt(i+1,j,bi,bj) - psiVelInt(i,j,bi,bj)) 
     &        + (EKE(i,j,1,bi,bj) + EKE(i-1,j,1,bi,bj))
     &        *(psiVelInt(i,j+1,bi,bj) - psiVelInt(i,j,bi,bj)))

            gLambdaAdv(i,j,bi,bj) = (-1/(2*R_low(i,j,bi,bj)
     &        *rA(i,j,bi,bj)))
     &        *((Lambda(i,j,1,bi,bj) + Lambda(i,j+1,1,bi,bj))
     &        *(psiVelInt(i+1,j+1,bi,bj) - psiVelInt(i,j+1,bi,bj)) 
     &        - (Lambda(i,j,1,bi,bj) + Lambda(i+1,j,1,bi,bj))
     &        *(psiVelInt(i+1,j+1,bi,bj) - psiVelInt(i+1,j,bi,bj))
     &        - (Lambda(i,j,1,bi,bj) + Lambda(i,j-1,1,bi,bj))
     &        *(psiVelInt(i+1,j,bi,bj) - psiVelInt(i,j,bi,bj)) 
     &        + (Lambda(i,j,1,bi,bj) + Lambda(i-1,j,1,bi,bj))
     &        *(psiVelInt(i,j+1,bi,bj) - psiVelInt(i,j,bi,bj)))

            gEKEr(i,j,bi,bj) = r_EKE(i,j,bi,bj)
     &        *EKE(i,j,1,bi,bj)

            gLambdar(i,j,bi,bj) = r_Lambda(i,j,bi,bj)
     &        *Lambda(i,j,1,bi,bj)

          ELSE
            gEKEAdv(i,j,bi,bj)=0. 
            gLambdaAdv(i,j,bi,bj)=0.
            gEKEr(i,j,bi,bj)=0.
            gLambdar(i,j,bi,bj)=0.

          END IF
            
        ENDDO
      ENDDO
C--------------------------------------------------------------
C     CALCULATE K AND LAMBDA DIFFUSION    
C--------------------------------------------------------------

      DO j=1,sNy
        DO i=1,sNx

C if one or both of R_low = 0 then product = 0
C if neither R_low are zero then product > 0 => ceiling >= 1
C min => get zero when one or both of R_low is zero and one when neither is zero.

          gEKEDiff(i,j,bi,bj) = (mu(i,j,bi,bj)/rA(i,j,bi,bj))
     &      *((EKE(i,j+1,1,bi,bj) - EKE(i,j,1,bi,bj))
     &      *min(ceiling(R_low(i,j+1,bi,bj)*R_low(i,j,bi,bj)),1)
     &      +(EKE(i+1,j,1,bi,bj) - EKE(i,j,1,bi,bj))
     &      *min(ceiling(R_low(i+1,j,bi,bj)*R_low(i,j,bi,bj)),1)
     &      -(EKE(i,j,bi,1,bj) - EKE(i,j-1,1,bi,bj))
     &      *min(ceiling(R_low(i,j,bi,bj)*R_low(i,j-1,bi,bj)),1)
     &      -(EKE(i,j,1,bi,bj) - EKE(i-1,j,1,bi,bj))
     &      *min(ceiling(R_low(i,j,bi,bj)*R_low(i-1,j,bi,bj)),1))

          gLambdaDiff(i,j,bi,bj) = (mu(i,j,bi,bj)/rA(i,j,bi,bj))
     &      *((Lambda(i,j+1,1,bi,bj) - Lambda(i,j,1,bi,bj))
     &      *min(ceiling(R_low(i,j+1,bi,bj)*R_low(i,j,bi,bj)),1)
     &      +(Lambda(i+1,j,1,bi,bj) - Lambda(i,j,1,bi,bj))
     &      *min(ceiling(R_low(i+1,j,bi,bj)*R_low(i,j,bi,bj)),1)
     &      -(Lambda(i,j,1,bi,bj) - Lambda(i,j-1,1,bi,bj))
     &      *min(ceiling(R_low(i,j,bi,bj)*R_low(i,j-1,bi,bj)),1)
     &      -(Lambda(i,j,1,bi,bj) - Lambda(i-1,j,1,bi,bj))
     &      *min(ceiling(R_low(i,j,bi,bj)*R_low(i-1,j,bi,bj)),1))
        ENDDO
      ENDDO

C--------------------------------------------------------------
C     CALCULATE ENERGY CONVERSION AND ENSTROPHY GENERATION
C--------------------------------------------------------------



C     INTERPOLATE K & Lambda TO YGXC
      DO j=1,sNy
        DO i=1,sNx

          EKE_YGXC(i,j,bi,bj) = (EKE(i,j+1,1,bi,bj) 
     &      + EKE(i,j,1,bi,bj))/2

          Lambda_YGXC(i,j,bi,bj) = (Lambda(i,j+1,1,bi,bj) 
     &      + Lambda(i,j,1,bi,bj))/2

          EKE_YCXG(i,j,bi,bj) = (EKE(i+1,j,1,bi,bj) 
     &      + EKE(i,j,1,bi,bj))/2

          Lambda_YCXG(i,j,bi,bj) = (Lambda(i+1,j,1,bi,bj) 
     &      + Lambda(i,j,1,bi,bj))/2
        ENDDO
      ENDDO

C     DEFINE Q (MOVE THIS TO OUTSIDE OF TIME LOOP)
C     Q DEFINED AT CELL CENTRE POINTS
      DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
          q_YCXC(i,j,bi,bj) = -fCori(i,j,bi,bj)/R_low(i,j,bi,bj)

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: i,j = ', i,j

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: q_YCXC(i,j) = ', q_YCXC(i,j,bi,bj)
        ENDDO
      ENDDO

C     DQDX AND DPSIDY AT YCXG AND  DQDY AND DPSIDX YGXC
      DO j=1,sNy
        DO i=1,sNx

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: i,j = ', i,j

          dqdx_YCXG(i,j,bi,bj) = (q_YCXC(i+1,j,bi,bj) 
     &      - q_YCXC(i,j,bi,bj))/dxC(i,j,bi,bj)

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: dqdx_YCXG(i,j) = ', dqdx_YCXG(i,j,bi,bj) 

          dpsidx_YGXC(i,j,bi,bj) = (psiVelInt(i+1,j,bi,bj) 
     &      - psiVelInt(i,j,bi,bj))/dxG(i,j,bi,bj)

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: dpsidx_YGXC(i,j) = ', dpsidx_YGXC(i,j,bi,bj) 

          dqdy_YGXC(i,j,bi,bj) = (q_YCXC(i,j+1,bi,bj) 
     &      - q_YCXC(i,j,bi,bj))/dyC(i,j,bi,bj)

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: dqdy_YGXC(i,j) = ', dqdy_YGXC(i,j,bi,bj) 

          dpsidy_YCXG(i,j,bi,bj) = (psiVelInt(i,j+1,bi,bj) 
     &      - psiVelInt(i,j,bi,bj))/dyG(i,j,bi,bj)

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: dpsidy_YCXG(i,j) = ', dpsidy_YCXG(i,j,bi,bj) 

        ENDDO
      ENDDO

C     4 POINT AVERAGE OF DQDX_YCXG TO YGXC AND DQDY_YGXC TO YCXG
      DO j=1,sNy
        DO i=1,sNx
          dqdx_YGXC(i,j,bi,bj) = (dqdx_YCXG(i,j,bi,bj) 
     &      + dqdx_YCXG(i,j+1,bi,bj) 
     &      + dqdx_YCXG(i+1,j,bi,bj) 
     &      + dqdx_YCXG(i+1,j+1,bi,bj))/4

         dqdy_YCXG(i,j,bi,bj) = (dqdy_YGXC(i,j,bi,bj) 
     &      + dqdy_YGXC(i,j+1,bi,bj) 
     &      + dqdy_YGXC(i+1,j,bi,bj) 
     &      + dqdy_YGXC(i+1,j+1,bi,bj))/4
        ENDDO
      ENDDO

      DO j=1,sNy
        DO i=1,sNx

C     MOG(GRAD(Q)), KAPPA_Q, BAR(Q'U'), BAR(Q'U')*DPSIDX (ENERGY CONVERSION) 
C     and BAR(Q'U')*DQDX (ENSTROPHY GENERATION)AT YGXC
          modGradq_YGXC(i,j,bi,bj) = sqrt(dqdx_YGXC(i,j,bi,bj)**2 
     &      + dqdy_YGXC(i+1,j,bi,bj)**2)
          kappa_q_YGXC(i,j,bi,bj) = (2*gamma_q(i,j,bi,bj)
     &      *sqrt(EKE_YGXC(i+1,j,bi,bj)*Lambda_YGXC(i+1,j,bi,bj)))
     &      /modGradq_YGXC(i,j,bi,bj)
     &      *min(ceiling(R_low(i,j,bi,bj)*R_low(i,j-1,bi,bj)),1)
          quEDDY_YGXC(i,j,bi,bj) = -kappa_q_YGXC(i,j,bi,bj)
     &      *dqdx_YGXC(i,j,bi,bj)
          energyConv_YGXC(i,j,bi,bj) = quEDDY_YGXC(i,j,bi,bj)
     &      *dpsidx_YGXC(i+1,j+1,bi,bj)
          enstrophyGen_YGXC(i,j,bi,bj) = quEDDY_YGXC(i,j,bi,bj)
     &      *dqdx_YGXC(i,j,bi,bj)

C     MOG(GRAD(Q)), KAPPA_Q, BAR(Q'V') AND BAR(Q'V')*DPSIDY (ENERGY CONVERSION) 
C     and BAR(Q'V')*DQDY (ENSTROPHY GENERATION)AT YCXG

          modGradq_YCXG(i,j,bi,bj) = sqrt(dqdx_YCXG(i,j+1,bi,bj)**2 
     &      + dqdy_YCXG(i,j,bi,bj)**2)
          kappa_q_YCXG(i,j,bi,bj) = (2*gamma_q(i,j,bi,bj)
     &      *sqrt(EKE_YCXG(i,j+1,bi,bj)*Lambda_YCXG(i,j+1,bi,bj)))
     &      /modGradq_YCXG(i,j,bi,bj)
     &      *min(ceiling(R_low(i,j,bi,bj)*R_low(i-1,j,bi,bj)),1)
          qvEDDY_YCXG(i,j,bi,bj) = -kappa_q_YCXG(i,j,bi,bj)
     &      *dqdY_YCXG(i,j,bi,bj)
          energyConv_YCXG(i,j,bi,bj) = qvEDDY_YCXG(i,j,bi,bj)
     &      *dpsidy_YCXG(i+1,j+1,bi,bj)
          enstrophyGen_YCXG(i,j,bi,bj) = qvEDDY_YCXG(i,j,bi,bj)
     &      *dqdy_YCXG(i,j,bi,bj)
        ENDDO
      ENDDO

C     ENERGY CONVERSION AND ENSTROPHY GENERATION AT YCXC
      DO j=1,sNy
        DO i=1,sNx
          gEKEConv(i,j,bi,bj) = (energyConv_YGXC(i,j,bi,bj) 
     &      + energyConv_YGXC(i,j+1,bi,bj) 
     &      + energyConv_YCXG(i,j,bi,bj) 
     &      + energyConv_YCXG(i+1,j,bi,bj))/2
          gLambdaGen(i,j,bi,bj) = (enstrophyGen_YGXC(i,j,bi,bj) 
     &      + enstrophyGen_YGXC(i,j+1,bi,bj) 
     &      + enstrophyGen_YCXG(i,j,bi,bj) 
     &      + enstrophyGen_YCXG(i+1,j,bi,bj))/2
        ENDDO
      ENDDO 

C--------------------------------------------------------------
C     CALCULATE TOTAL EKE AND LAMBDA TENDENCY  
C--------------------------------------------------------------

      DO j=1,sNy
        DO i=1,sNx
          gEKE(i,j,1,bi,bj) = -gEKEAdv(i,j,bi,bj) 
     &      + gEKEConv(i,j,bi,bj) 
     &      - gEKEr(i,j,bi,bj) 
     &      + gEKEDiff(i,j,bi,bj)

          gLambda(i,j,1,bi,bj) = -gLambdaAdv(i,j,bi,bj) 
     &      - gLambdaGen(i,j,bi,bj) 
     &      - gLambdar(i,j,bi,bj) 
     &      + gLambdaDiff(i,j,bi,bj)
        ENDDO
      ENDDO


C--------------------------------------------------------------
C     CALCULATE gEKE^(n+1/2) WITH AB3
C--------------------------------------------------------------

      CALL ADAMS_BASHFORTH3(
     I                        bi, bj, 1, 1,
     U                        gEKE(1-OLx,1-OLy,1,bi,bj), gEKENm,
     O                        gEKE_AB,
     I                        mom_StartAB, myIter, myThid )

      CALL ADAMS_BASHFORTH3(
     I                        bi, bj, 1, 1,
     U                        gLambda(1-OLx,1-OLy,1,bi,bj), gLambdaNm,
     O                        gLambda_AB,
     I                        mom_StartAB, myIter, myThid )



C--------------------------------------------------------------
C     STEP FORWARD EKE AND LAMBDA
C--------------------------------------------------------------

      DO j=jMin,jMax
        DO i=iMin,iMax
          EKE(i,j,1,bi,bj) = EKE(i,j,1,bi,bj)
     &     + deltaTMom*( gEKE(i,j,1,bi,bj) )

          Lambda(i,j,1,bi,bj) = Lambda(i,j,1,bi,bj)
     &     + deltaTMom*( gLambda(i,j,1,bi,bj) )
        ENDDO
      ENDDO

      WRITE(standardMessageUnit,*)
     &      'calc_EKE: EKE = ', EKE


      RETURN
      END

