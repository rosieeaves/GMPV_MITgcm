#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

C--  File CALC_EKE.F:
C--   Contents
C--   o CALC_EKE

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|
CBOP
C     !ROUTINE: CALC_EKE
C     !INTERFACE:
      SUBROUTINE CALC_EDDY_FLUX(
     I                     iMin,iMax,jMin,jMax,bi,bj,
     I                     myTime, myIter, myThid )
C     !DESCRIPTION: \bv
C     *==========================================================*
C     | S/R CALC_EDDY_FLUX
C     | o Calculates eddy fluxes by:
C     |     - Calculating EKE & Lambda advection
C     |     - Calculating EKE & Lambda diffusion
C     |     - Calculating EKE & Lambda dissipation
C     |     - Calculating EKE conversion and Lambda generation
C     |     - (Eddy fluxes calculated in above step)
C     |     - Step forward EKE and Lambda in time 
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

C     !LOCAL VARIABLES:
C     i,j       :: Loop counters
C     gEKEAdv    :: EKE tendency from advection
C     gEKEConv    :: EKE tendence from energy conversion
C     gEKEr    :: EKE tendency from linear drag
C     gEKEDiff    :: EKE tendency from diffusion
      INTEGER i, j
C     YCXC
      _RL  gEKE_AB(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  gLambda_AB(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  q_YCXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     YGXC
C      _RL  EKE_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  Lambda_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  dqdx_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  dqdy_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  dpsidx_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  modGradq_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  kappa_q_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  energyConv_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  enstrophyGen_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C     YCXG
C      _RL  EKE_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  Lambda_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  dqdx_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  dqdy_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  dpsidy_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  modGradq_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  kappa_q_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  energyConv_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  enstrophyGen_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)


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

          EKE_YGXC(i,j,bi,bj) = 0. _d 0
          Lambda_YGXC(i,j,bi,bj) = 0. _d 0
          dqdx_YGXC(i,j,bi,bj) = 0. _d 0
          dqdy_YGXC(i,j,bi,bj) = 0. _d 0
          dpsidx_YGXC(i,j,bi,bj) = 0. _d 0
          modGradq_YGXC(i,j,bi,bj) = 0. _d 0
          kappa_q_YGXC(i,j,bi,bj) = 0. _d 0
          energyConv_YGXC(i,j,bi,bj) = 0. _d 0
          enstrophyGen_YGXC(i,j,bi,bj) = 0. _d 0
          quEDDY_YGXC(i,j,bi,bj) = 0. _d 0

          EKE_YCXG(i,j,bi,bj) = 0. _d 0
          Lambda_YCXG(i,j,bi,bj) = 0. _d 0
          dqdx_YCXG(i,j,bi,bj) = 0. _d 0
          dqdy_YCXG(i,j,bi,bj) = 0. _d 0
          dpsidy_YCXG(i,j,bi,bj) = 0. _d 0
          modGradq_YCXG(i,j,bi,bj) = 0. _d 0
          kappa_q_YCXG(i,j,bi,bj) = 0. _d 0
          energyConv_YCXG(i,j,bi,bj) = 0. _d 0
          enstrophyGen_YCXG(i,j,bi,bj) = 0. _d 0
          qvEDDY_YCXG(i,j,bi,bj) = 0. _d 0

        ENDDO
      ENDDO 
C--------------------------------------------------------------
C     CALCULATE EKE & Lambda ADVECTION TERM AND LINEAR DRAG
C--------------------------------------------------------------
      WRITE(standardMessageUnit,*)
     &      'calc_EKE: EKE  = ', EKE

      DO j=1,sNy
        DO i=1,sNx
          IF (abs(R_low(i,j,bi,bj)).NE.0) THEN

C   H = abs(R_low)

            gEKEAdv(i,j,bi,bj) = (1/(2*abs(R_low(i,j,bi,bj))
     &        *rA(i,j,bi,bj)))
     &        *((EKE(i,j,1,bi,bj) + EKE(i,j+1,1,bi,bj))
     &        *(psiVelInt(i+1,j+1,bi,bj) - psiVelInt(i,j+1,bi,bj)) 
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j+1,bi,bj))),1)
     &        - (EKE(i,j,1,bi,bj) + EKE(i+1,j,1,bi,bj))
     &        *(psiVelInt(i+1,j+1,bi,bj) - psiVelInt(i+1,j,bi,bj))
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i+1,j,bi,bj))),1)
     &        - (EKE(i,j,1,bi,bj) + EKE(i,j-1,1,bi,bj))
     &        *(psiVelInt(i+1,j,bi,bj) - psiVelInt(i,j,bi,bj)) 
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j-1,bi,bj))),1)
     &        + (EKE(i,j,1,bi,bj) + EKE(i-1,j,1,bi,bj))
     &        *(psiVelInt(i,j+1,bi,bj) - psiVelInt(i,j,bi,bj)))
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i-1,j,bi,bj))),1)

            gLambdaAdv(i,j,bi,bj) = (1/(2*abs(R_low(i,j,bi,bj))
     &        *rA(i,j,bi,bj)))
     &        *((Lambda(i,j,1,bi,bj) + Lambda(i,j+1,1,bi,bj))
     &        *(psiVelInt(i+1,j+1,bi,bj) - psiVelInt(i,j+1,bi,bj)) 
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j+1,bi,bj))),1)
     &        - (Lambda(i,j,1,bi,bj) + Lambda(i+1,j,1,bi,bj))
     &        *(psiVelInt(i+1,j+1,bi,bj) - psiVelInt(i+1,j,bi,bj))
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i+1,j,bi,bj))),1)
     &        - (Lambda(i,j,1,bi,bj) + Lambda(i,j-1,1,bi,bj))
     &        *(psiVelInt(i+1,j,bi,bj) - psiVelInt(i,j,bi,bj)) 
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j-1,bi,bj))),1)
     &        + (Lambda(i,j,1,bi,bj) + Lambda(i-1,j,1,bi,bj))
     &        *(psiVelInt(i,j+1,bi,bj) - psiVelInt(i,j,bi,bj)))
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i-1,j,bi,bj))),1)

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

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: gEKEAdv(i,j)  = ', gEKEAdv(i,j,bi,bj)

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: gLambdaAdv(i,j)  = ', gLambdaAdv(i,j,bi,bj)

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: gEKEr(i,j)  = ', gEKEr(i,j,bi,bj)

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: gLambdar(i,j)  = ', gLambdar(i,j,bi,bj)
            
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
     &      *min(ceiling(abs(R_low(i,j+1,bi,bj)*R_low(i,j,bi,bj))),1)
     &      +(EKE(i+1,j,1,bi,bj) - EKE(i,j,1,bi,bj))
     &      *min(ceiling(abs(R_low(i+1,j,bi,bj)*R_low(i,j,bi,bj))),1)
     &      -(EKE(i,j,bi,1,bj) - EKE(i,j-1,1,bi,bj))
     &      *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j-1,bi,bj))),1)
     &      -(EKE(i,j,1,bi,bj) - EKE(i-1,j,1,bi,bj))
     &      *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i-1,j,bi,bj))),1))

          gLambdaDiff(i,j,bi,bj) = (mu(i,j,bi,bj)/rA(i,j,bi,bj))
     &      *((Lambda(i,j+1,1,bi,bj) - Lambda(i,j,1,bi,bj))
     &      *min(ceiling(abs(R_low(i,j+1,bi,bj)*R_low(i,j,bi,bj))),1)
     &      +(Lambda(i+1,j,1,bi,bj) - Lambda(i,j,1,bi,bj))
     &      *min(ceiling(abs(R_low(i+1,j,bi,bj)*R_low(i,j,bi,bj))),1)
     &      -(Lambda(i,j,1,bi,bj) - Lambda(i,j-1,1,bi,bj))
     &      *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j-1,bi,bj))),1)
     &      -(Lambda(i,j,1,bi,bj) - Lambda(i-1,j,1,bi,bj))
     &      *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i-1,j,bi,bj))),1))

           WRITE(standardMessageUnit,*)
     &      'calc_EKE: gEKEDiff(i,j)  = ', gEKEDiff(i,j,bi,bj)

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: gLambdaDiff(i,j)  = ', gLambdaDiff(i,j,bi,bj)
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
C     no more Infinites
      DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
          IF (R_low(i,j,bi,bj).NE.0) THEN
            q_YCXC(i,j,bi,bj) = fCori(i,j,bi,bj)/abs(R_low(i,j,bi,bj))

          ELSE
            q_YCXC(i,j,bi,bj) = 0.
          END IF

        ENDDO
      ENDDO

      DO j=1,sNy
        DO i=1,sNx

C     DQDX AND DPSIDY AT YCXG 

          dqdx_YCXG(i,j,bi,bj) = (q_YCXC(i+1,j,bi,bj) 
     &      - q_YCXC(i,j,bi,bj))/dxC(i,j,bi,bj)
     &      *min(ceiling(abs(R_low(i+1,j,bi,bj)*R_low(i,j,bi,bj))),1)

          dpsidx_YGXC(i,j,bi,bj) = (psiVelInt(i+1,j,bi,bj) 
     &      - psiVelInt(i,j,bi,bj))/dxG(i,j,bi,bj)

C     DQDY AND DPSIDX AT YGXC

          dqdy_YGXC(i,j,bi,bj) = (q_YCXC(i,j+1,bi,bj) 
     &      - q_YCXC(i,j,bi,bj))/dyC(i,j,bi,bj)
     &      *min(ceiling(abs(R_low(i,j+1,bi,bj)*R_low(i,j,bi,bj))),1)

          dpsidy_YCXG(i,j,bi,bj) = (psiVelInt(i,j+1,bi,bj) 
     &      - psiVelInt(i,j,bi,bj))/dyG(i,j,bi,bj)

        ENDDO
      ENDDO

C     4 POINT AVERAGE OF DQDX_YCXG TO YGXC AND DQDY_YGXC TO YCXG
C     no more Nans
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     need to check what happens at land boundaries
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
C     no more NaNs in first time step

          modGradq_YGXC(i,j,bi,bj) = max(sqrt(dqdx_YGXC(i,j,bi,bj)**2 
     &      + dqdy_YGXC(i,j,bi,bj)**2),1.E-16)

          kappa_q_YGXC(i,j,bi,bj) = (2*gamma_q(i,j,bi,bj)
     &      *sqrt(EKE_YGXC(i,j,bi,bj)*Lambda_YGXC(i,j,bi,bj)))
     &      /modGradq_YGXC(i,j,bi,bj)
     &      *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j-1,bi,bj))),1)
          
          quEDDY_YGXC(i,j,bi,bj) = -kappa_q_YGXC(i,j,bi,bj)
     &      *dqdx_YGXC(i,j,bi,bj)

          energyConv_YGXC(i,j,bi,bj) = quEDDY_YGXC(i,j,bi,bj)
     &      *dpsidx_YGXC(i,j,bi,bj)

          enstrophyGen_YGXC(i,j,bi,bj) = quEDDY_YGXC(i,j,bi,bj)
     &      *dqdx_YGXC(i,j,bi,bj)

C     MOG(GRAD(Q)), KAPPA_Q, BAR(Q'V') AND BAR(Q'V')*DPSIDY (ENERGY CONVERSION) 
C     and BAR(Q'V')*DQDY (ENSTROPHY GENERATION)AT YCXG

          modGradq_YCXG(i,j,bi,bj) = max(sqrt(dqdx_YCXG(i,j,bi,bj)**2 
     &      + dqdy_YCXG(i,j,bi,bj)**2),1.E-16)

          kappa_q_YCXG(i,j,bi,bj) = (2*gamma_q(i,j,bi,bj)
     &      *sqrt(EKE_YCXG(i,j,bi,bj)*Lambda_YCXG(i,j,bi,bj)))
     &      /modGradq_YCXG(i,j,bi,bj)
     &      *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i-1,j,bi,bj))),1)
          

          qvEDDY_YCXG(i,j,bi,bj) = -kappa_q_YCXG(i,j,bi,bj)
     &      *dqdy_YCXG(i,j,bi,bj)

          energyConv_YCXG(i,j,bi,bj) = qvEDDY_YCXG(i,j,bi,bj)
     &      *dpsidy_YCXG(i,j,bi,bj)

          enstrophyGen_YCXG(i,j,bi,bj) = qvEDDY_YCXG(i,j,bi,bj)
     &      *dqdy_YCXG(i,j,bi,bj)

        ENDDO
      ENDDO

C     ENERGY CONVERSION AND ENSTROPHY GENERATION AT YCXC
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C     need to check what happens at land boundaries
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

       DO j=1,sNy
        DO i=1,sNx
          WRITE(standardMessageUnit,*)
     &      'calc_EKE: AB gEKE(i,j)  = ', gEKE(i,j,1,bi,bj)

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: AB gLambda(i,j)  = ', gLambda(i,j,1,bi,bj)
        ENDDO
       ENDDO

C--------------------------------------------------------------
C     STEP FORWARD EKE AND LAMBDA
C--------------------------------------------------------------

      DO j=jMin,jMax
        DO i=iMin,iMax
          EKE(i,j,1,bi,bj) = EKE(i,j,1,bi,bj)
     &     + deltaTMom*( gEKE(i,j,1,bi,bj) )

          Lambda(i,j,1,bi,bj) = Lambda(i,j,1,bi,bj)
     &     + deltaTMom*( gLambda(i,j,1,bi,bj) )

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: EKE(i,j)  = ', EKE(i,j,1,bi,bj)

          WRITE(standardMessageUnit,*)
     &      'calc_EKE: Lambda(i,j)  = ', Lambda(i,j,1,bi,bj)
        ENDDO
      ENDDO


      RETURN
      END

