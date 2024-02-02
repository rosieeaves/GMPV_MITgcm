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
     I                     iMin,iMax,jMin,jMax, bi, bj,
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

C     !LOCAL VARIABLES:
C     i,j       :: Loop counters
C     gEKEAdv    :: EKE tendency from advection
C     gEKEConv    :: EKE tendence from energy conversion
C     gEKEr    :: EKE tendency from linear drag
C     gEKEDiff    :: EKE tendency from diffusion
      INTEGER i, j
C     YCXC
      _RL  gEKEAdv(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy) 
      _RL  gEKEConv(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  gEKEr(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  gEKEDiff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  gEKE_AB(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  gLambdaAdv(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  gLambdaGen(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  gLambdar(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  gLambdaDiff(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  gLambda_AB(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  energyConv_YCXC(1-OLx+1:sNx+OLx-2,
     &      1-OLy+1:sNy+OLy-2,nSx,nSy)
      _RL  enstrophyGen_YCXC(1-OLx+1:sNx+OLx-2,
     &      1-OLy+1:sNy+OLy-2,nSx,nSy)
      _RL  EKEDiff_YCXC(1-OLx+1:sNx+OLx-1,1-OLy+1:sNy+OLy-1,nSx,nSy)
      _RL  LambdaDiff_YCXC(1-OLx+1:sNx+OLx-1,
     &      1-OLy+1:sNy+OLy-1,nSx,nSy)
C     YGXC
      _RL  EKE_YGXC(1-OLx:sNx+OLx,1-OLy+1:sNy+OLy,nSx,nSy)
      _RL  Lambda_YGXC(1-OLx:sNx+OLx,1-OLy+1:sNy+OLy,nSx,nSy)
      _RL  dqdx_YGXC(1-OLx:sNx+OLx-1,1-OLy:sNy+OLy,nSx,nSy)
      _RL  dqdy_YGXC(1-OLx:sNx+OLx-1,1-OLy+1:sNy+OLy-1,nSx,nSy)
      _RL  dpsidx_YGXC(1-OLx:sNx+OLx-1,1-OLy:sNy+OLy,nSx,nSy)
      _RL  modGradq_YGXC(1-OLx:sNx+OLx-1,1-OLy+1:sNy+OLy-1,nSx,nSy)
      _RL  kappa_q_YGXC(1-OLx:sNx+OLx-1,1-OLy+1:sNy+OLy-1,nSx,nSy)
      _RL  energyConv_YGXC(1-OLx:sNx+OLx-1,1-OLy+1:sNy+OLy-1,nSx,nSy)
      _RL  enstrophyGen_YGXC(1-OLx:sNx+OLx-1,
     &      1-OLy+1:sNy+OLy-1,nSx,nSy)
C     YCXG
      _RL  EKE_YCXG(1-OLx+1:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  Lambda_YCXG(1-OLx+1:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  dqdx_YCXG(1-OLx+1:sNx+OLx-1,1-OLy:sNy+OLy-1,nSx,nSy)
      _RL  dqdy_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy-1,nSx,nSy)
      _RL  dpsidy_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy-1,nSx,nSy)
      _RL  modGradq_YCXG(1-OLx+1:sNx+OLx-1,1-OLy:sNy+OLy-1,nSx,nSy)
      _RL  kappa_q_YCXG(1-OLx+1:sNx+OLx-1,1-OLy:sNy+OLy-1,nSx,nSy)
      _RL  energyConv_YCXG(1-OLx+1:sNx+OLx-1,1-OLy:sNy+OLy-1,nSx,nSy)
      _RL  enstrophyGen_YCXG(1-OLx+1:sNx+OLx-1,
     &      1-OLy:sNy+OLy-1,nSx,nSy)
C     YGXG
      _RL  q_YGXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)


C    Initialize gEKE arrays with zero values
      DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
          gEKEAdv(i,j,bi,bj) = 0. _d 0
          gEKEConv(i,j,bi,bj) = 0. _d 0
          gEKEr(i,j,bi,bj) = 0. _d 0
          gEKEDiff(i,j,bi,bj) = 0. _d 0
        ENDDO
      ENDDO 

      WRITE(standardMessageUnit,*)
     &      'calc_EKE: psiVelInt = ', psiVelInt

      WRITE(standardMessageUnit,*)
     &      'calc_EKE: EKE = ', EKE

      WRITE(standardMessageUnit,*)
     &      'calc_EKE: R_low = ', R_low

      WRITE(standardMessageUnit,*)
     &      'calc_EKE: Ro_surf = ', Ro_surf

      WRITE(standardMessageUnit,*)
     &      'calc_EKE: hFacC = ', hFacC

      WRITE(standardMessageUnit,*)
     &      'calc_EKE: rA = ', rA

      WRITE(standardMessageUnit,*)
     &      'calc_EKE: topography = ', topography
C--------------------------------------------------------------
C     CALCULATE EKE & Lambda ADVECTION TERM AND LINEAR DRAG
C--------------------------------------------------------------

      DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx

          gEKEAdv(i,j,bi,bj) = (-1/(2*rA(i,j,bi,bj)))
     &      *((EKE(i,j,1,bi,bj) + EKE(i,j+1,1,bi,bj))
     &      *(psiVelInt(i+1,j+1,bi,bj) - psiVelInt(i,j+1,bi,bj)) 
     &      - (EKE(i,j,1,bi,bj) + EKE(i+1,j,1,bi,bj))
     &      *(psiVelInt(i+1,j+1,bi,bj) - psiVelInt(i+1,j,bi,bj))
     &      - (EKE(i,j,1,bi,bj) + EKE(i,j-1,1,bi,bj))
     &      *(psiVelInt(i+1,j,bi,bj) - psiVelInt(i,j,bi,bj)) 
     &      + (EKE(i,j,1,bi,bj) + EKE(i-1,j,1,bi,bj))
     &      *(psiVelInt(i,j+1,bi,bj) - psiVelInt(i,j,bi,bj)))

          gLambdaAdv(i,j,bi,bj) = (-1/(2*rA(i,j,bi,bj)))
     &      *((Lambda(i,j,1,bi,bj) + Lambda(i,j+1,1,bi,bj))
     &      *(psiVelInt(i+1,j+1,bi,bj) - psiVelInt(i,j+1,bi,bj)) 
     &      - (Lambda(i,j,1,bi,bj) + Lambda(i+1,j,1,bi,bj))
     &      *(psiVelInt(i+1,j+1,bi,bj) - psiVelInt(i+1,j,bi,bj))
     &      - (Lambda(i,j,1,bi,bj) + Lambda(i,j-1,1,bi,bj))
     &      *(psiVelInt(i+1,j,bi,bj) - psiVelInt(i,j,bi,bj)) 
     &      + (Lambda(i,j,1,bi,bj) + Lambda(i-1,j,1,bi,bj))
     &      *(psiVelInt(i,j+1,bi,bj) - psiVelInt(i,j,bi,bj)))

           gEKEr(i,j,bi,bj) = r_EKE*EKE(i,j,1,bi,bj)

           gLambdar(i,j,bi,bj) = r_Lambda*Lambda(i,j,1,bi,bj)
        ENDDO
      ENDDO
C--------------------------------------------------------------
C     CALCULATE K AND LAMBDA DIFFUSION    
C--------------------------------------------------------------

C      DO j=1-OLy+1,sNy+OLy-1
C        DO i=1-OLx+1,sNx+OLx-1
C          EKEDiff_YCXC(i,j,bi,bj) = (mu/rA)*(EKE(i,j+1,1,bi,bj) 
C     &      + EKE(i+1,j,1,bi,bj) 
C     &      + EKE(i,j-1,1,bi,bj) 
C     &      + EKE(i-1,j,1,bi,bj) 
C     &      - 4*EKE(i,j,1,bi,bj))

C          LambdaDiff_YCXC(i,j,bi,bj) = (mu/rA)*(Lambda(i,j+1,1,bi,bj) 
C     &      + Lambda(i+1,j,1,bi,bj) 
C     &      + Lambda(i,j-1,1,bi,bj) 
C     &      + Lambda(i-1,j,1,bi,bj) 
C     &      - 4*Lambda(i,j,1,bi,bj))
C        ENDDO
C      ENDDO

      WRITE(standardMessageUnit,*)
     &      'calc_EKE: gEKEAdv = ', gEKEAdv
C--------------------------------------------------------------
C     CALCULATE ENERGY CONVERSION AND ENSTROPHY GENERATION
C--------------------------------------------------------------

C     INTERPOLATE K & Lambda TO YGXC
      DO j=1-OLy+1,sNy+OLy
        DO i=1-OLx,sNx+OLx
          EKE_YGXC(i,j,bi,bj) = (EKE(i,j+1,1,bi,bj) 
     &      + EKE(i,j,1,bi,bj))/2
          Lambda_YGXC(i,j,bi,bj) = (Lambda(i,j+1,1,bi,bj) 
     &      + Lambda(i,j,1,bi,bj))/2
        ENDDO
      ENDDO
C     INTERPOLATE K & Lambda TO YCXG
      DO j=1-OLy,sNy+OLy
        DO i=1-OLx+1,sNx+OLx
          EKE_YCXG(i,j,bi,bj) = (EKE(i+1,j,1,bi,bj) 
     &      + EKE(i,j,1,bi,bj))/2
          Lambda_YCXG(i,j,bi,bj) = (Lambda(i+1,j,1,bi,bj) 
     &      + Lambda(i,j,1,bi,bj))/2
        ENDDO
      ENDDO

C     DEFINE Q (MOVE THIS TO OUTSIDE OF TIME LOOP AND ADD FACTOR 1/H IN)
      DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
          q_YGXG(i,j,bi,bj) = fCoriG(i,j,bi,bj)
        ENDDO
      ENDDO

C     DQDX AND DPSIDX AT YGXC
C      DO j=1-OLy,sNy+OLy
C        DO i=1-OLx,sNx+OLx-1
C          dqdx_YGXC(i,j,bi,bj) = (q_YGXG(i+1,j,bi,bj) 
C     &      - q_YGXG(i,j,bi,bj))/dxG
C          dpsidx_YGXC(i,j,bi,bj) = (psiVelInt(i+1,j,bi,bj) 
C     &      - psiVelInt(i,j,bi,bj))/dxG
C        ENDDO
C      ENDDO
C     DQDY AND DPSIDY AT YCXG
C      DO j=1-OLy,sNy+OLy-1
C        DO i=1-OLx,sNx+OLx
C          dqdy_YCXG(i,j,bi,bj) = (q_YGXG(i,j+1,bi,bj) 
C     &      - q_YGXG(i,j,bi,bj))/dyG
C          dpsidy_YCXG(i,j,bi,bj) = (psiVelInt(i,j+1,bi,bj) 
C     &      - psiVelInt(i,j,bi,bj))/dyG
C        ENDDO
C      ENDDO
C     4 POINT AVERAGE OF DQDX_YGXC TO YCXG
      DO j=1-OLy,sNy+OLy-1
        DO i=1-OLx+1,sNx+OLx-1
          dqdx_YCXG(i,j,bi,bj) = (dqdx_YGXC(i,j,bi,bj) 
     &      + dqdx_YGXC(i,j+1,bi,bj) 
     &      + dqdx_YGXC(i+1,j,bi,bj) 
     &      + dqdx_YGXC(i+1,j+1,bi,bj))/4
        ENDDO
      ENDDO
C     4 POINT AVERAGE OF DQDY TO YGXC
      DO j=1-OLy+1,sNy+OLy-1
        DO i=1-OLx,sNx+OLx-1
          dqdy_YGXC(i,j,bi,bj) = (dqdy_YGXC(i,j,bi,bj) 
     &      + dqdy_YGXC(i,j+1,bi,bj) 
     &      + dqdy_YGXC(i+1,j,bi,bj) 
     &      + dqdy_YGXC(i+1,j+1,bi,bj))/4
        ENDDO
      ENDDO
C     MOG(GRAD(Q)), KAPPA_Q, BAR(Q'U'), BAR(Q'U')*DPSIDX (ENERGY CONVERSION) 
C     and BAR(Q'U')*DQDX (ENSTROPHY GENERATION)AT YGXC
C      DO j=1-OLy+1,sNy+OLy-1
C        DO i=1-OLx,sNx+OLx-1
C          modGradq_YGXC(i,j,bi,bj) = sqrt(dqdx_YGXC(i,j+1,bi,bj)**2 
C     &      + dqdy_YGXC(i,j,bi,bj)**2)
C          kappa_q_YGXC(i,j,bi,bj) = (2*gamma_q
C     &      *sqrt(EKE_YGXC(i,j,bi,bj)*Lambda_YGXC(i,j,bi,bj)))
C     &      /modGradq_YGXC
C          quEDDY_YGXC(i,j,bi,bj) = -kappa_q_YGXC(i,j,bi,bj)
C     &      *dqdx_YGXC(i,j+1,bi,bj)
C          energyConv_YGXC(i,j,bi,bj) = quEDDY_YGXC(i,j,bi,bj)
C     &      *dpsidx_YGXC(i,j+1,bi,bj)
C          enstrophyGen_YGXC(i,j,bi,bj) = quEDDY_YGXC(i,j,bi,bj)
C     &      *dqdx_YGXC(i,j+1,bi,bj)
C        ENDDO
C      ENDDO
C     MOG(GRAD(Q)), KAPPA_Q, BAR(Q'V') AND BAR(Q'V')*DPSIDY (ENERGY CONVERSION) 
C     and BAR(Q'V')*DQDY (ENSTROPHY GENERATION)AT YCXG
C      DO j=1-OLy,sNy+OLy-1
C        DO i=1-OLx+1,sNx+OLx-1
C          modGradq_YCXG(i,j,bi,bj) = sqrt(dqdx_YCXG(i,j,bi,bj)**2 
C     &      + dqdy_YCXG(i+1,j,bi,bj)**2)
C          kappa_q_YCXG(i,j,bi,bj) = (2*gamma_q
C     &      *sqrt(EKE_YCXG(i,j,bi,bj)*Lambda_YCXG(i,j,bi,bj)))
C     &      /modGradq_YCXG
C          qvEDDY_YCXG(i,j,bi,bj) = -kappa_q_YCXG(i,j,bi,bj)
C     &      *dqdY_YCXG(i+1,j,bi,bj)
C          energyConv_YCXG(i,j,bi,bj) = qvEDDY_YCXG(i,j,bi,bj)
C     &      *dpsidy_YCXG(i+1,j,bi,bj)
C          enstrophyGen_YCXG(i,j,bi,bj) = qvEDDY_YCXG(i,j,bi,bj)
C     &      *dqdy_YCXG(i+1,j,bi,bj)
C        ENDDO
C      ENDDO

C     ENERGY CONVERSION AND ENSTROPHY GENERATION AT YCXC
      DO j=1-OLy+1,sNy+OLy-2
        DO i=1-OLx+1,sNx+OLx-2
          energyConv_YCXC(i,j,bi,bj) = (energyConv_YGXC(i+1,j,bi,bj) 
     &      + energyConv_YGXC(i+1,j+1,bi,bj) 
     &      + energyConv_YCXG(i,j+1,bi,bj) 
     &      + energyConv_YCXG(i+1,j+1,bi,bj))/2
          enstrophyGen_YCXC(i,j,bi,bj) = (enstrophyGen_YGXC(i+1,j,bi,bj) 
     &      + enstrophyGen_YGXC(i+1,j+1,bi,bj) 
     &      + enstrophyGen_YCXG(i,j+1,bi,bj) 
     &      + enstrophyGen_YCXG(i+1,j+1,bi,bj))/2
        ENDDO
      ENDDO 

C--------------------------------------------------------------
C     CALCULATE TOTAL EKE AND LAMBDA TENDENCY  
C--------------------------------------------------------------

      DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
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

      WRITE(standardMessageUnit,*)
     &      'calc_EKE: gEKE = ', gEKE

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

      WRITE(standardMessageUnit,*)
     &      'calc_EKE: deltaTMom = ', deltaTMom

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

