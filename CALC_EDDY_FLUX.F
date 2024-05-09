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
      INTEGER i, j, k
C     YCXC
      _RL  gEKE_AB(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
      _RL  gLambda_AB(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  q_YCXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
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
          dpsidx_YGXC(i,j,bi,bj) = 0. _d 0
          kappa_q_YGXC(i,j,bi,bj) = 0. _d 0
          energyConv_YGXC(i,j,bi,bj) = 0. _d 0
          enstrophyGen_YGXC(i,j,bi,bj) = 0. _d 0
          quEDDY_YGXC(i,j,bi,bj) = 0. _d 0

          EKE_YCXG(i,j,bi,bj) = 0. _d 0
          Lambda_YCXG(i,j,bi,bj) = 0. _d 0
          dpsidy_YCXG(i,j,bi,bj) = 0. _d 0
          kappa_q_YCXG(i,j,bi,bj) = 0. _d 0
          energyConv_YCXG(i,j,bi,bj) = 0. _d 0
          enstrophyGen_YCXG(i,j,bi,bj) = 0. _d 0
          qvEDDY_YCXG(i,j,bi,bj) = 0. _d 0

          uVelInt(i,j,bi,bj) = 0. _d 0
          vVelInt(i,j,bi,bj) = 0. _d 0

        ENDDO
      ENDDO 

C--------------------------------------------------------------
C     DEPTH INTEGRAL OF U AND V FOR ENERGY CONVERSION
C     AND ADVECTION
C--------------------------------------------------------------
      DO j=1,sNy
        DO i=1,sNx
          DO k=1,Nr
            uVelInt(i,j,bi,bj) = uVelInt(i,j,bi,bj) + 
     &          hFacW(i,j,k,bi,bj)*drF(k)*uVel(i,j,k,bi,bj)

            vVelInt(i,j,bi,bj) = vVelInt(i,j,bi,bj) + 
     &          hFacS(i,j,k,bi,bj)*drF(k)*vVel(i,j,k,bi,bj)
          ENDDO
        ENDDO
      ENDDO


      DO j=1,sNy
        DO i=1,sNx
          IF (abs(R_low(i,j,bi,bj)).NE.0) THEN

C--------------------------------------------------------------
C     CALCULATE EKE & Lambda ADVECTION TERM 
C--------------------------------------------------------------
C   H = abs(R_low)

            gEKEAdv(i,j,bi,bj) = (1/(4*abs(R_low(i,j,bi,bj))
     &        *rA(i,j,bi,bj)))

     &        *((EKE(i,j,1,bi,bj) + EKE(i,j+1,1,bi,bj))
     &        *vVelInt(i,j+1,bi,bj)*dxG(i,j+1,bi,bj)
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j+1,bi,bj))),1)

     &        + (EKE(i,j,1,bi,bj) + EKE(i+1,j,1,bi,bj))
     &        *uVelInt(i+1,j,bi,bj)*dyG(i+1,j,bi,bj)
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i+1,j,bi,bj))),1)

     &        - (EKE(i,j,1,bi,bj) + EKE(i,j-1,1,bi,bj))
     &        *vVelInt(i,j,bi,bj)*dxG(i,j,bi,bj)
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j-1,bi,bj))),1)

     &        - (EKE(i,j,1,bi,bj) + EKE(i-1,j,1,bi,bj))
     &        *uVelInt(i,j,bi,bj)*dyG(i,j,bi,bj)
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i-1,j,bi,bj))),1))

            gLambdaAdv(i,j,bi,bj) = (1/(4*abs(R_low(i,j,bi,bj))
     &        *rA(i,j,bi,bj)))

     &        *((Lambda(i,j,1,bi,bj) + Lambda(i,j+1,1,bi,bj))
     &        *vVelInt(i,j+1,bi,bj)*dxG(i,j+1,bi,bj)
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j+1,bi,bj))),1)

     &        + (Lambda(i,j,1,bi,bj) + Lambda(i+1,j,1,bi,bj))
     &        *uVelInt(i+1,j,bi,bj)*dyG(i+1,j,bi,bj)
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i+1,j,bi,bj))),1)

     &        - (Lambda(i,j,1,bi,bj) + Lambda(i,j-1,1,bi,bj))
     &        *vVelInt(i,j,bi,bj)*dxG(i,j,bi,bj)
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j-1,bi,bj))),1)

     &        - (Lambda(i,j,1,bi,bj) + Lambda(i-1,j,1,bi,bj))
     &        *uVelInt(i,j,bi,bj)*dyG(i,j,bi,bj)
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i-1,j,bi,bj))),1))



C--------------------------------------------------------------
C     CALCULATE EKE & Lambda LINEAR DRAG
C--------------------------------------------------------------

            gEKEr(i,j,bi,bj) = r_EKE(i,j,bi,bj)
     &        *EKE(i,j,1,bi,bj)

            gLambdar(i,j,bi,bj) = r_Lambda(i,j,bi,bj)
     &        *Lambda(i,j,1,bi,bj)

C--------------------------------------------------------------
C     CALCULATE K AND LAMBDA DIFFUSION    
C--------------------------------------------------------------

C if one or both of R_low = 0 then product = 0
C if neither R_low are zero then product > 0 => ceiling >= 1
C min => get zero when one or both of R_low is zero and one when neither is zero.

            gEKEDiff(i,j,bi,bj) = (mu(i,j,bi,bj)/(rA(i,j,bi,bj)
     &        *abs(R_low(i,j,bi,bj))))
     &        *((EKE(i,j+1,1,bi,bj)*abs(R_low(i,j+1,bi,bj))
     &        - EKE(i,j,1,bi,bj)*abs(R_low(i,j,bi,bj)))
     &        *min(ceiling(abs(R_low(i,j+1,bi,bj)*R_low(i,j,bi,bj))),1)
     &        +(EKE(i+1,j,1,bi,bj)*abs(R_low(i+1,j,bi,bj))
     &        - EKE(i,j,1,bi,bj)*abs(R_low(i,j,bi,bj)))
     &        *min(ceiling(abs(R_low(i+1,j,bi,bj)*R_low(i,j,bi,bj))),1)
     &        -(EKE(i,j,bi,1,bj)*abs(R_low(i,j,bi,bj))
     &        - EKE(i,j-1,1,bi,bj)*abs(R_low(i,j-1,bi,bj)))
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j-1,bi,bj))),1)
     &        -(EKE(i,j,1,bi,bj)*abs(R_low(i,j,bi,bj))
     &        - EKE(i-1,j,1,bi,bj)*abs(R_low(i-1,j,bi,bj)))
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i-1,j,bi,bj))),1))

            gLambdaDiff(i,j,bi,bj) = (mu(i,j,bi,bj)/(rA(i,j,bi,bj)
     &        *abs(R_low(i,j,bi,bj))))
     &        *((Lambda(i,j+1,1,bi,bj)*abs(R_low(i,j+1,bi,bj))
     &        - Lambda(i,j,1,bi,bj)*abs(R_low(i,j,bi,bj)))
     &        *min(ceiling(abs(R_low(i,j+1,bi,bj)*R_low(i,j,bi,bj))),1)
     &        +(Lambda(i+1,j,1,bi,bj)*abs(R_low(i+1,j,bi,bj))
     &        - Lambda(i,j,1,bi,bj)*abs(R_low(i,j,bi,bj)))
     &        *min(ceiling(abs(R_low(i+1,j,bi,bj)*R_low(i,j,bi,bj))),1)
     &        -(Lambda(i,j,1,bi,bj)*abs(R_low(i,j,bi,bj))
     &        - Lambda(i,j-1,1,bi,bj)*abs(R_low(i,j-1,bi,bj)))
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j-1,bi,bj))),1)
     &        -(Lambda(i,j,1,bi,bj)*abs(R_low(i,j,bi,bj))
     &        - Lambda(i-1,j,1,bi,bj)*abs(R_low(i-1,j,bi,bj)))
     &        *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i-1,j,bi,bj))),1))

          ELSE
            gEKEAdv(i,j,bi,bj)=0. 
            gLambdaAdv(i,j,bi,bj)=0.
            gEKEr(i,j,bi,bj)=0.
            gLambdar(i,j,bi,bj)=0.
            gEKEDiff(i,j,bi,bj)=0.
            gLambdaDiff(i,j,bi,bj)=0. 

          END IF


            
        ENDDO
      ENDDO



C--------------------------------------------------------------
C     CALCULATE ENERGY CONVERSION AND ENSTROPHY GENERATION
C--------------------------------------------------------------

      DO j=1,sNy
        DO i=1,sNx

C     INTERPOLATE K & Lambda TO YGXC
          EKE_YGXC(i,j,bi,bj) = max((EKE(i,j,1,bi,bj) 
     &      + EKE(i,j-1,1,bi,bj))/2,0.)

          Lambda_YGXC(i,j,bi,bj) = max((Lambda(i,j,1,bi,bj) 
     &      + Lambda(i,j-1,1,bi,bj))/2,0.)

          EKE_YCXG(i,j,bi,bj) = max((EKE(i,j,1,bi,bj) 
     &      + EKE(i-1,j,1,bi,bj))/2,0.)

          Lambda_YCXG(i,j,bi,bj) = max((Lambda(i,j,1,bi,bj) 
     &      + Lambda(i-1,j,1,bi,bj))/2,0.)

C     MOG(GRAD(Q)), KAPPA_Q, BAR(Q'U'), BAR(Q'U')*DPSIDX (ENERGY CONVERSION) 
C     and BAR(Q'U')*DQDX (ENSTROPHY GENERATION)AT YGXC
C     no more NaNs in first time step

          kappa_q_YGXC(i,j,bi,bj) = (2*gamma_q(i,j,bi,bj)
     &      *sqrt(EKE_YGXC(i,j,bi,bj)*Lambda_YGXC(i,j,bi,bj)))
     &      /modGradq_YGXC(i,j,bi,bj)
     &      *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i,j-1,bi,bj))),1)
          
          quEDDY_YGXC(i,j,bi,bj) = -kappa_q_YGXC(i,j,bi,bj)
     &      *dqdx_YGXC(i,j,bi,bj)

          energyConv_YGXC(i,j,bi,bj) = quEDDY_YGXC(i,j,bi,bj)
     &      *vVelInt(i,j,bi,bj)

          enstrophyGen_YGXC(i,j,bi,bj) = quEDDY_YGXC(i,j,bi,bj)
     &      *dqdx_YGXC(i,j,bi,bj)

C     MOG(GRAD(Q)), KAPPA_Q, BAR(Q'V') AND BAR(Q'V')*DPSIDY (ENERGY CONVERSION) 
C     and BAR(Q'V')*DQDY (ENSTROPHY GENERATION)AT YCXG

          kappa_q_YCXG(i,j,bi,bj) = (2*gamma_q(i,j,bi,bj)
     &      *sqrt(EKE_YCXG(i,j,bi,bj)*Lambda_YCXG(i,j,bi,bj)))
     &      /modGradq_YCXG(i,j,bi,bj)
     &      *min(ceiling(abs(R_low(i,j,bi,bj)*R_low(i-1,j,bi,bj))),1)

          qvEDDY_YCXG(i,j,bi,bj) = -kappa_q_YCXG(i,j,bi,bj)
     &      *dqdy_YCXG(i,j,bi,bj)

          energyConv_YCXG(i,j,bi,bj) = -qvEDDY_YCXG(i,j,bi,bj)
     &      *uVelInt(i,j,bi,bj)

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
C     CALCULATE EKE SOURCE TERM 
C--------------------------------------------------------------
C Calculate b
      DO k=1,Nr
        DO j=1,sNy
          DO i=1,sNx
            b(i,j,k,bi,bj) = gravity*(rhoInSitu(i,j,k,bi,bj) 
     &        - rhoConst)/rhoConst
          ENDDO
        ENDDO
      ENDDO
     
C calculate dbdx, dbdy, dbdz except at top and bottom
      DO k=2,Nr-1
        DO j=1,sNy
          DO i=1,sNx
            dbdx(i,j,k,bi,bj) = (b(i+1,j,k,bi,bj) - 
     &        b(i-1,j,k,bi,bj))/(2*dxC(i,j,bi,bj))
            dbdy(i,j,k,bi,bj) = (b(i,j+1,k,bi,bj) - 
     &        b(i,j-1,k,bi,bj))/(2*dyC(i,j,bi,bj))
            dbdz(i,j,k,bi,bj) = (b(i,j,k+1,bi,bj) - 
     &        b(i,j,k-1,bi,bj))/(drC(k) + drC(k+1))
          ENDDO
        ENDDO
      ENDDO
C calculate for k = 1 and k = Nr
      DO j=1,sNy
        DO i=1,sNx
          dbdx(i,j,1,bi,bj) = (b(i+1,j,1,bi,bj) - 
     &      b(i-1,j,1,bi,bj))/(2*dxC(i,j,bi,bj))
          dbdx(i,j,Nr,bi,bj) = (b(i+1,j,Nr,bi,bj) - 
     &      b(i-1,j,Nr,bi,bj))/(2*dxC(i,j,bi,bj))

          dbdy(i,j,1,bi,bj) = (b(i,j+1,1,bi,bj) - 
     &      b(i,j-1,1,bi,bj))/(2*dyC(i,j,bi,bj))
          dbdy(i,j,Nr,bi,bj) = (b(i,j+1,Nr,bi,bj) - 
     &      b(i,j-1,Nr,bi,bj))/(2*dyC(i,j,bi,bj))

          dbdz(i,j,1,bi,bj) = (b(i,j,2,bi,bj) - 
     &      b(i,j,1,bi,bj))/(drC(1))

          dbdz(i,j,Nr,bi,bj) = (b(i,j,Nr,bi,bj) - 
     &      b(i,j,Nr-1,bi,bj))/(drC(Nr))
        ENDDO
      ENDDO

C calculate |grad b| and flux = bar(b'w')
      DO k=1,Nr
        DO j=1,sNy
          DO i=1,sNx
            grad_b(i,j,k,bi,bj) = SQRT(dbdx(i,j,k,bi,bj)**2
     &        + dbdy(i,j,k,bi,bj)**2)
            flux(i,j,k,bi,bj) = (K*grad_b(i,j,k,bi,bj))/
     &        dbdz(i,j,k,bi,bj)
          ENDDO
        ENDDO
      ENDDO

C calculate depth integral of flux
      DO k=1,Nr
        DO j=1,sNy
          DO i=1,sNx
            flux_zInt(i,j,bi,bj) = flux_zInt(i,j,bi,bj) + 
     &        flux(i,j,k,bi,bj)*drF(k)*hFacC(i,j,k,bi,bj)
          ENDDO
        ENDDO
      ENDDO

C calculate EKE source term = depth average of flux term
      DO j=1,sNy
        DO i=1,sNx
          gEKESource(i,j,bi,bj) = flux_zInt(i,j,bi,bj)/
     &      abs(R_low(i,j,bi,bj))
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
     &      + gEKESource(i,j,bi,bj)

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
          EKE(i,j,1,bi,bj) = max(EKE(i,j,1,bi,bj)
     &     + deltaTMom*( gEKE(i,j,1,bi,bj) ),0.)

          Lambda(i,j,1,bi,bj) = max(Lambda(i,j,1,bi,bj)
     &     + deltaTMom*( gLambda(i,j,1,bi,bj) ),0.)

        ENDDO
      ENDDO


      RETURN
      END

