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
#include "GMREDI.h"

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
C    VARIABLES FOR SOURCE TERM
C      _RL  b(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL  dbdx(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL  dbdy(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL  dbdz(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL  grad_b(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL  flux(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr,nSx,nSy)
C      _RL  flux_zInt(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
C      _RL  kappa_GM(1-OLx:sNx+OLx,1-OLy:sNy+OLy)
      _RL  partition
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

C SECOND EKE SOURCE METHOD
       _RL  kGM_x_YCXGZG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr+1,nSx,nSy)
       _RL  kGM_y_YGXCZG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr+1,nSx,nSy)
       _RL  fluxX2_YCXGZG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr+1,nSx,nSy)
       _RL  fluxY2_YGXCZG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,Nr+1,nSx,nSy)
       _RL  fluxXInt2_YCXG(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
       _RL  fluxYInt2_YGXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
       _RL  fluxInt2_YCXC(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)
       _RL  gEKESource2(1-OLx:sNx+OLx,1-OLy:sNy+OLy,nSx,nSy)


C    Initialize arrays with zero values
      DO j=1-OLy,sNy+OLy
        DO i=1-OLx,sNx+OLx
          gEKEAdv(i,j,bi,bj) = 0. _d 0
          gEKEConv(i,j,bi,bj) = 0. _d 0
          gEKEr(i,j,bi,bj) = 0. _d 0
          gEKEDiff(i,j,bi,bj) = 0. _d 0
          gEKESource(i,j,bi,bj) = 0. _d 0
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

          fluxXInt2_YCXG(i,j,bi,bj) = 0. _d 0
          fluxYInt2_YGXC(i,j,bi,bj) = 0. _d 0
          fluxInt2_YCXC(i,j,bi,bj) = 0. _d 0
          gEKEsource2(i,j,bi,bj) = 0. _d 0

        ENDDO
      ENDDO 

      DO k=1,Nr
        DO j=1,sNy
          DO i=1,sNx
            b_ZC(i,j,k,bi,bj) = 0. _d 0
            dbdx_YCXGZC(i,j,k,bi,bj) = 0. _d 0
            dbdy_YGXCZC(i,j,k,bi,bj) = 0. _d 0
            kGM_x_YCXGZC(i,j,k,bi,bj) = 0. _d 0
            kGM_y_YGXCZC(i,j,k,bi,bj) = 0. _d 0
            numx_YCXGZC(i,j,k,bi,bj) = 0. _d 0
            numy_YGXCZC(i,j,k,bi,bj) = 0. _d 0
            numx_YCXGZG(i,j,k,bi,bj) = 0. _d 0
            numy_YGXCZG(i,j,k,bi,bj) = 0. _d 0
            dbdz_YCXCZG(i,j,k,bi,bj) = 0. _d 0
            dbdz_YCXGZG(i,j,k,bi,bj) = 0. _d 0
            dbdz_YGXCZG(i,j,k,bi,bj) = 0. _d 0
            flux_YCXGZG(i,j,k,bi,bj) = 0. _d 0
            flux_YGXCZG(i,j,k,bi,bj) = 0. _d 0

            kGM_x_YCXGZG(i,j,k,bi,bj) = 0. _d 0
            kGM_y_YGXCZG(i,j,k,bi,bj) = 0. _d 0
            fluxX2_YCXGZG(i,j,k,bi,bj) = 0. _d 0
            fluxY2_YGXCZG(i,j,k,bi,bj) = 0. _d 0
           ENDDO
        ENDDO
      ENDDO

      DO j=1,sNy
        DO i=1,sNx
          numx_YCXGZG(i,j,Nr+1,bi,bj) = 0. _d 0
          numy_YGXCZG(i,j,Nr+1,bi,bj) = 0. _d 0
          dbdz_YCXCZG(i,j,Nr+1,bi,bj) = 0. _d 0
          dbdz_YCXGZG(i,j,Nr+1,bi,bj) = 0. _d 0
          dbdz_YGXCZG(i,j,Nr+1,bi,bj) = 0. _d 0
          flux_YCXGZG(i,j,Nr+1,bi,bj) = 0. _d 0
          flux_YGXCZG(i,j,Nr+1,bi,bj) = 0. _d 0
          fluxInt_YCXG(i,j,bi,bj) = 0. _d 0
          fluxInt_YGXC(i,j,bi,bj) = 0. _d 0
          fluxInt_YCXC(i,j,bi,bj) = 0. _d 0
        ENDDO
      ENDDO

      partition = 1.

C--------------------------------l------------------------------
C     DEPTH INTEGRAL OF U AND V FOR ENERGY CONVERSION
C     AND ADVECTION
C--------------------------------------------------------------
      DO k=1,Nr
        DO j=1,sNy
          DO i=1,sNx
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
            b_ZC(i,j,k,bi,bj) = -gravity*
     &        (rhoInSitu(i,j,k,bi,bj))/rhoConst
          ENDDO
        ENDDO
      ENDDO

      DO k=1,Nr
        DO j=1,sNy+1
          DO i=1,sNx+1
C calculate dbdx, dbdy at cell faces at ZC
C set to zero at partial cells
            dbdx_YCXGZC(i,j,k,bi,bj) = ((b_ZC(i,j,k,bi,bj) 
     &        - b_ZC(i-1,j,k,bi,bj))/dxC(i,j,bi,bj))
     &        *floor(abs(hFacC(i,j,k,bi,bj)))
     &        *floor(abs(hFacC(i-1,j,k,bi,bj)))

            dbdy_YGXCZC(i,j,k,bi,bj) = ((b_ZC(i,j,k,bi,bj) 
     &        - b_ZC(i,j-1,k,bi,bj))/dyC(i,j,bi,bj))
     &        *floor(abs(hFacC(i,j,k,bi,bj)))
     &        *floor(abs(hFacC(i,j-1,k,bi,bj)))

C calculate kGM_x_ZC and kGM_y_ZC
C set to zero if either side is a partial cell
            kGM_x_YCXGZC(i,j,k,bi,bj) = GM_background_K
     &       *taperX(i,j,k,bi,bj)
     &       *floor(abs(hFacC(i,j,k,bi,bj)))
     &       *floor(abs(hFacC(i-1,j,k,bi,bj)))

            kGM_y_YGXCZC(i,j,k,bi,bj) = GM_background_K
     &       *taperY(i,j,k,bi,bj)
     &       *floor(abs(hFacC(i,j,k,bi,bj)))
     &       *floor(abs(hFacC(i,j-1,k,bi,bj)))
          ENDDO
        ENDDO
      ENDDO

      WRITE(standardMessageUnit,*)
     &      'dbdx_ZC:  = ', dbdx_YCXGZC
      WRITE(standardMessageUnit,*)
     &      'dbdy_ZC:  = ', dbdy_YGXCZC

      WRITE(standardMessageUnit,*)
     &      'dsigmadx:  = ', BT_dsigmadx

      WRITE(standardMessageUnit,*)
     &      'dsigmady:  = ', BT_dsigmady

      WRITE(standardMessageUnit,*)
     &      'kGM_x_YCXGZC:  = ', kGM_x_YCXGZC
      WRITE(standardMessageUnit,*)
     &      'kGM_y_YGXCZC:  = ', kGM_y_YGXCZC


C caluclate numerator of flux terms:
C multiply kGM with dbdx^2 and dbdy^2
      DO k=1,Nr
        DO j=1,sNy+1
          DO i=1,sNx+1
            numx_YCXGZC(i,j,k,bi,bj) = kGM_x_YCXGZC(i,j,k,bi,bj)*
     &        (dbdx_YCXGZC(i,j,k,bi,bj)**2)

            numy_YGXCZC(i,j,k,bi,bj) = kGM_y_YGXCZC(i,j,k,bi,bj)*
     &        (dbdy_YGXCZC(i,j,k,bi,bj)**2)
          ENDDO
        ENDDO
      ENDDO


C interpolate numerator of flux terms to ZG points
C don't calculate at top and bottom (should be ero there) 

      DO k=2,Nr
        DO j=1,sNy+1
          DO i=1,sNx+1
            numx_YCXGZG(i,j,k,bi,bj) = numx_YCXGZC(i,j,k-1,bi,bj) + 
     &        0.5*drF(k-1)*((numx_YCXGZC(i,j,k,bi,bj) 
     &        - numx_YCXGZC(i,j,k-1,bi,bj))/drC(k))

            numy_YGXCZG(i,j,k,bi,bj) = numy_YGXCZC(i,j,k-1,bi,bj) + 
     &        0.5*drF(k-1)*((numy_YGXCZC(i,j,k,bi,bj) 
     &        - numy_YGXCZC(i,j,k-1,bi,bj))/drC(k))
          ENDDO
        ENDDO
      ENDDO

C calculate dbdz at vertical cell interfaces
C don't calculate at top and bottom
C min(0.5,hFacC(k)) is for partial cells at bottom
C floor(abc(hFacC(i,j,k-1,bi,bj))) is for below seafloor
C dbdz_ZG will be zero below the seafloor

      DO k=2,Nr
        DO j=1,sNy+1
          DO i=1,sNx+1
            dbdz_YCXCZG(i,j,k,bi,bj) = max(
     &        floor(abs(hFacC(i,j,k-1,bi,bj)))*
     &        (b_ZC(i,j,k-1,bi,bj) - b_ZC(i,j,k,bi,bj))
     &        /(0.5*drF(k-1) + 
     &        0.5*hFacC(i,j,k,bi,bj)*drF(k)),
     &        1.E-7)
          ENDDO
        ENDDO
      ENDDO

      WRITE(standardMessageUnit,*)
     &      'dbdz_YCXC:  = ', dbdz_YCXCZG

C interpolate dbdz to YGXC and YCXG
C don't need to set to zero at lateral boundaries
C because dbdx and dbdy values will be zero there

      DO k=2,Nr
        DO j=1,sNy+1
          DO i=1,sNx+1
            dbdz_YCXGZG(i,j,k,bi,bj) = (dbdz_YCXCZG(i,j,k,bi,bj)
     &        + dbdz_YCXCZG(i-1,j,k,bi,bj))/2

            dbdz_YGXCZG(i,j,k,bi,bj) = (dbdz_YCXCZG(i,j,k,bi,bj)
     &        + dbdz_YCXCZG(i,j-1,k,bi,bj))/2
          ENDDO
        ENDDO
      ENDDO

      WRITE(standardMessageUnit,*)
     &      'dbdz_YCXG  = ', dbdz_YCXGZG
     
      WRITE(standardMessageUnit,*)
     &      'dsigmadzW:  = ', BT_dsigmadzW

      WRITE(standardMessageUnit,*)
     &      'dbdz_YGXC  = ', dbdz_YGXCZG

      WRITE(standardMessageUnit,*)
     &      'dsigmadzS:  = ', BT_dsigmadzS
C calculate x and y componenets of flux

      DO k=2,Nr
        DO j=1,sNy+1
          DO i=1,sNx+1
C           x-component 
            flux_YCXGZG(i,j,k,bi,bj) = numx_YCXGZG(i,j,k,bi,bj)/
     &        dbdz_YCXGZG(i,j,k,bi,bj)

C           y-component
            flux_YGXCZG(i,j,k,bi,bj) = numy_YGXCZG(i,j,k,bi,bj)/
     &        dbdz_YGXCZG(i,j,k,bi,bj)
            
          ENDDO
        ENDDO
      ENDDO

C calculate depth integral of x and y components of flux
      DO k=1,Nr+1
        DO j=1,sNy+1
          DO i=1,sNx+1
            fluxInt_YCXG(i,j,bi,bj) = 
     &        fluxInt_YCXG(i,j,bi,bj) + 
     &        flux_YCXGZG(i,j,k,bi,bj)*drC(k)

             fluxInt_YGXC(i,j,bi,bj) = 
     &        fluxInt_YGXC(i,j,bi,bj) + 
     &        flux_YGXCZG(i,j,k,bi,bj)*drC(k)
          ENDDO
        ENDDO
      ENDDO

C interpolate flux to YCXC

      DO j=1,sNy
          DO i=1,sNx
            fluxInt_YCXC(i,j,bi,bj) = (
     &        fluxInt_YCXG(i,j,bi,bj) + 
     &        fluxInt_YCXG(i+1,j,bi,bj) +
     &        fluxInt_YGXC(i,j,bi,bj) + 
     &        fluxInt_YGXC(i,j+1,bi,bj))/4
          ENDDO
        ENDDO
        
C calculate EKE source term = depth average of flux term
      DO j=1,sNy
        DO i=1,sNx
          IF (abs(R_low(i,j,bi,bj)).NE.0) THEN
            gEKESource(i,j,bi,bj) = 
     &        partition*fluxInt_YCXC(i,j,bi,bj)/
     &        abs(R_low(i,j,bi,bj))  
          ELSE 
            gEKESource(i,j,bi,bj) = 0.
          END IF
        ENDDO
      ENDDO

C-------------------------------------
C SECOND METHOD
C-------------------------------------

C interpolate kappa_GM to ZG
      DO k=2,Nr
        DO j=1,sNy+1
          DO i=1,sNx+1
            kGM_x_YCXGZG(i,j,k,bi,bj) = kGM_x_YCXGZC(i,j,k-1,bi,bj) + 
     &        0.5*drF(k-1)*((kGM_x_YCXGZC(i,j,k,bi,bj) 
     &        - kGM_x_YCXGZC(i,j,k-1,bi,bj))/drC(k))

            kGM_y_YGXCZG(i,j,k,bi,bj) = kGM_y_YGXCZC(i,j,k-1,bi,bj) + 
     &        0.5*drF(k-1)*((kGM_y_YGXCZC(i,j,k,bi,bj) 
     &        - kGM_y_YGXCZC(i,j,k-1,bi,bj))/drC(k))
          ENDDO
        ENDDO
      ENDDO

      WRITE(standardMessageUnit,*)
     &      'kGM_x_YCXGZG:  = ', kGM_x_YCXGZG

      WRITE(standardMessageUnit,*)
     &      'kGM_y_YGXCZG:  = ', kGM_y_YGXCZG

      DO k=2,Nr
        DO j=1,sNy+1
          DO i=1,sNx+1
C            x-component 
            fluxX2_YCXGZG(i,j,k,bi,bj) = kGM_x_YCXGZG(i,j,k,bi,bj)*
     &        (BT_dsigmadx(i,j,k,bi,bj)**2)/BT_dsigmadzW(i,j,k,bi,bj)

C           y-component
            fluxY2_YGXCZG(i,j,k,bi,bj) = kGM_y_YGXCZG(i,j,k,bi,bj)*
     &        (BT_dsigmady(i,j,k,bi,bj)**2)/BT_dsigmadzS(i,j,k,bi,bj)
          ENDDO
        ENDDO
      ENDDO

      WRITE(standardMessageUnit,*)
     &      'fluxX_YCXGZG:  = ', flux_YCXGZG
     
      WRITE(standardMessageUnit,*)
     &      'fluxY_YGXCZG:  = ', flux_YGXCZG

      WRITE(standardMessageUnit,*)
     &      'fluxX2_YCXGZG:  = ', fluxX2_YCXGZG
     
      WRITE(standardMessageUnit,*)
     &      'fluxY2_YGXCZG:  = ', fluxY2_YGXCZG

C calculate depth integral of x and y components of flux
      DO k=1,Nr+1
        DO j=1,sNy+1
          DO i=1,sNx+1
            fluxXInt2_YCXG(i,j,bi,bj) = 
     &        fluxXInt2_YCXG(i,j,bi,bj) + 
     &        fluxX2_YCXGZG(i,j,k,bi,bj)*drC(k)

             fluxYInt2_YGXC(i,j,bi,bj) = 
     &        fluxYInt2_YGXC(i,j,bi,bj) + 
     &        fluxY2_YGXCZG(i,j,k,bi,bj)*drC(k)
          ENDDO
        ENDDO
      ENDDO

      WRITE(standardMessageUnit,*)
     &      'fluxXInt_YCXG:  = ', fluxInt_YCXG
     
      WRITE(standardMessageUnit,*)
     &      'fluxYInt_YGXC:  = ', fluxInt_YGXC

      WRITE(standardMessageUnit,*)
     &      'fluxXInt2_YCXG:  = ', fluxXInt2_YCXG
     
      WRITE(standardMessageUnit,*)
     &      'fluxYInt2_YGXC:  = ', fluxYInt2_YGXC

C interpolate flux to YCXC

      DO j=1,sNy
        DO i=1,sNx
          fluxInt2_YCXC(i,j,bi,bj) = (
     &      fluxXInt2_YCXG(i,j,bi,bj) + 
     &      fluxXInt2_YCXG(i+1,j,bi,bj) +
     &      fluxYInt2_YGXC(i,j,bi,bj) + 
     &      fluxYInt2_YGXC(i,j+1,bi,bj))/4
        ENDDO
      ENDDO

      WRITE(standardMessageUnit,*)
     &      'fluxInt_YCXC:  = ', fluxInt_YCXC
     
      WRITE(standardMessageUnit,*)
     &      'fluxInt2_YCXC:  = ', fluxInt2_YCXC

C calculate EKE source term = depth average of flux term
      DO j=1,sNy
        DO i=1,sNx
          IF (abs(R_low(i,j,bi,bj)).NE.0) THEN
            gEKESource2(i,j,bi,bj) = 
     &        partition*fluxInt2_YCXC(i,j,bi,bj)/
     &        abs(R_low(i,j,bi,bj))  
          ELSE 
            gEKESource2(i,j,bi,bj) = 0.
          END IF
        ENDDO
      ENDDO

      WRITE(standardMessageUnit,*)
     &      'gEKEcSource:  = ', gEKEsource

      WRITE(standardMessageUnit,*)
     &      'gEKESource2:  = ', gEKEsource2

C-------------------------------------
C SECOND METHOD
C-------------------------------------

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

C--------------------------------------------------------------
C     CALL DIAGNOSTICS
C--------------------------------------------------------------

C       custom: parameterized EKE
      CALL DIAGNOSTICS_FILL(EKE,'EKE     ',1,1,2,bi,bj,myThid)
C       custom: parameterized EKE tendency
      CALL DIAGNOSTICS_FILL(gEKE,'gEKE    ',1,1,2,bi,bj,myThid)
C      custom: parameterized EKE conversion
      CALL DIAGNOSTICS_FILL(gEKEConv,'gEKEConv',1,1,2,bi,bj,myThid)
C       custom: parameterized EKE advection
      CALL DIAGNOSTICS_FILL(gEKEAdv,'gEKEAdv ',1,1,2,bi,bj,myThid)
C       custom: parameterized EKE drag
      CALL DIAGNOSTICS_FILL(gEKEr,'gEKEr   ',1,1,2,bi,bj,myThid)
C       custom: parameterized EKE diffusion
      CALL DIAGNOSTICS_FILL(gEKEDiff,'gEKEDiff',1,1,2,bi,bj,myThid)
C     SOURCE TERM
C       custom: parameterized EKE source term
      CALL DIAGNOSTICS_FILL(gEKESource,'gEKESour',1,1,2,bi,bj,myThid)
C       custom: parameterized Lambda
      CALL DIAGNOSTICS_FILL(Lambda,'Lambda  ',1,1,2,bi,bj,myThid)
C       custom: parameterized Lambda tendency
      CALL DIAGNOSTICS_FILL(gLambda,'gLambda ',1,1,2,bi,bj,myThid)
C       custom: parameterized Lambda generation
      CALL DIAGNOSTICS_FILL(gLambdaGen,'gLambdaGen',
     &     1,1,2,bi,bj,myThid)
C       custom: parameterized Lambda advection
      CALL DIAGNOSTICS_FILL(gLambdaAdv,'gLambdaAdv',
     &     1,1,2,bi,bj,myThid)
C       custom: parameterized Lambda drag
      CALL DIAGNOSTICS_FILL(gLambdar,'gLambdar',
     &     1,1,2,bi,bj,myThid)
C       custom: parameterized Lambda diffusion
      CALL DIAGNOSTICS_FILL(gLambdaDiff,'gLambdaDiff',
     &     1,1,2,bi,bj,myThid)
C       custom: parameterized bar(q'u')
      CALL DIAGNOSTICS_FILL(quEDDY_YGXC,'quEDDY  ',
     &     1,1,2,bi,bj,myThid)

C       custom: parameterized bar(q'v')
      CALL DIAGNOSTICS_FILL(qvEDDY_YCXG,'qvEDDY  ',
     &     1,1,2,bi,bj,myThid)
C       custom: fluxInt_YCXG
        CALL DIAGNOSTICS_FILL(fluxInt_YCXG,'fluxXInt',
     &     1,1,2,bi,bj,myThid)
C       custom: fluxInt_YGXC
        CALL DIAGNOSTICS_FILL(fluxInt_YGXC,'fluxYInt',
     &     1,1,2,bi,bj,myThid)
C       custom: fluxInt_YCXC
        CALL DIAGNOSTICS_FILL(fluxInt_YCXC,'fluxInt ',
     &     1,1,2,bi,bj,myThid)



C       custom: qYCXC
      CALL DIAGNOSTICS_FILL(q_YCXC,'qYCXC   ',1,1,2,bi,bj,myThid)
C       custom: dqdxYCXG
      CALL DIAGNOSTICS_FILL(dqdx_YCXG,'dqdxYCXG',1,1,2,bi,bj,myThid)
C       custom: dqdyYGXC
      CALL DIAGNOSTICS_FILL(dqdy_YGXC,'dqdyYGXC',1,1,2,bi,bj,myThid)
C       custom: dqdxYGXC
      CALL DIAGNOSTICS_FILL(dqdx_YGXC,'dqdxYGXC',1,1,2,bi,bj,myThid)
C       custom: dqdyYCXG
      CALL DIAGNOSTICS_FILL(dqdy_YCXG,'dqdyYCXG',1,1,2,bi,bj,myThid)
C       custom: modGradq_YCXG
      CALL DIAGNOSTICS_FILL(modGradq_YCXG,'modqYCXG',1,1,2,bi,bj,myThid)
C       custom: modGradq_YGXC
      CALL DIAGNOSTICS_FILL(modGradq_YGXC,'modqYGXC',1,1,2,bi,bj,myThid)



      

C     SOURCE TERM TESTING
      DO k = 1,Nr
C       custom: parameterized buoyancy
        CALL DIAGNOSTICS_FILL(b_ZC(:,:,k,:,:),'b       ',
     &     k,1,2,bi,bj,myThid)
C       custom: dbdx_YCXGZC
        CALL DIAGNOSTICS_FILL(dbdx_YCXGZC(:,:,k,:,:),'dbdxZC  ',
     &     k,1,2,bi,bj,myThid)
C       custom: dbdy_YGXCZC
        CALL DIAGNOSTICS_FILL(dbdy_YGXCZC(:,:,k,:,:),'dbdyZC  ',
     &     k,1,2,bi,bj,myThid)
C       custom: kGM_x_YCXGZC
        CALL DIAGNOSTICS_FILL(kGM_x_YCXGZC(:,:,k,:,:),'kGMxZC  ',
     &     k,1,2,bi,bj,myThid)
C       custom: kGM_y_YGXCZC
        CALL DIAGNOSTICS_FILL(kGM_y_YGXCZC(:,:,k,:,:),'kGMyZC  ',
     &     k,1,2,bi,bj,myThid)
C       custom: num_YCXGZC
        CALL DIAGNOSTICS_FILL(numx_YCXGZC(:,:,k,:,:),'numxZC  ',
     &     k,1,2,bi,bj,myThid)
C       custom: num_YGXCZC
        CALL DIAGNOSTICS_FILL(numy_YGXCZC(:,:,k,:,:),'numyZC  ',
     &     k,1,2,bi,bj,myThid)
C       custom: num_YCXGZG
        CALL DIAGNOSTICS_FILL(numx_YCXGZG(:,:,k,:,:),'numxZG  ',
     &     k,1,2,bi,bj,myThid)
C       custom: num_YGXCZG
        CALL DIAGNOSTICS_FILL(numy_YGXCZG(:,:,k,:,:),'numyZG  ',
     &     k,1,2,bi,bj,myThid)
C       custom: dbdz_YCXCZG
        CALL DIAGNOSTICS_FILL(dbdz_YCXCZG(:,:,k,:,:),'dbdzYCXC',
     &     k,1,2,bi,bj,myThid)
C       custom: dbdz_YCXGZG
        CALL DIAGNOSTICS_FILL(dbdz_YCXGZG(:,:,k,:,:),'dbdzYCXG',
     &     k,1,2,bi,bj,myThid)
C       custom: dbdz_YGXCZG
        CALL DIAGNOSTICS_FILL(dbdz_YGXCZG(:,:,k,:,:),'dbdzYGXC',
     &     k,1,2,bi,bj,myThid)
C       custom: flux_YCXGZG
        CALL DIAGNOSTICS_FILL(flux_YCXGZG(:,:,k,:,:),'fluxX   ',
     &     k,1,2,bi,bj,myThid)
C       custom: flux_YGXCZG
        CALL DIAGNOSTICS_FILL(flux_YGXCZG(:,:,k,:,:),'fluxY   ',
     &     k,1,2,bi,bj,myThid)
C       custom: GM taperX
        CALL DIAGNOSTICS_FILL(taperX(:,:,k,:,:),'taperX  ',
     &     k,1,2,bi,bj,myThid)
C       custom: GM taperY
        CALL DIAGNOSTICS_FILL(taperY(:,:,k,:,:),'taperY  ',
     &     k,1,2,bi,bj,myThid)
C       custom: isoneutral density
        CALL DIAGNOSTICS_FILL(BT_isoRho(:,:,k,:,:),'isoRho  ',
     &     k,1,2,bi,bj,myThid)

      ENDDO


      RETURN
      END

