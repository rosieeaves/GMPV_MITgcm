#include "MOM_COMMON_OPTIONS.h"

CBOP
C     !ROUTINE: MOM_DIAGNOSTICS_INIT
C     !INTERFACE:
      SUBROUTINE MOM_DIAGNOSTICS_INIT( myThid )

C     !DESCRIPTION:
C     Routine to initialize Momentum diagnostics
C      (from mom_common, mom_fluxform & mom_vecinv packages)

C     !USES:
      IMPLICIT NONE
C     === Global variables ===
#include "SIZE.h"
#include "EEPARAMS.h"
#include "PARAMS.h"

C     !INPUT/OUTPUT PARAMETERS:
C     === Routine arguments ===
C     myThid :: my Thread Id. number
      INTEGER myThid
CEOP

#ifdef ALLOW_DIAGNOSTICS
#ifdef ALLOW_MOM_COMMON
C     !LOCAL VARIABLES:
C     === Local variables ===
C     msgBuf      :: Informational/error message buffer
c     CHARACTER*(MAX_LEN_MBUF) msgBuf

      INTEGER       diagNum
      INTEGER       diagMate
      CHARACTER*8   diagName
      CHARACTER*16  diagCode
      CHARACTER*16  diagUnits
      CHARACTER*(80) diagTitle

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

C---  Add diagnostics to the (long) list of available diagnostics:

C--   add diagnostics for viscosity coefficients:
      diagName  = 'VISCAHZ '
      diagTitle = 'Harmonic Visc Coefficient (m2/s) (Zeta Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SZ      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VISCA4Z '
      diagTitle = 'Biharmonic Visc Coefficient (m4/s) (Zeta Pt)'
      diagUnits = 'm^4/s           '
      diagCode  = 'SZ      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VISCAHD '
      diagTitle = 'Harmonic Viscosity Coefficient (m2/s) (Div Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VISCA4D '
      diagTitle = 'Biharmonic Viscosity Coefficient (m4/s) (Div Pt)'
      diagUnits = 'm^4/s           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VISCAHW '
      diagTitle = 'Harmonic Viscosity Coefficient (m2/s) (W Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'WM      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VISCA4W '
      diagTitle = 'Biharmonic Viscosity Coefficient (m4/s) (W Pt)'
      diagUnits = 'm^4/s           '
      diagCode  = 'WM      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VAHZMAX '
      diagTitle = 'CFL-MAX Harm Visc Coefficient (m2/s) (Zeta Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SZ      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VA4ZMAX '
      diagTitle = 'CFL-MAX Biharm Visc Coefficient (m4/s) (Zeta Pt)'
      diagUnits = 'm^4/s           '
      diagCode  = 'SZ      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VAHDMAX '
      diagTitle = 'CFL-MAX Harm Visc Coefficient (m2/s) (Div Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VA4DMAX '
      diagTitle = 'CFL-MAX Biharm Visc Coefficient (m4/s) (Div Pt)'
      diagUnits = 'm^4/s           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VAHZMIN '
      diagTitle = 'RE-MIN Harm Visc Coefficient (m2/s) (Zeta Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SZ      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VA4ZMIN '
      diagTitle = 'RE-MIN Biharm Visc Coefficient (m4/s) (Zeta Pt)'
      diagUnits = 'm^4/s           '
      diagCode  = 'SZ      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VAHDMIN '
      diagTitle = 'RE-MIN Harm Visc Coefficient (m2/s) (Div Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VA4DMIN '
      diagTitle = 'RE-MIN Biharm Visc Coefficient (m4/s) (Div Pt)'
      diagUnits = 'm^4/s           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VAHZLTH '
      diagTitle = 'Leith Harm Visc Coefficient (m2/s) (Zeta Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SZ      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VA4ZLTH '
      diagTitle = 'Leith Biharm Visc Coefficient (m4/s) (Zeta Pt)'
      diagUnits = 'm^4/s           '
      diagCode  = 'SZ      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VAHDLTH '
      diagTitle = 'Leith Harm Visc Coefficient (m2/s) (Div Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VA4DLTH '
      diagTitle = 'Leith Biharm Visc Coefficient (m4/s) (Div Pt)'
      diagUnits = 'm^4/s           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VAHZLTHD'
      diagTitle = 'LeithD Harm Visc Coefficient (m2/s) (Zeta Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SZ      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VA4ZLTHD'
      diagTitle = 'LeithD Biharm Visc Coefficient (m4/s) (Zeta Pt)'
      diagUnits = 'm^4/s           '
      diagCode  = 'SZ      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VAHDLTHD'
      diagTitle = 'LeithD Harm Visc Coefficient (m2/s) (Div Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VA4DLTHD'
      diagTitle = 'LeithD Biharm Visc Coefficient (m4/s) (Div Pt)'
      diagUnits = 'm^4/s           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#ifdef ALLOW_LEITH_QG
      diagName  = 'VAHZLTHQ'
      diagTitle = 'LeithQG Harm Visc Coefficient (m2/s) (Zeta Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SZ      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VAHDLTHQ'
      diagTitle = 'LeithQG Harm Visc Coefficient (m2/s) (Div Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif

      diagName  = 'VAHZSMAG'
      diagTitle = 'Smagorinsky Harm Visc Coefficient (m2/s) (Zeta Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SZ      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VA4ZSMAG'
      diagTitle = 'Smagorinsky Biharm Visc Coeff. (m4/s) (Zeta Pt)'
      diagUnits = 'm^4/s           '
      diagCode  = 'SZ      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VAHDSMAG'
      diagTitle = 'Smagorinsky Harm Visc Coefficient (m2/s) (Div Pt)'
      diagUnits = 'm^2/s           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'VA4DSMAG'
      diagTitle = 'Smagorinsky Biharm Visc Coeff. (m4/s) (Div Pt)'
      diagUnits = 'm^4/s           '
      diagCode  = 'SM      MR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C--   add diagnostics of simple derivative quantities:
      diagName  = 'momKE   '
      diagTitle = 'Kinetic Energy (in momentum Eq.)'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2/s^2         '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

#ifdef ALLOW_MOM_VECINV
      diagName  = 'momHDiv '
      diagTitle = 'Horizontal Divergence (in momentum Eq.)'
      diagCode  = 'SMR     MR      '
      diagUnits = 's^-1            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'momVort3'
      diagTitle = '3rd component (vertical) of Vorticity'
      diagCode  = 'SZR     MR      '
      diagUnits = 's^-1            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: relative vorticity squared

      diagName  = 'VORTSQ  '
      diagTitle = '3rd component (vertical) of Vorticity squared'
      diagCode  = 'SZR     MR      '
      diagUnits = 'ms^-2            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized EKE

      diagName  = 'EKE     '
      diagTitle = 'parameterized EKE'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-2            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized gEKE

      diagName  = 'gEKE    '
      diagTitle = 'parameterized EKE tendency'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized gEKEConv

      diagName  = 'gEKEConv'
      diagTitle = 'parameterized EKE conversion'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized gEKEAdv

      diagName  = 'gEKEAdv '
      diagTitle = 'parameterized EKE advection'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized gEKEDiff

      diagName  = 'gEKEDiff'
      diagTitle = 'parameterized EKE diffusion'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized gEKEr

      diagName  = 'gEKEr   '
      diagTitle = 'parameterized EKE linear drag'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized Lambda

      diagName  = 'Lambda  '
      diagTitle = 'parameterized Lambda'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-2            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized gLambda

      diagName  = 'gLambda '
      diagTitle = 'parameterized Lambda tendency'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized gLambdaGen

      diagName  = 'gLambdaGen'
      diagTitle = 'parameterized Lambda generation'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized gLambdaAdv

      diagName  = 'gLambdaAdv'
      diagTitle = 'parameterized Lambda advection'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized gLambdaDiff

      diagName  = 'gLambdaDiff'
      diagTitle = 'parameterized Lambda diffusion'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized gLambdar

      diagName  = 'gLambdar'
      diagTitle = 'parameterized Lambda linear drag'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized bar(q'u')

      diagName  = 'quEDDY  '
      diagTitle = 'quEDDY'
      diagCode  = 'SMR     MR      '
      diagUnits = 's^-2            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C     custom: parameterized bar(q'v')

      diagName  = 'qvEDDY  '
      diagTitle = 'qvEDDY'
      diagCode  = 'SMR     MR      '
      diagUnits = 's^-2            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )


C --------------------------------------------------------------------
C SOURCE TERM
C --------------------------------------------------------------------

C     custom: parameterized gEKESource

      diagName  = 'gEKESour'
      diagTitle = 'parameterized EKE source term'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C --------------------------------------------------------------------
C SOURCE TERM TESTING
C --------------------------------------------------------------------

      diagName  = 'b'
      diagTitle = 'BUOYANCY'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'dbdxYCXG'
      diagTitle = 'dbdx_YCXG'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'dbdxYCXC'
      diagTitle = 'dbdx_YCXC'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'dbdyYGXC'
      diagTitle = 'dbdy_YGXC'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'dbdyYCXC'
      diagTitle = 'dbdy_YCXC'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'dbdz_ZG'
      diagTitle = 'dbdz_ZG'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'kappa_GM'
      diagTitle = 'kappa_GM'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'gradb_ZC'
      diagTitle = 'grad_b_ZC'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'gradb_ZG'
      diagTitle = 'grad_b_ZG'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'flux_ZG'
      diagTitle = 'flux_ZG'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'fluxzInt'
      diagTitle = 'flux_zInt'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'taperX  '
      diagTitle = 'taperX'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'taperY  '
      diagTitle = 'taperY'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C --------------------------------------------------------------------
C TESTING
C --------------------------------------------------------------------

      diagName  = 'EKE_YGXC'
      diagTitle = 'EKE_YGXC'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'L_YGXC  '
      diagTitle = 'L_YGXC'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'dqdxYGXC'
      diagTitle = 'dqdxYGXC'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'dqdyYGXC'
      diagTitle = 'dqdyYGXC'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'dpsidx  '
      diagTitle = 'dpsidx'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'modqYGXC'
      diagTitle = 'modGradq_YGXC'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'kqYGXC  '
      diagTitle = 'kappa_q_YGXC'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'ECYGXC  '
      diagTitle = 'energyConv_YGXC'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EGYGXC  '
      diagTitle = 'enstrophyGen_YGXC'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EKE_YCXG'
      diagTitle = 'EKE_YCXG'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'L_YCXG  '
      diagTitle = 'L_YCXG'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'dqdxYCXG'
      diagTitle = 'dqdxYCXG'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'dqdyYCXG'
      diagTitle = 'dqdyYCXG'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'dpsidy  '
      diagTitle = 'dpsidy_YCXG'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'modqYCXG'
      diagTitle = 'modGradq_YCXG'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'kqYCXG  '
      diagTitle = 'kappa_q_YCXG'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'ECYCXG  '
      diagTitle = 'energyConv_YCXG'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'EGYCXG  '
      diagTitle = 'enstrophyGen_YCXG'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'q_YCXC  '
      diagTitle = 'q_YCXC'
      diagCode  = 'SMR     MR      '
      diagUnits = 'm^-2s^-3            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C --------------------------------------------------------------------
C TESTING
C --------------------------------------------------------------------

      diagName  = 'Strain  '
      diagTitle = 'Horizontal Strain of Horizontal Velocities'
      diagCode  = 'SZR     MR      '
      diagUnits = 's^-1            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'Tension '
      diagTitle = 'Horizontal Tension of Horizontal Velocities'
      diagCode  = 'SMR     MR      '
      diagUnits = 's^-1            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* ALLOW_MOM_VECINV */

      diagName  = 'Stretch '
      diagTitle = 'Vortex stretching from QG Leith dynamic viscosity'
      diagCode  = 'SM      MR      '
      diagUnits = 's^-1            '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C--   add diagnostics of tendencies from each terms
      diagUnits = 'm/s^2           '

      diagName  = 'USidDrag'
      diagTitle = 'U momentum tendency from Side Drag'
      diagCode  = 'UUR     MR      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'VSidDrag'
      diagTitle = 'V momentum tendency from Side Drag'
      diagCode  = 'VVR     MR      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

      diagName  = 'Um_Diss '
      diagTitle = 'U momentum tendency from Dissipation'
     &          //' (Explicit part)'
      diagCode  = 'UUR     MR      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'Vm_Diss '
      diagTitle = 'V momentum tendency from Dissipation'
     &          //' (Explicit part)'
      diagCode  = 'VVR     MR      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

      diagName  = 'Um_ImplD'
      diagTitle = 'U momentum tendency from Dissipation'
     &          //' (Implicit part)'
      diagCode  = 'UUR     MR      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'Vm_ImplD'
      diagTitle = 'V momentum tendency from Dissipation'
     &          //' (Implicit part)'
      diagCode  = 'VVR     MR      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

      diagName  = 'Um_Advec'
      diagTitle = 'U momentum tendency from Advection terms'
      diagCode  = 'UUR     MR      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'Vm_Advec'
      diagTitle = 'V momentum tendency from Advection terms'
      diagCode  = 'VVR     MR      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

      diagName  = 'Um_Cori '
      diagTitle = 'U momentum tendency from Coriolis term'
      diagCode  = 'UUR     MR      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'Vm_Cori '
      diagTitle = 'V momentum tendency from Coriolis term'
      diagCode  = 'VVR     MR      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

      diagName  = 'Um_dPhiX'
      diagTitle = 'U momentum tendency from Pressure/Potential grad'
      diagCode  = 'UUR     MR      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'Vm_dPhiY'
      diagTitle = 'V momentum tendency from Pressure/Potential grad'
      diagCode  = 'VVR     MR      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

      diagName  = 'Um_Ext  '
      diagTitle = 'U momentum tendency from external forcing'
      diagCode  = 'UUR     MR      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'Vm_Ext  '
      diagTitle = 'V momentum tendency from external forcing'
      diagCode  = 'VVR     MR      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

#ifdef ALLOW_MOM_VECINV
      diagName  = 'Um_AdvZ3'
      diagTitle = 'U momentum tendency from Vorticity Advection'
      diagCode  = 'UUR     MR      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'Vm_AdvZ3'
      diagTitle = 'V momentum tendency from Vorticity Advection'
      diagCode  = 'VVR     MR      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

      diagName  = 'Um_AdvRe'
      diagTitle = 'U momentum tendency from vertical Advection'
     &          //' (Explicit part)'
      diagCode  = 'UUR     MR      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'Vm_AdvRe'
      diagTitle = 'V momentum tendency from vertical Advection'
     &          //' (Explicit part)'
      diagCode  = 'VVR     MR      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

#endif /* ALLOW_MOM_VECINV */

#ifdef ALLOW_NONHYDROSTATIC
C-    vertical momentum tendencies
      IF ( usingPCoords ) THEN
        diagUnits = 'Pa/s^2          '
      ELSE
        diagUnits = 'm/s^2           '
      ENDIF
      diagName  = 'Wm_Diss '
      diagTitle = 'W momentum tendency from Dissipation'
      diagCode  = 'WMr     LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'Wm_Advec'
      diagTitle = 'W momentum tendency from Advection terms'
      diagCode  = 'WMr     LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

      diagName  = 'WSidDrag'
      diagTitle = 'Vertical momentum tendency from Side Drag'
      diagCode  = 'WMr     LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* ALLOW_NONHYDROSTATIC */

C--   add diagnostics of advective & viscous flux :
C-    bottom fluxes:
      diagName  = 'botTauX '
      diagTitle = 'zonal bottom stress, >0 increases uVel'
      diagUnits = 'N/m^2           '
      diagCode  = 'UU      U1      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'botTauY '
      diagTitle = 'meridional bottom stress, >0 increases vVel'
      diagUnits = 'N/m^2           '
      diagCode  = 'VV      U1      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )

      IF ( usingPCoords ) THEN
        diagUnits = 'Pa.m^3/s^2      '
      ELSE
        diagUnits = 'm^4/s^2         '
      ENDIF

C-     Advective flux of Zonal momentum:
#ifdef ALLOW_MOM_FLUXFORM
      diagName  = 'ADVx_Um '
      diagTitle = 'Zonal      Advective Flux of U momentum'
      diagCode  = 'UM      MR      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'ADVy_Um '
      diagTitle = 'Meridional Advective Flux of U momentum'
      diagCode  = 'VZ      MR      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'ADVrE_Um'
      diagTitle = 'Vertical   Advective Flux of U momentum'
     &          //' (Explicit part)'
      diagCode  = 'WU      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* ALLOW_MOM_FLUXFORM */
c     diagName  = 'ADVrI_Um'
c     diagTitle = 'Vertical   Advective Flux of U momentum'
c    &          //' (Implicit part)'
c     diagCode  = 'WU      LR      '
c     CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c    I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C-    Advective flux of Meridional momentum:
#ifdef ALLOW_MOM_FLUXFORM
      diagName  = 'ADVx_Vm '
      diagTitle = 'Zonal      Advective Flux of V momentum'
      diagCode  = 'UZ      MR      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'ADVy_Vm '
      diagTitle = 'Meridional Advective Flux of V momentum'
      diagCode  = 'VM      MR      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'ADVrE_Vm'
      diagTitle = 'Vertical   Advective Flux of V momentum'
     &          //' (Explicit part)'
      diagCode  = 'WV      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* ALLOW_MOM_FLUXFORM */
c     diagName  = 'ADVrI_Vm'
c     diagTitle = 'Vertical   Advective Flux of V momentum'
c    &          //' (Implicit part)'
c     diagCode  = 'WV      LR      '
c     CALL DIAGNOSTICS_ADDTOLIST( diagNum,
c    I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C-    Viscous flux of Zonal momentum :
#ifdef ALLOW_MOM_FLUXFORM
      diagName  = 'VISCx_Um'
      diagTitle = 'Zonal      Viscous Flux of U momentum'
      diagCode  = 'UM      MR      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'VISCy_Um'
      diagTitle = 'Meridional Viscous Flux of U momentum'
      diagCode  = 'VZ      MR      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'VISrE_Um'
      diagTitle = 'Vertical   Viscous Flux of U momentum'
     &          //' (Explicit part)'
      diagCode  = 'WU      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* ALLOW_MOM_FLUXFORM */
      diagName  = 'VISrI_Um'
      diagTitle = 'Vertical   Viscous Flux of U momentum'
     &          //' (Implicit part)'
      diagCode  = 'WU      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C-    Viscous flux of Meridional momentum :
#ifdef ALLOW_MOM_FLUXFORM
      diagName  = 'VISCx_Vm'
      diagTitle = 'Zonal      Viscous Flux of V momentum'
      diagCode  = 'UZ      MR      '
      diagMate  = diagNum + 2
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'VISCy_Vm'
      diagTitle = 'Meridional Viscous Flux of V momentum'
      diagCode  = 'VM      MR      '
      diagMate  = diagNum
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I   diagName, diagCode, diagUnits, diagTitle, diagMate, myThid )
      diagName  = 'VISrE_Vm'
      diagTitle = 'Vertical   Viscous Flux of V momentum'
     &          //' (Explicit part)'
      diagCode  = 'WV      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )
#endif /* ALLOW_MOM_FLUXFORM */
      diagName  = 'VISrI_Vm'
      diagTitle = 'Vertical   Viscous Flux of V momentum'
     &          //' (Implicit part)'
      diagCode  = 'WV      LR      '
      CALL DIAGNOSTICS_ADDTOLIST( diagNum,
     I          diagName, diagCode, diagUnits, diagTitle, 0, myThid )

C---+----1----+----2----+----3----+----4----+----5----+----6----+----7-|--+----|

#endif /* ALLOW_MOM_COMMON */
#endif /* ALLOW_DIAGNOSTICS */

      RETURN
      END
