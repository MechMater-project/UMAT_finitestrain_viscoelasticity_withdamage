!----------------------------------------------------------------------
!    UMAT FOR FINITE VISCOELASTIC THEORY WITH DAMAGE 
!    Based on Alan Jason Correa implementation of Reese and Govindjee
!    https://github.com/thealanjason/umat_finite_viscoelasticity
!    WITH W(I1,I2)+U(J) HYPERELASTIC COMPRESSIBLE MODEL
!    IMPLEMENTED FOR SOLID PROPELLANT MECHANICAL FEATURES - FLORIAN GOUHIER
!    https://github.com/MechMater-project/UMAT_finitestrain_viscoelasticity_withdamage 
!    This is the version for non hybrid elements
!    If use please refer to the howtoUse document to cite the
!    appropriate references.
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see https://www.gnu.org/licenses/
! ---------------------------------------------------------------------
!     INPUTS
! ---------------------------------------------------------------------
! 
!     NSTATV    - NUMBER OF INTERNAL VARIABLE TO BE SOLVED LOCALLY (Depvar)
!     DFGRD0    - DEFORMATION GRADIENT AT THE BEGINNING OF THE STEP (AT N-1)
!     DFGRD1    - DEFORMATION GRADIENT AT THE CURRENT STWEP (AT N)
!     STRESS    - TOTAL CAUCHY STRESS TENSOR (VOIGT NOTATION)
!     DDSDDE    - TOTAL ABAQUS TANGENT STIFNESS MATRIX
!     TEMP      - TEMPERATURE AT THE START OF THE INCREMENT
!     DTEMP     - INCREMENT OF TEMPERATURE
!     KSTEP     - STEP NUMBER
!     KINC      - INCREMENT NUMBER
!     DTIME     - TIME INCREMENT
!     TIME(1)   - VALUE OF STEP TIME AT THE BEGINNING OF THE CURRENT INCREMENT
!     TIME(2)   - VALUE OF TOTAL TIME AT THE BEGINNING OF THE CURRENT INCREMENT
! 
!     PROPS(1)  - TEMPREF
!     PROPS(2)  - CWLF1
!     PROPS(3)  - CWLF2
!     PROPS(4)  - DAMCHOICE
!     PROPS(5)  - AEXPD
!     PROPS(6)  - BEXPD
!     PROPS(7)  - AEXPH
!     PROPS(8)  - BEXPH
!     PROPS(9)  - W_PH
!     PROPS(10) - P_SAT
!     PROPS(11) - MU
!     PROPS(12) - KELAS
!     PROPS(13) - MUVIS
!     PROPS(14) - KVIS
!     PROPS(15) - ETADEV
!     PROPS(16) - ETAVOL
! 
!----------------------------------------------------------------------
      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)
! 
      INCLUDE 'ABA_PARAM.INC'
! 
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),
     3 DFGRD0(3, 3), DFGRD1(3, 3), TIME(2)
! 
! ---------------------------------------------------------------------
!     LOCAL ARRAYS - NEQ PART
! ---------------------------------------------------------------------
! 
!     BeOLD          - ELASTIC LEFT CAUCHY-GREEN DEFORMATION TENSOR (AT N-1)
!     BeOLDINV       - INVERSE(BeOLD) 
!     CiOLD          - VISCOUS RIGHT CAUCHY-GREEN DEFORMATION TENSOR (AT N-1)
!     CiOLDINV       - INVERSE(CiOLD)
!     BeTR           - ELASTIC TRIAL LEFT CAUCHY-GREEN TENSOR
!     BeTR_          - COPY OF BeTR FOR EIGEN VALUE & VECTOR DECOMPOSITION 
!     PVBeTR         - EIGEN VALUES OF BeTR
!     PDBeTR         - EIGEN VECTORS OF BeTR
!     EPSeTR         - TRIAL NEQ PRINCIPAL(EIGEN) STRAIN
!     EPSe           - ELASTIC LOGARITHMIC PRINCIPAL(EIGEN) STRAIN 
!     DELEPSe        - DELTA EPSe FOR NEWTON RAPHSON'S ITERATION 
!     Be             - ELASTIC LEFT CAUCHY-GREEN TENSOR (AT N)
!     PVBe           - EIGEN VALUES OF Be 
!     PVBeBAR        - EIGEN VALUES OF ELASTIC ISOCHORIC LEFT-CAUCHY GREEN TENSOR
!     RESVEC         - RESIDUAL VECTOR FOR NEWTON RAPHSON'S ITERATION
!     KMAT           - K MATRIX FOR NEWTON RAPHSON'S ITERATION
!     KINV           - INVERSE OF K MATRIX
!     TAUNEQ         - NEQ KIRCHOFF STRESS TENSOR
!     TAUNEQD        - DEVIATORIC PART OF NEQ KIRCHOFF STRESS TENSOR
!     TAUNEQH        - HYDROSTATIC PART OF NEQ KIRCHOFF STRESS TENSOR
!     PVTAUNEQ       - EIGEN VALUES OF TAUNEQ
!     PVTAUNEQD      - EIGEN VALUES OF TAUNEQD
!     PVTAUNEQH      - EIGEN VALUES OF TAUNEQH
!     PDTAUNEQ       - EIGEN VECTORS OF TAUNEQ
!     SIGMANEQ       - NEQ CAUCHY STRESS TENSOR
!     DPVTAUDEPSe    - D(PVTAUNEQ) / D(EPSe)
!     DPVDEVTAUDEPSe - D(PVTAUNEQD) / D(EPSe)
!     DPVHYDTAUDEPSe - D(PVTAUNEQH) / D(EPSe)
!     CALG           - NEQ ALGORITHMIC TANGENT MODULUS
!     L4NEQ          - NEQ INTERMEDIATE MATERIAL TANGENT STIFFNESS MODULUS
!     C4NEQ          - NEQ SPATIAL TANGENT STIFFNESS MODULUS
!     C4NEQJ         - NEQ SPATIAL TANGENT STIFFNESS MODULUS (FE)
!     CNEQ           - NEQ SPATIAL TANGENT STIFFNESS MODULUS (VOIGT)
! 
! ---------------------------------------------------------------------
! 
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: BeOLD
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: BeOLDINV
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: CiOLD
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: CiOLDINV
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: BeTR
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: BeTR_
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3) :: PVBeTR
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: PDBeTR
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3) :: EPSeTR
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3) :: EPSe
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3) :: DELEPSe
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: Be
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3) :: PVBe
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3) :: PVBeBAR     
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3) :: RESVEC
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: KMAT
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: KINV
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: TAUNEQ
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: TAUNEQD
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: TAUNEQH 
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3) :: PVTAUNEQ
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3) :: PVTAUNEQD
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) :: PVTAUNEQH
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: PDTAUNEQ
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3):: SIGMANEQ
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: DPVTAUDEPSe
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: DPVDEVTAUDEPSe
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) :: DPVHYDTAUDEPSe
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3) :: CALG
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3,3,3) :: L4NEQ
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3,3,3) :: C4NEQ
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3,3,3,3) :: C4NEQJ
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,6,6) :: CNEQ
! 
! ---------------------------------------------------------------------
!     LOCAL ARRAYS - EQ PART
! ---------------------------------------------------------------------
! 
!     PVTAUEQ       - EIGEN VALUES OF EQ KIRCHOFF STRESS TENSOR
!     PVTAUEQD      - EIGEN VALUES OF DEVIATORIC PART OF EQ KIRCHOFF STRESS TENSOR
!     PVTAUEQH      - EIGEN VALUES OF HYDROSTATIC PART OF EQ KIRCHOFF STRESS TENSOR
!     SIGMAEQ       - EQ CAUCHY STRESS TENSOR
!     SIGMAEQD      - DEVIATORIC PART OF EQ CAUCHY STRESS TENSOR
!     SIGMAEQH      - HYDROSTATIC PART OF EQ CAUCHY STRESS TENSOR
!     CAB           - CAB FOR EQ SPATIAL TANGENT STIFFNESS MODULUS
!     CABDEV        - DEVIATORIC CAB FOR EQ SPATIAL TANGENT STIFFNESS MODULUS
!     CABHYD        - HYDROSTATIC CAB FOR EQ SPATIAL TANGENT STIFFNESS MODULUS
!     GAB           - GAB FOR EQ SPATIAL TANGENT STIFFNESS MODULUS
!     C4EQMAT       - COEFFICIENTS OF EQ SPATIAL TANGENT STIFFNESS MODULUS 
!     C4EQ          - EQ SPATIAL TANGENT STIFFNESS MODULUS
!     C4EQJ         - EQ SPATIAL TANGENT STIFFNESS MODULUS (FE)
!     CEQ           - EQ SPATIAL TANGENT STIFFNESS MODULUS (VOIGT)
! 
! ---------------------------------------------------------------------
!  
      DOUBLE PRECISION, DIMENSION(3) :: PVTAUEQ
      DOUBLE PRECISION, DIMENSION(3) :: PVTAUEQD
      DOUBLE PRECISION               :: PVTAUEQH
      DOUBLE PRECISION, DIMENSION(3,3) :: SIGMAEQ
      DOUBLE PRECISION, DIMENSION(3,3) :: SIGMAEQD
      DOUBLE PRECISION, DIMENSION(3,3) :: SIGMAEQH
      DOUBLE PRECISION, DIMENSION(3,3) :: CAB
      DOUBLE PRECISION :: CABHYD
      DOUBLE PRECISION, DIMENSION(3,3) :: GAB
      DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C4EQMAT
      DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C4EQ
      DOUBLE PRECISION, DIMENSION(3,3,3,3) :: C4EQJ
      DOUBLE PRECISION, DIMENSION(6,6) :: CEQ
! 
! ---------------------------------------------------------------------
!     LOCAL ARRAYS - TOTAL
! ---------------------------------------------------------------------
! 
!     BTOT          - LEFT CAUCHY-GREEN TENSOR
!     BTOTOLD       - OLD LEFT CAUCHY-GREEN TENSOR
!     BTOT_         - COPY OF BTOT FOR EIGEN VALUE/VECTOR DECOMPOSITION
!     BTOTOOLD_     - COPY OF BTOTOLD FOR EIGEN VALUE/VECTOR DECOMPOSITION
!     BTOTV         - LEFT CAUCHY-GREEN TENSOR (VOIGT)
!     PVBTOT        - EIGEN VALUES OF BTOT
!     PVBTOTOLD     - EIGEN VALUES OF BTOTOLD
!     PDBTOT        - EIGEN VECTORS OF BTOT
!     PDBTOTOLD     - EIGEN VECTORS OF BTOTOLD
!     PVBBAR        - EIGEN VALUES OF ISOCHORIC LEFT CAUCHY-GREEN TENSOR
!     PVBBAROLD     - EIGEN VALUES OF OLD ISOCHORIC LEFT CAUCHY-GREEN TENSOR
!     STRESSTOT     - TOTAL CAUCHY STRESS
!     PRESSURE      - HYDROSTATIC PRESSURE
!     IDT2          - 2ND ORDER IDENTITY TENSOR
! 
! ---------------------------------------------------------------------
! 
      DOUBLE PRECISION, DIMENSION(3,3) :: BTOT
      DOUBLE PRECISION, DIMENSION(3,3) :: BTOTOLD
      DOUBLE PRECISION, DIMENSION(3,3) :: BTOT_
      DOUBLE PRECISION, DIMENSION(3,3) :: BTOTOLD_
      DOUBLE PRECISION, DIMENSION(6) :: BTOTV
      DOUBLE PRECISION, DIMENSION(3) :: PVBTOT
      DOUBLE PRECISION, DIMENSION(3) :: PVBTOTOLD
      DOUBLE PRECISION, DIMENSION(3,3) :: PDBTOT
      DOUBLE PRECISION, DIMENSION(3,3) :: PDBTOTOLD
      DOUBLE PRECISION, DIMENSION(3) :: PVBBAR
      DOUBLE PRECISION, DIMENSION(3) :: PVBBAROLD
      DOUBLE PRECISION, DIMENSION(3,3) :: STRESSTOT
      DOUBLE PRECISION, DIMENSION(3,3) :: IDT2      
!      
! ---------------------------------------------------------------------
!     LOCAL VARIABLES - WLF LAW
! --------------------------------------------------------------------- 
! 
!     TEMP       - TEMPERATURE
!     TEMPREF    - REFERENCE TEMPERATURE FOR WLF LAW
!     CWLF1      - 1RST CONSTANT C1 FROM WLF LAW
!     CWLF2      - 2ND CONSTANT C2 FROM WLF LAW
!     LOGaT      - SHIFT FACTOR FROM WLF LAW
! 
! ---------------------------------------------------------------------  
! 
      DOUBLE PRECISION TEMPREF
      DOUBLE PRECISION CWLF1
      DOUBLE PRECISION CWLF2
      DOUBLE PRECISION LOGaT  
! 
! ---------------------------------------------------------------------
!     LOCAL VARIABLES - DAMAGE LAW
! --------------------------------------------------------------------- 
! 
!     AEXPD       - 1RST PARAMETER A_EXP FROM THE EXPONENTIAL DEVIATORIC DAMAGE LAW
!     BEXPD       - 2ND PARAMETER B_EXP FROM THE EXPONENTIAL DEVIATORIC DAMAGE LAW
!     AEXPH       - 1RST PARAMETER A_EXP FROM THE EXPONENTIAL HYDROSTATIC DAMAGE LAW
!     BEXPH       - 2ND PARAMETER B_EXP FROM THE EXPONENTIAL HYDROSTATIC DAMAGE LAW
!     DAMCHOICE   - CHOICE OF DAMAGE VARIABLE
!     DAMVARINIT  - INITIAL DAMAGE VARIABLE VALUE
! 
!     DAMVARD_OLD - DEVIATORIC DAMAGE VARIABLE VALUE (AT N-1)
!     DAMVARD0    - DEVIATORIC DAMAGE VARIABLE VALUE WITHOUT PRESSURE (AT N-1)
!     DAMVARD     - DEVIATORIC DAMAGE VARIABLE VALUE (AT N)
!     DAMVARD1    - DEVIATORIC DAMAGE VARIABLE VALUE WITHOUT PRESSURE (AT N)
!     DDAMVARDDT  - TIME DERIVATIVE OF DEVIATORIC DAMAGE VARIABLE VALUE (AT N)
!     DAMAGED_OLD - DEVIATORIC DAMAGE FUNCTION VALUE (AT N-1)
!     DAMAGED     - DEVIATORIC DAMAGE FUNCTION VALUE (AT N)
!     DDAMAGEDDT  - TIME DERIVATIVE OF DEVIATORIC DAMAGE FUNCTION VALUE (AT N)
! 
!     DAMVARH_OLD - HYDROSTATIC DAMAGE VARIABLE VALUE (AT N-1)
!     DAMVARH     - HYDROSTATIC DAMAGE VARIABLE VALUE (AT N)
!     DDAMVARHDT  - TIME DERIVATIVE OF HYDROSTATIC DAMAGE VARIABLE VALUE (AT N)
!     DAMAGEH_OLD - HYDROSTATIC DAMAGE FUNCTION VALUE (AT N-1)
!     DAMAGEH     - HYDROSTATIC DAMAGE FUNCTION VALUE (AT N)
!     DDAMAGEHDT  - TIME DERIVATIVE OF HYDROSTATIC DAMAGE FUNCTION VALUE (AT N)
! 
! ---------------------------------------------------------------------       
! 
      INTEGER DAMCHOICE   
      DOUBLE PRECISION DAMVARINIT
      
      DOUBLE PRECISION AEXPD
      DOUBLE PRECISION BEXPD
      DOUBLE PRECISION DAMVARD_OLD
      DOUBLE PRECISION DAMVARD0
      DOUBLE PRECISION DAMVARD2
      DOUBLE PRECISION DAMVARD
      DOUBLE PRECISION DAMVARD1
      DOUBLE PRECISION DAMVARD3
      DOUBLE PRECISION DDAMVARDDT 
      DOUBLE PRECISION DAMAGED_OLD
      DOUBLE PRECISION DAMAGED
      DOUBLE PRECISION DDAMAGEDDT
      
      DOUBLE PRECISION AEXPH
      DOUBLE PRECISION BEXPH
      DOUBLE PRECISION DAMVARH_OLD
      DOUBLE PRECISION DAMVARH
      DOUBLE PRECISION DDAMVARHDT
      DOUBLE PRECISION DAMAGEH_OLD 
      DOUBLE PRECISION DAMAGEH
      DOUBLE PRECISION DDAMAGEHDT
! 
! ---------------------------------------------------------------------
!     LOCAL VARIABLES - HYDROSTATIC PRESSURE CONSIDERATIONS
! ---------------------------------------------------------------------
! 
!     PRESSURE       - HYDROSTATIC PRESSURE = (1/3) tr(SIGMA)
!     W_PH           - 1RST CONSTANT omega_Ph FROM DAMAGE LAW WITH PRESSURE
!     P_SAT          - 2ND CONSTANT P_sat FROM DAMAGE LAW WITH PRESSURE
! 
! ---------------------------------------------------------------------  
!  
      DOUBLE PRECISION PRESSURE
      DOUBLE PRECISION W_PH
      DOUBLE PRECISION P_SAT
! 
! ---------------------------------------------------------------------
!     LOCAL VARIABLES/ARRAYS - ISOTROPIC STRAIN ENERGY INVARIANTS
! ---------------------------------------------------------------------
! 
!     ALPHANEQ       - ALPHA INVARIANT FROM NEQ ISOTROPIC STRAIN ENERGY DENSITY
!     DALPHANEQDI1   - D(ALPHANEQ) / D(I1NEQ)
!     DALPHANEQDI2   - D(ALPHANEQ) / D(I2NEQ)
!     BETANEQ        - BETA INVARIANT FROM NEQ ISOTROPIC STRAIN ENERGY DENSITY
!     DBETANEQDI1    - D(BETANEQ) / D(I1NEQ)
!     DBETANEQDI2    - D(BETANEQ) / D(I2NEQ)
!     DI1eDEPSe      - D(I1e) / D(EPSe)
!     DI2eDEPSe      - D(I2e) / D(EPSe)
!     ALPHAEQ        - ALPHA INVARIANT FROM EQ ISOTROPIC STRAIN ENERGY DENSITY
!     DALPHAEQDI1    - D(ALPHAEQ) / D(I1EQ)
!     DALPHAEQDI2    - D(ALPHAEQ) / D(I2EQ)
!     BETAEQ         - BETA INVARIANT FROM EQ ISOTROPIC STRAIN ENERGY DENSITY
!     DBETAEQDI1     - D(BETAEQ) / D(I1EQ)
!     DBETAEQDI2     - D(BETAEQ) / D(I2EQ)
!     DI1DEPS        - D(I1) / D(EPS)
!     DI2DEPS        - D(I2) / D(EPS)
!  
! ---------------------------------------------------------------------  
!    
      DOUBLE PRECISION ALPHAEQ
      DOUBLE PRECISION DALPHAEQDI1
      DOUBLE PRECISION DALPHAEQDI2
      DOUBLE PRECISION BETAEQ
      DOUBLE PRECISION BETAEQ1
      DOUBLE PRECISION BETAEQ2
      DOUBLE PRECISION, DIMENSION(3) :: DI1DEPS
      DOUBLE PRECISION, DIMENSION(3) :: DI2DEPS
      
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) ::  ALPHANEQ
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) ::  DALPHANEQDI1
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) ::  DALPHANEQDI2
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) ::  BETANEQ
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) ::  DBETANEQDI1
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) ::  DBETANEQDI2
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3) :: DI1eDEPSe
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6,3) :: DI2eDEPSe
! 
! ---------------------------------------------------------------------
!     LOCAL VARIABLES
! ---------------------------------------------------------------------
! 
!     MU         - EQ SHEAR MODULUS FROM NEOHOOKEAN STRAIN ENERGY
!     KELAS      - EQ BULK MODULUS 
!     DET        - JACOBIAN OF TOTAL DEFORMATION GRADIENT
!     DETOLD     - JACOBIAN OF OLD TOTAL DEFORMATION GRADIENT
!     MUVIS      - NEQ SHEAR MODULUS FROM NEOHOOKEAN STRAIN ENERGY
!     KVIS       - NEQ BULK MODULUS  
!     ETADEV     - DEVIATORIC VISCOSITY 
!     ETAVOL     - VOLUMETRIC VISCOSITY   
!     Je         - JACOBIAN OF NEQ DEFORMATION GRADIENT
!     NORMRES    - NORM OF RESIDUAL VECTOR
!     RESTOL     - TOLERANCE FOR NEWTON RHAPSON CONVERGENCE
!     EPS        - TOLERANCE FOR EQUIVALENCE OF FLOAT VALUES   
!     ITER       - LOCAL ITERATION NUMBER
!  
! ---------------------------------------------------------------------
! 
      DOUBLE PRECISION MU
      DOUBLE PRECISION KELAS
      
      DOUBLE PRECISION DET
      DOUBLE PRECISION DETOLD
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) :: MUVIS
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) :: KVIS
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) :: ETADEV
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) :: ETAVOL
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) :: Je
      DOUBLE PRECISION, DIMENSION((NSTATV - 4)/6) :: NORMRES
 
      INTEGER ITER, I, J, K, L, M, N, N_B, N_BRANCHES
      LOGICAL OK_FLAG
      DOUBLE PRECISION :: RESTOL = 1.0D-12
      DOUBLE PRECISION :: EPS = 1.0D-12   
! 
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0,
     1 FIVE=5.D0, SIX=6.D0, EIGHT=8.D0, NINE=9.D0, TEN=10.D0)
!     
! ---------------------------------------------------------------------
!     STEP 0 - DEFINING THE NUMBER OF MAXWELL VISCOELASTIC BRANCHES
! ---------------------------------------------------------------------
! 
      N_BRANCHES = (NSTATV - FOUR) / SIX
! 
! ---------------------------------------------------------------------
!     STEP 1 - COMPUTING NECESSARY VARIABLES AND CONSTANTS
! ---------------------------------------------------------------------
! 
!     DEFINE SECOND ORDER IDENTITY TENSOR
      IDT2 = ZERO
      DO I = 1,3
        IDT2(I,I) = ONE
      END DO
! 
!     CALCULATE JACOBIAN OF DEFORMATION GRADIENT
      DET = DFGRD1(1,1)*DFGRD1(2,2)*DFGRD1(3,3)
     1     -DFGRD1(1,2)*DFGRD1(2,1)*DFGRD1(3,3)
      IF(NSHR.EQ.3) THEN
        DET=DET+DFGRD1(1,2)*DFGRD1(2,3)*DFGRD1(3,1)
     1         +DFGRD1(1,3)*DFGRD1(3,2)*DFGRD1(2,1)
     2         -DFGRD1(1,3)*DFGRD1(3,1)*DFGRD1(2,2)
     3         -DFGRD1(2,3)*DFGRD1(3,2)*DFGRD1(1,1)
      END IF
!  
!     ASSIGN PARAMETERS OF THE WLF LAW
      TEMPREF     = PROPS(1)
      CWLF1       = PROPS(2)
      CWLF2       = PROPS(3)
      LOGaT       = (CWLF1*(TEMPREF - TEMP))/(CWLF2 + (TEMP - TEMPREF))
! 
!     ASSIGN PARAMETERS OF THE DAMAGE EXPONENTIAL LAW
      DAMCHOICE   = PROPS(4)
      AEXPD       = PROPS(5)
      BEXPD       = PROPS(6)
      AEXPH       = PROPS(7)
      BEXPH       = PROPS(8)
! 
!     DEFINING DAMAGE INTERNAL STATE VARIABLES FOR FIRST ITERATION
      IF (KSTEP.EQ.1 .AND. KINC.EQ.1) THEN
            STATEV(NSTATV-3) = ZERO
            STATEV(NSTATV-2) = ZERO
            STATEV(NSTATV-1) = ZERO
            STATEV(NSTATV)   = ZERO
      END IF
! 
!     GET DAMAGE INTERNAL STATE VARIABLES FROM PREVIOUS ITERATION (N-1)
      DAMVARD_OLD = STATEV(NSTATV-3)
      DAMVARH_OLD = STATEV(NSTATV-2)
      DAMAGED_OLD = STATEV(NSTATV-1)
      DAMAGEH_OLD = STATEV(NSTATV)
! 
!     CALCULATE HYDROSTATIC PRESSURE FROM PREVIOUS ITERATION (N-1) 
      PRESSURE = - (ONE/THREE)*(STRESS(1)+STRESS(2)+STRESS(3))
! 
!     ASSIGN PARAMETERS OF THE DAMAGE LAW WITH PRESSURE CONSIDERATIONS
      W_PH        = PROPS(9)
      P_SAT       = PROPS(10)
!  
!     ASSIGN ELASTIC PROPERTIES TO CONSTANTS
      MU          = PROPS(11)
      KELAS       = PROPS(12)
! 
!     VISCOELASTIC LOOP
      DO N_B = 1, N_BRANCHES
!       DEFINING VISCOELASTIC INTERNAL STATE VARIABLES FOR FIRST ITERATION
        IF (KSTEP.EQ.1 .AND. KINC.EQ.1) THEN
              STATEV(1 + 6*(N_B - 1)) = ONE
              STATEV(2 + 6*(N_B - 1)) = ONE
              STATEV(3 + 6*(N_B - 1)) = ONE
              STATEV(4 + 6*(N_B - 1)) = ZERO
              STATEV(5 + 6*(N_B - 1)) = ZERO
              STATEV(6 + 6*(N_B - 1)) = ZERO
        END IF
!     
!       GET VISCOELASTIC INTERNAL STATE VARIABLES FROM PREVIOUS ITERATION (N-1)
        BeOLD(N_B,1,1) = STATEV(1 + 6*(N_B - 1))
        BeOLD(N_B,1,2) = STATEV(4 + 6*(N_B - 1))
        BeOLD(N_B,1,3) = STATEV(5 + 6*(N_B - 1))
        BeOLD(N_B,2,1) = STATEV(4 + 6*(N_B - 1))
        BeOLD(N_B,2,2) = STATEV(2 + 6*(N_B - 1))
        BeOLD(N_B,2,3) = STATEV(6 + 6*(N_B - 1))
        BeOLD(N_B,3,1) = STATEV(5 + 6*(N_B - 1))
        BeOLD(N_B,3,2) = STATEV(6 + 6*(N_B - 1))
        BeOLD(N_B,3,3) = STATEV(3 + 6*(N_B - 1))
! 
!       ASSIGN VISCOELASTIC PROPERTIES TO CONSTANTS CONSIDERING TEMPERATURE
        MUVIS(N_B)        = PROPS(13 + 4*(N_B - 1))
        KVIS(N_B)         = PROPS(14 + 4*(N_B - 1))
        ETADEV(N_B)       = PROPS(15 + 4*(N_B - 1))*(TEN**LOGaT)
        ETAVOL(N_B)       = PROPS(16 + 4*(N_B - 1))*(TEN**LOGaT)
      END DO
! 
! ---------------------------------------------------------------------
!     STEP 2 - CALCULATION OF TRIAL LEFT CAUCHY-GREEN TENSOR
! ---------------------------------------------------------------------
!     VISCOELASTIC LOOP
      DO N_B = 1, N_BRANCHES
        CALL M33INV(BeOLD(N_B, :, :), BeOLDINV(N_B, :, :), OK_FLAG)
        CiOLD(N_B, :, :) = MATMUL(TRANSPOSE(DFGRD0), MATMUL(BeOLDINV(N_B, :, :), DFGRD0))
        CALL M33INV(CiOLD(N_B, :, :), CiOLDINV(N_B, :, :), OK_FLAG)
        BeTR(N_B, :, :)  = MATMUL(DFGRD1, MATMUL(CiOLDINV(N_B, :, :), TRANSPOSE(DFGRD1)))
      END DO
! 
! ---------------------------------------------------------------------
!     STEP 3 - LOCAL NEWTON RHAPSON FOR NEQ EVOLUTION EQUATION
! ---------------------------------------------------------------------
!     IN:
!           BeTR
!           DTIME
!           MUVIS, KVIS, ETADEV, ETAVOL
!
!     OUT:
!           Be, TAUNEQ
!           PVBeTR, PVTAU, PDTAU
!           CALG
! ---------------------------------------------------------------------
!     VISCOELASTIC LOOP
      DO N_B = 1, N_BRANCHES
!       GET EIGEN VALUES AND EIGEN VECTORS OF BeTR
        DO I = 1,3
          DO J = 1,3
            BeTR_(N_B,I,J) = BeTR(N_B,I,J)
          END DO
        END DO
        CALL DSYEVJ3(BeTR_(N_B, :, :), PDBeTR(N_B, :, :), PVBeTR(N_B, :))
        PDBeTR(N_B, :, :) = TRANSPOSE(PDBeTR(N_B, :, :))
! 
!       CALCULATE EPSeTR
        DO I = 1,3
          EPSeTR(N_B,I)=ONE/TWO*LOG(PVBeTR(N_B,I))
        END DO
! 
!       PERFORM LOCAL ITERATION
!       1. Initialize Iteration Variables
        DO I = 1,3
          PVBe(N_B,I) = PVBeTR(N_B,I)
          EPSe(N_B,I) = EPSeTR(N_B,I)
        END DO
!       2. Newton Rhapsons Iteration with MAXITER=200
        DO ITER = 1, 200 
!         - Calculating Jacobian
          Je(N_B) = (PVBe(N_B,1)*PVBe(N_B,2)*PVBe(N_B,3))**(ONE/TWO)
!         - Calculating Principal Values of BeBAR
          PVBeBAR(N_B,1) = Je(N_B)**(-TWO/THREE)*PVBe(N_B,1)
          PVBeBAR(N_B,2) = Je(N_B)**(-TWO/THREE)*PVBe(N_B,2)
          PVBeBAR(N_B,3) = Je(N_B)**(-TWO/THREE)*PVBe(N_B,3)
! _________________________________________________________
!         - Modification required when changing Wneq(Ce) 
!         - Define invariants for the deviatoric NEQ deviatoric strain energy
! _________________________________________________________
          ALPHANEQ(N_B)     = MUVIS(N_B)/TWO
          DALPHANEQDI1(N_B) = ZERO
          DALPHANEQDI2(N_B) = ZERO
          BETANEQ(N_B)      = ZERO
          DBETANEQDI1(N_B)  = ZERO
          DBETANEQDI2(N_B)  = ZERO
!         - Calculating Principal Values of Deviatoric Kirchoff Stress
          PVTAUNEQD(N_B,1) = ALPHANEQ(N_B) *
     1                     ((FOUR/THREE)*PVBeBAR(N_B,1)
     2                      -(TWO/THREE)*PVBeBAR(N_B,2)
     3                      -(TWO/THREE)*PVBeBAR(N_B,3))
     4                     + BETANEQ(N_B) *
     5                     ((FOUR/THREE)*(PVBeBAR(N_B,1)**TWO)
     6                      -(TWO/THREE)*(PVBeBAR(N_B,2)**TWO)
     7                      -(TWO/THREE)*(PVBeBAR(N_B,3)**TWO))
          PVTAUNEQD(N_B,2) = ALPHANEQ(N_B) *
     1                     ((FOUR/THREE)*PVBeBAR(N_B,2)
     2                      -(TWO/THREE)*PVBeBAR(N_B,3)
     3                      -(TWO/THREE)*PVBeBAR(N_B,1))
     4                     + BETANEQ(N_B) *
     5                     ((FOUR/THREE)*(PVBeBAR(N_B,2)**TWO)
     6                      -(TWO/THREE)*(PVBeBAR(N_B,3)**TWO)
     7                      -(TWO/THREE)*(PVBeBAR(N_B,1)**TWO))
          PVTAUNEQD(N_B,3) = ALPHANEQ(N_B) *
     1                     ((FOUR/THREE)*PVBeBAR(N_B,3)
     2                      -(TWO/THREE)*PVBeBAR(N_B,1)
     3                      -(TWO/THREE)*PVBeBAR(N_B,2))
     4                     + BETANEQ(N_B) *
     5                     ((FOUR/THREE)*(PVBeBAR(N_B,3)**TWO)
     6                      -(TWO/THREE)*(PVBeBAR(N_B,1)**TWO)
     7                      -(TWO/THREE)*(PVBeBAR(N_B,2)**TWO))
! _________________________________________________________
!         - Modification required when changing U(J) 
!         - Calculating Principal Values of Hydrostatic Kirchoff Stress
! _________________________________________________________
          PVTAUNEQH(N_B) = KVIS(N_B)/TWO*(Je(N_B)*Je(N_B)-ONE)
!         - Calculating Residual Vector
          DO N = 1, 3
            RESVEC(N_B,N) = EPSe(N_B,N) - EPSeTR(N_B,N) 
     1                    + DTIME * 
     2		            (ONE/(TWO*ETADEV(N_B))*PVTAUNEQD(N_B,N)
     2                      + ONE/(THREE*ETAVOL(N_B))*PVTAUNEQH(N_B))
          END DO
!         - Calculating Norm of Residual Vector
          NORMRES(N_B) = (RESVEC(N_B,1)**TWO 
     1                  + RESVEC(N_B,2)**TWO
     2                  + RESVEC(N_B,3)**TWO)**(ONE/TWO)
!         - Terminate if Norm of Residual Vector is less than Tolerance
          IF ((ITER.GT.1).AND.ABS(NORMRES(N_B)).LT.RESTOL) THEN
              EXIT
          END IF
!         - Calculating d(I1e)/dEPSe and d(I2e)/dEPSe
          DI1eDEPSe(N_B,1) = (FOUR/THREE)*PVBeBAR(N_B,1)
     1                     - (TWO/THREE)*PVBeBAR(N_B,2)
     2                     - (TWO/THREE)*PVBeBAR(N_B,3)
          DI2eDEPSe(N_B,1) = (TWO/THREE)*PVBeBAR(N_B,1) *
     1                     (PVBeBAR(N_B,2) + PVBeBAR(N_B,3))
     2                     - (FOUR/THREE)*PVBeBAR(N_B,2)*PVBeBAR(N_B,3)
          DI1eDEPSe(N_B,2) = (FOUR/THREE)*PVBeBAR(N_B,2)
     1                     - (TWO/THREE)*PVBeBAR(N_B,3)
     2                     - (TWO/THREE)*PVBeBAR(N_B,1)
          DI2eDEPSe(N_B,2) = (TWO/THREE)*PVBeBAR(N_B,2) *
     1                     (PVBeBAR(N_B,3) + PVBeBAR(N_B,1))
     2                     - (FOUR/THREE)*PVBeBAR(N_B,3)*PVBeBAR(N_B,1)
          DI1eDEPSe(N_B,3) = (FOUR/THREE)*PVBeBAR(N_B,3)
     1                     - (TWO/THREE)*PVBeBAR(N_B,1)
     2                     - (TWO/THREE)*PVBeBAR(N_B,2)
          DI2eDEPSe(N_B,3) =(TWO/THREE)*PVBeBAR(N_B,3) *
     1                     (PVBeBAR(N_B,1) + PVBeBAR(N_B,2))
     2                     - (FOUR/THREE)*PVBeBAR(N_B,1)*PVBeBAR(N_B,2)
!         - Calculating d(PVTAUNEQD)/dEPSe
          DPVDEVTAUDEPSe(N_B,1,1) = (DALPHANEQDI1(N_B)*DI1eDEPSe(N_B,1)
     1                    + DALPHANEQDI2(N_B)*DI2eDEPSe(N_B,1)) *
     2                    ((FOUR/THREE)*PVBeBAR(N_B,1)
     2                     - (TWO/THREE)*PVBeBAR(N_B,2)
     2                     - (TWO/THREE)*PVBeBAR(N_B,3))
     3                    + ALPHANEQ(N_B) *
     3                      (((FOUR*FOUR)/NINE)*PVBeBAR(N_B,1)
     3                       + (FOUR/NINE)*PVBeBAR(N_B,2)
     3                       + (FOUR/NINE)*PVBeBAR(N_B,3))
     4                    + (DBETANEQDI1(N_B)*DI1eDEPSe(N_B,1) 
     4                       + DBETANEQDI2(N_B)*DI2eDEPSe(N_B,1)) *
     5                      ((FOUR/THREE)*(PVBeBAR(N_B,1)**TWO)
     5                       - (TWO/THREE)*(PVBeBAR(N_B,2)**TWO)
     5                       - (TWO/THREE)*(PVBeBAR(N_B,3)**TWO))
     6                    + BETANEQ(N_B) *
     6                      (((FOUR*EIGHT)/NINE)*(PVBeBAR(N_B,1)**TWO)
     6                       + (EIGHT/NINE)*(PVBeBAR(N_B,2)**TWO)
     6                       + (EIGHT/NINE)*(PVBeBAR(N_B,3)**TWO) )
          DPVDEVTAUDEPSe(N_B,2,2) = (DALPHANEQDI1(N_B)*DI1eDEPSe(N_B,2) 
     1                    + DALPHANEQDI2(N_B)*DI2eDEPSe(N_B,2)) *
     2                    ((FOUR/THREE)*PVBeBAR(N_B,2)
     2                     - (TWO/THREE)*PVBeBAR(N_B,3)
     2                     - (TWO/THREE)*PVBeBAR(N_B,1))
     3                    + ALPHANEQ(N_B) *
     3                      (((FOUR*FOUR)/NINE)*PVBeBAR(N_B,2)
     3                       + (FOUR/NINE)*PVBeBAR(N_B,3)
     3                       + (FOUR/NINE)*PVBeBAR(N_B,1))
     4                    + (DBETANEQDI1(N_B)*DI1eDEPSe(N_B,2) 
     4                       + DBETANEQDI2(N_B)*DI2eDEPSe(N_B,2)) *
     5                      ((FOUR/THREE)*(PVBeBAR(N_B,2)**TWO)
     5                       - (TWO/THREE)*(PVBeBAR(N_B,3)**TWO)
     5                       - (TWO/THREE)*(PVBeBAR(N_B,1)**TWO))
     6                    + BETANEQ(N_B) *
     6                      (((FOUR*EIGHT)/NINE)*(PVBeBAR(N_B,2)**TWO)
     6                       + (EIGHT/NINE)*(PVBeBAR(N_B,3)**TWO)
     6                       + (EIGHT/NINE)*(PVBeBAR(N_B,1)**TWO))
          DPVDEVTAUDEPSe(N_B,3,3) = (DALPHANEQDI1(N_B)*DI1eDEPSe(N_B,3) 
     1                    + DALPHANEQDI2(N_B)*DI2eDEPSe(N_B,3)) *
     2                    ((FOUR/THREE)*PVBeBAR(N_B,3)
     2                     - (TWO/THREE)*PVBeBAR(N_B,1)
     2                     - (TWO/THREE)*PVBeBAR(N_B,2))
     3                    + ALPHANEQ(N_B) *
     3                      (((FOUR*FOUR)/NINE)*PVBeBAR(N_B,3)
     3                       + (FOUR/NINE)*PVBeBAR(N_B,1)
     3                       + (FOUR/NINE)*PVBeBAR(N_B,2))
     4                    + (DBETANEQDI1(N_B)*DI1eDEPSe(N_B,3) 
     4                       + DBETANEQDI2(N_B)*DI2eDEPSe(N_B,3)) *
     5                      ((FOUR/THREE)*(PVBeBAR(N_B,3)**TWO)
     5                       - (TWO/THREE)*(PVBeBAR(N_B,1)**TWO)
     5                       - (TWO/THREE)*(PVBeBAR(N_B,2)**TWO) )
     6                    + BETANEQ(N_B) *
     6                      (((FOUR*EIGHT)/NINE)*(PVBeBAR(N_B,3)**TWO)
     6                       + (EIGHT/NINE)*(PVBeBAR(N_B,1)**TWO)
     6                       + (EIGHT/NINE)*(PVBeBAR(N_B,2)**TWO))
          DPVDEVTAUDEPSe(N_B,1,2) = (DALPHANEQDI1(N_B)*DI1eDEPSe(N_B,2) 
     1                    + DALPHANEQDI2(N_B)*DI2eDEPSe(N_B,2)) *
     2                    ((FOUR/THREE)*PVBeBAR(N_B,1)
     2                     - (TWO/THREE)*PVBeBAR(N_B,2)
     2                     - (TWO/THREE)*PVBeBAR(N_B,3) )
     3                    + ALPHANEQ(N_B) *
     3                      ((FOUR/NINE)*PVBeBAR(N_B,3)
     3                       - (EIGHT/NINE)*PVBeBAR(N_B,1)
     3                       - (EIGHT/NINE)*PVBeBAR(N_B,2))
     4                    + (DBETANEQDI1(N_B)*DI1eDEPSe(N_B,2) 
     4                       + DBETANEQDI2(N_B)*DI2eDEPSe(N_B,2)) *
     5                      ((FOUR/THREE)*(PVBeBAR(N_B,1)**TWO)
     5                       - (TWO/THREE)*(PVBeBAR(N_B,2)**TWO)
     5                       - (TWO/THREE)*(PVBeBAR(N_B,3)**TWO))
     6                    + BETANEQ(N_B) *
     6                      ((EIGHT/NINE)*(PVBeBAR(N_B,3)**TWO)
     6                       - ((TWO*EIGHT)/NINE)*(PVBeBAR(N_B,1)**TWO)
     6                       - ((TWO*EIGHT)/NINE)*(PVBeBAR(N_B,2)**TWO))
          DPVDEVTAUDEPSe(N_B,1,3) = (DALPHANEQDI1(N_B)*DI1eDEPSe(N_B,3) 
     1                    + DALPHANEQDI2(N_B)*DI2eDEPSe(N_B,3)) *
     2                    ((FOUR/THREE)*PVBeBAR(N_B,1)
     2                     - (TWO/THREE)*PVBeBAR(N_B,2)
     2                     - (TWO/THREE)*PVBeBAR(N_B,3) )
     3                    + ALPHANEQ(N_B) *
     3                      ((FOUR/NINE)*PVBeBAR(N_B,2)
     3                       - (EIGHT/NINE)*PVBeBAR(N_B,3)
     3                       - (EIGHT/NINE)*PVBeBAR(N_B,1))
     4                    + (DBETANEQDI1(N_B)*DI1eDEPSe(N_B,3) 
     4                       + DBETANEQDI2(N_B)*DI2eDEPSe(N_B,3)) *
     5                      ((FOUR/THREE)*(PVBeBAR(N_B,1)**TWO)
     5                       - (TWO/THREE)*(PVBeBAR(N_B,2)**TWO)
     5                       - (TWO/THREE)*(PVBeBAR(N_B,3)**TWO))
     6                    + BETANEQ(N_B) *
     6                      ((EIGHT/NINE)*(PVBeBAR(N_B,2)**TWO)
     6                       - ((TWO*EIGHT)/NINE)*(PVBeBAR(N_B,3)**TWO)
     6                       - ((TWO*EIGHT)/NINE)*(PVBeBAR(N_B,1)**TWO))
          DPVDEVTAUDEPSe(N_B,2,3) = (DALPHANEQDI1(N_B)*DI1eDEPSe(N_B,3) 
     1                    + DALPHANEQDI2(N_B)*DI2eDEPSe(N_B,3)) *
     2                    ((FOUR/THREE)*PVBeBAR(N_B,2)
     2                     - (TWO/THREE)*PVBeBAR(N_B,3)
     2                     - (TWO/THREE)*PVBeBAR(N_B,1) )
     3                    + ALPHANEQ(N_B) *
     3                      ((FOUR/NINE)*PVBeBAR(N_B,1)
     3                       - (EIGHT/NINE)*PVBeBAR(N_B,2)
     3                       - (EIGHT/NINE)*PVBeBAR(N_B,3))
     4                    + (DBETANEQDI1(N_B)*DI1eDEPSe(N_B,3) 
     4                       + DBETANEQDI2(N_B)*DI2eDEPSe(N_B,3)) *
     5                      ((FOUR/THREE)*(PVBeBAR(N_B,2)**TWO)
     5                       - (TWO/THREE)*(PVBeBAR(N_B,3)**TWO)
     5                       - (TWO/THREE)*(PVBeBAR(N_B,1)**TWO) )
     6                    + BETANEQ(N_B) *
     6                      ((EIGHT/NINE)*(PVBeBAR(N_B,1)**TWO)
     6                       - ((TWO*EIGHT)/NINE)*(PVBeBAR(N_B,2)**TWO)
     6                       - ((TWO*EIGHT)/NINE)*(PVBeBAR(N_B,3)**TWO))
          DPVDEVTAUDEPSe(N_B,2,1) = DPVDEVTAUDEPSe(N_B,1,2)
          DPVDEVTAUDEPSe(N_B,3,1) = DPVDEVTAUDEPSe(N_B,1,3)
          DPVDEVTAUDEPSe(N_B,3,2) = DPVDEVTAUDEPSe(N_B,2,3)
! _________________________________________________________
!         - Modification required when changing U(J) 
!         - Calculating d(PVTAUNEQD)/dEPSe
! _________________________________________________________
          DPVHYDTAUDEPSe(N_B) = KVIS(N_B)*Je(N_B)*Je(N_B)
!
!         - Calculating K Matrix for Newton Rhapson
!
          DO I = 1,3
            DO J = 1,3
              KMAT(N_B,I,J) = IDT2(I,J)
     1                   + DTIME * 
     2		           ((ONE/(TWO*ETADEV(N_B)))*DPVDEVTAUDEPSe(N_B,I,J) 
     2                    + (ONE/(THREE*ETAVOL(N_B))*DPVHYDTAUDEPSe(N_B))) 
            END DO
          END DO
!         - Calculating KINV
          CALL M33INV(KMAT(N_B, :, :), KINV(N_B, :, :), OK_FLAG)
!         - Calculating NEQ Strain Increment
          DELEPSe(N_B,1) = -(KINV(N_B,1,1)*RESVEC(N_B,1)
     1                      +KINV(N_B,1,2)*RESVEC(N_B,2)
     2                      +KINV(N_B,1,3)*RESVEC(N_B,3))
                      
          DELEPSe(N_B,2) = -(KINV(N_B,2,1)*RESVEC(N_B,1)
     1                      +KINV(N_B,2,2)*RESVEC(N_B,2)
     2                      +KINV(N_B,2,3)*RESVEC(N_B,3))
         
          DELEPSe(N_B,3) = -(KINV(N_B,3,1)*RESVEC(N_B,1)
     1                      +KINV(N_B,3,2)*RESVEC(N_B,2)
     2                      +KINV(N_B,3,3)*RESVEC(N_B,3))
!         - Calculating Updated NEQ Strain
          DO I = 1,3
          EPSe(N_B,I) = EPSe(N_B,I) + DELEPSe(N_B,I)
          END DO
!         - Updating Eigen Value of Be
          DO I = 1,3
            PVBe(N_B,I) = EXP(TWO*EPSe(N_B,I))
          END DO
        END DO
!       3. Non-Convergence Test for Newton Rhapson
!        IF (ITER.GT.200) THEN
!              WRITE(*,*) "LOCAL ITERATION DID NOT CONVERGE" 
!        END IF
! 
!       CALCULATE REQUIRED QUANTITIES
!       1. Be
        DO I=1,3
          DO J=1,3
              Be(N_B,J,I) = PVBe(N_B,1)*PDBeTR(N_B,1,J)*PDBeTR(N_B,1,I)
     1                    + PVBe(N_B,2)*PDBeTR(N_B,2,J)*PDBeTR(N_B,2,I)
     2                    + PVBe(N_B,3)*PDBeTR(N_B,3,J)*PDBeTR(N_B,3,I)
          END DO
        END DO
! 
!       2. PVTAUNEQ
        DO I=1,3
              PVTAUNEQ(N_B,I) = PVTAUNEQD(N_B,I) + PVTAUNEQH(N_B)
        END DO
! 
!       3. PDTAUNEQ
        DO I=1,3
          DO J=1,3
            PDTAUNEQ(N_B,J,I) = PDBeTR(N_B,J,I)
          END DO
        END DO
! 
!       4. CALCULATING DEVIATORIC, HYDROSTATIC AND TOTAL NEQ KIRCHHOFF STRESS TENSOR
        DO I=1,3
          DO J=1,3
            TAUNEQD(N_B,J,I) = PVTAUNEQD(N_B,1)*PDTAUNEQ(N_B,1,J)*PDTAUNEQ(N_B,1,I)
     1                       + PVTAUNEQD(N_B,2)*PDTAUNEQ(N_B,2,J)*PDTAUNEQ(N_B,2,I)
     2                       + PVTAUNEQD(N_B,3)*PDTAUNEQ(N_B,3,J)*PDTAUNEQ(N_B,3,I)
            TAUNEQH(N_B,J,I) = PVTAUNEQH(N_B)*PDTAUNEQ(N_B,1,J)*PDTAUNEQ(N_B,1,I)
     1                       + PVTAUNEQH(N_B)*PDTAUNEQ(N_B,2,J)*PDTAUNEQ(N_B,2,I)
     2                       + PVTAUNEQH(N_B)*PDTAUNEQ(N_B,3,J)*PDTAUNEQ(N_B,3,I)
            TAUNEQ(N_B,J,I) = TAUNEQD(N_B,J,I) + TAUNEQH(N_B,J,I)
          END DO
        END DO
! 
!       5. CALG
!       - Calculating Calculating d(PVTAUNEQ)/dEPSe
        DO I=1,3
          DO J=1,3
            DPVTAUDEPSe(N_B,J,I) = DPVDEVTAUDEPSe(N_B,J,I) + DPVHYDTAUDEPSe(N_B)
          END DO
        END DO
!       - Calculating CALG
        CALG(N_B, :, :) = MATMUL(DPVTAUDEPSe(N_B, :, :), KINV(N_B, :, :))
      END DO
!     
! ---------------------------------------------------------------------
!     STEP 4 - UPDATE INTERNAL STATE VARIABLE - Be
! ---------------------------------------------------------------------
!     SET INTERNAL STATE-VARIABLE FROM FINAL LOCAL ITERATION
      DO N_B = 1, N_BRANCHES
        STATEV(1 + 6*(N_B - 1)) = Be(N_B,1,1)
        STATEV(2 + 6*(N_B - 1)) = Be(N_B,2,2)
        STATEV(3 + 6*(N_B - 1)) = Be(N_B,3,3)
        STATEV(4 + 6*(N_B - 1)) = Be(N_B,1,2)
        STATEV(5 + 6*(N_B - 1)) = Be(N_B,1,3)
        STATEV(6 + 6*(N_B - 1)) = Be(N_B,2,3)
      END DO
!
! ---------------------------------------------------------------------
!     STEP 5 - CALCULATION OF NEQ CAUCHY STRESS & ELASTICITY TENSOR
! ---------------------------------------------------------------------
!     INITIALISATION FOR THE NEQ CAUCHY STRESS TENSOR
      SIGMANEQ = ZERO
!     INITIALISATION FOR THE MATERIAL TENSOR (L4NEQ) IN INTERMEDIATE CONFIGURATION 
!     EXPRESSED IN EIGEN BASIS 
      L4NEQ = ZERO
!     INITIALISATION FOR THE MATERIAL TENSOR (C4NEQ) IN CURRENT CONFIGURATION
!     EXPRESSED IN STANDARD BASIS
      C4NEQ = ZERO
!     INITIALISATION FOR THE JAUMANN RATE FORM OF MATERIAL TENSOR
      C4NEQJ = ZERO
!     INITIALISATION FOR THE GENERATE ABAQUS NEQ TANGENT STIFFNESS MATRIX (VOIGT NOTATION)
      CNEQ = ZERO
!     
!     VISCOELASTIC LOOP
      DO N_B = 1, N_BRANCHES
!       CALCULATE NEQ CAUCHY STRESS TENSOR
        DO J=1,3
          DO K=1,3
            SIGMANEQ(N_B,J,K) = TAUNEQ(N_B,J,K)/DET
          END DO
        END DO
! 
!       CALCULATE MATERIAL TENSOR (L4NEQ) IN INTERMEDIATE CONFIGURATION 
!       EXPRESSED IN EIGEN BASIS 
!       1. CALCULATE COEFFICIENTS FOR AABB
        DO J=1,3
          DO K=1,3
            L4NEQ(N_B,J,J,K,K) = (CALG(N_B,J,K) - PVTAUNEQ(N_B,J)*TWO*IDT2(J,K))
     1                           /(PVBeTR(N_B,J)*PVBeTR(N_B,K))
          END DO
        END DO 
!       2. CALCULATE COEFFICIENTS FOR ABAB / ABBA
        DO J=1,3
          DO K=1,3
            IF (J.EQ.K) THEN
              CYCLE
            END IF
!           - In case of Equal PVBeTR, Use Alternate Form
            IF (ABS(PVBeTR(N_B,J)-PVBeTR(N_B,K)).LT.EPS) THEN
              L4NEQ(N_B,J,K,J,K) = (L4NEQ(N_B,J,J,J,J)-L4NEQ(N_B,J,J,K,K))/TWO
            ELSE
              L4NEQ(N_B,J,K,J,K) = ((PVTAUNEQ(N_B,K)/PVBeTR(N_B,K))-(PVTAUNEQ(N_B,J)/PVBeTR(N_B,J)))
     1                             /(PVBeTR(N_B,K)-PVBeTR(N_B,J))
            END IF
            L4NEQ(N_B,J,K,K,J) = L4NEQ(N_B,J,K,J,K)
          END DO
        END DO
!       CALCULATE MATERIAL TENSOR (C4NEQ) IN CURRENT CONFIGURATION
!       EXPRESSED IN STANDARD BASIS
        DO J=1,3
          DO K=1,3
            DO L=1,3
              DO M=1,3
                C4NEQ(N_B,J,K,L,M) = 
     1          L4NEQ(N_B,1,1,1,1)*PVBeTR(N_B,1)*PVBeTR(N_B,1)
     1                            *PDBeTR(N_B,1,J)*PDBeTR(N_B,1,K)
     1                            *PDBeTR(N_B,1,L)*PDBeTR(N_B,1,M)
     1         +L4NEQ(N_B,1,1,2,2)*PVBeTR(N_B,1)*PVBeTR(N_B,2)
     1                            *PDBeTR(N_B,1,J)*PDBeTR(N_B,1,K)
     1                            *PDBeTR(N_B,2,L)*PDBeTR(N_B,2,M)
     1         +L4NEQ(N_B,1,1,3,3)*PVBeTR(N_B,1)*PVBeTR(N_B,3)
     1                            *PDBeTR(N_B,1,J)*PDBeTR(N_B,1,K)
     1                            *PDBeTR(N_B,3,L)*PDBeTR(N_B,3,M)
     1         +L4NEQ(N_B,2,2,1,1)*PVBeTR(N_B,2)*PVBeTR(N_B,1)
     1                            *PDBeTR(N_B,2,J)*PDBeTR(N_B,2,K)
     1                            *PDBeTR(N_B,1,L)*PDBeTR(N_B,1,M)
     1         +L4NEQ(N_B,2,2,2,2)*PVBeTR(N_B,2)*PVBeTR(N_B,2)
     1                            *PDBeTR(N_B,2,J)*PDBeTR(N_B,2,K)
     1                            *PDBeTR(N_B,2,L)*PDBeTR(N_B,2,M)
     1         +L4NEQ(N_B,2,2,3,3)*PVBeTR(N_B,2)*PVBeTR(N_B,3)
     1                            *PDBeTR(N_B,2,J)*PDBeTR(N_B,2,K)
     1                            *PDBeTR(N_B,3,L)*PDBeTR(N_B,3,M)
     1         +L4NEQ(N_B,3,3,1,1)*PVBeTR(N_B,3)*PVBeTR(N_B,1)
     1                            *PDBeTR(N_B,3,J)*PDBeTR(N_B,3,K)
     1                            *PDBeTR(N_B,1,L)*PDBeTR(N_B,1,M)
     1         +L4NEQ(N_B,3,3,2,2)*PVBeTR(N_B,3)*PVBeTR(N_B,2)
     1                            *PDBeTR(N_B,3,J)*PDBeTR(N_B,3,K)
     1                            *PDBeTR(N_B,2,L)*PDBeTR(N_B,2,M)
     1         +L4NEQ(N_B,3,3,3,3)*PVBeTR(N_B,3)*PVBeTR(N_B,3)
     1                            *PDBeTR(N_B,3,J)*PDBeTR(N_B,3,K)
     1                            *PDBeTR(N_B,3,L)*PDBeTR(N_B,3,M)
     1         +L4NEQ(N_B,1,2,1,2)*PVBeTR(N_B,1)*PVBeTR(N_B,2)
     1                            *PDBeTR(N_B,1,J)*PDBeTR(N_B,2,K)
     1                            *PDBeTR(N_B,1,L)*PDBeTR(N_B,2,M)
     1         +L4NEQ(N_B,1,2,2,1)*PVBeTR(N_B,1)*PVBeTR(N_B,2)
     1                            *PDBeTR(N_B,1,J)*PDBeTR(N_B,2,K)
     1                            *PDBeTR(N_B,2,L)*PDBeTR(N_B,1,M)
     1         +L4NEQ(N_B,1,3,1,3)*PVBeTR(N_B,1)*PVBeTR(N_B,3)
     1                            *PDBeTR(N_B,1,J)*PDBeTR(N_B,3,K)
     1                            *PDBeTR(N_B,1,L)*PDBeTR(N_B,3,M)
     1         +L4NEQ(N_B,1,3,3,1)*PVBeTR(N_B,1)*PVBeTR(N_B,3)
     1                            *PDBeTR(N_B,1,J)*PDBeTR(N_B,3,K)
     1                            *PDBeTR(N_B,3,L)*PDBeTR(N_B,1,M)
     1         +L4NEQ(N_B,2,1,2,1)*PVBeTR(N_B,2)*PVBeTR(N_B,1)
     1                            *PDBeTR(N_B,2,J)*PDBeTR(N_B,1,K)
     1                            *PDBeTR(N_B,2,L)*PDBeTR(N_B,1,M)
     1         +L4NEQ(N_B,2,1,1,2)*PVBeTR(N_B,2)*PVBeTR(N_B,1)
     1                            *PDBeTR(N_B,2,J)*PDBeTR(N_B,1,K)
     1                            *PDBeTR(N_B,1,L)*PDBeTR(N_B,2,M)
     1         +L4NEQ(N_B,2,3,2,3)*PVBeTR(N_B,2)*PVBeTR(N_B,3)
     1                            *PDBeTR(N_B,2,J)*PDBeTR(N_B,3,K)
     1                            *PDBeTR(N_B,2,L)*PDBeTR(N_B,3,M)
     1         +L4NEQ(N_B,2,3,3,2)*PVBeTR(N_B,2)*PVBeTR(N_B,3)
     1                            *PDBeTR(N_B,2,J)*PDBeTR(N_B,3,K)
     1                            *PDBeTR(N_B,3,L)*PDBeTR(N_B,2,M)
     1         +L4NEQ(N_B,3,1,3,1)*PVBeTR(N_B,3)*PVBeTR(N_B,1)
     1                            *PDBeTR(N_B,3,J)*PDBeTR(N_B,1,K)
     1                            *PDBeTR(N_B,3,L)*PDBeTR(N_B,1,M)
     1         +L4NEQ(N_B,3,1,1,3)*PVBeTR(N_B,3)*PVBeTR(N_B,1)
     1                            *PDBeTR(N_B,3,J)*PDBeTR(N_B,1,K)
     1                            *PDBeTR(N_B,1,L)*PDBeTR(N_B,3,M)
     1         +L4NEQ(N_B,3,2,3,2)*PVBeTR(N_B,3)*PVBeTR(N_B,2)
     1                            *PDBeTR(N_B,3,J)*PDBeTR(N_B,2,K)
     1                            *PDBeTR(N_B,3,L)*PDBeTR(N_B,2,M)
     1         +L4NEQ(N_B,3,2,2,3)*PVBeTR(N_B,3)*PVBeTR(N_B,2)
     1                            *PDBeTR(N_B,3,J)*PDBeTR(N_B,2,K)
     1                            *PDBeTR(N_B,2,L)*PDBeTR(N_B,3,M)
              END DO
            END DO
          END DO
        END DO
! 
!       CALCULATE JAUMANN RATE FORM OF MATERIAL TENSOR
        DO I=1,3
          DO J=1,3
            DO K=1,3
              DO L=1,3
                C4NEQJ(N_B,I,J,K,L) = C4NEQ(N_B,I,J,K,L)/DET
     1                              + (1 - DAMAGED_OLD) * (IDT2(I,K)*TAUNEQ(N_B,J,L)
     1                               +IDT2(I,L)*TAUNEQ(N_B,J,K)
     1                               +IDT2(J,K)*TAUNEQ(N_B,I,L)
     1                               +IDT2(J,L)*TAUNEQ(N_B,I,K))/(TWO*DET)        
              END DO
            END DO
          END DO
        END DO
! 
!       GENERATE ABAQUS NEQ TANGENT STIFFNESS MATRIX (VOIGT NOTATION)
        CNEQ(N_B,1,1) = C4NEQJ(N_B,1,1,1,1)
        CNEQ(N_B,1,2) = C4NEQJ(N_B,1,1,2,2)
        CNEQ(N_B,1,3) = C4NEQJ(N_B,1,1,3,3)
        CNEQ(N_B,1,4) = C4NEQJ(N_B,1,1,1,2)
        CNEQ(N_B,1,5) = C4NEQJ(N_B,1,1,1,3)
        CNEQ(N_B,1,6) = C4NEQJ(N_B,1,1,2,3)
        CNEQ(N_B,2,1) = CNEQ(N_B,1,2)
        CNEQ(N_B,2,2) = C4NEQJ(N_B,2,2,2,2)
        CNEQ(N_B,2,3) = C4NEQJ(N_B,2,2,3,3)
        CNEQ(N_B,2,4) = C4NEQJ(N_B,2,2,1,2)
        CNEQ(N_B,2,5) = C4NEQJ(N_B,2,2,1,3)
        CNEQ(N_B,2,6) = C4NEQJ(N_B,2,2,2,3)
        CNEQ(N_B,3,1) = CNEQ(N_B,1,3)
        CNEQ(N_B,3,2) = CNEQ(N_B,2,3)
        CNEQ(N_B,3,3) = C4NEQJ(N_B,3,3,3,3)
        CNEQ(N_B,3,4) = C4NEQJ(N_B,3,3,1,2)
        CNEQ(N_B,3,5) = C4NEQJ(N_B,3,3,1,3)
        CNEQ(N_B,3,6) = C4NEQJ(N_B,3,3,2,3)
        CNEQ(N_B,4,1) = CNEQ(N_B,1,4)
        CNEQ(N_B,4,2) = CNEQ(N_B,2,4)
        CNEQ(N_B,4,3) = CNEQ(N_B,3,4)
        CNEQ(N_B,4,4) = C4NEQJ(N_B,1,2,1,2)
        CNEQ(N_B,4,5) = C4NEQJ(N_B,1,2,1,3)
        CNEQ(N_B,4,6) = C4NEQJ(N_B,1,2,2,3)
        CNEQ(N_B,5,1) = CNEQ(N_B,1,5)
        CNEQ(N_B,5,2) = CNEQ(N_B,2,5)
        CNEQ(N_B,5,3) = CNEQ(N_B,3,5)
        CNEQ(N_B,5,4) = CNEQ(N_B,4,5)
        CNEQ(N_B,5,5) = C4NEQJ(N_B,1,3,1,3)
        CNEQ(N_B,5,6) = C4NEQJ(N_B,1,3,2,3)
        CNEQ(N_B,6,1) = CNEQ(N_B,1,6)
        CNEQ(N_B,6,2) = CNEQ(N_B,2,6)
        CNEQ(N_B,6,3) = CNEQ(N_B,3,6)
        CNEQ(N_B,6,4) = CNEQ(N_B,4,6)
        CNEQ(N_B,6,5) = CNEQ(N_B,5,6)
        CNEQ(N_B,6,6) = C4NEQJ(N_B,2,3,2,3)
      END DO
! 
! ---------------------------------------------------------------------
!     STEP 6 - CALCULATION OF EQ CAUCHY STRESS & ELASTICITY TENSOR
! ---------------------------------------------------------------------
! 
!     CALCULATE LEFT CAUCHY-GREEN DEFORMATION TENSOR
      BTOT = MATMUL(DFGRD1, TRANSPOSE(DFGRD1))
! 
!     SETUP VOIGT NOTATION FOR BTOT 
      BTOTV(1) = BTOT(1,1)
      BTOTV(2) = BTOT(2,2)
      BTOTV(3) = BTOT(3,3)
      BTOTV(4) = BTOT(1,2)+BTOT(2,1)
      BTOTV(5) = BTOT(1,3)+BTOT(3,1)
      BTOTV(6) = BTOT(2,3)+BTOT(3,2)
! 
!     GET EIGEN VALUES AND EIGEN VECTORS OF BTOT
      DO I = 1,3
        DO J = 1,3
          BTOT_(I,J) = BTOT(I,J)
        END DO
      END DO     
      CALL DSYEVJ3(BTOT_, PDBTOT, PVBTOT)
      PDBTOT = TRANSPOSE(PDBTOT)
! 
!     CALCULATE EIGEN VALUES OF ISOCHORIC LEFT CAUCHY-GREEN
!     DEFORMATION TENSOR
      DO J=1,3
        PVBBAR(J) = DET**(-TWO/THREE)*PVBTOT(J)
      END DO
! _________________________________________________________
!         - Modification required when changing Weq(C) 
!     DEFINE INVARIANTS FOR THE DEVIATORIC EQ STRAIN ENERGY
! _________________________________________________________
      ALPHAEQ     = MU/TWO
      DALPHAEQDI1 = ZERO
      DALPHAEQDI2 = ZERO
      BETAEQ      = ZERO
      DBETAEQDI1  = ZERO
      DBETAEQDI2  = ZERO
!
!     CALCULATING PRINCIPAL VALUES OF DEVIATORIC KIRCHOFF STRESS
      PVTAUEQD(1) = ALPHAEQ * ((FOUR/THREE)*PVBBAR(1)
     1              - (TWO/THREE)*PVBBAR(2) 
     2              - (TWO/THREE)*PVBBAR(3))
     3            + BETAEQ * ((FOUR/THREE)*(PVBBAR(1)**TWO)
     4              - (TWO/THREE)*(PVBBAR(2)**TWO) 
     5              - (TWO/THREE)*(PVBBAR(3)**TWO))
      PVTAUEQD(2) = ALPHAEQ * ((FOUR/THREE)*PVBBAR(2)
     1              - (TWO/THREE)*PVBBAR(3) 
     2              - (TWO/THREE)*PVBBAR(1))
     3            + BETAEQ * ((FOUR/THREE)*(PVBBAR(2)**TWO)
     4              - (TWO/THREE)*(PVBBAR(3)**TWO) 
     5              - (TWO/THREE)*(PVBBAR(1)**TWO))
      PVTAUEQD(3) = ALPHAEQ * ((FOUR/THREE)*PVBBAR(3)
     1              - (TWO/THREE)*PVBBAR(1) 
     2              - (TWO/THREE)*PVBBAR(2))
     3            + BETAEQ * ((FOUR/THREE)*(PVBBAR(3)**TWO)
     4              - (TWO/THREE)*(PVBBAR(1)**TWO) 
     5              - (TWO/THREE)*(PVBBAR(2)**TWO))
! 
! _________________________________________________________
!         - Modification required when changing U(J) 
!     CALCULATING PRINCIPAL VALUES OF HYDROSTATIC KIRCHOFF STRESS
! _________________________________________________________
      PVTAUEQH = KELAS/TWO*(DET*DET - ONE)
! 
!     CALCULATING PRINCIPAL VALUES OF TOTAL KIRCHOFF STRESS
      PVTAUEQ(1) = PVTAUEQD(1) + PVTAUEQH
      PVTAUEQ(2) = PVTAUEQD(2) + PVTAUEQH
      PVTAUEQ(3) = PVTAUEQD(3) + PVTAUEQH
! 
!     CALCULATING DEVIATORIC, HYDROSTATIC AND TOTAL EQ CAUCHY STRESS TENSOR
      DO J=1,3
        DO K=1,3
          SIGMAEQD(K,J) = (PVTAUEQD(1)*PDBTOT(1,K)*PDBTOT(1,J)
     1                    +PVTAUEQD(2)*PDBTOT(2,K)*PDBTOT(2,J)
     2                    +PVTAUEQD(3)*PDBTOT(3,K)*PDBTOT(3,J))/DET
     	  SIGMAEQH(K,J) = (PVTAUEQH*PDBTOT(1,K)*PDBTOT(1,J)
     1                    +PVTAUEQH*PDBTOT(2,K)*PDBTOT(2,J)
     2                    +PVTAUEQH*PDBTOT(3,K)*PDBTOT(3,J))/DET
          SIGMAEQ(K,J) = SIGMAEQD(K,J) + SIGMAEQH(K,J)
        END DO
      END DO
! 
!     Calculating d(I1)/dEPS and d(I2)/dEPS
          DI1DEPS(1) = (FOUR/THREE)*PVBBAR(1)
     1                - (TWO/THREE)*PVBBAR(2)
     2                - (TWO/THREE)*PVBBAR(3)
          DI2DEPS(1) = (TWO/THREE)*PVBBAR(1) * 
     1                (PVBBAR(2) + PVBBAR(3))
     2              - (FOUR/THREE)*PVBBAR(2)*PVBBAR(3)
          DI1DEPS(2) = (FOUR/THREE)*PVBBAR(2)
     1                - (TWO/THREE)*PVBBAR(3)
     2                - (TWO/THREE)*PVBBAR(1)
          DI2DEPS(2) = (TWO/THREE)*PVBBAR(2) * 
     1                (PVBBAR(3) + PVBBAR(1))
     2              - (FOUR/THREE)*PVBBAR(3)*PVBBAR(1)
          DI1DEPS(3) = (FOUR/THREE)*PVBBAR(3)
     1                - (TWO/THREE)*PVBBAR(1)
     2                - (TWO/THREE)*PVBBAR(2)
          DI2DEPS(3) = (TWO/THREE)*PVBBAR(3) * 
     1                (PVBBAR(1) + PVBBAR(2))
     2              - (FOUR/THREE)*PVBBAR(1)*PVBBAR(2)
! 
! _________________________________________________________
!         - Modification required when changing U(J) 
!     CALCULATING HYDROSTATIC PART OF CAB
! _________________________________________________________
      CABHYD = KELAS*DET*DET 
! 
!     CALCULATING TOTAL CAB
      CAB(1,1) = (DALPHAEQDI1*DI1DEPS(1) + DALPHAEQDI2*DI2DEPS(1)) *
     1         ((FOUR/THREE)*PVBBAR(1)
     1          - (TWO/THREE)*(PVBBAR(2) + PVBBAR(3)))
     2         + ALPHAEQ *
     2           (((FOUR*FOUR)/NINE)*PVBBAR(1)
     2            + (FOUR/NINE)*(PVBBAR(2) + PVBBAR(3)))
     3         + (DBETAEQDI1*DI1DEPS(1) + DBETAEQDI2*DI2DEPS(1)) *
     4         ((FOUR/THREE)*(PVBBAR(1)**TWO)
     4          - (TWO/THREE)*((PVBBAR(2)**TWO) + (PVBBAR(3)**TWO)))
     5         + BETAEQ *
     5           (((FOUR*EIGHT)/NINE)*(PVBBAR(1)**TWO)
     5            + (EIGHT/NINE)*((PVBBAR(2)**TWO) + (PVBBAR(3)**TWO)))
     6         + CABHYD
      CAB(2,2) = (DALPHAEQDI1*DI1DEPS(2) + DALPHAEQDI2*DI2DEPS(2)) *
     1         ((FOUR/THREE)*PVBBAR(2)
     1          - (TWO/THREE)*(PVBBAR(3) + PVBBAR(1)))
     2         + ALPHAEQ * 
     2           (((FOUR*FOUR)/NINE)*PVBBAR(2)
     2            + (FOUR/NINE)*(PVBBAR(3) + PVBBAR(1)))
     3         + (DBETAEQDI1*DI1DEPS(2) + DBETAEQDI2*DI2DEPS(2)) *
     4         ((FOUR/THREE)*(PVBBAR(2)**TWO)
     4          - (TWO/THREE)*((PVBBAR(3)**TWO) + (PVBBAR(1)**TWO)))
     5         + BETAEQ * 
     5           (((FOUR*EIGHT)/NINE)*(PVBBAR(2)**TWO)
     5            + (EIGHT/NINE)*((PVBBAR(3)**TWO) + (PVBBAR(1)**TWO)))
     6         + CABHYD
      CAB(3,3) = (DALPHAEQDI1*DI1DEPS(3) + DALPHAEQDI2*DI2DEPS(3)) *
     1         ((FOUR/THREE)*PVBBAR(3)
     1          - (TWO/THREE)*(PVBBAR(1) + PVBBAR(2)))
     2         + ALPHAEQ * 
     2           (((FOUR*FOUR)/NINE)*PVBBAR(3)
     2            + (FOUR/NINE)*(PVBBAR(1) + PVBBAR(2)))
     3         + (DBETAEQDI1*DI1DEPS(3) + DBETAEQDI2*DI2DEPS(3)) *
     4         ((FOUR/THREE)*(PVBBAR(3)**TWO)
     4          - (TWO/THREE)*((PVBBAR(1)**TWO) + (PVBBAR(2)**TWO)))
     5         + BETAEQ *
     5           (((FOUR*EIGHT)/NINE)*(PVBBAR(3)**TWO)
     5            + (EIGHT/NINE)*((PVBBAR(1)**TWO) + (PVBBAR(2)**TWO)))
     6         + CABHYD
      CAB(1,2) = (DALPHAEQDI1*DI1DEPS(2) + DALPHAEQDI2*DI2DEPS(2)) *
     1         ((FOUR/THREE)*PVBBAR(1)
     1          - (TWO/THREE)*(PVBBAR(2) + PVBBAR(3)))
     2         + ALPHAEQ *
     2           ((FOUR/NINE)*PVBBAR(3)
     2            - (EIGHT/NINE)*(PVBBAR(1) + PVBBAR(2)) )
     3         + (DBETAEQDI1*DI1DEPS(2) + DBETAEQDI2*DI2DEPS(2)) *
     4         ((FOUR/THREE)*(PVBBAR(1)**TWO)
     4          - (TWO/THREE)*((PVBBAR(2)**TWO) + (PVBBAR(3)**TWO)))
     5         + BETAEQ *
     5           ((EIGHT/NINE)*(PVBBAR(3)**TWO)
     5            - ((TWO*EIGHT)/NINE)*((PVBBAR(1)**TWO) + (PVBBAR(2)**TWO)))
     6         + CABHYD
      CAB(1,3) = (DALPHAEQDI1*DI1DEPS(3) + DALPHAEQDI2*DI2DEPS(3)) *
     1         ((FOUR/THREE)*PVBBAR(1)
     1          - (TWO/THREE)*(PVBBAR(2) + PVBBAR(3)) )
     2         + ALPHAEQ * 
     2           ((FOUR/NINE)*PVBBAR(2)
     2            - (EIGHT/NINE)*(PVBBAR(3) + PVBBAR(1)) )
     3         + (DBETAEQDI1*DI1DEPS(3) + DBETAEQDI2*DI2DEPS(3)) *
     4         ((FOUR/THREE)*(PVBBAR(1)**TWO)
     4          - (TWO/THREE)*((PVBBAR(2)**TWO) + (PVBBAR(3)**TWO)) )
     5         + BETAEQ *
     5           ((EIGHT/NINE)*(PVBBAR(2)**TWO)
     5            - ((TWO*EIGHT)/NINE)*((PVBBAR(3)**TWO) + (PVBBAR(1)**TWO)) )
     6         + CABHYD
      CAB(2,3) = (DALPHAEQDI1*DI1DEPS(3) + DALPHAEQDI2*DI2DEPS(3)) *
     1         ((FOUR/THREE)*PVBBAR(2)
     1          - (TWO/THREE)*(PVBBAR(3) + PVBBAR(1)) )
     2         + ALPHAEQ * 
     2           ((FOUR/NINE)*PVBBAR(1)
     2            - (EIGHT/NINE)*(PVBBAR(2) + PVBBAR(3)) )
     3         + (DBETAEQDI1*DI1DEPS(3) + DBETAEQDI2*DI2DEPS(3)) *
     4         ((FOUR/THREE)*(PVBBAR(2)**TWO)
     4          - (TWO/THREE)*((PVBBAR(3)**TWO) + (PVBBAR(1)**TWO)) )
     5         + BETAEQ * 
     5           ((EIGHT/NINE)*(PVBBAR(1)**TWO) 
     5            - ((TWO*EIGHT)/NINE)*((PVBBAR(2)**TWO) + (PVBBAR(3)**TWO)) )
     6         + CABHYD
      CAB(2,1) = CAB(1,2)
      CAB(3,1) = CAB(1,3)
      CAB(3,2) = CAB(2,3)
! 
!     CALCULATING GAB 
      DO J=1,3
        DO K=1,3
          IF(J.EQ.K) THEN
            GAB(J,K) = FOUR * ALPHAEQ * PVBBAR(J)
     1               + EIGHT * BETAEQ * (PVBBAR(J)**TWO)
          ELSE
            IF (ABS(PVBTOT(J)-PVBTOT(K)).LT.EPS) THEN
	      GAB(J,K) = (DALPHAEQDI1*DI1DEPS(J)+DALPHAEQDI2*DI2DEPS(J)) *
     1                   (TWO*PVBBAR(J) - TWO*PVBBAR(K))
     2                 + ALPHAEQ * 
     2                   ((EIGHT/THREE)*PVBBAR(J) 
     2                   + (FOUR/THREE)*PVBBAR(K))
     3                 + (DBETAEQDI1*DI1DEPS(J)+DBETAEQDI2*DI2DEPS(J)) *
     3                   (TWO*(PVBBAR(J)**TWO) - TWO*(PVBBAR(K)**TWO))
     4                 + BETAEQ * 
     4                   (((TWO*EIGHT)/THREE)*(PVBBAR(J)**TWO) 
     4                   + (EIGHT/THREE)*(PVBBAR(K)**TWO))
            ELSE
              GAB(J,K) = (PVTAUEQ(J)*PVBBAR(K) - PVTAUEQ(K)*PVBBAR(J))
     1                 / (PVBBAR(J) - PVBBAR(K))     
            END IF
          END IF
        END DO
      END DO
! 
!     CALCULATE MATERIAL TENSOR COEFFICIENTS (C4EQMAT)
!     EXPRESSED IN EIGEN BASIS 
      C4EQMAT = ZERO
!     1. CALCULATE COEFFICIENTS FOR AABB
      DO J=1,3
        DO K=1,3
          C4EQMAT(J,J,K,K) = (CAB(J,K) - TWO*PVTAUEQ(J)*IDT2(J,K))
     1                       /DET
        END DO
      END DO 
!     2. CALCULATE COEFFICIENTS FOR ABAB / ABBA
      DO J=1,3
        DO K=1,3
          IF (J.EQ.K) THEN
            CYCLE
          END IF
          C4EQMAT(J,K,J,K) = GAB(J,K) / (TWO*DET)
          C4EQMAT(J,K,K,J) = C4EQMAT(J,K,J,K)
        END DO
      END DO
!     
!     CALCULATE MATERIAL TENSOR (C4EQ) IN CURRENT CONFIGURATION
!     EXPRESSED IN STANDARD BASIS
      C4EQ = ZERO
      DO J=1,3
        DO K=1,3
          DO L=1,3
            DO M=1,3
              C4EQ(J,K,L,M) = 
     1          C4EQMAT(1,1,1,1)*PDBTOT(1,J)*PDBTOT(1,K)
     1                          *PDBTOT(1,L)*PDBTOT(1,M)
     1         +C4EQMAT(1,1,2,2)*PDBTOT(1,J)*PDBTOT(1,K)
     1                          *PDBTOT(2,L)*PDBTOT(2,M)
     1         +C4EQMAT(1,1,3,3)*PDBTOT(1,J)*PDBTOT(1,K)
     1                          *PDBTOT(3,L)*PDBTOT(3,M)
     1         +C4EQMAT(2,2,1,1)*PDBTOT(2,J)*PDBTOT(2,K)
     1                          *PDBTOT(1,L)*PDBTOT(1,M)
     1         +C4EQMAT(2,2,2,2)*PDBTOT(2,J)*PDBTOT(2,K)
     1                          *PDBTOT(2,L)*PDBTOT(2,M)
     1         +C4EQMAT(2,2,3,3)*PDBTOT(2,J)*PDBTOT(2,K)
     1                          *PDBTOT(3,L)*PDBTOT(3,M)
     1         +C4EQMAT(3,3,1,1)*PDBTOT(3,J)*PDBTOT(3,K)
     1                          *PDBTOT(1,L)*PDBTOT(1,M)
     1         +C4EQMAT(3,3,2,2)*PDBTOT(3,J)*PDBTOT(3,K)
     1                          *PDBTOT(2,L)*PDBTOT(2,M)
     1         +C4EQMAT(3,3,3,3)*PDBTOT(3,J)*PDBTOT(3,K)
     1                          *PDBTOT(3,L)*PDBTOT(3,M)
     1         +C4EQMAT(1,2,1,2)*PDBTOT(1,J)*PDBTOT(2,K)
     1                          *PDBTOT(1,L)*PDBTOT(2,M)
     1         +C4EQMAT(1,2,2,1)*PDBTOT(1,J)*PDBTOT(2,K)
     1                          *PDBTOT(2,L)*PDBTOT(1,M)
     1         +C4EQMAT(1,3,1,3)*PDBTOT(1,J)*PDBTOT(3,K)
     1                          *PDBTOT(1,L)*PDBTOT(3,M)
     1         +C4EQMAT(1,3,3,1)*PDBTOT(1,J)*PDBTOT(3,K)
     1                          *PDBTOT(3,L)*PDBTOT(1,M)
     1         +C4EQMAT(2,1,2,1)*PDBTOT(2,J)*PDBTOT(1,K)
     1                          *PDBTOT(2,L)*PDBTOT(1,M)
     1         +C4EQMAT(2,1,1,2)*PDBTOT(2,J)*PDBTOT(1,K)
     1                          *PDBTOT(1,L)*PDBTOT(2,M)
     1         +C4EQMAT(2,3,2,3)*PDBTOT(2,J)*PDBTOT(3,K)
     1                          *PDBTOT(2,L)*PDBTOT(3,M)
     1         +C4EQMAT(2,3,3,2)*PDBTOT(2,J)*PDBTOT(3,K)
     1                          *PDBTOT(3,L)*PDBTOT(2,M)
     1         +C4EQMAT(3,1,3,1)*PDBTOT(3,J)*PDBTOT(1,K)
     1                          *PDBTOT(3,L)*PDBTOT(1,M)
     1         +C4EQMAT(3,1,1,3)*PDBTOT(3,J)*PDBTOT(1,K)
     1                          *PDBTOT(1,L)*PDBTOT(3,M)
     1         +C4EQMAT(3,2,3,2)*PDBTOT(3,J)*PDBTOT(2,K)
     1                          *PDBTOT(3,L)*PDBTOT(2,M)
     1         +C4EQMAT(3,2,2,3)*PDBTOT(3,J)*PDBTOT(2,K)
     1                          *PDBTOT(2,L)*PDBTOT(3,M)
            END DO
          END DO
        END DO
      END DO
! 
!     CALCLUATE JAUMANN RATE FORM OF MATERIAL TENSOR
      C4EQJ = ZERO
      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
              C4EQJ(I,J,K,L) = C4EQ(I,J,K,L)/DET
     1                        + (1 - DAMAGED_OLD) * (IDT2(I,K)*SIGMAEQ(J,L)
     1                         +IDT2(I,L)*SIGMAEQ(J,K)
     1                         +IDT2(J,K)*SIGMAEQ(I,L)
     1                         +IDT2(J,L)*SIGMAEQ(I,K))/(TWO)        
            END DO
          END DO
        END DO
      END DO
! 
!     GENERATE ABAQUS EQ TANGENT STIFFNESS MATRIX (VOIGT NOTATION)
      CEQ = ZERO
      CEQ(1,1) = C4EQJ(1,1,1,1)
      CEQ(1,2) = C4EQJ(1,1,2,2)
      CEQ(1,3) = C4EQJ(1,1,3,3)
      CEQ(1,4) = C4EQJ(1,1,1,2)
      CEQ(1,5) = C4EQJ(1,1,1,3)
      CEQ(1,6) = C4EQJ(1,1,2,3)
      CEQ(2,1) = CEQ(1,2)
      CEQ(2,2) = C4EQJ(2,2,2,2)
      CEQ(2,3) = C4EQJ(2,2,3,3)
      CEQ(2,4) = C4EQJ(2,2,1,2)
      CEQ(2,5) = C4EQJ(2,2,1,3)
      CEQ(2,6) = C4EQJ(2,2,2,3)
      CEQ(3,1) = CEQ(1,3)
      CEQ(3,2) = CEQ(2,3)
      CEQ(3,3) = C4EQJ(3,3,3,3)
      CEQ(3,4) = C4EQJ(3,3,1,2)
      CEQ(3,5) = C4EQJ(3,3,1,3)
      CEQ(3,6) = C4EQJ(3,3,2,3)
      CEQ(4,1) = CEQ(1,4)
      CEQ(4,2) = CEQ(2,4)
      CEQ(4,3) = CEQ(3,4)
      CEQ(4,4) = C4EQJ(1,2,1,2)
      CEQ(4,5) = C4EQJ(1,2,1,3)
      CEQ(4,6) = C4EQJ(1,2,2,3)
      CEQ(5,1) = CEQ(1,5)
      CEQ(5,2) = CEQ(2,5)
      CEQ(5,3) = CEQ(3,5)
      CEQ(5,4) = CEQ(4,5)
      CEQ(5,5) = C4EQJ(1,3,1,3)
      CEQ(5,6) = C4EQJ(1,3,2,3)
      CEQ(6,1) = CEQ(1,6)
      CEQ(6,2) = CEQ(2,6)
      CEQ(6,3) = CEQ(3,6)
      CEQ(6,4) = CEQ(4,6)
      CEQ(6,5) = CEQ(5,6)
      CEQ(6,6) = C4EQJ(2,3,2,3)
! 
! ---------------------------------------------------------------------
!     STEP 7 - UPDATE OF THE DAMAGE VALUES
! ---------------------------------------------------------------------
!     INITIALISATION FOR THE CURRENT DAMAGE VARIABLES AND DERIVATIVES
      DAMVARD0   = ZERO
      DAMVARD1   = ZERO
      DAMVARD2   = ZERO
      DAMVARD3   = ZERO
      DAMVARD    = ZERO
      DDAMVARDDT = ZERO
      DAMAGED    = ZERO
      DAMVARH    = ZERO
      DDAMVARHDT = ZERO
      DAMAGEH    = ZERO
! 
!     CALCULATE JACOBIAN OF OLD DEFORMATION GRADIENT
      DETOLD = DFGRD0(1,1)*DFGRD0(2,2)*DFGRD0(3,3)
     1        -DFGRD0(1,2)*DFGRD0(2,1)*DFGRD0(3,3)
      IF(NSHR.EQ.3) THEN
        DETOLD=DETOLD+DFGRD0(1,2)*DFGRD0(2,3)*DFGRD0(3,1)
     1         +DFGRD0(1,3)*DFGRD0(3,2)*DFGRD0(2,1)
     2         -DFGRD0(1,3)*DFGRD0(3,1)*DFGRD0(2,2)
     3         -DFGRD0(2,3)*DFGRD0(3,2)*DFGRD0(1,1)
      END IF
! 
!     CALCULATE OLD LEFT CAUCHY-GREEN DEFORMATION TENSOR
      BTOTOLD = MATMUL(DFGRD0, TRANSPOSE(DFGRD0))
! 
!     GET EIGEN VALUES AND EIGEN VECTORS OF BTOTOLD
      DO I = 1,3
        DO J = 1,3
          BTOTOLD_(I,J) = BTOTOLD(I,J)
        END DO
      END DO     
      CALL DSYEVJ3(BTOTOLD_, PDBTOTOLD, PVBTOTOLD)
! 
!     CALCULATE EIGEN VALUES OF ISOCHORIC LEFT CAUCHY-GREEN
!     DEFORMATION TENSOR
      DO J=1,3
        PVBBAROLD(J) = DETOLD**(-TWO/THREE)*PVBTOTOLD(J)
      END DO
! 
!     CALCULATE DEVIATORIC DAMAGE INTERNAL STATE VARIABLE
      CALL DAMAGEVAR(DAMCHOICE, PVBBAR, DAMVARD1, DAMVARINIT)
      CALL DAMAGEVAR(DAMCHOICE, PVBBAROLD, DAMVARD0, DAMVARINIT)
!     CALCULATE TIME DERIVATIVE OF DEVIATORIC DAMAGE VARIABLE
      DDAMVARDDT = ( ( DAMVARD1 - DAMVARD0 ) / DTIME ) *
     1             ( ONE - W_PH * ( ONE - EXP(-PRESSURE/P_SAT) ) )
!     CALCULATE DEVIATORIC DAMAGE VARIABLE
      DAMVARD = DTIME * DDAMVARDDT + DAMVARD_OLD
!     CHECK DEVIATORIC DAMAGE INTERNAL STATE VARIABLE EVOLUTION
      IF (DAMVARD.GT.DAMVARD_OLD) THEN
! 	UPDATE OF DEVIATORIC DAMAGE VARIABLE (STATEV)
        STATEV(NSTATV-3) = DAMVARD
! 	CALCULATE TIME DERIVATIVE OF DEVIATORIC DAMAGE FUNCTION
	DDAMAGEDDT = AEXPD * BEXPD * ((STATEV(NSTATV-3) - DAMVARINIT)**(AEXPD - ONE)) *
     1             EXP( - BEXPD*( ( STATEV(NSTATV-3) - DAMVARINIT )**AEXPD ) ) *
     2             DDAMVARDDT * ( ONE - W_PH * ( ONE - EXP(-PRESSURE/P_SAT) ) )
! 	CALCULATE DEVIATORIC DAMAGE FUNCTION 
!	DAMAGED = DTIME * DDAMAGEDDT + DAMAGED_OLD
	DAMAGED = ONE - EXP( - BEXPD*( ( STATEV(NSTATV-3) - DAMVARINIT )**AEXPD ) )
      ELSE
! 	UPDATE OF DEVIATORIC DAMAGE VARIABLE (STATEV-3)
        STATEV(NSTATV-3) = DAMVARD_OLD
! 	CALCULATE DEVIATORIC DAMAGE FUNCTION 
	DAMAGED = DAMAGED_OLD
      END IF
!     UPDATE OF DEVIATORIC DAMAGE FUNCTION (STATEV-1)
      STATEV(NSTATV-1) = DAMAGED
! 
!     CALCULATE TIME DERIVATIVE OF DEVIATORIC DAMAGE VARIABLE
      DDAMVARHDT = ( ( DAMVARD1 - DAMVARD0 ) / DTIME ) *
     1             ( ONE - (NINE/TEN) * ( ONE - EXP(-PRESSURE/P_SAT) ) )
!     CALCULATE HYDROSTATIC DAMAGE INTERNAL STATE VARIABLE
      DAMVARH = DTIME * DDAMVARHDT + DAMVARH_OLD
      IF (DAMVARH.GT.DAMVARH_OLD) THEN
! 	UPDATE OF DEVIATORIC DAMAGE VARIABLE (STATEV)
        STATEV(NSTATV-2) = DAMVARH
!     	CALCULATE HYDROSTATIC DAMAGE FUNCTION
      	DAMAGEH = ONE - EXP( - BEXPH*( ( STATEV(NSTATV-2) - DAMVARINIT )**AEXPH ) )
      ELSE
! 	UPDATE OF HYDROSTATIC DAMAGE VARIABLE (STATEV)
        STATEV(NSTATV-2) = DAMVARH_OLD
! 	CALCULATE HYDROSTATIC DAMAGE FUNCTION 
	DAMAGEH = DAMAGEH_OLD
      END IF
!     UPDATE OF HYDROSTATIC DAMAGE FUNCTION (STATEV)
      STATEV(NSTATV) = DAMAGEH
! 
! ---------------------------------------------------------------------
!     STEP 8 - CALCULATION OF TOTAL CAUCHY STRESS & ELASTICITY TENSOR
! ---------------------------------------------------------------------
!     CALCULATE TOTAL CAUCHY STRESS TENSOR
      STRESSTOT = ZERO
      DO J=1,3
        DO K=1,3
          STRESSTOT(J,K) = (ONE-DAMAGED)*SIGMAEQD(J,K)
     1                   + (ONE-DAMAGEH)*SIGMAEQH(J,K)
!         VISCOELASTIC LOOP
          DO N_B = 1, N_BRANCHES
            STRESSTOT(J,K) = STRESSTOT(J,K) + (ONE-DAMAGED)*TAUNEQD(N_B,J,K)/DET
     1                     + (ONE-DAMAGEH)*TAUNEQH(N_B,J,K)/DET
          END DO
        END DO
      END DO

!     VOIGT NOTATION FOR CAUCHY STRESS
      STRESS(1) = STRESSTOT(1,1)
      STRESS(2) = STRESSTOT(2,2)
      STRESS(3) = STRESSTOT(3,3)
      STRESS(4) = STRESSTOT(1,2)
      STRESS(5) = STRESSTOT(1,3)
      STRESS(6) = STRESSTOT(2,3)
!  
!     CALCULATE TOTAL ABAQUS TANGENT STIFNESS MATRIX
      DO J=1,6
        DO K=1,6
          DDSDDE(J,K) = CEQ(J,K)
!         VISCOELASTIC LOOP
          DO N_B = 1, N_BRANCHES
            DDSDDE(J,K) = DDSDDE(J,K) + CNEQ(N_B,J,K)
          END DO
        END DO
      END DO
 
      RETURN
      END SUBROUTINE UMAT

      
! ---------------------------------------------------------------------
!     DAMAGE SUBROUTINES
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
!     DVAR  -  Compute the choosen damage variable.
! 
!     DAMCHOICE  -  Input  - Choice of the damage variable
!     PVBBAR     -  Input  - Eigen values of the BBAR tensor
!     DAMVAR     -  Output  - Value of the calculated damage variable
! ---------------------------------------------------------------------
      SUBROUTINE DAMAGEVAR (DAMCHOICE, PVBBAR, DAMVAR, DAMVARINIT)
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: DAMCHOICE
            DOUBLE PRECISION, DIMENSION(3), INTENT(IN) :: PVBBAR
            DOUBLE PRECISION, INTENT(OUT) :: DAMVAR, DAMVARINIT
            
            INTEGER J
            DOUBLE PRECISION :: I1, I2
            DOUBLE PRECISION :: ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, FOUR=4.D0
            DOUBLE PRECISION :: FIVE=5.D0, SIX=6.D0, NINE=9.D0
!            
!           INITIALISATION FOR USEFUL QUANTITIES        
            I1 = ZERO
            I2 = ZERO 
! 
!           DAMAGEVAR 0: NO DAMAGE
            IF (DAMCHOICE == 0) THEN
              DAMVAR = ZERO
              DAMVARINIT = ZERO
              RETURN    
!           DAMAGEVAR 1: sqrt(I1/3) - 1             
            ELSEIF (DAMCHOICE == 1) THEN
              DO J=1,3
                I1 = I1 + PVBBAR(J)
              END DO
              DAMVAR = SQRT( (I1/THREE) ) - ONE
              DAMVARINIT = ZERO
              RETURN
!           DAMAGEVAR 2: magnitude = sqrt(I1**2 - 2*I2)                
            ELSEIF (DAMCHOICE == 2) THEN
              DO J=1,3
                DAMVAR = DAMVAR + PVBBAR(J)**TWO
              END DO
              DAMVAR = SQRT(DAMVAR)
              DAMVARINIT = SQRT(THREE)
              RETURN
!           DAMAGEVAR 3: Igamma = (1/6)*sqrt( 2*(I1**2) - 6*I2 )                   
            ELSEIF (DAMCHOICE == 3) THEN
              DO J=1,3
                I1 = I1 + PVBBAR(J)
                DAMVAR = DAMVAR + PVBBAR(J)**TWO
              END DO
              DAMVAR = SQRT( TWO*(I1**TWO) - THREE*( (I1**TWO) - DAMVAR ) )
              DAMVAR = (ONE/SIX)*DAMVAR
              DAMVARINIT = ZERO
              RETURN
!           DAMAGEVAR 4: heq = sqrt( (2/3)*(h1**2 + h2**2 + h3**2) )  
            ELSEIF (DAMCHOICE == 4) THEN
              DO J=1,3
                DAMVAR = DAMVAR + (TWO/THREE) * ( ( (ONE/TWO) * LOG(PVBBAR(J)) )**TWO )
              END DO
              DAMVAR = SQRT(DAMVAR)
              DAMVARINIT = ZERO
              RETURN                            
            END IF

      END SUBROUTINE DAMAGEVAR

! ---------------------------------------------------------------------
!     UTILITY SUBROUTINES
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
!     M33INV  -  Compute the Inverse of a 3x3 Matrix.
!     SOURCE:  David G. Simpson - NASA Goddard Space Flight Center
! 
!     A       - Input  - 3x3 Matrix to be Inverted
!     AINV    - Output - 3x3 Inverse of Matrix A
!     OK_FLAG - .TRUE. If A is non-singular and A is Inverted
! ---------------------------------------------------------------------
      SUBROUTINE M33INV (A, AINV, OK_FLAG)
            IMPLICIT NONE
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
            DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
            LOGICAL, INTENT(OUT) :: OK_FLAG

            DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-16
            DOUBLE PRECISION :: DET
            DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR
! 
! 
            DET = A(1,1) * A(2,2) * A(3,3)  
     1          - A(1,1) * A(2,3) * A(3,2)  
     2          - A(1,2) * A(2,1) * A(3,3)  
     3          + A(1,2) * A(2,3) * A(3,1)  
     4          + A(1,3) * A(2,1) * A(3,2)  
     5          - A(1,3) * A(2,2) * A(3,1)
! 
            IF (ABS(DET) .LE. EPS) THEN
            AINV = 0.0D0
            OK_FLAG = .FALSE.
            RETURN
            END IF
! 
            COFACTOR(1,1) = +(A(2,2) * A(3,3) - A(2,3) * A(3,2))
            COFACTOR(1,2) = -(A(2,1) * A(3,3) - A(2,3) * A(3,1))
            COFACTOR(1,3) = +(A(2,1) * A(3,2) - A(2,2) * A(3,1))
            COFACTOR(2,1) = -(A(1,2) * A(3,3) - A(1,3) * A(3,2))
            COFACTOR(2,2) = +(A(1,1) * A(3,3) - A(1,3) * A(3,1))
            COFACTOR(2,3) = -(A(1,1) * A(3,2) - A(1,2) * A(3,1))
            COFACTOR(3,1) = +(A(1,2) * A(2,3) - A(1,3) * A(2,2))
            COFACTOR(3,2) = -(A(1,1) * A(2,3) - A(1,3) * A(2,1))
            COFACTOR(3,3) = +(A(1,1) * A(2,2) - A(1,2) * A(2,1))
! 
            AINV = TRANSPOSE(COFACTOR) / DET
! 
            OK_FLAG = .TRUE.
            RETURN
! 
      END SUBROUTINE M33INV

* ----------------------------------------------------------------------------
      SUBROUTINE DSYEVJ3(A, Q, W)
* ----------------------------------------------------------------------------
* From Kopp (2008) doi :10.1142/S0129183108012303
* Calculates the eigenvalues and normalized eigenvectors of a symmetric 3x3
* matrix A using the Jacobi algorithm.
* The upper triangular part of A is destroyed during the calculation,
* the diagonal elements are read but not destroyed, and the lower
* triangular elements are not referenced at all.
* ----------------------------------------------------------------------------
* Parameters:
*   A: The symmetric input matrix
*   Q: Storage buffer for eigenvectors
*   W: Storage buffer for eigenvalues
* ----------------------------------------------------------------------------
*     .. Arguments ..
      DOUBLE PRECISION A(3,3)
      DOUBLE PRECISION Q(3,3)
      DOUBLE PRECISION W(3)

*     .. Parameters ..
      INTEGER          N
      PARAMETER        ( N = 3 )
    
*     .. Local Variables ..
      DOUBLE PRECISION SD, SO
      DOUBLE PRECISION S, C, T
      DOUBLE PRECISION G, H, Z, THETA
      DOUBLE PRECISION THRESH
      INTEGER          I, X, Y, R

*     Initialize Q to the identitity matrix
*     --- This loop can be omitted if only the eigenvalues are desired ---
      DO 10 X = 1, N
        Q(X,X) = 1.0D0
        DO 11, Y = 1, X-1
          Q(X, Y) = 0.0D0
          Q(Y, X) = 0.0D0
   11   CONTINUE
   10 CONTINUE

*     Initialize W to diag(A)
      DO 20 X = 1, N
        W(X) = A(X, X)
   20 CONTINUE

*     Calculate SQR(tr(A))  
      SD = 0.0D0
      DO 30 X = 1, N
        SD = SD + ABS(W(X))
   30 CONTINUE
      SD = SD**2
 
*     Main iteration loop
      DO 40 I = 1, 50
*       Test for convergence
        SO = 0.0D0
        DO 50 X = 1, N
          DO 51 Y = X+1, N
            SO = SO + ABS(A(X, Y))
   51     CONTINUE
   50   CONTINUE
        IF (SO .EQ. 0.0D0) THEN
          RETURN
        END IF

        IF (I .LT. 4) THEN
          THRESH = 0.2D0 * SO / N**2
        ELSE
          THRESH = 0.0D0
        END IF

*       Do sweep
        DO 60 X = 1, N
          DO 61 Y = X+1, N
            G = 100.0D0 * ( ABS(A(X, Y)) )
            IF ( I .GT. 4 .AND. ABS(W(X)) + G .EQ. ABS(W(X))
     $                    .AND. ABS(W(Y)) + G .EQ. ABS(W(Y)) ) THEN
              A(X, Y) = 0.0D0
            ELSE IF (ABS(A(X, Y)) .GT. THRESH) THEN
*             Calculate Jacobi transformation
              H = W(Y) - W(X)
              IF ( ABS(H) + G .EQ. ABS(H) ) THEN
                T = A(X, Y) / H
              ELSE
                THETA = 0.5D0 * H / A(X, Y)
                IF (THETA .LT. 0.0D0) THEN
                  T = -1.0D0 / (SQRT(1.0D0 + THETA**2) - THETA)
                ELSE
                  T = 1.0D0 / (SQRT(1.0D0 + THETA**2) + THETA)
                END IF
              END IF

              C = 1.0D0 / SQRT( 1.0D0 + T**2 )
              S = T * C
              Z = T * A(X, Y)
              
*             Apply Jacobi transformation
              A(X, Y) = 0.0D0
              W(X)    = W(X) - Z
              W(Y)    = W(Y) + Z
              DO 70 R = 1, X-1
                T       = A(R, X)
                A(R, X) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   70         CONTINUE
              DO 80, R = X+1, Y-1
                T       = A(X, R)
                A(X, R) = C * T - S * A(R, Y)
                A(R, Y) = S * T + C * A(R, Y)
   80         CONTINUE
              DO 90, R = Y+1, N
                T       = A(X, R)
                A(X, R) = C * T - S * A(Y, R)
                A(Y, R) = S * T + C * A(Y, R)
   90         CONTINUE

*             Update eigenvectors
*             --- This loop can be omitted if only the eigenvalues are desired ---
              DO 100, R = 1, N
                T       = Q(R, X)
                Q(R, X) = C * T - S * Q(R, Y)
                Q(R, Y) = S * T + C * Q(R, Y)
  100         CONTINUE
            END IF
   61     CONTINUE
   60   CONTINUE
   40 CONTINUE

      PRINT *, "DSYEVJ3: No convergence."
            
      END SUBROUTINE
* End of subroutine DSYEVJ3

