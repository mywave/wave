MODULE WAM_INTERFACE_MODULE

! ---------------------------------------------------------------------------- !
!          THIS MODULE COLLECTS PROCEDURES, WHICH ARE USED IN THE WAM MODEL    !
!          TO COMPUTE PARAMETERS FROM SPECTRA AND TO INTERPOLATE SPECTRA.      !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     MODULE VARIABLES.                                                        !
!     -----------------                                                        !

USE WAM_FRE_DIR_MODULE, ONLY: ML, KL, CO, TFAK, FR, DFIM, SINTH, COSTH,        &
&                             DFIMOFR, DFFR, DFFR2, DELTH,                     &
&                             MO_TAIL, MM1_TAIL, MP1_TAIL, MP2_TAIL

USE WAM_GENERAL_MODULE, ONLY: G, PI, ZPI

IMPLICIT NONE

PRIVATE

REAL      :: EMIN = 1.0E-12    !! REPLACES THE INTRINSIC TINY

! ---------------------------------------------------------------------------- !
!                                                                              !
!     A.  GENERIC INTERFACES (THIS MODULE CONTAINS THE PROCEDURES).            !
!                                                                              !
! ---------------------------------------------------------------------------- !

INTERFACE FEMEAN                           !! MEAN FREQUENCY AND WAVE NUMBER.
   MODULE PROCEDURE FEMEAN
END INTERFACE
PUBLIC FEMEAN

INTERFACE INTSPEC                          !! INTERPOLATES SPECTRA.
   MODULE PROCEDURE INTSPEC
END INTERFACE
PUBLIC INTSPEC

INTERFACE MEAN_DIRECTION                   !! MEAN DIRECTION AND SPREAD
   MODULE PROCEDURE MEAN_DIRECTION_1       !! SCALAR VERSION
   MODULE PROCEDURE MEAN_DIRECTION_B       !! VECTOR VERSION
END INTERFACE
PUBLIC MEAN_DIRECTION

INTERFACE PEAK_PERIOD                      !! COMPUTATES PEAK PERIOD.
   MODULE  PROCEDURE PEAK_PERIOD_1         !! SCALAR VERSION
   MODULE  PROCEDURE PEAK_PERIOD_B         !! VECTOR VERSION
END INTERFACE
PUBLIC PEAK_PERIOD

INTERFACE ROTSPEC                          !! ROTATE A SPECTRUM.
   MODULE PROCEDURE ROTSPEC
END INTERFACE
PUBLIC ROTSPEC

INTERFACE STRSPEC                          !! STRETCH A SPECTRUM.
   MODULE PROCEDURE STRSPEC
END INTERFACE
PUBLIC STRSPEC

INTERFACE STOKES_DRIFT                     !! COMPUTATES STOKES DRIFT.
   MODULE PROCEDURE STOKES_DRIFT
END INTERFACE
PUBLIC STOKES_DRIFT

INTERFACE TM1_TM2_PERIODS                  !! COMPUTATES TM1 AND/OR TM2 PERIODS.
   MODULE  PROCEDURE TM1_TM2_PERIODS_1     !! SCALAR VERSION
   MODULE  PROCEDURE TM1_TM2_PERIODS_B     !! VECTOR VERSION
END INTERFACE
PUBLIC TM1_TM2_PERIODS

INTERFACE TOTAL_ENERGY                     !! COMPUTES TOTAL ENERGY.
   MODULE  PROCEDURE TOTAL_ENERGY_1        !! SCALAR VERSION
   MODULE  PROCEDURE TOTAL_ENERGY_B        !! VECTOR VERSION
END INTERFACE
PUBLIC TOTAL_ENERGY

INTERFACE COS2_SPR                         !! COSINE SQUARE SPREAD.
   MODULE  PROCEDURE COS2_SPR_1            !! SCALAR VERSION
   MODULE  PROCEDURE COS2_SPR_B            !! VECTOR VERSION
END INTERFACE
PUBLIC COS2_SPR

INTERFACE WM1_WM2_WAVENUMBER               !! WM1 AND/OR WM2 WAVENUMBER.
   MODULE  PROCEDURE WM1_WM2_WAVENUMBER_1  !! SCALAR VERSION
   MODULE  PROCEDURE WM1_WM2_WAVENUMBER_B  !! VECTOR VERSION
END INTERFACE
PUBLIC WM1_WM2_WAVENUMBER

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

CONTAINS

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE FEMEAN (F, EMEAN, FM, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   FEMEAN - COMPUTATION OF MEAN FREQUENCY.                                    !
!                                                                              !
!     S.D. HASSELMANN                                                          !
!     OPTIMIZED BY : L. ZAMBRESKY AND H. GUENTHER                              !
!     H. GUNTHER     GKSS         DECEMBER 2001    FT90                        !
!                                                                              !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE MEAN FREQUENCY AT EACH GRID POINT.                             !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)            :: F(:,:,:)      !! BLOCK OF SPECTRA.
REAL,    INTENT(IN)            :: EMEAN(:)      !! TOTAL ENERGY.
REAL,    INTENT(OUT)           :: FM   (:)      !! MEAN FREQUENCY.
LOGICAL, INTENT(IN),  OPTIONAL :: MASK (:,:,:)  !! INTERATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL     :: TEMP2(SIZE(F,1),SIZE(F,3))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !

IF (PRESENT(MASK)) THEN
   TEMP2 = SUM(F,DIM=2, MASK=MASK)
ELSE
   TEMP2 = SUM(F,DIM=2)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES.                                           !
!        ---------------------------                                           !

FM = MATMUL(TEMP2,DFIMOFR)            !! INTEGRATE OVER FREQUENCIES.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. ADD TAIL ENERGY.                                                      !
!        ----------------                                                      !

FM = FM + MM1_TAIL*TEMP2(:,SIZE(F,3)) !! ADD TAIL.
FM = EMEAN/MAX(FM,EMIN)               !! NORMALIZE WITH TOTAL ENERGY.
 
END SUBROUTINE FEMEAN

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE INTSPEC (DEL12, DEL1L,  F1, FMEAN1, EMEAN1, THETM1,                 &
&                   F2, FMEAN2, EMEAN2, THETM2, FL, FMEAN, EMEAN, THETM )

! ---------------------------------------------------------------------------- !
!                                                                              !
!   INTSPEC  -  INTERPOLATION OF SPECTRA.                                      !
!                                                                              !
!     SUSANNE HASSELMANN  MPI        JUNE 1990.                                !
!     H. GUNTHER          GKSS/ECMWF JAN. 1991   MODIFIED FOR CYCLE_4          !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       INTERPOLATION OF SPECTRA.                                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       ROTATE SPECTRA ACCORDING TO MEAN OF MEAN ANGLES, TRANSFORM             !
!       FREQUENCIES ACCORDING TO MEAN OF MEAN FREQUENCIES ,ADJUST ENERGY       !
!       ACCORDCING TO MEAN OF TOTAL ENERGY AND INTERPOLATE RESULTING           !
!       SPECTRA.                                                               !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       K.HASSELMANN, 1990,                                                    !
!          INTERPOLATION OF WAVE SPECTRA. WAM NOTE 6/6/90.                     !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)  :: DEL12         !! DISTANCE SPECTRUM 2 - SPECTRUM 1.
REAL,    INTENT(IN)  :: DEL1L         !! DISTANCE SPECTRUM L - SPECTRUM 1.
REAL,    INTENT(IN)  :: F1(:,:)       !! SPECTRUM 1.
REAL,    INTENT(IN)  :: FMEAN1        !! MEAN FREQUENCY OF F1.
REAL,    INTENT(IN)  :: EMEAN1        !! MEAN ENERGY OF F1.
REAL,    INTENT(IN)  :: THETM1        !! MEAN DIRECTION OF F1.
REAL,    INTENT(IN)  :: F2(:,:)       !! SPECTRUM 2.
REAL,    INTENT(IN)  :: FMEAN2        !! MEAN FREQUENCY OF F2.
REAL,    INTENT(IN)  :: EMEAN2        !! MEAN ENERGY OF F2.
REAL,    INTENT(IN)  :: THETM2        !! MEAN DIRECTION OF F2.
REAL,    INTENT(OUT) :: FL(: ,:)      !! INTEPOLATED SPECTRUM.
REAL,    INTENT(OUT) :: FMEAN         !! INTEPOLATED MEAN FREQUENCY.
REAL,    INTENT(OUT) :: EMEAN         !! INTEPOLATED MEAN ENERGY.
REAL,    INTENT(OUT) :: THETM         !! INTEPOLATED MEAN DIRECTION.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL :: GW1, GW2
REAL :: F_L(SIZE(F1,1),SIZE(F1,2)), F_R(SIZE(F1,1),SIZE(F1,2))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTERPOLATION WEIGHTS.                                                !
!        ----------------------                                                !

GW2 = DEL1L/DEL12
GW1 = 1. - GW2

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. LEFT WEIGHT OR ENERGY OF LEFT SPECTRUM IS ZERO.                       !
!        -----------------------------------------------                       !

IF (ABS(GW1).LT.EPSILON(1.) .OR. EMEAN1.LT.EPSILON(1.)) THEN
   FL = GW2*F2
   EMEAN = GW2*EMEAN2
   FMEAN = GW2*FMEAN2
   THETM = GW2*THETM2
   RETURN
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. RIGHT WEIGHT OR ENERGY OF RIGHT SPECTRUM IS ZERO.                     !
!        -------------------------------------------------                     !

IF (ABS(GW2).LT.EPSILON(1.) .OR. EMEAN2.LT.EPSILON(1.)) THEN
   FL = GW1*F1
   EMEAN = GW1*EMEAN1
   FMEAN = GW1*FMEAN1
   THETM = GW1*THETM1
   RETURN
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. ENERGY AND WEIGHTS OF BOTH SPECTRA ARE GT ZERO.                       !
!        -----------------------------------------------                       !

!     3.1 INTERPOLATE MEAN VALUES.                                             !

EMEAN = GW1*EMEAN1+GW2*EMEAN2
FMEAN = GW1*FMEAN1+GW2*FMEAN2
THETM = ATAN2 (GW1*SIN(THETM1)+GW2*SIN(THETM2),GW1*COS(THETM1)+GW2*COS(THETM2))

!     3.2 ADJUST LEFT SPECTRUM TO MEAN VALUES.                                 !

CALL ROTSPEC (F1, FL, THETM-THETM1)        !! ROTATE.
CALL STRSPEC (FL, F_L, FMEAN1/FMEAN)       !! STRETCH.
GW1 = GW1*EMEAN/EMEAN1                     !! ADJUST ENERGY.

!    3.3 ADJUST RIGHT SPECTRUM TO MEAN VALUES.                                 !

CALL ROTSPEC (F2, FL, THETM-THETM2)        !! ROTATE.
CALL STRSPEC (FL, F_R, FMEAN2/FMEAN)       !! STRETCH.
GW2 = GW2*EMEAN/EMEAN2                     !! ADJUST ENERGY.

!      3.4 LINEAR INTERPOLATION TO NEW SPECTRA.                                !

FL = GW1*F_L + GW2*F_R

END SUBROUTINE INTSPEC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MEAN_DIRECTION_B (F3, THQ, SPREAD, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!  MEAN_DIRECTION_B - COMPUTATION OF MEAN WAVE DIRECTION FOR BLOCK.            !
!                                                                              !
!     S.D. HASSELMANN                                                          !
!     OPTIMIZED BY L. ZAMBRESKY                                                !
!     MODIFIED FOR K-MODEL BY C.SCHNEGGENBURGER                                !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO COMPUTE MEAN WAVE DIRECTION FROM ENERGY DENSITY AT EACH GRID POINT. !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGRATION OF SPECTRUM TIMES SIN AND COS OVER DIRECTION.              !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLE                                                       !

REAL,    INTENT(IN)            :: F3(:,:,:)   !! BLOCK OF DENSITY SPECTRA.
REAL,    INTENT(OUT), OPTIONAL :: THQ(:)      !! MEAN DIRECTION [RAD].
REAL,    INTENT(OUT), OPTIONAL :: SPREAD(:)   !! MEAN SPREAD [RAD].
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLE                                                           !

INTEGER         :: K, M
REAL            :: SI(1:SIZE(F3,1)), CI(1:SIZE(F3,1)), TEMP(1:SIZE(F3,1))
REAL*8          :: TEMP_DBL(1:SIZE(F3,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIALISE SIN AND COS ARRAYS.                                        !
!        ------------------------------                                        !

SI = 0.
CI = 0.

IF (PRESENT(SPREAD)) SPREAD = 0.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.                            !
!        ------------------------------------------                            !

IF (PRESENT(MASK)) THEN
   DO K = 1,SIZE(F3,2)
      TEMP = 0.
      DO M = 1,SIZE(F3,3)
         WHERE (MASK(:,K,M)) TEMP = TEMP + F3(:,K,M)*DFIM(M)
      END DO
      SI = SI + SINTH(K)*TEMP
      CI = CI + COSTH(K)*TEMP
      IF (PRESENT(SPREAD)) SPREAD = SPREAD + TEMP
   END DO
ELSE
   DO K = 1,SIZE(F3,2)
      TEMP_DBL = MATMUL(DBLE(F3(:,K,:)),DBLE(DFIM))
      TEMP = TEMP_DBL
      SI = SI + SINTH(K)*TEMP
      CI = CI + COSTH(K)*TEMP
      IF (PRESENT(SPREAD)) SPREAD = SPREAD + TEMP
   END DO
END IF
! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COMPUTE MEAN DIRECTION.                                               !
!        -----------------------                                               !

IF (PRESENT(THQ)) THEN
   WHERE (CI.EQ.0.) CI = 0.1E-30
   THQ = ATAN2(SI,CI)
   WHERE (THQ.LT.0.) THQ = THQ + ZPI
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     4. COMPUTE MEAN SPREAD.                                                  !
!        --------------------                                                  !

IF (PRESENT(SPREAD)) THEN
   WHERE (ABS(CI) .LT. 0.1E-15) CI = SIGN(0.1E-15,CI)
   WHERE (ABS(SI) .LT. 0.1E-15) SI = SIGN(0.1E-15,SI)
   WHERE (ABS(SPREAD) .LT. 0.1E-15) SPREAD = SIGN(0.1E-15,SPREAD)
   SPREAD = 2.*(1.-SQRT(SI**2 + CI**2)/SPREAD)
   WHERE (SPREAD.LE.0)
      SPREAD = TINY(1.)
   ELSEWHERE
      SPREAD = SQRT(SPREAD)
   END WHERE
END IF

END SUBROUTINE MEAN_DIRECTION_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE MEAN_DIRECTION_1 (F3, THQ, SPREAD)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   MEAN_DIRECTION_1 - COMPUTATION OF MEAN WAVE DIRECTION ONE SPECTRUM.        !
!                                                                              !
!     S.D. HASSELMANN                                                          !
!     OPTIMIZED BY L. ZAMBRESKY                                                !
!     MODIFIED FOR K-MODEL BY C.SCHNEGGENBURGER                                !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO COMPUTE MEAN WAVE DIRECTION FROM ONE SPECTRUM.                      !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGRATION OF SPECTRUM TIMES SIN AND COS OVER DIRECTION.              !
!                                                                              !
!     EXTERNALS.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLE                                                       !

REAL, INTENT(IN)            :: F3(:,:)  !! DENSITY SPECTRUM.
REAL, INTENT(OUT), OPTIONAL :: THQ      !! MEAN DIRECTION [RAD].
REAL, INTENT(OUT), OPTIONAL :: SPREAD   !! MEAN SPREAD [RAD].

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLE                                                           !

REAL     :: SI, CI
REAL     :: TEMP(1:SIZE(F3,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.                            !
!        ------------------------------------------                            !

TEMP = MATMUL(F3, DFIM)
SI = DOT_PRODUCT(TEMP,SINTH)
CI = DOT_PRODUCT(TEMP,COSTH)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN DIRECTION.                                               !
!        -----------------------                                               !

IF (PRESENT(THQ)) THEN
   IF (CI.EQ.0.) CI = 0.1E-30
   THQ = ATAN2(SI,CI)
   IF (THQ.LT.0.) THQ = THQ + ZPI
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. COMPUTE MEAN SPREAD.                                                  !
!        --------------------                                                  !

IF (PRESENT(SPREAD)) THEN
   SPREAD = SUM(TEMP)
   IF (ABS(CI) .LT. 0.1E-15) CI = SIGN(0.1E-15,CI)
   IF (ABS(SI) .LT. 0.1E-15) SI = SIGN(0.1E-15,SI)
   IF (ABS(SPREAD) .LT. 0.1E-15) SPREAD = SIGN(0.1E-15,SPREAD)
   SPREAD = 2.*(1.-SQRT(SI**2 + CI**2)/SPREAD)
   IF (SPREAD.LE.0) THEN
      SPREAD = TINY(1.)
   ELSE
      SPREAD = SQRT(SPREAD)
   END IF
END IF

END SUBROUTINE MEAN_DIRECTION_1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PEAK_PERIOD_B (F, PEAKP, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    PEAK_PERIOD_B - COMPUTATES PEAK PERIOD (VECTOR VERSION).                  !
!                                                                              !
!     H. GUNTHER      ECMWF            DECEMBER 1989                           !
!     (CODE REMOVED FROM SUB. FEMEAN)                                          !
!     H. GUNTHER      GKSS            FEBRUARY 2002  CHANGED TO PERIOD.        !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE PEAK PERIOD AT EACH GRID POINT.                                !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE FREQUENCY INDEX OF THE 1-D SPECTRA ARE COMPUTED AND                !
!       CONVERTED TO PERIODS.                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)            :: F(:,:,:)    !! BLOCK OF SPECTRA.
REAL,    INTENT(OUT)           :: PEAKP(:)    !! PEAK PERIODS.
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER  :: IJ
INTEGER  :: IPEAK(SIZE(F,1))
REAL     :: EED1D(SIZE(F,1),SIZE(F,3))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTE 1-D SPECTRUM (WITHOUT DELTA THETA).                           !
!        -------------------------------------------                           !

IF (PRESENT(MASK)) THEN
   EED1D = SUM(F, DIM=2, MASK=MASK)
ELSE
   EED1D = SUM(F, DIM=2)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. DEFINE PEAK INDEX.                                                    !
!        ------------------                                                    !

DO IJ = 1,SIZE(F,1)
   IPEAK(IJ:IJ) = MAXLOC(EED1D(IJ,:))
END DO

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. CALCULATE PEAK PERIOD FROM PEAK INDEX.                                !
!        --------------------------------------                                !

PEAKP = 1./FR(IPEAK)
WHERE (IPEAK.EQ.1) PEAKP = 1.

END SUBROUTINE PEAK_PERIOD_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE PEAK_PERIOD_1 (F, PEAKP, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!    PEAK_PERIOD_1 - COMPUTATES PEAK PERIOD (SCALAR VERSION).                  !
!                                                                              !
!     H. GUNTHER      ECMWF            DECEMBER 1989                           !
!     (CODE REMOVED FROM SUB. FEMEAN)                                          !
!     H. GUNTHER      GKSS            FEBRUARY 2002  CHANGED TO PERIOD.        !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE PEAK PERIOD.                                                   !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       THE FREQUENCY INDEX OF THE 1-D SPECTRUM IS COMPUTED AND                !
!       CONVERTED TO PERIOD.                                                   !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)            :: F(:,:)     !! SPECTRUM.
REAL,    INTENT(OUT)           :: PEAKP      !! PEAK PERIOD.
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:)  !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER  :: IPEAK(1:1)
REAL     :: EED1D(SIZE(F,2))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COMPUTE 1-D SPECTRUM (WITHOUT DELTA THETA).                           !
!        -------------------------------------------                           !

IF (PRESENT(MASK)) THEN
   EED1D = SUM(F, DIM=1, MASK=MASK)
ELSE
   EED1D = SUM(F, DIM=1)
END IF


! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. DEFINE PEAK INDEX.                                                    !
!        ------------------                                                    !

IPEAK(1:1) = MAXLOC(EED1D)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. CALCULATE PEAK PERIOD FROM PEAK INDEX.                                !
!        --------------------------------------                                !

PEAKP = 1./FR(IPEAK(1))
IF (IPEAK(1).EQ.1) PEAKP = 1.

END SUBROUTINE PEAK_PERIOD_1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE ROTSPEC (F_IN, F_OUT, RTHET)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   ROTSPEC - ROUTINE TO ROTATE THE SPECTRUM.                                  !
!                                                                              !
!     EVA BAUER      MPI  HAMBURG    MAY 1990.                                 !
!     H. GUNTHER          GKSS/ECMWF JAN. 1991   MODIFIED FOR CYCLE_4          !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!       TO ROTATE THE SPECTRUM.                                                !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLE.                                                      !
!     -------------------                                                      !

REAL,    INTENT(IN)  :: F_IN(:,:)      !! SPECTRUM TO BE ROTATED.
REAL,    INTENT(OUT) :: F_OUT(:,:)     !! ROTATED SPECTRUM.
REAL,    INTENT(IN) ::  RTHET          !! TURNING ANGLE [RAD], CLOCKWISE.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER         :: INC
REAL            :: ADIF, BDIF

! ---------------------------------------------------------------------------- !

ADIF = RTHET * REAL(SIZE(F_IN,1)) / ZPI
INC = -FLOOR(ADIF)
ADIF = ADIF + REAL(INC)
BDIF = 1. - ADIF

F_OUT = BDIF * CSHIFT(F_IN, INC) + ADIF * CSHIFT(F_IN, INC-1)

END SUBROUTINE ROTSPEC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE STRSPEC (F_IN, F_OUT, GAMMA)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   STRSPEC - ROUTINE TO STRETCH A SPECTRUM.                                   !
!                                                                              !
!      EVA BAUER      MPI  HAMBURG    MAY 1990.                                !
!      H. GUNTHER     GKSS/ECMWF      JAN 1991  MODIFIED FOR CYCLE_4.          !
!      H. GUNTHER     GKSS            JAN 2002  FT90.                          !
!                                               ERROR FOR SHIFT TO HIGER       !
!                                               FREQUENCIES CORRECTED.         !
!                                                                              !
!     PURPOSE.                                                                 !
!     -------                                                                  !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCES.                                                              !
!     -----------                                                              !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)  :: F_IN (:,:)    !! INPUT SPECTRUM.
REAL,    INTENT(OUT) :: F_OUT(:,:)    !! OUTPUT SPECTRUM.
REAL,    INTENT(IN)  :: GAMMA         !! STRETCHING PARAMETER.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

INTEGER         :: ML, INC
REAL            :: ADIF, BDIF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIALIZATION.                                                       !
!        ---------------                                                       !

IF (GAMMA.EQ.1.0) THEN
   F_OUT = F_IN
   RETURN
END IF

F_OUT = 0.0
ML = SIZE(F_IN,2)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. DETERMINE ACROSS HOW MANY FREQUENCY BINS THE  STRETCHING IS ACTING    !
!        AND THE INTERPOLATION WEIGHTS.                                        !
!        --------------------------------------------------------------------  !

INC = FLOOR(LOG10(GAMMA)/LOG10(CO))
IF (ABS(INC).GE.ML-1) RETURN      !! ENERGY IS SHIFTED OUT OF FREQUENCY RANGE

ADIF = (CO -GAMMA*CO**(-INC))/(CO-1.)
BDIF = 1. - ADIF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. STRECH SPECTRUM.                                                      !
!        ----------------                                                      !

IF (INC.GE.0) THEN

!     3.1 SHIFT TO LOWER FREQUENCIES.                                          !

   F_OUT(:,1:ML-INC-1) = ADIF*F_IN(:,1+INC:ML-1) + BDIF*F_IN(:,2+INC:ML)
   F_OUT(:,ML-INC)     = ADIF*F_IN(:,ML)
ELSE

!      3.2 SHIFT TO HIGHER FREQUENCIES.                                        !

   F_OUT(:,1-INC:ML) = ADIF*F_IN(:,1:ML+INC) + BDIF*F_IN(:,2:ML+INC+1)
   F_OUT(:,-INC)     = BDIF*F_IN(:,1)
END IF

END SUBROUTINE STRSPEC

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE STOKES_DRIFT (F3, UST, VST, IN)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   STOKES_DRIFT - COMPUTES STOKES DRIFT FROM SPECTRA FOR DEEP WATER.          !
!                                                                              !
!     M. REISTAD AND O SAETRA     DNMI     AUGUST 1997                         !
!     H. GUNTHER    GKSS          NOVEMBER 2010  FT90.                         !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       STOKES DRIFT FIELDS ARE CREATED FROM SPECTRUM FIELDS.                  !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!        METHOD PROPOSED BY Alastair Jenkins                                   !
!                                                                              !
!     REFERENCES.                                                              !
!      -----------                                                             !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !
    
REAL,              INTENT(IN)  :: F3(:,:,:)  !! BLOCK OF SPECTRA.
REAL,              INTENT(OUT) :: UST(:)     !! U COMPONENTS OF STOCKES DRIFT. 
REAL,              INTENT(OUT) :: VST(:)     !! V COMPONENTS OF STOCKES DRIFT.
INTEGER, OPTIONAL, INTENT(IN)  :: IN(:)      !! SHALLOW WATER TABLE INDEX

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: FAK_DEEP = 16.*PI**3/G
REAL, PARAMETER :: FAK_SHALLOW = 4.*PI
REAL            :: FAK, TAILFAC
REAL            :: FAKT(1:SIZE(F3,1))
REAL            :: SI(1:SIZE(F3,1)), CI(1:SIZE(F3,1))
INTEGER         :: IJ, K, M

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INITIAL.                                                              !
!        --------                                                              !

UST(:) = 0.
VST(:) = 0.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES AND DIRECTIONS.                            !
!        ------------------------------------------                            !

IF (PRESENT(IN)) THEN     !! SHALLOW WATER
   FREQ_LOOP1: DO M = 1,ML
      FAKT(:)  = FAK_SHALLOW * FR(M) * TFAK(IN(:),M) * DFIM(M)
      DO IJ=1,SIZE(F3,1)   
         SI(IJ) = F3(IJ,1,M)*SINTH(1)
         CI(IJ) = F3(IJ,1,M)*COSTH(1)
      END DO
      DIR_LOOP1: DO K=2,KL
         DO IJ=1,SIZE(F3,1)   
            SI(IJ) = SI(IJ) + F3(IJ,K,M)*SINTH(K)
            CI(IJ) = CI(IJ) + F3(IJ,K,M)*COSTH(K)
         END DO
      END DO DIR_LOOP1

      DO IJ=1,SIZE(F3,1)   
         SI(IJ) = FAKT(IJ) * SI(IJ) 
         CI(IJ) = FAKT(IJ) * CI(IJ)
         UST(IJ) = UST(IJ) + SI(IJ)
         VST(IJ) = VST(IJ) + CI(IJ)
      END DO
   END DO FREQ_LOOP1

ELSE                        !! DEEP WATER

   FREQ_LOOP2: DO M = 1,ML
      FAK   = FAK_DEEP * FR(M)**3 * DFIM(M)
      DO IJ=1,SIZE(F3,1)   
         SI(IJ) = F3(IJ,1,M)*SINTH(1)
         CI(IJ) = F3(IJ,1,M)*COSTH(1)
      END DO
      DIR_LOOP2: DO K=2,KL
         DO IJ=1,SIZE(F3,1)   
            SI(IJ) = SI(IJ) + F3(IJ,K,M)*SINTH(K)
            CI(IJ) = CI(IJ) + F3(IJ,K,M)*COSTH(K)
         END DO
      END DO DIR_LOOP2

      DO IJ=1,SIZE(F3,1)   
         SI(IJ) = FAK * SI(IJ) 
         CI(IJ) = FAK * CI(IJ)
         UST(IJ) = UST(IJ) + SI(IJ)
         VST(IJ) = VST(IJ) + CI(IJ)
      END DO
   END DO FREQ_LOOP2
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. ADD CONTRIBUTION FROM TAIL.                                           !
!        ---------------------------                                           !

tailfac = dfim(ml)/delth
TAILFAC = FR(ML)**2 / (tailfac *(FR(ML)+0.5*tailfac))

UST(:) = UST(:) + TAILFAC*SI(:)
VST(:) = VST(:) + TAILFAC*CI(:)

END SUBROUTINE STOKES_DRIFT

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TM1_TM2_PERIODS_B (F, EMEAN, TM1, TM2, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   TM1_TM2_PERIODS_B - COMPUTES TM1 AND/OR TM2 PERIODS (VECTOR VESION).       !
!                                                                              !
!     C.SCHNEGGENBURGER 08/97.                                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE TM1 AND TM2 PERIODS.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGARATION OF SPECTRA AND ADDING OF TAIL FACTORS.                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN )           :: F(:,:,:)    !! BLOCK OF SPECTRA.
REAL,    INTENT(IN )           :: EMEAN(:)    !! TOTAL ENERGY [M*M].
REAL,    INTENT(OUT), OPTIONAL :: TM1(:)      !! TM1 PERIOD [S].
REAL,    INTENT(OUT), OPTIONAL :: TM2(:)      !! TM2 PERIOD [S].
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL    :: TEMP(SIZE(F,1),SIZE(F,3))
REAL*8  :: TM1_DBL(1:SIZE(F,1)),TM2_DBL(1:SIZE(F,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !

IF (PRESENT(MASK)) THEN
   TEMP = SUM(F, DIM=2, MASK=MASK)
ELSE
   TEMP = SUM(F, DIM=2)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TM1 PERIOD.                                                           !
!        -----------                                                           !

IF (PRESENT(TM1)) THEN
   TM1_DBL = MATMUL(DBLE(TEMP), DBLE(DFFR))       !! FREQUENCY INTEGRATION.
   TM1 = TM1_DBL

   TM1  = TM1 + MP1_TAIL * TEMP(:,SIZE(F,3))

   WHERE (EMEAN.GT.EMIN)                          !! NORMALIZE WITH ENERGY.
      TM1 = EMEAN/TM1   
   ELSEWHERE
      TM1 = 1.
   END WHERE
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. TM2 PERIOD.                                                           !
!        -----------                                                           !

IF (PRESENT(TM2)) THEN
   TM2_DBL = MATMUL(DBLE(TEMP), DBLE(DFFR2))      !! FREQUENCY INTEGRATION.
   TM2 = TM2_DBL

   TM2    = TM2 + MP2_TAIL*TEMP(:,SIZE(F,3))      !! ADD TAIL CORRECTION.

   WHERE (EMEAN.GT.EMIN)                          !! NORMALIZE WITH ENERGY.
      TM2 = SQRT(EMEAN/TM2)   
   ELSEWHERE
      TM2 = 1.
   END WHERE
END IF

END SUBROUTINE TM1_TM2_PERIODS_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TM1_TM2_PERIODS_1 (F, EMEAN, TM1, TM2, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   TM1_TM2_PERIODS_1 - COMPUTES TM1 AND/OR TM2 PERIODS (SCALAR VESION).       !
!                                                                              !
!     C.SCHNEGGENBURGER 08/97.                                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE TM1 AND TM2 PERIODS.                                           !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGARATION OF SPECTRA AND ADDING OF TAIL FACTORS.                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN )           :: F(:,:)    !! BLOCK OF SPECTRA.
REAL,    INTENT(IN )           :: EMEAN     !! TOTAL ENERGY [M*M].
REAL,    INTENT(OUT), OPTIONAL :: TM1       !! TM1 PERIOD [S].
REAL,    INTENT(OUT), OPTIONAL :: TM2       !! TM2 PERIOD [S].
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL    :: TEMP(SIZE(F,2))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !

IF (PRESENT(MASK)) THEN
   TEMP = SUM(F, DIM=1, MASK=MASK)
ELSE
   TEMP = SUM(F, DIM=1)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. TM1 PERIOD.                                                           !
!        -----------                                                           !

IF (PRESENT(TM1)) THEN
   TM1 = DOT_PRODUCT(TEMP, DFFR)              !! FREQUENCY INTEGRATION.

   TM1    = TM1 + MP1_TAIL*TEMP(SIZE(F,2))    !! ADD TAIL CORRECTION.
   TM1    = MAX(TM1, EMIN)

   IF (EMEAN.GT.EMIN) THEN                    !! NORMALIZE WITH TOTAL ENERGY.
      TM1 = EMEAN/TM1
   ELSE
      TM1 = 1.
   END IF
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. TM2 PERIOD.                                                           !
!        -----------                                                           !

IF (PRESENT(TM2)) THEN
   TM2 = DOT_PRODUCT(TEMP, DFFR2)              !! FREQUENCY INTEGRATION.

   TM2    = TM2 + MP2_TAIL*TEMP(SIZE(F,2))     !! ADD TAIL CORRECTION.
   TM2    = MAX(TM2, EMIN)

   IF (EMEAN.GT.EMIN) THEN                     !! NORMALIZE WITH TOTAL ENERGY.
      TM2 = SQRT(EMEAN/TM2)
   ELSE
      TM2 = 1.
   END IF
END IF

END SUBROUTINE TM1_TM2_PERIODS_1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TOTAL_ENERGY_B (F3, EMEAN, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   TOTAL_ENERGY_B - COMPUTES TOTAL ENERGY (VECTOR VERSION).                   !
!                                                                              !
!     S.D. HASSELMANN.                                                         !
!     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER                               !
!     H. GUENTHER     GKSS   DECEMBER 2001  FT90                               !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.                            !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGRATION OVER DIRECTION AND FREQUENCY. A TAIL CORRECTION IS ADDED.  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)            :: F3(:,:,:)    !! BLOCK OF SPECTRA.
REAL,    INTENT(OUT)           :: EMEAN(:)     !! TOTAL ENERGY.
LOGICAL, INTENT(IN), OPTIONAL  :: MASK(:,:,:)  !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL     :: TEMP(SIZE(F3,1),SIZE(F3,3))
REAL*8   :: EMEAN_DBL(1:SIZE(F3,1))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTION (WITHOUT DELTH).                             !
!        -----------------------------------------                             !

IF (PRESENT(MASK)) THEN
   TEMP = SUM(F3, DIM=2, MASK=MASK)
ELSE
   TEMP = SUM(F3, DIM=2)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES.                                           !
!        ---------------------------                                           !

EMEAN_DBL = MATMUL(DBLE(TEMP),DBLE(DFIM))
EMEAN = EMEAN_DBL

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. ADD TAIL ENERGY.                                                      !
!        ----------------                                                      !

EMEAN = EMEAN + MO_TAIL*TEMP(:,SIZE(F3,3))
EMEAN = MAX(EMEAN, EMIN)

END SUBROUTINE TOTAL_ENERGY_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE TOTAL_ENERGY_1 (F3, EMEAN, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   TOTAL_ENERGY_1 - COMPUTES TOTAL ENERGY (SCALAR VERSION).                   !
!                                                                              !
!     S.D. HASSELMANN.                                                         !
!     OPTIMIZED BY: L. ZAMBRESKY AND H. GUENTHER                               !
!     H. GUENTHER     GKSS   DECEMBER 2001  FT90                               !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       TO COMPUTE TOTAL ENERGY AT EACH GRID POINT.                            !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGRATION OVER DIRECTION AND FREQUENCY. A TAIL CORRECTION IS ADDED.  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN)           :: F3(:,:)    !! SPECTRUM.
REAL,    INTENT(OUT)          :: EMEAN      !! TOTAL ENERGY.
LOGICAL, INTENT(IN), OPTIONAL :: MASK(:,:)  !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL     :: TEMP(SIZE(F3,2))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTION (WITHOUT DELTH).                             !
!        -----------------------------------------                             !

IF (PRESENT(MASK)) THEN
   TEMP = SUM(F3 ,DIM=1,MASK=MASK)
ELSE
   TEMP = SUM(F3, DIM=1)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. INTEGRATE OVER FREQUENCIES.                                           !
!        ---------------------------                                           !

EMEAN = DOT_PRODUCT(TEMP,DFIM)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     3. ADD TAIL ENERGY.                                                      !
!        ----------------                                                      !

EMEAN = EMEAN + MO_TAIL*TEMP(SIZE(F3,2))
EMEAN = MAX(EMEAN, EMIN)

END SUBROUTINE TOTAL_ENERGY_1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE COS2_SPR_1 (TH, THES, ST)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     COS2_SPR - ROUTINE TO COMPUTE SPREADING FACTOR (SCALAR VERSION).         !
!                                                                              !
!     SUSANNE HASSELMANN  JULY 1986.                                           !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTATION OF COS**2 SPREADING FUNCTION.                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN) :: TH(:)       !! DIRECTIONS.
REAL,    INTENT(IN) :: THES        !! MEAN WAVE DIRECTION.
REAL,    INTENT(OUT) :: ST(:)      !! SPREADING FUNCTION.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: ZDP=2./PI

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COSINE SQUARE SPREAD.                                                         !
!     ------------------------                                                         !

ST(:) = MAX(0. ,COS(TH(:)-THES))
ST = ZDP*ST**2
WHERE (ST.LT.0.1E-08) ST = 0.

END SUBROUTINE COS2_SPR_1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE COS2_SPR_B (TH, THES, ST)

! ---------------------------------------------------------------------------- !
!                                                                              !
!     COS2_SPR - ROUTINE TO COMPUTE SPREADING FACTOR (VECTOR VERSION).         !
!                                                                              !
!     SUSANNE HASSELMANN  JULY 1986.                                           !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTATION OF COS**2 SPREADING FUNCTION.                              !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN) :: TH(:)       !! DIRECTIONS.
REAL,    INTENT(IN) :: THES(:)     !! MEAN WAVE DIRECTIONS.
REAL,    INTENT(OUT) :: ST(:,:)    !! SPREADING FUNCTION.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL, PARAMETER :: ZDP=2./PI
INTEGER         :: K

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. COSINE SQUARE SPREAD.                                                         !
!     ------------------------                                                         !

DO K = 1,SIZE(TH)
   ST(:,K) = MAX(0. ,COS(TH(K)-THES(:)))
END DO
ST = ZDP*ST**2
WHERE (ST.LT.0.1E-08) ST = 0.

END SUBROUTINE COS2_SPR_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WM1_WM2_WAVENUMBER_B (F, EMEAN, WM1, WM2, IN, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   WM1_WM2_WAVENUMBER_B - COMPUTES WM1 AND/OR WM2 WAVENUMBERS (VECTOR VESION).!
!                                                                              !
!     C.SCHNEGGENBURGER 08/97.                                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE  WM1 AND/OR WM2 WAVENUMBERS                                    !
!          WM1 IS SQRT(1/K)*F INTGRATION                                       !
!          WM2 IS SQRT(K)*F INTGRATION                                         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGARATION OF SPECTRA AND ADDING OF TAIL FACTORS.                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN )           :: F(:,:,:)    !! BLOCK OF SPECTRA.
REAL,    INTENT(IN )           :: EMEAN(:)    !! TOTAL ENERGY [M*M].
REAL,    INTENT(OUT), OPTIONAL :: WM1(:)      !! WM1 WAVENUMBER [M].
REAL,    INTENT(OUT), OPTIONAL :: WM2(:)      !! WM2 WAVENUMBER [M].
INTEGER, INTENT(IN),  OPTIONAL :: IN   (:)    !! DEPTH TABLE INDEX.
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL    :: DEL2
REAL    :: TEMP(SIZE(F,1),SIZE(F,3))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !

IF (PRESENT(MASK)) THEN
   TEMP = SUM(F, DIM=2, MASK=MASK)
ELSE
   TEMP = SUM(F, DIM=2)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN WAVE NUMBER WM1.                                         !
!        -----------------------------                                         !

IF (PRESENT(WM1)) THEN
   DEL2 = SQRT(G)/ZPI
   IF (PRESENT(IN)) THEN
!AB
      WM1 = MATMUL(TEMP/(SQRT(TFAK(IN,:))+0.0000001), DFIM)   !! INTEGRATE OVER FREQUENCY.
   ELSE
      WM1 = MATMUL(TEMP, DEL2*DFIMOFR)
   END IF

   WM1 = WM1 +  MM1_TAIL*DEL2*TEMP(:,SIZE(F,3))   !! ADD TAIL.

   WHERE (EMEAN.GT.EMIN) 
      WM1 = (EMEAN/WM1)**2                        !! NORMALIZE WITH ENERGY.
   ELSEWHERE
      WM1 = 1.
   END WHERE   
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN WAVE NUMBER WM2.                                         !
!        -----------------------------                                         !

IF (PRESENT(WM2)) THEN
   DEL2 = ZPI/SQRT(G)
   IF (PRESENT(IN)) THEN
      WM2 = MATMUL(TEMP*SQRT(TFAK(IN,:)), DFIM)   !! INTEGRATE OVER FREQUENCY.
   ELSE
      WM2 = DEL2*MATMUL(TEMP, DFFR)
   END IF

   WM2 = WM2 + MP1_TAIL*DEL2*TEMP(:,SIZE(F,3))   !! ADD TAIL.
 
   WHERE (EMEAN.GT.EMIN) 
      WM2 = (WM2/EMEAN)**2                        !! NORMALIZE WITH ENERGY.
   ELSEWHERE
      WM2 = 1.
   END WHERE   
END IF

END SUBROUTINE WM1_WM2_WAVENUMBER_B

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

SUBROUTINE WM1_WM2_WAVENUMBER_1 (F, EMEAN, WM1, WM2, IN, MASK)

! ---------------------------------------------------------------------------- !
!                                                                              !
!   WM1_WM2_WAVENUMBER_1 - COMPUTES WM1 AND/OR WM2 WAVENUMBER (SCALAR VESION). !
!                                                                              !
!     C.SCHNEGGENBURGER 08/97.                                                 !
!                                                                              !
!     PURPOSE.                                                                 !
!     --------                                                                 !
!                                                                              !
!       COMPUTE  WM1 AND/OR WM2 WAVENUMBERS                                    !
!          WM1 IS SQRT(1/K)*F INTGRATION                                       !
!          WM2 IS SQRT(K)*F INTGRATION                                         !
!                                                                              !
!     METHOD.                                                                  !
!     -------                                                                  !
!                                                                              !
!       INTEGARATION OF SPECTRA AND ADDING OF TAIL FACTORS.                    !
!                                                                              !
!     REFERENCE.                                                               !
!     ----------                                                               !
!                                                                              !
!       NONE.                                                                  !
!                                                                              !
! ---------------------------------------------------------------------------- !
!                                                                              !
!     INTERFACE VARIABLES.                                                     !
!     --------------------                                                     !

REAL,    INTENT(IN )           :: F(:,:)    !! SPECTRUM.
REAL,    INTENT(IN )           :: EMEAN     !! TOTAL ENERGY [M*M].
REAL,    INTENT(OUT), OPTIONAL :: WM1       !! WM1 WAVENUMBER [M].
REAL,    INTENT(OUT), OPTIONAL :: WM2       !! WM2 WAVENUMBER [M].
INTEGER, INTENT(IN),  OPTIONAL :: IN        !! DEPTH TABLE INDEX.
LOGICAL, INTENT(IN),  OPTIONAL :: MASK(:,:) !! INTEGRATION MASK.

! ---------------------------------------------------------------------------- !
!                                                                              !
!     LOCAL VARIABLES.                                                         !
!     ----------------                                                         !

REAL    :: DEL2
REAL    :: TEMP(SIZE(F,2))

! ---------------------------------------------------------------------------- !
!                                                                              !
!     1. INTEGRATE OVER DIRECTIONS (WITHOUT DELTH).                            !
!        ------------------------------------------                            !

IF (PRESENT(MASK)) THEN
   TEMP = SUM(F, DIM=1, MASK=MASK)
ELSE
   TEMP = SUM(F, DIM=1)
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN WAVE NUMBER WM1.                                         !
!        -----------------------------                                         !

IF (PRESENT(WM1)) THEN
   DEL2 = SQRT(G)/ZPI
   IF (PRESENT(IN)) THEN
      WM1 = DOT_PRODUCT(TEMP/SQRT(TFAK(IN,:)), DFIM)  !! FREQUENCY INTEGRATION.
   ELSE
      WM1 = DOT_PRODUCT(TEMP, DEL2*DFIMOFR)
   END IF

   WM1 = WM1 + MM1_TAIL*DEL2*TEMP(SIZE(F,2))          !! ADD TAIL.
   
   IF (EMEAN.GT.EMIN) THEN
      WM1 = (EMEAN/WM1)**2                            !! NORMALIZE WITH ENERGY.
   ELSE
      WM1 = 1.
   END IF  
END IF

! ---------------------------------------------------------------------------- !
!                                                                              !
!     2. COMPUTE MEAN WAVE NUMBER WM2.                                         !
!        -----------------------------                                         !

IF (PRESENT(WM2)) THEN
   DEL2 = ZPI/SQRT(G)
   IF (PRESENT(IN)) THEN
      WM2 = DOT_PRODUCT(TEMP*SQRT(TFAK(IN,:)), DFIM)  !! FREQUENCY INTEGRATION.
   ELSE
      WM2 = DEL2*DOT_PRODUCT(TEMP, DFFR)
   END IF

   WM2 = WM2 + MP1_TAIL*DEL2*TEMP(SIZE(F,2))          !! ADD TAIL.

   IF (EMEAN.GT.EMIN) THEN
      WM2 = (WM2/EMEAN)**2                            !! NORMALIZE WITH ENERGY.
   ELSE
      WM2 = 1.
   END IF  
END IF

END SUBROUTINE WM1_WM2_WAVENUMBER_1

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

END MODULE WAM_INTERFACE_MODULE
