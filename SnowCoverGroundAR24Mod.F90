module SnowCoverGroundAR24Mod

!!! Compute ground snow cover fraction based on Niu and Yang (2007, JGR) scheme

  use Machine
  use NoahmpVarType
  use ConstantDefineMod

  implicit none

contains

  subroutine SnowCoverGroundAR24(noahmp)

! ------------------------ Code history -----------------------------------
! Original Noah-MP subroutine: None (embedded in ENERGY subroutine)
! Original code: Guo-Yue Niu and Noah-MP team (Niu et al. 2011)
! Refactered code: C. He, P. Valayamkunnath, & refactor team (He et al. 2023)
! Code updated: Abolafia-Rosenzweig and He and NCAR's Noah-MP team (Abolafia-Rosenzweig et al., 2024)
! -------------------------------------------------------------------------

    implicit none

    type(noahmp_type), intent(inout) :: noahmp

! local variable
    real(kind=kind_noahmp)           :: SnowDensBulk   ! bulk density of snow [Kg/m3]
    real(kind=kind_noahmp)           :: MeltFac        ! melting factor for snow cover frac
    real(kind=kind_noahmp)           :: SnowMeltFac        ! snowmelt m parameter (scale-dependent)
    real(kind=kind_noahmp)           :: SnowCoverFac        ! snow cover factor [scfac] (scale-dependent)
! --------------------------------------------------------------------
    associate(                                                     &
              SnowDepth      => noahmp%water%state%SnowDepth      ,& ! in,  snow depth [m]
              SnowWaterEquiv => noahmp%water%state%SnowWaterEquiv ,& ! in,  snow water equivalent [mm]
              GridSize       => noahmp%config%domain%GridSize     ,& ! in,    noahmp model grid spacing [m]
              SnowCoverFrac  => noahmp%water%state%SnowCoverFrac   & ! out, snow cover fraction
             )
! ----------------------------------------------------------------------

    SnowCoverFrac = 0.0
    if ( SnowDepth > 0.0 ) then
         !calculate SCF parameters as a function of grid size:
         SnowMeltFac = 0.9713 + tanh(0.7436*(GridSize/1000));
         SnowCoverFac = 0.0062*sinh(0.0555*(GridSize/1000))+ 0.0555;
         !using scale-dependent parameters, employ the Niu-Yang 07 SCF soluiton
         SnowDensBulk  = SnowWaterEquiv / SnowDepth
         MeltFac       = (SnowDensBulk / 100.0)**SnowMeltFac
         SnowCoverFrac = tanh( SnowDepth /(SnowCoverFac * MeltFac)) ! C.He: bring hard-coded 2.5*z0 to MPTABLE
    endif

    end associate

  end subroutine SnowCoverGroundAR24

end module SnowCoverGroundAR24Mod
