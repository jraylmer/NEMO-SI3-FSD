MODULE par_ice
   !!======================================================================
   !!                        ***  par_ice  ***
   !! Ice :   set the ice parameters
   !!======================================================================
   !! History :  4.x  !  2023     (Rousset)  Original code
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   USE par_kind          ! kind parameters

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !!                   shared namelist parameters
   !!----------------------------------------------------------------------
   !                                     !!** ice-generic parameters namelist (nampar) **
   INTEGER           , PUBLIC ::   jpl              !: number of ice  categories
   INTEGER           , PUBLIC ::   nlay_i           !: number of ice  layers
   INTEGER           , PUBLIC ::   nlay_s           !: number of snow layers
   LOGICAL           , PUBLIC ::   ln_virtual_itd   !: virtual ITD mono-category parameterization (T) or not (F)
   LOGICAL           , PUBLIC ::   ln_icedyn        !: flag for ice dynamics (T) or not (F)
   LOGICAL           , PUBLIC ::   ln_icethd        !: flag for ice thermo   (T) or not (F)
   REAL(wp)          , PUBLIC ::   rn_amax_n        !: maximum ice concentration Northern hemisphere
   REAL(wp)          , PUBLIC ::   rn_amax_s        !: maximum ice concentration Southern hemisphere
   CHARACTER(len=256), PUBLIC ::   cn_icerst_in     !: suffix of ice restart name (input)
   CHARACTER(len=256), PUBLIC ::   cn_icerst_out    !: suffix of ice restart name (output)
   CHARACTER(len=256), PUBLIC ::   cn_icerst_indir  !: ice restart input directory
   CHARACTER(len=256), PUBLIC ::   cn_icerst_outdir !: ice restart output directory

   !                                     !! ** namelist (namini) **
   LOGICAL, PUBLIC  ::   ln_iceini        !: Ice initialization or not
   INTEGER, PUBLIC  ::   nn_iceini_file   !: Ice initialization:
   !                                      !        0 = Initialise sea ice based on SSTs
   !                                      !        1 = Initialise sea ice from single category netcdf file
   !                                      !        2 = Initialise sea ice from multi category restart file

   !                                     !!** ice-itd namelist (namitd) **
   REAL(wp), PUBLIC ::   rn_himin         !: minimum ice thickness
   !
   !                                     !!** ice-surface boundary conditions namelist (namsbc) **
   REAL(wp), PUBLIC ::   rn_snwblow       !: coef. for partitioning of snowfall between leads and sea ice
                                          ! -- icethd_zdf and icealb -- !
   INTEGER , PUBLIC ::   nn_snwfra        !: calculate the fraction of ice covered by snow
   !                                      !   = 0  fraction = 1 (if snow) or 0 (if no snow)
   !                                      !   = 1  fraction = 1-exp(-0.2*rhos*hsnw) [MetO formulation]
   !                                      !   = 2  fraction = hsnw / (hsnw+0.02)    [CICE formulation]
                                          ! -- icesbc -- !
   REAL(wp), PUBLIC ::   rn_Cd_io         !: drag coefficient for oceanic stress
                                          ! -- icethd_zdf -- !
   LOGICAL , PUBLIC ::   ln_cndflx        !: use conduction flux as surface boundary condition (instead of qsr and qns)
   LOGICAL , PUBLIC ::   ln_cndemulate    !: emulate conduction flux (if not provided)
   INTEGER , PUBLIC ::   nn_qtrice        !: Solar flux transmitted thru the surface scattering layer:
   !                                      !   = 0  Grenfell and Maykut 1977 (depends on cloudiness and is 0 when there is snow)
   !                                      !   = 1  Lebrun 2019 (equals 0.3 anytime with different melting/dry snw conductivities)

   !                                     !!** namelist (namthd) **
   LOGICAL , PUBLIC ::   ln_icedH         ! activate ice thickness change from growing/melting (T) or not (F)
   LOGICAL , PUBLIC ::   ln_icedA         ! activate lateral melting param. (T) or not (F)
   LOGICAL , PUBLIC ::   ln_icedO         ! activate ice growth in open-water (T) or not (F)
   LOGICAL , PUBLIC ::   ln_leadhfx       ! heat in the leads is used to melt sea-ice before warming the ocean
   !
   !
   !                                     !!** ice-vertical diffusion namelist (namthd_zdf) **
   !                                      ! Conduction flux as surface forcing or not
   INTEGER , PUBLIC, PARAMETER ::   np_cnd_OFF = 0  !: no forcing from conduction flux (ice thermodynamics forced via qsr and qns)
   INTEGER , PUBLIC, PARAMETER ::   np_cnd_ON  = 1  !: forcing from conduction flux (SM0L) (compute qcn and qsr_tr via sbcblk.F90 or sbccpl.F90)
   INTEGER , PUBLIC, PARAMETER ::   np_cnd_EMU = 2  !: emulate conduction flux via icethd_zdf.F90 (BL99) (1st round compute qcn and qsr_tr, 2nd round use it)
   LOGICAL , PUBLIC ::   ln_cndi_U64      !: thermal conductivity: Untersteiner (1964)
   LOGICAL , PUBLIC ::   ln_cndi_P07      !: thermal conductivity: Pringle et al (2007)
   REAL(wp), PUBLIC ::   rn_kappa_i       !: coef. for the extinction of radiation in sea ice, Grenfell et al. (2006) [1/m]
   REAL(wp), PUBLIC ::   rn_kappa_s       !: coef. for the extinction of radiation in snw (nn_qtrice=0) [1/m]
   REAL(wp), PUBLIC ::   rn_kappa_smlt    !: coef. for the extinction of radiation in melt snw (nn_qtrice=1) [1/m]
   REAL(wp), PUBLIC ::   rn_kappa_sdry    !: coef. for the extinction of radiation in dry  snw (nn_qtrice=1) [1/m]
   LOGICAL , PUBLIC ::   ln_zdf_chkcvg    !: check convergence of heat diffusion scheme

   !                                     !!** ice-salinity namelist (namthd_sal) **
   INTEGER , PUBLIC ::   nn_icesal        !: salinity configuration used in the model
   !                                      ! 1 - constant salinity in both space and time
   !                                      ! 2 - prognostic salinity (s(z,t))
   !                                      ! 3 - salinity profile, constant in time
   !                                      ! 4 - flushing and gravity drainage
   REAL(wp), PUBLIC ::   rn_icesal        !: bulk salinity (ppt) in case of constant salinity
   REAL(wp), PUBLIC ::   rn_sinew         !: fraction of sss that is kept in new ice
   REAL(wp), PUBLIC ::   rn_simin         !: minimum ice salinity [PSU]
   LOGICAL , PUBLIC ::   ln_sal_chk       !: sanity checks for salt drainage and flushing
   INTEGER , PUBLIC ::   nn_liquidus      !: formulation of liquidus
   !                                        1 = linear liquidus
   !                                        2 = Vancopenolle et al (2019) formulation
   !                                        3 = Weast formulation (used in RJW2014)

   !                                     !!** ice-ponds namelist (namthd_pnd)
   LOGICAL , PUBLIC ::   ln_pnd           !: Melt ponds (T) or not (F)
   LOGICAL , PUBLIC ::   ln_pnd_TOPO      !: Topographic Melt ponds scheme (Flocco et al 2007, 2010)
   LOGICAL , PUBLIC ::   ln_pnd_LEV       !: Simple melt pond scheme
   LOGICAL , PUBLIC ::   ln_pnd_CST       !: Melt ponds scheme with constant fraction and depth
   LOGICAL,  PUBLIC ::   ln_pnd_lids      !: Allow ponds to have frozen lids
   LOGICAL,  PUBLIC ::   ln_pnd_rain      !: rain added to melt ponds
   LOGICAL , PUBLIC ::   ln_pnd_alb       !: melt ponds affect albedo
   REAL(wp), PUBLIC ::   rn_pnd_hl_min    !: pond lid thickness below which full pond area used in albedo calculation
   REAL(wp), PUBLIC ::   rn_pnd_hl_max    !: pond lid thickness above which ponds disappear from albedo calculation
   !
   !                                     !!** ice-diagnostics namelist (namdia) **
   LOGICAL , PUBLIC ::   ln_icediachk     !: flag for ice diag (T) or not (F)
   REAL(wp), PUBLIC ::   rn_icechk_cel    !: rate of ice spuriously gained/lost (at any gridcell)
   REAL(wp), PUBLIC ::   rn_icechk_glo    !: rate of ice spuriously gained/lost (globally)
   LOGICAL , PUBLIC ::   ln_icectl        !: flag for sea-ice points output (T) or not (F)
   INTEGER , PUBLIC ::   iiceprt          !: debug i-point
   INTEGER , PUBLIC ::   jiceprt          !: debug j-point
   REAL(wp), PUBLIC ::   rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft !: conservation diagnostics
   
   !                                     !!** ice-floe size distribution namelist (namfsd)
   LOGICAL , PUBLIC ::   ln_fsd           !: Use FSD (T) or not (F)
   INTEGER , PUBLIC ::   nn_nfsd          !: Number of floe-size categories
   REAL(wp), PUBLIC ::   rn_floeshape     !: Floe shape parameter (dimensionless)
   LOGICAL , PUBLIC ::   ln_fsd_brit      !: Activate FSD brittle fractue scheme or not

   !                                     !!** wave-ice interactions namelist (namwav)
   LOGICAL , PUBLIC ::   ln_ice_wav       !: Activate wave ice interactions (T) or not (F)
   LOGICAL , PUBLIC ::   ln_ice_wav_attn  !: Activate wave-in-ice attenuation scheme (Roach et al. 2018)
   LOGICAL , PUBLIC ::   ln_ice_wav_spec  !: Read full wave spectrum from model/forcing

   !!----------------------------------------------------------------------
   !!                   shared other parameters
   !!----------------------------------------------------------------------
   INTEGER , PUBLIC ::   kt_ice           !: iteration number
   REAL(wp), PUBLIC ::   rDt_ice          !: ice time step
   REAL(wp), PUBLIC ::   r1_Dt_ice        !: = 1. / rDt_ice
   REAL(wp), PUBLIC ::   r1_nlay_i        !: 1 / nlay_i
   REAL(wp), PUBLIC ::   r1_nlay_s        !: 1 / nlay_s
   REAL(wp), PUBLIC, PARAMETER ::   ppepsi02 = 1.e-02_wp   !: Small number
   REAL(wp), PUBLIC, PARAMETER ::   epsi06 = 1.e-06_wp  !: small number
   REAL(wp), PUBLIC, PARAMETER ::   epsi10 = 1.e-10_wp  !: small number
   REAL(wp), PUBLIC, PARAMETER ::   epsi20 = 1.e-20_wp  !: small number

   
#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE par_ice
