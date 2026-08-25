MODULE icewav
   !!======================================================================
   !!                       ***  MODULE icewav ***
   !!   sea-ice : ocean wave-ice interactions
   !!======================================================================
   !! History :  5.0  !  2025     (J.R. Aylmer)         Original code based
   !!                                                   on CPOM-CICE and
   !!                                                   CICE/Icepack
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3' :                                     SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_wav_newice  : calculate floe size of new ice from local wave properties
   !!   ice_wav_attn    : calculate attenuated wave spectrum under ice
   !!   ice_wav_calc    : calculate wave properties under ice
   !!   ice_wav_frac    : wave breakup of ice, main routine called by ice_stp
   !!   wav_frac_z16    : wave fracture scheme (Zhang et al. 2016)
   !!   wav_frac_y24a   : wave fracture scheme (Yang et al. 2024 method A)
   !!   wav_frac_y24b   : wave fracture scheme (Yang et al. 2024 method B)
   !!   wav_frac_ht15   : wave fracture scheme (Horvat and Tziperman 2015)
   !!   wav_spec_bret() : Bretschneider wave spectrum formula
   !!   wav_spec_rayl() : Rayleigh spectrum expressed as function of frequency
   !!   ice_wav_init    : initialisation of wave-ice interaction module
   !!----------------------------------------------------------------------

   USE par_oce           ! ocean parameters
   USE dom_oce           ! ocean space and time domain

   USE par_ice           ! SI3 parameters
   USE phycst , ONLY :   rpi, grav, rhoi
   USE sbc_oce, ONLY :   wndm, ln_wave, ln_wave_spec, nn_nwfreq            ! SBC module
   USE sbcwave           ! SBC wave variables
   USE ice               ! sea-ice: variables
   USE icefsd , ONLY :   a_ifsd, nf_newice, floe_sl, floe_sc, floe_su, floe_ds     ! floe size distribution parameters/variables
   USE icefsd , ONLY :   ice_fsd_tstep, ice_fsd_cor, ice_fsd_dia, floe_size_dist   ! floe size distribution functions/routines

   USE in_out_manager    ! I/O manager (needed for lwm and lwp logicals)
   USE iom               ! I/O manager library (needed for iom_put)
   USE lib_mpp           ! MPP library (needed for read_nml_substitute.h90 and mpi_comm_oce)
   USE lbclnk            ! lateral boundary conditions (or mpp links)
   USE timing            ! Timing

   IMPLICIT NONE
   PRIVATE

   INTERFACE wav_merge_glo
      !!----------------------------------------------------------------------
      !!                  ***  INTERFACE wav_merge_glo  ***
      !!----------------------------------------------------------------------
      !! ** Purpose :   Merge global-domain arrays from all processors for wave attenuation scheme
      !!
      !! ** Method  :   Call mpp_sum for general input array where each processor has computed the
      !!                corresponding sub-domain portion of its copy of the global domain array,
      !!                with all other values being zero. Since all the latter values are calculated
      !!                on the other processors, summing the arrays from all processors and
      !!                distributing the result back to all processors (i.e., MPI_ALLREDUCE with
      !!                MPI_SUM operation, which is eventually called through this interface) gives
      !!                each processor the correct, completely filled global domain.
      !!
      !!                NOTE: therefore, all sub-domain calculations (in each processor) of global
      !!                arrays must NOT include halo cells in loop otherwise they are duplicated!
      !!----------------------------------------------------------------------
      MODULE PROCEDURE wav_merge_glo_2d
      MODULE PROCEDURE wav_merge_glo_3d
   END INTERFACE

   PUBLIC ::   ice_wav_newice   ! routine called by ice_thd_do
   PUBLIC ::   ice_wav_attn     ! routine called by ice_stp
   PUBLIC ::   ice_wav_frac     ! routine called by ice_stp
   PUBLIC ::   ice_wav_init     ! routine called by ice_init

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)   ::   x1d         ! 1D subdomain for wave fracture in HT15 scheme
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   Bfrac_uni   ! uniform fracture redistributor (Z16 and Y24a schemes)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wgtQ_y24b   ! weight factor (Y24b scheme; Q term)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   wgtB_y24b   ! weight factor (Y24b scheme; B term)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   stmer       ! meridional distance across T cells (m)

   INTEGER , PARAMETER ::   jpfrac_z16  = 1         ! option for nn_frac_scheme -> Zhang et al. (2016) scheme
   INTEGER , PARAMETER ::   jpfrac_y24a = 2         ! option for nn_frac_scheme -> Yang et al. (2024) scheme A
   INTEGER , PARAMETER ::   jpfrac_y24b = 3         ! option for nn_frac_scheme -> Yang et al. (2024) scheme B
   INTEGER , PARAMETER ::   jpfrac_ht15 = 4         ! option for nn_frac_scheme -> Horvat and Tziperman (2015) scheme

   LOGICAL ::   l_attn_calc_spec   ! whether spectrum needs to be calculated in subroutine ice_wav_attn
   LOGICAL ::   l_frac_calc_spec   ! whether spectrum needs to be calculated in subroutine ice_wav_frac

   ! General thresholds on wave parameters to compute wave spectra and/or trigger breakup events, where needed:
   ! (routines can use individual values but helpful to have some global 'defaults' for consistency):
   REAL(wp), PARAMETER ::   minhsw      = .01_wp    ! threshold minimum sig. wave height (m)
   REAL(wp), PARAMETER ::   minwpf      = .001_wp   ! threshold minimum peak frequency (Hz)
   REAL(wp), PARAMETER ::   minwmp      = .01_wp    ! threshold minimum wave mean period (s)
   REAL(wp), PARAMETER ::   minwspec    = epsi06    ! threshold minimum max. of wave spectrum (m2.Hz-1)

   ! Global-domain arrays needed for attenuation (ice_wav_attn) -- which also must be 'global' in module scope:
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   glamt_glo   ! T-grid longitude (degrees east)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   gphit_glo   ! T-grid latitude (degrees north)

   !                                     !!* namelist (namwav) *
   REAL(wp)        ::   rn_attn_lam_tol   !: Attenuation scheme longitude tolerance for identifying meridians (degrees E)
   !                                      !  Terms in attenuation coefficient, quadratic approx. for ln[a(T,h)] from HT15:
   REAL(wp)        ::   rn_attn_c0        !:    constant term
   REAL(wp)        ::   rn_attn_ch        !:    factor of ice thickness, h
   REAL(wp)        ::   rn_attn_ct        !:    factor of wave period, T
   REAL(wp)        ::   rn_attn_ch2       !:    factor of h^2
   REAL(wp)        ::   rn_attn_ct2       !:    factor of T^2
   REAL(wp)        ::   rn_attn_cht       !:    factor of h*T
   REAL(wp)        ::   rn_attn_tun       !: Overall tuning factor on a(h,T) [NOT logarithm of a]
   INTEGER         ::   nn_frac_scheme    !: Selection of wave-ice fracture scheme
   REAL(wp)        ::   rn_ice_wav_ecri   !: Critical strain at which ice fractures due to waves (dimensionless)
   LOGICAL         ::   ln_z16_const      !: Use constant participation factor (--> recover Zhang et al. 20*15* scheme)
   REAL(wp)        ::   rn_z16_cb         !: Constant participation factor if ln_z16_const
   REAL(wp)        ::   rn_z16_k          !: Z16 scheme dimensionless parameter 'k'
   REAL(wp)        ::   rn_z16_a          !: Z16 scheme dimensionless parameter 'a'
   REAL(wp)        ::   rn_z16_b          !: Z16 scheme dimensionless parameter 'b'
   REAL(wp)        ::   rn_z16_hc         !: Z16 scheme cutoff ice thickness (m)
   INTEGER         ::   nn_z16_nf         !: Z16 scheme domain size (number of surrounding grid cells) for fetch parameter
   REAL(wp)        ::   rn_y24a_cw        !: Parameter in Y24A fracture scheme
   REAL(wp)        ::   rn_y24a_alpha     !: Parameter in Y24A fracture scheme
   INTEGER         ::   nn_ht15_nx1d      !: Size of 1D subdomain for wave fracture calculation (HT15 only)
   REAL(wp)        ::   rn_ht15_dx1d      !: Increment of 1D subdomain for wave fracture calculation (m; HT15 only)
   INTEGER         ::   nn_ht15_rmin      !: Radius of smallest floes affected by wave fracture in units of rn_ht15_dx1d
   LOGICAL         ::   ln_ht15_rand      !: Use random phases for SSH in HT15 wave fracture calculation
   REAL(wp)        ::   rn_ht15_tol       !: Max. fractional difference between iterations for HT15 to 'converge' (ln_ht15_rand=T)
   INTEGER         ::   nn_ht15_maxiter   !: Max. iterations to convergence for HT15 when ln_ht15_rand=T
   LOGICAL         ::   ln_ht15_warn      !: Write warnings when HT15 scheme does not converge
   !
   ! Note: other namwav parameters declared in par_ice as they are needed by module
   ! ----- sbcwave which cannot access this module (would create circular dependency)

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"

CONTAINS

   SUBROUTINE wav_merge_glo_2d( paglo )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE wav_merge_glo_2d ***
      !!----------------------------------------------------------------------
      !! ** Purpose :   Merge 2D global domain array from all processors.
      !! ** Method  :   Call mpp_sum for general input 2D array.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(inout)   ::  paglo   ! input global domain array
      REAL(wp), DIMENSION(:)  , ALLOCATABLE     ::  zwork   ! flattened array for mpp
      INTEGER                                   ::  isize   ! number of values in paglo
      INTEGER                                   ::  ierr    ! allocate status return value
      !!----------------------------------------------------------------------
#if ! defined key_mpi_off

      isize = SIZE(paglo)

      ALLOCATE( zwork(isize), STAT=ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'wav_merge_glo_2d: unable to allocate array' )

      zwork(:) = RESHAPE( paglo(:,:), (/ isize /) )

      CALL mpp_sum( 'icewav', zwork, mpi_comm_oce )

      paglo(:,:) = RESHAPE( zwork, SHAPE(paglo) )

      IF( ALLOCATED(zwork) )   DEALLOCATE( zwork )   ! necessary?

#endif
   END SUBROUTINE wav_merge_glo_2d


   SUBROUTINE wav_merge_glo_3d( paglo )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE wav_merge_glo_3d ***
      !!----------------------------------------------------------------------
      !! ** Purpose :   Merge 3D global domain array from all processors.
      !! ** Method  :   Call mpp_sum for general input 3D array.
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout)   ::  paglo   ! input global domain array
      REAL(wp), DIMENSION(:)    , ALLOCATABLE     ::  zwork   ! flattened array for mpp
      INTEGER                                     ::  isize   ! number of values in paglo
      INTEGER                                     ::  ierr    ! allocate status return value
      !!----------------------------------------------------------------------
#if ! defined key_mpi_off

      isize = SIZE(paglo)

      ALLOCATE( zwork(isize), STAT=ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'wav_merge_glo_3d: unable to allocate array' )

      zwork(:) = RESHAPE( paglo(:,:,:), (/ isize /) )

      CALL mpp_sum( 'icewav', zwork, mpi_comm_oce )

      paglo(:,:,:) = RESHAPE( zwork, SHAPE(paglo) )

      IF( ALLOCATED(zwork) )   DEALLOCATE( zwork )   ! necessary?

#endif
   END SUBROUTINE wav_merge_glo_3d


   SUBROUTINE ice_wav_newice( phsw, pwpf, kcat )
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE ice_wav_newice ***
      !!
      !! ** Purpose :   Calculate floe size (category) of new ice grown in
      !!                open water dependent on local wave conditions
      !!
      !! ** Method  :   Maximum floe diameter is determined from:
      !!
      !!                   s_max = sqrt[ (2*C2*Lp^2) / (Wa*g*rhoi*pi^3) ]
      !!
      !!                where:
      !!
      !!                   C2   : tensile stress mode parameter (kg.m-1.s-2)
      !!                   Lp   : wavelength of peak wave energy density (m)
      !!                   Wa   : wave amplitude (m)
      !!                   g    : acceleration due to gravity (m.s-2)
      !!                   rhoi : density of sea ice (kg.m-3)
      !!
      !!                Following Roach et al. (2019), Wa is taken to be half
      !!                of the local significant wave height Hs, Lp is computed
      !!                from the local peak frequency fp using the dispersion
      !!                relation for deep water surface gravity waves:
      !!
      !!                   Lp = g / (2*pi*fp^2),
      !!
      !!                and C2 is a hard-coded constant, hence:
      !!
      !!                   s_max = [1/(pi^2*fp^2)] * sqrt[ (C2*g) / (pi*rhoi*Hs) ]
      !!
      !! ** Inputs  :   phsw : local significant wave height (m)
      !!                pwpf : local wave frequency of peak energy (Hz)
      !!
      !! ** Outputs :   kcat : integer index of floe size distribution category
      !!                       to which new ice formation is added
      !!
      !! ** Callers :   ice_thd_do -> [ice_wav_newice]
      !!
      !! ** Note    :   if fracture scheme is Z16, which does not use wave data, the calculation
      !!                is bypassed and instead returns the default new floe size category
      !!
      !! ** References
      !!    ----------
      !!    Roach, L. A., Bitz, C. M., Horvat, C., & Dean, S. M. (2019).
      !!              Advances in modeling interactions between sea ice and ocean surface waves.
      !!              Journal of Advances in Modeling Earth Systems, 11, 4167-4181.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), INTENT(in)    ::   phsw   ! local significant wave height (m)
      REAL(wp), INTENT(in)    ::   pwpf   ! local wave peak frequency (Hz)
      INTEGER , INTENT(inout) ::   kcat   ! FSD category index for new ice growth
      !
      REAL(wp), PARAMETER ::   zC2 = .167_wp   ! tensile mode stress parameter (kg.m-1.s-2)
      REAL(wp)            ::   zsmax           ! maximum new ice floe size (m)
      INTEGER             ::   jf              ! dummy loop index
      !
      !!-------------------------------------------------------------------

      ! Do not need as strict condition on hsw or wpf as in global module variables minhsw, minwpf
      ! (just need to make sure we avoid division by zero):
      IF( (nn_frac_scheme == jpfrac_z16) .OR. (phsw < epsi06) .OR. (pwpf < epsi06) ) THEN
         !
         kcat = nf_newice   ! no waves present => set to default new floe size category
         !
      ELSE
         !
         ! --- Calculate maximum floe size from wave conditions, r_max:
         zsmax = SQRT( zC2 * grav / (rpi * rhoi * phsw) ) / (rpi**2 * pwpf**2)
         !
         ! --- Find FSD category that r_max belongs to:
         !
         IF( zsmax < floe_sl(1) ) THEN
            ! Smaller than smallest 'resolved' floe size => put in smallest cat. anyway:
            kcat = 1
            !
         ELSE
            ! Find FSD category that s_max belongs to. Note that if s_max exceeds upper
            ! limit of largest floe size category, it goes into that category anyway:
            !
            DO jf = nn_nfsd, 1, -1
               IF( zsmax > floe_sl(jf) ) THEN
                  kcat = jf
                  EXIT   ! found kcat => stop iterating
               ENDIF
            ENDDO
            !
         ENDIF
         !
      ENDIF

   END SUBROUTINE ice_wav_newice


   SUBROUTINE ice_wav_attn( kt )
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE ice_wav_attn ***
      !!-------------------------------------------------------------------
      !!
      !! ** Purpose :   Calculate attenuated wave spectrum and derived wave properties in sea ice
      !!                grid cells using wave data from non-sea ice covered grid cells.
      !!
      !! ** Method  :   For each ice-covered grid cell -- the "target" grid cell -- locate
      !!                the nearest equatorward open-ocean grid cell along meridians that also
      !!                contains wave data -- the "source" grid cell. Waves are assumed to
      !!                propagate from the source to target grid cell being attenuated along the
      !!                way according to the ice properties in each grid cell encountered. This
      !!                'waves-in-ice' attenuation scheme follows the method of Roach et al. (2018).
      !!
      !!                Attenuation is a function of wave frequency, f, mean ice thickness, h,
      !!                and the number of floes, N. Specifically, when waves with energy spectrum
      !!                E0(f) travel to a neighbouring grid cell, the energy spectrum E1(f) at the
      !!                neighbouring grid cell is then given by:
      !!
      !!                   E1(f) = E0(f) * exp[ -N1 * a(h1,f) ]
      !!
      !!                where a(h,f) is an attenuation coefficient set to a quadratic fit to a
      !!                (more complex) theoretical formula, given by Horvat and Tziperman (2015):
      !!
      !!                   a(h,f) = c0 + ch*h + cT*T + ch2*h^2 + cT2*T^2 + chT*h*T
      !!
      !!                with T = 1/f is the wave component period and the various c are constants.
      !!                As waves travel encountering multiple grid cells k with different ice
      !!                properties from the source to the target grid cell, the resulting spectrum
      !!                is calculated as:
      !!
      !!                   E_target(f) = E_source(f) * exp[ -SUM_k{ Nk * a(hk,f) } ]
      !!
      !!                This routine proceeds as follows:
      !!
      !!                1.   Calculate the attenuation exponents, Nk * a(hk,f), for all ice covered grid
      !!                     cells, and (if necessary, according to module flag l_attn_calc_spec) the
      !!                     wave spectrum from significant wave height/peak frequency inputs at open
      !!                     ocean (source) grid cells.
      !!
      !!                2.   Copy these terms into global-domain, global (i.e., module scope) arrays. This
      !!                     is necessary because for any given target grid cell the source grid cell may
      !!                     be in a different computational subdomain, so all values must be available to
      !!                     all subdomains before step 3. A call to lib_mpp routine mppsync (a wrapper for
      !!                     MPP_BARRIER) ensures this is possible. Similar global arrays are required for
      !!                     the grid cell longitude and latitudes, but these are computed in advance, in
      !!                     subroutine ice_wav_init.
      !!
      !!                3.   For each target grid cell (in the computational subdomain), search for the
      !!                     source grid cell according to the criteria above, which is applied using
      !!                     masks on the global coordinate arrays computed/described in step 2.
      !!
      !!                     If a source is found, calculate the attenuated wave spectrum at the target
      !!                     grid cell according to the above equation.
      !!
      !!                4.   Calculate wave properties at the target cell (significant wave height, etc.)
      !!                     for the actual arrays defining such quantities declared in sbc_wave.F90.
      !!
      !! ** Callers :   ice_stp --> [ice_wav_attn]
      !! ** Calls   :               [ice_wav_attn] --> mppsync
      !!                                           --> lbc_lnk
      !!                                           --> wav_spec_bret
      !!                                           --> ice_wav_calc
      !!
      !! ** Notes   :   This routine is only called when ln_ice_wav_attn=T. It either uses the wave spectrum
      !!                in the nearest open ocean if available (ln_wave_spec=T, ln_ice_wav_spec=T), otherwise
      !!                it estimates the spectrum using the Bretschneider formula on the significant wave height
      !!                and peak frequency wave field inputs at the same location.
      !!
      !! ** References
      !!    ----------
      !!    Horvat, C., & Tziperman, E. (2015).
      !!              A prognostic model of the sea-ice floe size and thickness distribution.
      !!              The Cryosphere, 9, 2119-2134.
      !!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018).
      !!              An emergent sea ice floe size disribution in a global coupled ocean-sea ice model.
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!
      !!-------------------------------------------------------------------
      !
      INTEGER , INTENT(in)               ::   kt            ! ocean time step
      !
      REAL(wp), DIMENSION(nn_nwfreq)     ::   zattxp        ! local damping exponent (number of floes x attenuation coefficient)
      REAL(wp)                           ::   zdmean        ! grid cell mean floe diameter (m)
      REAL(wp)                           ::   znfloes       ! number of floes
      REAL(wp)                           ::   zhi           ! mean ice thickness (m)
      REAL(wp)                           ::   zloga         ! natural logarithm of attenuation coefficient
      REAL(wp)                           ::   zattxp_tot    ! total damping exponent along trajectory
      !
      REAL(wp), DIMENSION(jpiglo,jpjglo,nn_nwfreq) ::   zwspec_glo   ! GLOBAL DOMAIN ARRAY: wave spectrum
      REAL(wp), DIMENSION(jpiglo,jpjglo,nn_nwfreq) ::   zattxp_glo   ! GLOBAL DOMAIN ARRAY: local damping exponents (all grid cells)
      !
      LOGICAL , DIMENSION(jpiglo,jpjglo) ::   ll_mask_ice   ! GLOBAL DOMAIN ARRAY: ice-present  mask
      LOGICAL , DIMENSION(jpiglo,jpjglo) ::   ll_mask_mer   ! GLOBAL DOMAIN ARRAY: meridional   mask
      LOGICAL , DIMENSION(jpiglo,jpjglo) ::   ll_mask_eqt   ! GLOBAL DOMAIN ARRAY: equatorward  mask
      LOGICAL , DIMENSION(jpiglo,jpjglo) ::   ll_mask_pol   ! GLOBAL DOMAIN ARRAY: poleward     mask
      LOGICAL , DIMENSION(jpiglo,jpjglo) ::   ll_mask_tot   ! GLOBAL DOMAIN ARRAY: total        mask
      !
      INTEGER                            ::   ierr                                 ! allocate status return value
      INTEGER                            ::   ji, jj, jf, jl, jw, ji_glo, jj_glo   ! dummy loop indices
      INTEGER , DIMENSION(2)             ::   isource                              ! indices of source grid cells (in global domain)
      !
      REAL(wp), PARAMETER ::   zf_noice = -1._wp       ! dummy flag value < 0        for ice-free ocean grid cell
      REAL(wp), PARAMETER ::   zf_land  = -2._wp       ! dummy flag value < zf_noice for land grid cell
      !
      !!-------------------------------------------------------------------

      ! Control:
      IF( ln_timing )   CALL timing_start('ice_wav_attn')

      zwspec_glo(:,:,:) = 0._wp
      zattxp_glo(:,:,:) = 0._wp

      ! =================================================== !
      ! 1  Calculate attenuation exponents -- all ice cells !
      ! =================================================== !
      DO_2D(0, 0, 0, 0)
         !
         ! Use a less-strict threshold than wave fracture routine here as
         ! wave-dependent growth needs the derived wave fields (hsw, wpf)
         !
         IF( at_i(ji,jj) > epsi06 ) THEN

            ! Mean floe diameter:
            zdmean = 0._wp
            DO jf = 1, nn_nfsd
               DO jl = 1, jpl
                  zdmean = zdmean + floe_sc(jf) * a_i(ji,jj,jl) * a_ifsd(ji,jj,jf,jl)
               ENDDO
            ENDDO

            ! Number of floes per unit distance encountered ~ sea ice conc. / zdmean
            ! => number of floes encountered by waves travelling meridionally across grid cell:
            IF( zdmean > 0._wp ) THEN
               znfloes = stmer(ji,jj) * at_i(ji,jj) / zdmean   ! stmer(ji,jj) = meridional distance across T cells
            ELSE
               znfloes = 0._wp
            ENDIF

            ! Grid cell mean ice thickness:
            zhi = vt_i(ji,jj) / at_i(ji,jj)

            DO jw = 1, nn_nwfreq
               !
               ! Calculate the (natural) logarithm of the attenuation coefficient for
               ! this period (= 1/frequency) using the quadratic approx. given in HT15:
               !
               zloga = rn_attn_c0 + rn_attn_ch  * zhi             + rn_attn_ct  / wfreq(jw)      &
                  &               + rn_attn_ch2 * zhi**2          + rn_attn_ct2 / wfreq(jw)**2   &
                  &               + rn_attn_cht * zhi / wfreq(jw)
               !
               ! Save the exponent of the damping factor:
               zattxp(jw) = znfloes * rn_attn_tun * EXP(zloga)
               !
            ENDDO

            ! Wave spectrum: force to be zero under ice (we will either calculate it in step 3,
            ! or there will be no 'source' wave data in which case this needs to be left as zero):
            wspec(ji,jj,:) = 0._wp

         ELSEIF( tmask(ji,jj,1) < 1._wp ) THEN
            zattxp(:) = zf_land      ! set land flag
            wspec(ji,jj,:) = 0._wp   ! set no wave spectrum data (not used anyway)
         ELSE
            zattxp(:) = zf_noice     ! set no-ice flag
            !
            ! Set wave spectrum for this open-ocean grid cell if needed:
            !    if reading spectrum from file/model, then l_attn_calc_spec=F => do nothing
            !    otherwise, calculate from hsw and wpf, but only if hsw is large enough:
            IF( l_attn_calc_spec ) CALL wav_spec_bret( hsw(ji,jj), wpf(ji,jj), wspec(ji,jj,:) )
            !
         ENDIF

         ! ======================================================== !
         ! 2  Transfer data into global domain, global scope arrays !
         ! ======================================================== !
         ji_glo = mig(ji,nn_hls)  ! computational domain indices
         jj_glo = mjg(jj,nn_hls)  ! ---> global domain indices

         ! Fill necessary global arrays:
         zwspec_glo(ji_glo,jj_glo,:) = wspec (ji,jj,:)
         zattxp_glo(ji_glo,jj_glo,:) = zattxp(:)

      END_2D

      CALL wav_merge_glo( zwspec_glo )   ! merge global arrays to/from all processors
      CALL wav_merge_glo( zattxp_glo )   !

      ! Compute logical mask on the global domain identifying ice-covered grid
      ! cells (T) or not (F), deduced from zattxp_glo, for use in loop below:
      ll_mask_ice(:,:) = ( zattxp_glo(:,:,1) > zf_noice )

      ! =========================================== !
      ! 3  Locate source grid cells for each target !
      ! =========================================== !
      DO_2D(0, 0, 0, 0)
         IF( at_i(ji,jj) > epsi06 ) THEN

            ! Compute logical mask on the global domain that identifies a
            ! meridian within a tolerance +/- rn_attn_lam_tol:
            ll_mask_mer(:,:) =       (glamt_glo(:,:) >= (glamt(ji,jj) - rn_attn_lam_tol))   &
               &               .AND. (glamt_glo(:,:) <= (glamt(ji,jj) + rn_attn_lam_tol))

            ! Locate the source grid cell (nearest open-ocean grid cell).
            ! This is the most poleward latitude along the meridian defined by
            ! mask l_mask_mer that is not ice covered, is more equatorward than
            ! the target grid cell latitude, has wave data, and is not land.
            !
            ! Compute logical mask on the global domain that identifies points
            ! equatorward of the target grid cell (depends on hemisphere):
            IF( gphit(ji,jj) > 0._wp ) THEN
               !
               ! == Target grid cell is in the northern hemisphere == !
               !
               ll_mask_eqt(:,:) = gphit_glo(:,:) < gphit(ji,jj)
               !
               ! Identify (global array) indices of the source grid cell:
               isource = MAXLOC( gphit_glo(:,:),   &
                  &              MASK=ll_mask_mer .AND. (.NOT. ll_mask_ice) .AND. ll_mask_eqt)
               !
            ELSE
               !
               ! == Target grid cell is in the southern hemisphere == !
               !
               ll_mask_eqt(:,:) = gphit_glo(:,:) > gphit(ji,jj)
               !
               ! Identify (global array) indices of the source grid cell:
               isource = MINLOC( gphit_glo(:,:),   &
                  &              MASK=ll_mask_mer .AND. (.NOT. ll_mask_ice) .AND. ll_mask_eqt)
               !
            ENDIF

            IF( (SUM(isource) > 1) .AND. (zattxp_glo(isource(1),isource(2),1) > zf_land) ) THEN
               !
               ! Found source grid cell: compute net attenuation exponent at the target
               ! grid cell. Need additional mask for grid cells poleward of source:
               IF( gphit(ji,jj) > 0._wp ) THEN
                  ll_mask_pol(:,:) = (gphit_glo(:,:) >= gphit_glo(isource(1),isource(2)))
               ELSE
                  ll_mask_pol(:,:) = (gphit_glo(:,:) <= gphit_glo(isource(1),isource(2)))
               ENDIF
               !
               ! Final mask (all conditions) of grid cells to integrate attenuation:
               ll_mask_tot(:,:) =       ll_mask_mer(:,:) .AND. ll_mask_ice(:,:)   &
                  &               .AND. ll_mask_eqt(:,:) .AND. ll_mask_pol(:,:)
               !
               DO jw = 1, nn_nwfreq
                  ! Total attenuation exponent for this frequency:
                  zattxp_tot = SUM(zattxp_glo(:,:,jw), MASK=ll_mask_tot)
                  !
                  ! Attenuated wave spectrum at target grid cell:
                  wspec(ji,jj,jw) = zwspec_glo(isource(1),isource(2),jw) * EXP(-zattxp_tot)
               ENDDO
               !
            ENDIF
            !
            ! =========================================================== !
            ! 4  Calculate attenuated wave properties at target grid cell !
            ! =========================================================== !
            ! Important to do this regardless of whether source grid cell found:
            !    if there was, wspec will have just been updated as attenuated value
            !       otherwise, wspec will be 0 here (set in step 1 of this routine)
            ! Either way, wave fields hsw, wpf, and wmp need to be (re)calculated under ice
            !
            CALL ice_wav_calc( wspec(ji,jj,:), hsw(ji,jj), wpf(ji,jj), wmp(ji,jj) )
            !
         ENDIF
      END_2D

      CALL lbc_lnk('ice_wav_attn', hsw(:,:)    , 'T', 1._wp, wpf(:,:), 'T', 1._wp, wmp(:,:), 'T', 1._wp)
      CALL lbc_lnk('ice_wav_attn', wspec(:,:,:), 'T', 1._wp)

      IF( ln_timing )   CALL timing_stop('ice_wav_attn')

   END SUBROUTINE ice_wav_attn


   SUBROUTINE ice_wav_calc( pwspec, phsw, pwpf, pwmp )
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE ice_wav_calc ***
      !!-------------------------------------------------------------------
      !!
      !! ** Purpose :   Calculate wave properties from the wave spectrum
      !!
      !! ** Method  :   Significant wave height, Hs, is given by:
      !!
      !!                   Hs = 4 * sqrt{ int[ E(f)df ] }
      !!
      !!                where E(f) is the wave energy spectrum as a function of
      !!                frequency, f. Peak frequency, fp, satisfies:
      !!
      !!                   E(fp) = max[ E(f) ]
      !!
      !!                The wave mean period is given by the ratio of zeroth to the
      !!                first moments of the wave spectrum:
      !!
      !!                   Tm = int[ E(f)df ] / int[ fE(f)df ]
      !!
      !!                See, e.g., WMO (2018).
      !!
      !! ** Inputs  :   pwspec : wave energy spectrum at one location,
      !!                         as a function of frequency (m2.Hz-1)
      !!
      !! ** Outputs :   phsw   : significant wave height (m)
      !!                pwpf   : wave peak frequency (Hz)
      !!                pwmp   : wave mean period (s)
      !!
      !! ** Callers :   ice_wav_attn --> [ice_wav_calc]
      !!
      !! ** References
      !!    ----------
      !!    World Meteorological Organization (WMO), 2018.
      !!              Guide to Wave Analysis and Forecasting.
      !!              2018 ed., Geneva, Switzerland.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(nn_nwfreq), INTENT(in)    ::   pwspec   ! wave energy spectrum (m2.Hz-1)
      REAL(wp),                       INTENT(inout) ::   phsw     ! significant wave height (m)
      REAL(wp),                       INTENT(inout) ::   pwpf     ! wave peak frequency (Hz)
      REAL(wp),                       INTENT(inout) ::   pwmp     ! wave mean period (s)
      !
      INTEGER  ::   imax    ! index of maximum in pwspec
      REAL(wp) ::   zm0     ! zeroth-moment of the wave spectrum (m2)
      REAL(wp) ::   zm1     ! first-moment of the wave spectrum (m2.s-1)
      !
      !!-------------------------------------------------------------------

      IF( MAXVAL(pwspec(:)) >= minwspec ) THEN
         ! Moments of the wave spectrum:
         zm0 = SUM(            pwspec(:) * wdfreq(:) )
         zm1 = SUM( wfreq(:) * pwspec(:) * wdfreq(:) )

         ! Significant wave height:
         phsw = 4._wp * SQRT( zm0 )

         ! Wave peak frequency:
         imax = MAXLOC(pwspec(:), DIM=1)
         pwpf = wfreq(imax)

         ! Wave mean period:
         pwmp = zm0 / zm1
      ELSE
         phsw = 0._wp
         pwmp = 0._wp
         pwpf = 0._wp
      ENDIF

   END SUBROUTINE ice_wav_calc


   SUBROUTINE ice_wav_frac( kt )
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE ice_wav_frac  ***
      !!
      !! ** Purpose :   Evolve the floe size distribution (FSD) subject to fracture by ocean waves
      !!
      !! ** Method  :   1.   Calculate wave energy spectrum, E(f), if needed. This is done
      !!                     when ln_wave_spec=F OR ln_ice_wav_spec=F, in which case E(f) is
      !!                     estimated as a Bretschneider spectrum from local (under ice)
      !!                     significant wave height and peak frequency (from file or wave model).
      !!
      !!                     If in addition ln_ice_wav_attn=T, then this step is skipped as E(f)
      !!                     under ice will already have been calculated in subroutine ice_wav_attn.
      !!
      !!                     If ln_wave_spec=T AND ln_ice_wav_spec=T, this step is also skipped
      !!                     as these options mean the wave spectrum is being read in from a file
      !!                     or (more likely) from a wave model.
      !!
      !!                2.   Calculate the source terms Q(s) and B(s',s) in the equation for the
      !!                     tendency of the floe size distribution, f(s), due to wave fracture:
      !!
      !!                        df(s)/dt = -Q(s)*f(s) + int [ B(s',s)*Q(s')*f(s') ds' ]
      !!
      !!                     The first term on the left-hand side represents loss of floes of size s
      !!                     due to them fracturing into smaller floes, and the second term represents
      !!                     gain of floes of size s due to fracturing of floes of size s'.
      !!
      !!                     Q(s) is the probability that floes of size s are fractured by waves, per
      !!                     unit time, and B(s',s) is a redistribution function quantifying how floes
      !!                     of size s' are transferred to size s; specifically, B(s',s)ds is the fraction
      !!                     of the original area of floes in the size interval [s', s'+ds'] transferred
      !!                     into the size interval [s,s+ds]. Hence, B(s',s) is normalised for each s'
      !!                     such that: int[ B(s',s)ds ] = 1.
      !!
      !!                     This general formulation from Zhang et al. (2015) allows different wave
      !!                     fracture schemes to be used via different choices of Q(s) and B(s',s).
      !!                     These are done by the local subroutines wav_frac_*
      !!
      !!                     In the discretised implementation, f(s) is replaced by L(s) which
      !!                     corresponds to model variables a_ifsd(ji,jj,jf,jl) / floe_dr(jf), while
      !!                     f(s')ds' in the integral is replaced directly with a_ifsd(ji,jj,:,jl).
      !!                     Hence there is an additional factor of floe_ds(jf) that is multiplied
      !!                     through and combined with the B(s',s) for which the local variable zBfrac
      !!                     then corresponds to B(s',s)ds (strictly speaking, B integrated over the
      !!                     floe size category interval, representing the total area fraction of
      !!                     fractured floes transferred into that category).
      !!
      !!                3.   Evolve the FSD using adaptive time stepping (Horvat and Tziperman, 2017).
      !!                     Additional time-step restrictions apply when evolving f(s), so a smaller step
      !!                     is (possibly) required, calculated in subroutine ice_fsd_tstep in module icefsd.
      !!
      !! ** Callers :   ice_stp -> [ice_wav_frac]
      !! ** Calls   :              [ice_wav_frac] -> wav_frac_{z16,y24a,y24b,ht15}
      !!                           [ice_wav_frac] -> ice_fsd_tstep
      !!                           [ice_wav_frac] -> wav_spec_bret()   (ln_ice_wav_spec = F AND ln_ice_wav_attn = F)
      !!
      !! ** Notes   :   ** On the computation of Bfrac and Qfrac for new fracture schemes
      !!
      !!                   Qfrac(jf) represents the floe-category mean fracture probability per unit time.
      !!
      !!                   In principle (depending on the assumptions of the chosen scheme), Bfrac is non-zero
      !!                   along the diagonal because it is possible for floes to fracture into sizes in the same
      !!                   category range. See ice_wav_init for explanation in terms of the uniform redistributor.
      !!
      !!                   A special case is the smallest floe size category, Bfrac(1,1), **which must always be 1**.
      !!                   This is not enforced by this routine, but it ensures floes in the smallest category that
      !!                   fracture [if Qfrac(1) /= 0] have a loss term always cancelled by its own gain term.
      !!
      !! ** References
      !!    ----------
      !!    Horvat, C., & Tziperman, E. (2015).
      !!              A prognostic model of the sea-ice floe size and thickness distribution.
      !!              The Cryosphere, 9, 2119-2134.
      !!    Horvat, C., & Tziperman, E. (2017).
      !!              The evolution of scaling laws in the sea ice floe size distribution.
      !!              Journal of Geophysical Research: Oceans, 122(9), 7630-7650.
      !!    Zhang, J., Schweiger, A., Steele, M., & Stern, H. (2015).
      !!              Sea ice floe size distribution in the marginal ice zone: Theory and numerical experiments.
      !!              Journal of Geophysical Research: Oceans, 120(5), 3484-3498.
      !!
      !!-------------------------------------------------------------------
      !
      INTEGER , INTENT(in)                    ::   kt                   ! ocean time step
      !
      REAL(wp), DIMENSION(A2D(0),nn_nfsd,jpl) ::   za_ifsdb             ! a_ifsd before fracture (for diagnostics)
      !
      REAL(wp), DIMENSION(nn_nfsd,nn_nfsd)    ::   zBfrac               ! fracture redistribution function B(s',s)ds
      REAL(wp), DIMENSION(nn_nfsd)            ::   zQfrac               ! fracture probability function (s-1)
      REAL(wp), DIMENSION(nn_nfsd)            ::   za_ifsd_tend         ! tendency of FSD due to wave fracture
      REAL(wp)                                ::   zh_i                 ! mean ice thickness
      REAL(wp)                                ::   zfsd_res             ! correction term for area conservation
      REAL(wp)                                ::   zt_elapsed           ! time elapsed during adaptive time stepping (units: s)
      INTEGER                                 ::   isubt                ! number of iterations used in adaptive time stepping
      INTEGER                                 ::   ji, jj, jl, jf       ! dummy loop indices
      !
      !                                                                 ! --  Z16 scheme only -- !
      REAL(wp)                                ::   z1mfetch, znumcell   ! for local fetch parameter calculation
      INTEGER                                 ::   jx, jy               ! more dummy loop indices
      !
      !                                                                 ! -- HT15 scheme only -- !
      INTEGER                                 ::   zstat, znumfrac      ! to track how often convergence not reached
      !
      REAL(wp), PARAMETER                     ::   zat_i_min = .01_wp   ! minimum concentration for fracture to occur
      !
      !!-------------------------------------------------------------------

      ! Control:
      IF( ln_timing )   CALL timing_start('ice_wav_frac')

      z1mfetch = 0._wp   ! initialise
      znumcell = 0._wp

      IF( kt == nit000 ) THEN   ! at first time-step
         ! Compute constant weight terms for Y24B fracture scheme:
         IF( nn_frac_scheme == jpfrac_y24b ) CALL y24b_weights
      ENDIF

      za_ifsdb(A2D(0),:,:) = a_ifsd(A2D(0),:,:)   ! save a_ifsd before fracture for tendency diagnostics

      zstat    = 0   ! counter, number of grid points for which convergence NOT reached (HT15 only)
      znumfrac = 0   ! counter, number of grid points for which HT15 scheme is calculated

      !-----------------!
      ! Begin main loop !
      !-----------------!
      DO_2D( 0, 0, 0, 0 )
         !
         ! Do not calculate fracture for total ice concentration below threshold:
         !
         IF( at_i(ji,jj) > zat_i_min ) THEN
            !
            zQfrac(:)   = 0._wp   ! reset from previous iterations
            zBfrac(:,:) = 0._wp   ! (done in scheme subroutines too but repeat here to be safe)
            !
            ! (1) Calculate wave spectrum, if needed. Condition depends on combination of various
            ! namelist flags; the net condition is saved in module variable l_frac_calc_spec
            !
            IF( l_frac_calc_spec ) CALL wav_spec_bret( hsw(ji,jj), wpf(ji,jj), wspec(ji,jj,:) )

            ! (2) Calculate source terms for the wave fracture equation
            !     This depends on the fracture scheme selected
            !
            ! Terms in common among > 1 schemes:
            zh_i = vt_i(ji,jj) / at_i(ji,jj)   ! mean ice thickness
            !
            SELECT CASE( nn_frac_scheme )
               CASE( jpfrac_z16  )
                  !
                  IF( .NOT. ln_z16_const ) THEN
                     ! Calculate 1 minus the fetch parameter for Z16 participation factor:
                     z1mfetch = 0._wp   ! 1 - f0, where f0 = mean open water fraction
                     znumcell = 0._wp   ! number of grid cells averaged over
                     !
                     ! Open water fraction is averaged over local + N surrounding grid cells
                     ! [(2N+1)x(2N+1) grid cells, where N is namelist par. nn_z16_nf
                     !
                     ! We do not necessarily average over all (2N+1)^2 cells, as some may be
                     ! land or not exist if we are close to the model domain edge
                     DO jx = MAX(1,ji-nn_z16_nf), MIN(jpi, ji+nn_z16_nf)
                        DO jy = MAX(1,jj-nn_z16_nf), MIN(jpj, jj+nn_z16_nf)
                           IF( tmask(jx,jy,1) > 0._wp ) THEN   ! => not a land cell
                              z1mfetch = z1mfetch + SUM(a_i(jx,jy,:))
                              znumcell = znumcell + 1._wp
                           ENDIF
                        ENDDO
                     ENDDO
                     IF( znumcell > 0._wp ) z1mfetch = z1mfetch / znumcell
                  ENDIF
                  !
                  CALL wav_frac_z16( wndm(ji,jj) , zh_i , z1mfetch                    ,   &
                     &               floe_size_dist( a_ifsd(ji,jj,:,:), a_i(ji,jj,:) ),   &
                     &               zQfrac(:)   , zBfrac(:,:)                            )
                  !
               CASE( jpfrac_y24a )
                  !
                  CALL wav_frac_y24a( hsw(ji,jj), wmp(ji,jj), zh_i, zQfrac(:), zBfrac(:,:) )
                  !
               CASE( jpfrac_y24b )
                  !
                  CALL wav_frac_y24b( wmp(ji,jj), wspec(ji,jj,:), zh_i, zQfrac(:), zBfrac(:,:) )
                  !
               CASE( jpfrac_ht15 )
                  !
                  ! Do not do this calculation if local wave spectrum is too weak
                  ! (note: other schemes are much cheaper so similar checks not needed):
                  IF( MAXVAL( wspec(ji,jj,:) ) >= minwspec ) THEN
                     znumfrac = znumfrac + 1   ! add to number of grid cells undergoing fracture
                     CALL wav_frac_ht15( wspec(ji,jj,:), zh_i, zQfrac(:), zBfrac(:,:), zstat )
                  ELSE
                     zQfrac(:)   = 0._wp
                     zBfrac(:,:) = 0._wp
                  ENDIF
                  !
               CASE DEFAULT
                  !
                  zQfrac(:)   = 0._wp
                  zBfrac(:,:) = 0._wp
                  !
            END SELECT

            ! Proceed only if some fractures can occur, implied by non-zero zQfrac.
            ! Fracturing quantified by zQfrac is applied to each ice thickness category
            ! if possible (enough ice to begin with), in proportion to its concentration
            !
            IF( MAXVAL(zQfrac(1:)) > epsi10 ) THEN
               !--------------------------------------!
               ! Begin sub-loop: thickness categories !
               !--------------------------------------!
               DO jl = 1, jpl
                  !
                  CALL ice_fsd_cor( a_ifsd(ji,jj,:,jl) )   ! necessary? Small value corrections/renormalisation
                  !
                  ! Conditions for wave fracture to occur in this ITD category:
                  !  - ice concentration (in this category) cannot be too small
                  !  - cannot have all ice in smallest floe size category
                  !  - total FSD cannot be zero (needed?)
                  !
                  IF(          (a_i(ji,jj,jl) > epsi06)                     &
                     &   .AND. (a_ifsd(ji,jj,1,jl) < 1._wp)                 &
                     &   .AND. (SUM(a_ifsd(ji,jj,:,jl)) >  epsi10) ) THEN
                     !
                     ! (3) Evolve the FSD with adaptive time stepping
                     !
                     ! Initialise:
                     zt_elapsed = 0._wp   ! time elapsed during adaptive time stepping
                     isubt      = 0       ! number of sub time steps taken
                     !
                     DO WHILE (zt_elapsed < rDt_ice)
                        !
                        ! Exit loop if all ice already in smallest floe size category:
                        IF( a_ifsd(ji,jj,1,jl) >= 1._wp - epsi10 ) EXIT
                        !
                        ! Calculate FSD tendency due to wave fracture:
                        DO jf = 1, nn_nfsd
                           za_ifsd_tend(jf) = SUM( zBfrac(:,jf) * zQfrac(:) * a_ifsd(ji,jj,:,jl) )   &
                              &               - zQfrac(jf) * a_ifsd(ji,jj,jf,jl)
                        ENDDO
                        !
                        WHERE( ABS(za_ifsd_tend) < epsi10 ) za_ifsd_tend = 0._wp
                        !
                        ! Evolve a_ifsd over maximum stable time step and increase zt_elapsed accordingly:
                        CALL ice_fsd_tstep('ice_wav_frac', a_ifsd(ji,jj,:,jl), za_ifsd_tend(:), zt_elapsed, isubt)
                        !
                     ENDDO ! while loop
                     !
                     ! === Corrections and normalisation ===
                     !
                     ! The implementation of the wave fracture equation may lead to an FSD
                     ! in this category that no longer sums to exactly 1; it must, by
                     ! definition (a_ifsd is the 'per ITD category' FSD), and if it does
                     ! not, that represents loss or gain of sea ice area.
                     !
                     ! Since wave fracture physically cannot/should not lead to loss of sea
                     ! ice area (at least, not directly), correct for any residual FSD here
                     ! by adding it back to the smallest floe size category (if some area
                     ! is lost) or by taking it away from the largest floe size category
                     ! that has at least that residual available, if some area has been gained.
                     !
                     zfsd_res = SUM(a_ifsd(ji,jj,:,jl)) - 1._wp
                     !
                     IF( zfsd_res <= 0._wp ) THEN
                        a_ifsd(ji,jj,1,jl) = a_ifsd(ji,jj,1,jl) + ABS(zfsd_res)
                     ELSE
                        DO jf = nn_nfsd, 1, -1
                           IF( a_ifsd(ji,jj,jf,jl) > zfsd_res) THEN
                              a_ifsd(ji,jj,jf,jl) = a_ifsd(ji,jj,jf,jl) - ABS(zfsd_res)
                              EXIT
                           ENDIF
                        ENDDO
                     ENDIF
                     !
                     CALL ice_fsd_cor( a_ifsd(ji,jj,:,jl) )   ! small value corrections/re-normalisation
                     !
                  ENDIF ! category jl can fracture
               ENDDO ! -- sub-loop (ice thickness categories)
            ENDIF ! ----- MAXVAL(zQfrac) > 0
         ENDIF ! -------- at_i(ji,jj) > zat_i_min
      END_2D ! ---------- main loop

      ! Write warnings if HT15 scheme did not converge for some grid cells
      IF( nn_frac_scheme == jpfrac_ht15 ) THEN
         IF ( ln_ht15_warn .AND. (zstat > 0) ) THEN
            ! Use WRITE directly as ctl_warn only writes on the rank 0 processor
            ! This means some warnings go into stdout, but not sure how else to do it...
            WRITE(numout,'(A,I5,A,I5,A,I10,A,I4)')                                   &
               &   'ice_wav_frac: HT15 scheme did not converge in', zstat, ' of ',   &
               &   znumfrac, ' grid cells, at kt = ', kt, ', rank ', mpprank
            WRITE(numout,*) '             Consider increasing nn_ht15_maxiter or rn_ht15_tol'
            !
            ! ... call ctl_warn anyway with general message, to append warning count:
            CALL ctl_warn('ice_wav_frac: HT15 scheme did not converge everywhere -- check stdout AND ocean.output')
            !
         ENDIF
      ENDIF

      ! Write FSD tendency diagnostics due to wave fractue:
      CALL ice_fsd_dia( 'wav', za_ifsdb, a_ifsd(A2D(0),:,:), a_i(A2D(0),:), a_i(A2D(0),:) )

      ! Control:
      IF( ln_timing )   CALL timing_stop('ice_wav_frac')

   END SUBROUTINE ice_wav_frac


   SUBROUTINE wav_frac_z16( puatm, ph_i, p1mfetch, pfsd, pQfrac, pBfrac )
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE wav_frac_z16 ***
      !!
      !! ** Purpose :   Calculate the probability and redistribution functions for the
      !!                equation evolving the FSD due to wave fracture for the scheme of
      !!                Zhang et al. (2016).
      !!
      !! ** Method  :   Fracture scheme that does *not* use wave forcing data (file or coupled model).
      !!                The probability of a floe of size s undergoing fracture, Q(s), is given by:
      !!
      !!                   Q(s) = MAX[ 0, 1 - int( f(s')ds' )/cb ]
      !!
      !!                where f(s') is the floe size distribution (integrated over thickness) and cb
      !!                is called a 'participation factor' representing the area fraction of ice that
      !!                could participate in fracture. In Zhang et al. (2015) this is set to a constant
      !!                and this behaviour can be achieved by setting namelist parameter ln_z16_const=T
      !!                and rn_z16_cb as required. In Zhang et al. (2016), this scheme was upgraded to
      !!                work in a GCM setting and cb is instead calculated from local wind speed U and
      !!                ice properties, as follows, which is the behaviour when ln_z16_const=F:
      !!
      !!                   c_b = [kU/MAX(hi,hc)] * EXP[ -a(1 - f0) - b(1 - savg/smax) ] * dt
      !!
      !!                where k, a, and b are dimensionless constants (set in namelist), hi is mean ice
      !!                thickness, hc is a cutoff thickness (namelist), f0 is open water fraction,
      !!                savg is mean floe size, smax is the largest resolved floe size, and dt is the model
      !!                timestep. The open water fraction is calculated as an average over the
      !!                (2N+1) x (2N+1) grid cells centred on the local grid cell (N = nn_z16_nf).
      !!
      !!                The redistribution function is determined by assuming any floe of size s'
      !!                that is undergoing wave fracture is equally likely to fracture into any
      !!                other floe size s < s'. This 'unifom redistributor' is calculated in
      !!                subroutine ice_wav_init as it is a constant.
      !!
      !! ** Inputs  :   puatm                   :   local wind speed (m/s)
      !!                ph_i                    :   local mean ice thickness (m)
      !!                p1mfetch                :   1 minus local fetch parameter
      !!                pfsd                    :   floe size distribution, f(s)ds
      !!
      !! ** Outputs :   pQfrac(nn_nfsd)         :   fracture probability function (s-1)
      !!                pBfrac(nn_nfsd,nn_nfsd) :   fracture redistribution function, B(s',s)ds
      !!                                            (note: first  index corresponds to original floe size s',
      !!                                                   second index corresponds to fractured floe size s)
      !!
      !! ** Callers :   ice_wav_frac --> [wav_frac_z16]
      !!
      !! ** Notes   :   A key distinction of this scheme (beyond it not using wave data) is that it
      !!                technically includes any wind-driven fracture process. In the MIZ, this is
      !!                dominantly wave mediated, but in the pack ice it is mediated by deformation
      !!                and internal stresses. So, with this scheme, there may be a need to adjust
      !!                the brittle fracture (routine ice_fsd_brit) parameters to compensate.
      !!
      !! ** References
      !!    ----------
      !!    Zhang, J., Schweiger, A., Steele, M., & Stern, H. (2015).
      !!              Sea ice floe size distribution in the marginal ice zone: Theory and numerical experiments.
      !!              Journal of Geophysical Research: Oceans, 120(5), 3484-3498.
      !!    Zhang, J., Stern, H., Hwang, B., Schweiger, A., Steele, M., Stark, M., & Graber, H.C.
      !!              Modeling the seasonal evolution of the Arctic sea ice floe size distribution.
      !!              Elementa, 4(000126)
      !!-------------------------------------------------------------------
      !
      REAL(wp)                            , INTENT(in)    ::   puatm     ! local near-surface wind speed (m/s)
      REAL(wp)                            , INTENT(in)    ::   ph_i      ! local mean sea ice thickness (m)
      REAL(wp)                            , INTENT(in)    ::   p1mfetch  ! 1 minus fetch parameter
      REAL(wp), DIMENSION(nn_nfsd)        , INTENT(in)    ::   pfsd      ! floe size distribution, f(s)ds
      REAL(wp), DIMENSION(nn_nfsd)        , INTENT(inout) ::   pQfrac    ! wave fracture probability function (s-1)
      REAL(wp), DIMENSION(nn_nfsd,nn_nfsd), INTENT(inout) ::   pBfrac    ! wave fracture redistribution function, B(s',s)ds
      !
      REAL(wp)                                            ::   zsavg     ! mean floe size
      REAL(wp)                                            ::   zc_b      ! participation factor
      INTEGER                                             ::   jf, jl    ! dummy loop index
      !
      !!-------------------------------------------------------------------

      ! --- Calculate participation factor, zc_b --- !
      !
      IF( ln_z16_const ) THEN
         zc_b = rn_z16_cb   ! constant participation factor (Zhang et al. 2015, Eq. 15)
      ELSE
         ! Variable participation factor (as in Zhang et al. 2016, Eq. 5)
         !
         ! Mean floe size:
         zsavg = 0._wp
         DO jf = 1, nn_nfsd
            zsavg = zsavg + floe_sc(jf) * pfsd(jf)
         ENDDO
         !
         zc_b = EXP( -rn_z16_a * p1mfetch - rn_z16_b * (1._wp - zsavg / floe_sc(nn_nfsd)) )
         zc_b = zc_b * rn_z16_k * puatm * rDt_ice / MAX( ph_i, rn_z16_hc )
      ENDIF

      ! --- Calculate source terms --- !
      DO jf = 1, nn_nfsd
         pQfrac(jf) = MAX( 0._wp, 1._wp - SUM(pfsd(jf:nn_nfsd)) / zc_b ) * r1_Dt_ice
      ENDDO

      pBfrac(:,:) = Bfrac_uni(:,:)   ! uniform redistribution

   END SUBROUTINE wav_frac_z16


   SUBROUTINE wav_frac_y24a( phsw, pwmp, ph_i, pQfrac, pBfrac )
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE wav_frac_y24a ***
      !!
      !! ** Purpose :   Calculate the probability and redistribution functions for the
      !!                equation evolving the FSD due to wave fracture for the first scheme
      !!                of Yang et al. (2024), called scheme A here.
      !!
      !! ** Method  :   Assume waves break ice if the strain e exceeds a critical value,
      !!                where the strain is calculated from input wave statistics:
      !!
      !!                   e = 5 * pi^4 * h_i * H_s / (2 * g^2 * T_m^4)
      !!
      !!                where h_i is mean ice thickness, H_s is significant wave height, g is
      !!                gravity, and T_m is the mean wave period. This formula is as in Yang et al.
      !!                (2024, Eqs. 15-16) except mean wavelength has been converted to T_m using
      !!                the deep water surface gravity wave dispersion relation, lambda = g*T^2/(2*pi).
      !!
      !!                If e >= e_crit (namelist parameter rn_ice_wav_ecri), the probability of floes
      !!                fracturing is given by a simple exponential profile:
      !!
      !!                   Q(s) = c_w * exp[ -alpha * (1 - s/smax)]
      !!
      !!                where c_w is a constant (sets maximum probability of the largest floes
      !!                fracturing), alpha is a constant (determines strength of probability
      !!                reduction at smaller floe sizes), and smax is the largest resolved floe size.
      !!                This formula differs slightly from Eq. (13) in Yang et al. (2024) which
      !!                is presumed to have a typesetting mistake in the exponent.
      !!
      !!                The redistribution function is determined by assuming any floe of size s'
      !!                that is undergoing wave fracture is equally likely to fracture into any
      !!                other floe size s < s'. This 'uniform redistributor' is calculated in
      !!                subroutine ice_wav_init as it is a constant.
      !!
      !! ** Inputs  :   phsw                    :   local significant wave height (m)
      !!                pwmp                    :   local wave mean period (s)
      !!                ph_i                    :   local (grid cell) mean sea ice thickness (m)
      !!
      !! ** Outputs :   pQfrac(nn_nfsd)         :   fracture probability function (s-1)
      !!                pBfrac(nn_nfsd,nn_nfsd) :   fracture redistribution function, B(s',s)ds
      !!                                            (note: first  index corresponds to original floe size s',
      !!                                                   second index corresponds to fractured floe size s)
      !!
      !! ** Callers :   ice_wav_frac --> [wav_frac_y24a]
      !!
      !! ** References
      !!    ----------
      !!    Yang, C.-Y., Liu, J., & Chen, D. (2024).
      !!              Understanding the influence of ocean waves on Arctic sea ice simulation:
      !!              a modelling study with an atmosphere-ocean-wave-sea ice coupled model.
      !!              The Cryosphere, 18(3), 1215-1239.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp)                            , INTENT(in)    ::   phsw     ! grid cell significant wave height (m)
      REAL(wp)                            , INTENT(in)    ::   pwmp     ! grid cell wave mean period (s)
      REAL(wp)                            , INTENT(in)    ::   ph_i     ! grid cell mean ice thickness (m)
      REAL(wp), DIMENSION(nn_nfsd)        , INTENT(inout) ::   pQfrac   ! wave fracture probability function (s-1)
      REAL(wp), DIMENSION(nn_nfsd,nn_nfsd), INTENT(inout) ::   pBfrac   ! wave fracture redistribution function, B(s',s)ds
      !
      INTEGER             ::   jf                ! dummy loop indices
      REAL(wp)            ::   zstrain           ! ice strain due to waves
      !
      !!-------------------------------------------------------------------

      pQfrac(:) = 0._wp   ! reset/initialise, and return value if critical strain not exceeded

      IF( (phsw >= minhsw) .AND. (pwmp >= minwmp) ) THEN                 ! <-- sufficient wave presence
         zstrain = 2.5_wp * rpi**4 * ph_i * phsw / (grav**2 * pwmp**4)   ! <-- strain experienced by ice
         IF( zstrain >= rn_ice_wav_ecri ) THEN
            DO jf = 1, nn_nfsd
               pQfrac(jf) = rn_y24a_cw * r1_Dt_ice   &
                  &                    * EXP( -rn_y24a_alpha * (1._wp - floe_sc(jf) / floe_sc(nn_nfsd)) )
            ENDDO
         ENDIF
      ENDIF

      pBfrac(:,:) = Bfrac_uni(:,:)   ! uniform redistribution

   END SUBROUTINE wav_frac_y24a


   SUBROUTINE y24b_weights
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE y24b_weights  ***
      !!
      !! ** Purpose :   Calculate 'weight' factors for Y24b wave fracture scheme.
      !!                ** This routine should only be called once as these are constants **
      !!
      !! ** Method  :   'Weights' are a function of spectral class and floe size category and
      !!                represent the sub-category--class fraction that satisfy the imposed
      !!                conditions on wave fracture, which in Y24b are functions of wavelength
      !!                and floe size only. In the actual scheme, local wave conditions are
      !!                used to determine whether fracture actually could occur as a function
      !!                of wavelength which, combined with these weights, determines the wave
      !!                fracture source terms (see routines ice_wav_frac and wav_frac_y24b).
      !!
      !!                The formulae/conditions used below are not straightforward to explain
      !!                in text... wgtQ_y24b represents the area fraction of the box ds by dL
      !!                (floe size and wavelength category widths) satisfying s <= L, and
      !!                wgtB_y24b is the fraction of the interval dL for which L/2 is contained
      !!                by each floe size category. In most cases (well, for wgtQ_y24B) these are
      !!                0 or 1, where the conditions are satisfied by none or all sub-category
      !!                values of s and L. The remaining value are between 0 and 1 and the exact
      !!                formula to calculate them depends on exactly how the category bounds
      !!                geometrically line up with the relevant condition (s <= L or s' = L/2).
      !!
      !! ** Outputs :   Allocates and calculates the private module variables: wgt{Q,B}_y24b
      !!
      !!-------------------------------------------------------------------
      !
      INTEGER ::   jf, jw, ierr   ! dummy loop indices and allocate status return value
      !
      !!-------------------------------------------------------------------

      ALLOCATE( wgtQ_y24b(nn_nwfreq,nn_nfsd), wgtB_y24b(nn_nwfreq,nn_nfsd), STAT=ierr )
      IF( ierr /= 0 ) CALL ctl_stop('icewav: unable to allocate wgt{Q,B}_y24b array(s)')

      wgtQ_y24b(:,:) = 0._wp   ! initialise
      wgtB_y24b(:,:) = 0._wp   ! default value (most of this array will = 0 anyway)

      DO jw = 1, nn_nwfreq
         DO jf = 1, nn_nfsd
            !
            ! --- Compute weights for Q(s) --- !
            !
            IF( wlam_u(jw) <= floe_sl(jf) ) THEN
               ! => all wavelengths in class jw can fracture any floes in category jf
               wgtQ_y24b(jw,jf) = 1._wp
               !
            ELSEIF( wlam_l(jw) >= floe_su(jf) ) THEN
               ! => no wavelengths in class jw can fracture any floes in category jf
               wgtQ_y24b(jw,jf) = 0._wp
               !
            ELSE
               ! => only some wavelength/floe size combinations sub-spectral class/
               !    sub-floe size category can participate in fracture, so we have to work
               !    out the fraction that can. Hard to explain where the following four expressions
               !    come from in text... sketch out the various cases of a box of width/height
               !    (ds,dL) overlapping the region L < s, and get the area fraction of overlap...
               !
               IF( wlam_u(jw) >= floe_su(jf) ) THEN
                  IF( wlam_l(jw) < floe_sl(jf) ) THEN
                     wgtQ_y24b(jw,jf) = (floe_sl(jf) - wlam_l(jw) + .5_wp * floe_ds(jf)) / wdlam(jw)
                  ELSE
                     wgtQ_y24b(jw,jf) = .5_wp * (floe_su(jf) - wlam_l(jw))**2 / (floe_ds(jf) * wdlam(jw))
                  ENDIF
               ELSE   ! floe_su(jf) > wlam_u(jw)
                  IF( wlam_l(jw) < floe_sl(jf) ) THEN
                     wgtQ_y24b(jw,jf) = 1._wp - .5_wp * (wlam_u(jw) - floe_sl(jf))**2 / (floe_ds(jf) * wdlam(jw))
                  ELSE
                     wgtQ_y24b(jw,jf) = (floe_su(jf) - wlam_u(jw) + .5_wp * wdlam(jw) ) / floe_ds(jf)
                  ENDIF
               ENDIF
            ENDIF
            !
            ! --- Compute weights for B(r,r') --- !
            ! Weight wgtB_y24b(jw,jf) represents fraction of wavelength interval jw for which
            ! s' = L/2 is contained within the floe size category interval jf. Again these are
            ! hard to explain, just have to sketch out the possible cases...
            !
            ! If first joint condition below is not satisfied, it means none of the range of
            ! fracture sizes s' = L/2 for this spectral class jw fit into this floe size
            ! category jf, so we do nothing (array initialised to 0 already)
            !
            IF( (wlam_u(jw) > 2._wp * floe_sl(jf)) .AND. (wlam_l(jw) < 2._wp * floe_su(jf)) ) THEN
               IF( wlam_u(jw) >= 2._wp * floe_su(jf) ) THEN
                  IF( wlam_l(jw) < 2._wp * floe_sl(jf) ) THEN
                     wgtB_y24b(jw,jf) = 2._wp * floe_ds(jf) / wdlam(jw)
                  ELSE
                     wgtB_y24b(jw,jf) = (2._wp * floe_su(jf) - wlam_l(jw)) / wdlam(jw)
                  ENDIF
               ELSE   ! 2*floe_su(jf) > wlam_u(jw)
                  IF( wlam_l(jw) < 2._wp * floe_sl(jf) ) THEN
                     wgtB_y24b(jw,jf) = (wlam_u(jw) - 2._wp * floe_sl(jf)) / wdlam(jw)
                  ELSE
                     wgtB_y24b(jw,jf) = 1._wp
                  ENDIF
               ENDIF
            ENDIF
            !
         ENDDO   ! jf (floe size category)
      ENDDO   ! -- jw (spectral class)

   END SUBROUTINE y24b_weights


   SUBROUTINE wav_frac_y24b( pwmp, pWspec, ph_i, pQfrac, pBfrac )
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE wav_frac_y24b ***
      !!
      !! ** Purpose :   Calculate the probability and redistribution functions for the equation
      !!                evolving the FSD due to wave fracture based on the second ('semi empirical')
      !!                fracture scheme of Yang et al. (2024), called scheme B here.
      !!
      !! ** Method  :   This scheme asserts that each spectral component of the local wave field
      !!                can be considered individually to determine whether it fractures ice or not
      !!                under a strain failure criterion as function of frequency:
      !!
      !!                   e(f) =  2 * pi^2 * A(f) * h_i / L(f)^2
      !!
      !!                where A(f) and L(f) are the wave component amplitude and wavelength:
      !!
      !!                   A(f) = sqrt( 2 * E(f) * df(f) )
      !!                   L(f) = g / (2 * pi * f^2)
      !!
      !!                where E(f) is the local wave energy spectrum and df is the spectral bandwidth.
      !!                Note L(f) is calculated already as a module constant. The above e(f) formula
      !!                is the maximum of the strain (h_i/2)y"(x) for a plane wave
      !!                y(x) = A(f)sin(kx + wt) with k = 2pi/L(f).
      !!
      !!                Assuming a Rayleigh distribution of wave frequencies, P(f)df, determines the
      !!                probability, q(L)dL, of each spectral class leading to fracture:
      !!
      !!                   q(L)dL = / P(f)df   e(f) >= e_crit
      !!                            \ 0        e(f) <  e_crit
      !!
      !!                [NB. Rayleigh spectrum is usually expressed in terms of wave period T = 1/f,
      !!                but function wav_spec_rayl() returns it in terms of f (to simplify the code
      !!                by avoiding having to define coordinate arrays for T). Here, it would make
      !!                more sense to express it in terms of wavelength (see conditions on floe
      !!                size below), but since we are always going to integrate over the spectral
      !!                class, we can just go from P(f)df = P(T)dT = P(L)dL directly.]
      !!
      !!                This scheme also makes assumptions relating wavelength L to floe size (diameter) s:
      !!
      !!                   1) Wavelength L can only break a floe of diameter s if L < s
      !!                   2) Fractured floe diameter s is half the triggering wavelength, so s = L/2
      !!
      !!                From this follows:
      !!
      !!                   Q(s) = int[ q(L') ] dL'   where the integral is from L' = 0 to s
      !!
      !!                             int[ q(L') * Theta(s' - L') * delta(2s - L') ] dL'
      !!                   B(s',s) = --------------------------------------------------
      !!                                     int[ q(L') * Theta(s' - L') ] dL'
      !!
      !!                where Theta(x) = {1 if x >= 0; 0 otherwise}, delta(x) = {1 if x = 0; 0 otherwise},
      !!                and the integral is over all wavelengths L. In the code we just check the conditions
      !!                are met explicitly for each combination of s', s, and L rather than using the last
      !!                formal equation. Those are quantified by 'weights' variables, computed once in routine
      !!                y24b_weights, as the conditions (1-2) between each spectral class/floe size
      !!                category, which are 0, 1, or somewhere between (accounting for partial fulfillment
      !!                of the conditions sub-category).
      !!
      !! ** Inputs  :   pwmp                    :   local wave mean period (s)
      !!                pWspec(nn_nwfreq)       :   local wave spectrum (spectral energy density; m2.Hz-1)
      !!                ph_i                    :   local (grid cell) mean sea ice thickness (m)
      !!
      !! ** Outputs :   pQfrac(nn_nfsd)         :   fracture probability function (s-1)
      !!                pBfrac(nn_nfsd,nn_nfsd) :   fracture redistribution function, B(s',s)ds
      !!                                            (note: first  index corresponds to original floe size s',
      !!                                                   second index corresponds to fractured floe size s)
      !!
      !! ** Callers :   ice_wav_frac --> [wav_frac_y24b]
      !!
      !! ** References
      !!    ----------
      !!    Yang, C.-Y., Liu, J., & Chen, D. (2024).
      !!              Understanding the influence of ocean waves on Arctic sea ice simulation:
      !!              a modelling study with an atmosphere-ocean-wave-sea ice coupled model.
      !!              The Cryosphere, 18(3), 1215-1239.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp)                            , INTENT(in)    ::   pwmp     ! grid cell wave mean period (s)
      REAL(wp)                            , INTENT(in)    ::   ph_i     ! grid cell mean ice thickness (m)
      REAL(wp), DIMENSION(nn_nwfreq)      , INTENT(in)    ::   pWspec   ! local wave spectral energy density (m2.Hz-1)
      REAL(wp), DIMENSION(nn_nfsd)        , INTENT(inout) ::   pQfrac   ! wave fracture probability function (s-1)
      REAL(wp), DIMENSION(nn_nfsd,nn_nfsd), INTENT(inout) ::   pBfrac   ! wave fracture redistribution function, B(s',s)ds
      !
      INTEGER                        ::   jf1, jf2, jw   ! dummy loop indices
      REAL(wp), DIMENSION(nn_nwfreq) ::   zamp           ! spectral amplitudes (m)
      REAL(wp), DIMENSION(nn_nwfreq) ::   zstrain        ! ice strain due to waves per spectral class
      REAL(wp), DIMENSION(nn_nwfreq) ::   zprayl         ! Rayleigh spectrum (Hz-1)
      REAL(wp), DIMENSION(nn_nwfreq) ::   zprob          ! probability of fracturing per spectral class
      !
      !!-------------------------------------------------------------------

      zamp    = SQRT( 2._wp * pWspec(:) * wdfreq(:) )            ! spectral amplitudes (m)
      zstrain = 2._wp * rpi**2 * ph_i * zamp(:) / (wlam(:)**2)   ! ice strain per spectral class

      CALL wav_spec_rayl( pwmp, zprayl(:) )                      ! calculate local Rayleigh spectrum (Hz-1)

      ! Probability of waves of each frequency class resulting in fracture:
      zprob(:) = 0._wp
      WHERE( zstrain >= rn_ice_wav_ecri ) zprob = zprayl * wdfreq   ! = q(f)df = q(L)dL

      pQfrac(:)   = 0._wp
      pBfrac(:,:) = 0._wp

      DO jf1 = 1, nn_nfsd
         DO jw = 1, nn_nwfreq
            !
            ! Update Q(s) <==> pQfrac(jf1) (note dL is already implicitly multipled into zprob):
            !
            ! Weight factor wgtQ_y24b is fraction of spectral class jw/floe size category jf
            ! for which the condition for fracture (wavelength L < s') is satisfied:
            pQfrac(jf1) = pQfrac(jf1) + wgtQ_y24b(jw,jf1) * zprob(jw)
            !
            ! Update B(s,r)dr <== > pBfrac(jf1,jf2):
            !
            ! Weight factor wgtB_y24b is the fraction of spectral class jw width for which the
            ! fracture size (s = L/2) is contained within the transferred floe size category jf2:
            DO jf2 = 1, nn_nfsd
                pBfrac(jf1,jf2) = pBfrac(jf1,jf2) +   wgtQ_y24b(jw,jf1) * zprob(jw)      &
                   &                                * wgtB_y24b(jw,jf2) * floe_ds(jf2)
            ENDDO
         ENDDO
         !
         ! Normalise probability rate to time step (units -> s-1):
         pQfrac(jf1) = pQfrac(jf1) * r1_Dt_ice
         !
         ! Ensure B(s',s)ds is normalised (integrates to 1):
         IF( SUM(pBfrac(jf1,:)) > 0._wp ) pBfrac(jf1,:) = pBfrac(jf1,:) / SUM(pBfrac(jf1,:))
         !
      ENDDO

   END SUBROUTINE wav_frac_y24b


   SUBROUTINE wav_frac_ht15( pWspec, ph_i, pQfrac, pBfrac, pstat )
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE wav_frac_ht15 ***
      !!
      !! ** Purpose :   Calculate the probability and redistribution functions for the
      !!                equation evolving the FSD due to wave fracture for the scheme
      !!                of Horvat and Tziperman (2015) and Roach et al. (2018)
      !!
      !! ** Method  :   This scheme calculates a distribution of fractured ice lengths from
      !!                the local sea surface height (SSH) field, n(x), defined along a 1D
      !!                sub-gridscale domain and computed from the local wave spectrum. This
      !!                is defined up to an arbitrary phase per spectral component. That is either
      !!                set to a constant (default) or, if ln_ht15_rand=T, random phases are used
      !!                and the calculations that follow are repeated until the 'fracture
      !!                distribution' (see below) converges to within a specified tolerance.
      !!
      !!                Sea ice is subject to strain due to flexure by the varying SSH
      !!                associated with the local wave field. In this subroutine, ice of
      !!                thickness h_i is assumed/conceptualised to cover the whole 1D sub-
      !!                domain for which the input sea surface height is defined. Ice fractures
      !!                at locations where the strain, e, exceeds a critical threshold, e_crit:
      !!
      !!                   e = 0.5 * h_i * |d^2 n/dx^2| >= e_crit
      !!
      !!                For each extremum in SSH, the nearest neighbouring extrema either side are
      !!                located. The strain is then calculated across such triplets of extrema in
      !!                SSH that are either {min., max., min.} or {max., min., max.} using a finite
      !!                differencing approximation across the triplet. Ice breaks at the central
      !!                extremum if e >= e_crit. The distances between all such breaking points along
      !!                the 1D domain (x) determines the lengths of fractured ice. The number
      !!                distribution of these lengths is then binned into the FSD category bins,
      !!                and the result is called the 'fracture distribution', W(s), such that
      !!                W(s)ds is the number of fracture lengths in the interval [s, s+ds]. From this,
      !!                the probability function in the wave fracture equation is:
      !!
      !!                   Q(s) = int[ s'W(s') ds' ] / D
      !!
      !!                where the integral limits are from the smallest floe size up to s, and the
      !!                redistribution function:
      !!
      !!                   B(s',s) = sW(s) / int[ xW(x) dx ]   for s < s'
      !!                           = 0                         otherwise
      !!
      !!                where the integral limits are from the smallest floe size up to s'.
      !!                Note that int[ B(s',s)ds ] = 1 as required (see subroutine ice_wav_frac).
      !!
      !!                W(s) satisfies int[ xW(x) dx ] = D, where the integral is over all floe
      !!                sizes, since the sum of all fracture lengths must equal the domain size, D.
      !!                So, Q(s) is the fraction of all fracture lengths smaller than s. This
      !!                normalisation factor cancels in the expression for B so there are no explicit
      !!                factors of D.
      !!
      !! ** Inputs  :   pWspec(nn_nwfreq)       :   local wave spectrum (spectral energy density; m2.Hz-1)
      !!                ph_i                    :   local (grid cell) mean sea ice thickness (m)
      !!                pstat                   :   integer, adds 1 to this if convergence is not reached
      !!
      !! ** Outputs :   pQfrac(nn_nfsd)         :   fracture probability function (s-1)
      !!                pBfrac(nn_nfsd,nn_nfsd) :   fracture redistribution function, B(s',s)ds
      !!                                            (note: first  index corresponds to original floe size s',
      !!                                                   second index corresponds to fractured floe size s)
      !!
      !! ** Callers :   ice_wav_frac --> [wav_frac_ht15]
      !!
      !! ** References
      !!    ----------
      !!    Horvat, C., & Tziperman, E. (2015).
      !!              A prognostic model of the sea-ice floe size and thickness distribution.
      !!              The Cryosphere, 9, 2119-2134.
      !!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018).
      !!              An emergent sea ice floe size disribution in a global coupled ocean-sea ice model.
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(nn_nwfreq)      , INTENT(in)    ::   pWspec   ! local wave spectral energy density (m2.Hz-1)
      REAL(wp)                            , INTENT(in)    ::   ph_i     ! grid cell mean ice thickness (m)
      REAL(wp), DIMENSION(nn_nfsd)        , INTENT(inout) ::   pQfrac   ! wave fracture probability function (s-1)
      REAL(wp), DIMENSION(nn_nfsd,nn_nfsd), INTENT(inout) ::   pBfrac   ! wave fracture redistribution function, B(s',s)ds
      INTEGER                             , INTENT(inout) ::   pstat    ! counter for whether convergence reached or not
      !
      INTEGER                              ::   jx, jy, jf              ! dummy loop indices
      INTEGER                              ::   jiter                   ! iteration counter for convergence
      INTEGER                              ::   iloop                   ! number of convergence loop iterations (=1 unless ln_ht15_rand=T)
      INTEGER                              ::   ixlo, ixhi              ! indices of x1d to locate extrema
      INTEGER                              ::   ixfrac                  ! number of fracture points along x1d
      LOGICAL , DIMENSION(nn_ht15_nx1d)    ::   llmin, llmax, llext     ! sea surface height is a min / is a max / is an extrema
      REAL(wp), DIMENSION(nn_nwfreq)       ::   zphi                    ! phase of wave spectrum components (rad)
      REAL(wp), DIMENSION(nn_ht15_nx1d)    ::   zssh                    ! sea surface height along x1d (m)
      REAL(wp)                             ::   zdx, zdxlo, zdxhi       ! distances between x1d points in finite difference computation (m)
      REAL(wp)                             ::   zstrain                 ! strain experienced by sea ice due to wave field
      REAL(wp), DIMENSION(nn_ht15_nx1d)    ::   zxfrac                  ! distances to points along x1d at which ice fractures
      REAL(wp)                             ::   zfrac_s                 ! floe size of a piece of fractured ice (m)
      INTEGER , DIMENSION(nn_nfsd)         ::   iWfrac                  ! fracture distribution as counts in each floe size category
      REAL(wp), DIMENSION(nn_nfsd)         ::   zsWfrac_b, zsWfrac      ! fracture distribution (multiplied by sds; dimensionless;
      !                                                                 !    one extra with _b for saving previous loop iteration)
      REAL(wp), DIMENSION(nn_nfsd)         ::   zsWfrac_err             ! fractional differences in zsWfrac between successive iterations
      !
      !!-------------------------------------------------------------------

      ! Control:
      IF( ln_timing )   CALL timing_start('wav_frac_ht15')

      ! Initialisation (some done at start of while loop below):
      llmin(:)       = .FALSE.
      llmax(:)       = .FALSE.
      iWfrac(:)      = 0         ! fracture length counts per FSD category (accumulates with each loop iteration)
      zsWfrac(:)     = 0._wp     ! fracture distribution multiply by s*ds (updated with each loop iteration)
      zsWfrac_b(:)   = 0._wp     ! zsWfrac at 'before'/previous iteration
      zsWfrac_err(:) = 1._wp     ! zsWfrac fractional error for convergence check (initially 1 => first check fails)
      pQfrac(:)      = 0._wp     ! probability term (returned)
      pBfrac(:,:)    = 0._wp     ! redistributor term (returned)

      ! While loop below is for iterating until convergence IF random phases are being used
      ! Otherwise (constant/uniform spectral phases), only one iteration is needed:
      IF( ln_ht15_rand ) THEN   ;   iloop = nn_ht15_maxiter
      ELSE                      ;   iloop = 1
      ENDIF
      !
      jiter = 0   ! iteration counter (*DO NOT CHANGE* start value: used in normalisation)
      !
      DO WHILE ( (jiter < iloop) .AND. (MAXVAL(zsWfrac_err) > rn_ht15_tol) )

         ! Initialise or reset for new iteration:
         ! *DO NOT RESET iWfrac* -- this accumulates values from all iterations
         ixfrac    = 1
         zxfrac(:) = 0._wp
         jiter     = jiter + 1   ! update now as used in normalisation of zrWfrac later

         ! Spectral phases:
         IF( ln_ht15_rand ) THEN
            ! Phase for each frequency is a random value between [0, 2pi)
            !
            ! NOTE: there are modules 'sto*' with various code structures for using random
            ! numbers and setting/reading the seed. However it does not seem well developed and
            ! only activates for stochastic EOS (it has a framework for other processes --
            ! e.g., ice strength 'p-star' parameter -- but these do not seem to do anything).
            !
            ! So, in the future this could potentially be added/interfaced with that and would
            ! (could?) make this reproducible (by fixing the seed)... for now, this leads to
            ! non bitwise reproducible behaviour:
            ! ------------------------------------------------------------------------------- !
            CALL RANDOM_NUMBER( zphi(:) )     ! random values, uniform distribution on [0, 1  )
            zphi(:) = zphi(:) * 2._wp * rpi   !                         scale to be on [0, 2pi)
            ! ------------------------------------------------------------------------------- !
         ELSE
            ! Phase is an arbitrary constant, same for all spectral classes:
            zphi(:) = rpi
         ENDIF

         ! Calculate sea surface height along 1D subdomain (x1d):
         zssh(:) = 0._wp
         DO jf = 1, nn_nwfreq
            zssh(:) = zssh(:) + SQRT( 2._wp * pWspec(jf) * wdfreq(jf) ) * COS( zphi(jf) + wknum(jf) * x1d(:) )
         ENDDO

         ! Find local extrema in sea surface height, defined to be minima or maxima over
         ! a 'moving window' of (2*nn_ht15_rmin + 1) points in the 1D subdomain x1d:
         !
         DO jx = 1 + nn_ht15_rmin, nn_ht15_nx1d - nn_ht15_rmin
            llmax(jx) = ( MAXLOC( zssh(jx-nn_ht15_rmin:jx+nn_ht15_rmin), DIM=1 ) == (nn_ht15_rmin + 1) )
            llmin(jx) = ( MINLOC( zssh(jx-nn_ht15_rmin:jx+nn_ht15_rmin), DIM=1 ) == (nn_ht15_rmin + 1) )
            llext(jx) = (llmin(jx) .OR. llmax(jx))
         ENDDO

         ! Loop over all points again to identify series of three consecutive, alternating
         ! extrema {min., max., min.} or {max., min., max.}, from which calculate strain
         ! and hence determine whether ice fractures there or not.
         !
         ! Loop start/end indices correspond to first/last possible index that could
         ! possibly be at the centre of a triplet.
         !
         DO jx = 2 + nn_ht15_rmin, nn_ht15_nx1d - nn_ht15_rmin - 1
            !
            ! Reset values for next loop iteration. Note: re-using local integer variables
            ! ixlo and ixhi from above; now they are the indices of x1d corresponding to the
            ! nearest extrema on either side of the current extrema being considered):
            !
            ixlo = 0
            ixhi = 0
            !
            IF( llext(jx) ) THEN
               !
               ! Identify nearest extrema on the left [such that x1d(ixlo) < x1d(jx)]:
               !
               DO jy = jx-1, 1, -1
                  IF( llext(jy) ) THEN
                     ixlo = jy
                     EXIT
                  ENDIF
               ENDDO
               !
               ! Identify nearest extrema on the right [such that x1d(jx) < x1d(ixhi)]:
               !
               DO jy = jx+1, nn_ht15_nx1d
                  IF( llext(jy) ) THEN
                     ixhi = jy
                     EXIT
                  ENDIF
               ENDDO
               !
               ! If we have a series of three extrema, with the central one being current jx, then both
               ! ixlo and and ixhi will have changed from 0. If they are alternating {max., min., max.}
               ! or {min., max., min.}, then calculate strain at x1d(jx) and determine if ice fractures
               ! there. If it does, append x1d(jx) to zxfrac array and increment ixfrac.
               !
               IF( (ixlo > 0) .AND. (ixhi > 0) ) THEN
                  !
                  IF(       ( llmax(ixlo) .AND. llmin(jx) .AND. llmax(ixhi) )   &
                     & .OR. ( llmin(ixlo) .AND. llmax(jx) .AND. llmin(ixhi) )    ) THEN
                     !
                     ! Calculate second derivative of SSH w.r.t. x at index jx
                     !
                     ! Centred finite difference for second derivative, forward and backward differences
                     ! on the first/inner derivatives, using the points x1d(ixlo) < x1d(jx) < x1d(ixhi).
                     ! Some simplifying algebra results in the calculation below, where also multiplying
                     ! by half of the mean ice thickness gives the strain at x1d(jx).
                     !
                     zdxlo = x1d(jx  ) - x1d(ixlo)
                     zdx   = x1d(ixhi) - x1d(ixlo)
                     zdxhi = x1d(ixhi) - x1d(jx  )
                     !
                     ! Note: zdx* are all strictly > 0, since ixlo <= jx - 1 and ixhi >= jx + 1
                     !
                     zstrain = ABS( .5_wp * ph_i * ( zssh(ixhi) * zdxlo - zssh(jx) * zdx + zssh(ixlo) * zdxhi)   &
                        &                          / ( zdxlo * zdx * zdxhi ) )
                     !
                     ! Only need to know whether this strain exceeds the critical strain
                     ! If it does, save it as a fracture point in array zxfrac:
                     IF( zstrain >= rn_ice_wav_ecri ) THEN
                        zxfrac(ixfrac) = x1d(jx)
                        ixfrac = ixfrac + 1
                     ENDIF
                     !
                  ENDIF ! pssh(jx)  is at the centre of a triplet of alternating min/max
               ENDIF ! -- pssh(jx)  is at the centre of a triplet of extrema
            ENDIF ! ----- llext(jx) [pssh(jx) is an extrema]
         ENDDO ! -------- jx        [loop of x1d points]

         ! Now have locations of strain points, zxfrac(1:ixfrac). The distances between such points, when
         ! converted to radii and binned into floe size categories, gives the fracture histogram.
         !
         ! In loop above, index ixfrac is used to populate zxfrac(:), so now, ixfrac - 1 = number of
         ! fracture points. So, if:
         !
         !    ixfrac == 1    ==> no fractures at all
         !    ixfrac == 2    ==> 1 fracture point
         !    ixfrac >= 3    ==> at least 2 fracture points
         !
         ! We only compute fracture lengths between fracture points; the end points, x = x1d(1) = 0 and
         ! x = x1d(nn_ice_wav_nx1d), do not count as fracture points (because we can never calculate the
         ! strain there). Therefore, only proceed if ixfrac is at least 3.
         !
         IF( ixfrac >= 3 ) THEN
            DO jx = 2, ixfrac - 1
               !
               zfrac_s = zxfrac(jx) - zxfrac(jx-1)   ! size (diameter) of fractured ice
               !
               ! Populate appropriate floe size category in iWfrac(:)
               !
               ! Note that iWfrac is NOT reset to 0 between outer while loop iterations, so that it
               ! can be updated until 'convergence' if using the random phase option
               !
               DO jf = 1, nn_nfsd - 1
                  IF( zfrac_s < floe_su(jf) ) THEN
                     iWfrac(jf) = iWfrac(jf) + 1
                     EXIT
                  ENDIF
               ENDDO
               !
               ! Separate check for largest fractures (even if it exceeds upper bound of largest
               ! floe size category, it goes into that category anyway; note similar for very small
               ! fractures accounted for in above loop anyway):
               IF( zfrac_s >= floe_sl(nn_nfsd) ) iWfrac(nn_nfsd) = iWfrac(nn_nfsd) + 1
               !
            ENDDO   ! jx           [loop of fracture points]
         ENDIF   ! -- jxfrac >= 3  [at least 2 fracture points]

         ! Compute sW(s)ds and normalise w.r.t. 1D subdomain size (i.e., divide by D <==> zsWfrac(:)
         !
         ! (note: 'normalisation' necessary at this stage even though W(s) is not a normalised
         ! distribution in the theory, for convergence checks)
         !
         ! For one iteration, zsWfrac is normalised by D. However, iWfrac is accumulating over
         ! multiple iterations -- effectively, multiple 'instances' of the 1D sub-domain. So, we
         ! need to normalise by D *multiplied* by the iteration count (hence why jiter needs to
         ! start at 0 then be incremented by 1 at the top of the while loop):
         !
         zsWfrac(:) = floe_sc(:) *  REAL(iWfrac(:), KIND=wp)   &
            &                    / (REAL(jiter    , KIND=wp) * (x1d(nn_ht15_nx1d) - x1d(1)))

         ! Fractional difference of current zsWfrac from previous iteration
         ! Breaks while loop if/when this is below tolerance value for all categories
         !
         ! First set to 0 if both before (zsWfrac_b) and current zsWfrac are ~0; otherwise, 1
         ! Then in latter case calculate fractional difference as along as zrWfrac_b /= 0
         !
         WHERE( (zsWfrac_b <= epsi10) .AND. (zsWfrac <= epsi10) )   ;   zsWfrac_err = 0._wp
         ELSEWHERE                                                  ;   zsWfrac_err = 1._wp
         ENDWHERE
         !
         WHERE( zsWfrac_b > epsi10 )   zsWfrac_err = ABS( zsWfrac - zsWfrac_b ) / zsWfrac_b

         zsWfrac_b(:) = zsWfrac(:)   ! save for next iteration

      ENDDO   ! while loop

      ! Add to counter if convergence not reached. This is (optionally) reported in ice_wav_frac:
      IF( ln_ht15_rand .AND. (jiter == iloop) )   pstat = pstat + 1

      ! Calculate the probability (pQfrac) and redistribution (pBfrac) functions
      ! from the fracture distribution [zsWfrac, which corresponds to sW(s)ds/D]:
      DO jf = 1, nn_nfsd
         !
         ! Q(s) = int[ xW(x)dx ] / D   for x < s:
         !
         ! Factor of 1/D already accounted for in calculation of zsWfrac
         !
         ! Floes can fracture into same category, so need to include up to jf, but weight it by 0.5
         ! to account that not all fracture sizes represented by sW(s)ds could have resulted from
         ! fractures of floes in the same range [s,s+ds] (this is assuming a uniform distribution of
         ! initial/fractured floes sub-category; unlike schemes Z16 and Y24* it is not possible to
         ! do anything more accurately here)
         !
         ! Normalise also by time step (units -> s-1):
         !
         pQfrac(jf) = ( SUM(zsWfrac(1:jf-1)) + .5_wp * zsWfrac(jf) ) * r1_Dt_ice

         ! B(s',s)ds = sW(s)ds / int[ xW(x)dx ] for s < s', x < s'
         ! [zsWfrac includes factor 1/D, but this cancels in this expression
         ! and normalisation below, so doesn't matter]:
         !
         pBfrac(jf,1:jf-1) =         zsWfrac(1:jf-1)
         pBfrac(jf,jf    ) = .5_wp * zsWfrac(jf    )   ! again, account for floes fracturing to same cat.

         ! Normalise B (effectively account for denominator in equation):
         IF( SUM(pBfrac(jf,:)) > 0._wp ) pBfrac(jf,:) = pBfrac(jf,:) / SUM(pBfrac(jf,:))
         !
      ENDDO

      ! Control:
      IF( ln_timing )   CALL timing_stop('wav_frac_ht15')

   END SUBROUTINE wav_frac_ht15


   SUBROUTINE wav_spec_bret( phsw_l, pwpf_l, pwspec_l )
      !!-------------------------------------------------------------------
      !!                *** ROUTINE wav_spec_bret ***
      !!-------------------------------------------------------------------
      !!
      !! ** Purpose :   Estimate local wave energy spectrum from local wave properties
      !!                calculated as Bretschneider spectrum
      !!
      !! ** Method  :   SB(f) = (5/16) * Hs^2 * (fp^4 / f^5) * EXP[ -(5/4)*(fp/f)^4 ]
      !!
      !!                where Hs is significant wave height (m), fp is the frequency
      !!                of peak wave energy (Hz) and f is frequency (Hz). SB(f) is the
      !!                Bretschneider wave energy spectrum (a.k.a., power spectral density)
      !!                and has units of m2.Hz-1 (i.e., m2.s).
      !!
      !! ** Inputs  :   phsw_l  :   local significant wave height (m)
      !!                pwpf_l  :   local peak frequency (Hz)
      !!
      !! ** Outputs :   pwspec_l:   local wave energy spectrum (Bretschneider; m2.Hz-1)
      !!
      !! ** Notes   :   The approach of using the Bretschneider formula to estimate the wave spectrum
      !!                for wave-ice interactions follows Horvat and Tziperman (2015) and Roach et al. (2018)
      !!                although the functional form here differs. The above matches Williams et al. (2013, Eq. 21)
      !!                and Bateson et al. (2020, Eq. 5), just converting from their expressions for SB(w) in terms
      !!                of angular frequency, w = 2*pi*f, using SB(w)dw = SB(w)(dw/df)df = 2*pi*SB(w)df,
      !!                hence SB(f) = 2*pi*SB(w). This form can be traced back to Bretschneider (1959).
      !!
      !! ** Callers :   ice_wav_frac  --> [wav_spec_bret]   (ln_ice_wav_spec=F AND ln_ice_wav_attn=F)
      !!                wav_attn_spec --> [wav_spec_bret]   (ln_ice_wav_spec=F AND ln_ice_wav_attn=T)
      !!
      !! ** References
      !!    ----------
      !!    Bateson, A. W., Feltham, D. L., Schroeder, D., Hosekova, L., Ridley, J. K., & Aksenov, Y. (2020).
      !!              Impact of sea ice floe size distribution on seasonal fragmentation and melt of Arctic sea ice.
      !!              The Cryosphere, 14, 403-428.
      !!    Bretschneider, C. L. (1959).
      !!              Wave variability and wave spectra for wind-generated gravity waves.
      !!              Technical Memorandum No. 118, Beach Erosion Board, U.S. Army Corps of Engineers, Washington, DC, USA.
      !!    Horvat, C., & Tziperman, E. (2015).
      !!              A prognostic model of the sea-ice floe size and thickness distribution.
      !!              The Cryosphere, 9, 2119-2134.
      !!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018).
      !!              An emergent sea ice floe size disribution in a global coupled ocean-sea ice model.
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!    Williams, T. D., Bennetts, L. G., Squire, V. A., Dumon, D. & Bertino, L. (2013).
      !!              Wave-ice interactions in the marginal ice zone. Part 1: Theoretical foundations.
      !!              Ocean Modelling, 71, 81-91.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp)                      , INTENT(in)    ::   phsw_l     ! local significant wave height (m)
      REAL(wp)                      , INTENT(in)    ::   pwpf_l     ! local peak frequency (Hz)
      REAL(wp), DIMENSION(nn_nwfreq), INTENT(inout) ::   pwspec_l   ! local wave energy spectrum (m2.Hz-1)
      !
      !!-------------------------------------------------------------------

      IF( (phsw_l >= minhsw) .AND. (pwpf_l >= minwpf) ) THEN
         pwspec_l(:) = .3125_wp * phsw_l**2 * pwpf_l**4 * EXP( -1.25_wp * ( pwpf_l / wfreq(:) )**4 ) / wfreq(:)**5
      ELSE
         pwspec_l(:) = 0._wp
      ENDIF

   END SUBROUTINE wav_spec_bret


   SUBROUTINE wav_spec_rayl( pwmp_l, prayl_l )
      !!-------------------------------------------------------------------
      !!                *** ROUTINE wav_spec_rayl ***
      !!-------------------------------------------------------------------
      !!
      !! ** Purpose :   Calculate local Rayleigh spectrum P(f)
      !!
      !! ** Method  :   P(f) = 2.7 * exp[ -0.675 / (T_m * f)^4 ] / (T_m^4 * f^5)
      !!
      !!                where f is frequency and T_m is mean wave period. The equation
      !!                is normalised so that that the integral of P over all f is 1.
      !!
      !! ** Inputs  :   pwmp_l  : local wave mean period
      !! ** Outputs :   prayl_l : local Rayleigh spectrum of frequencies, P(f) (Hz-1)
      !!
      !! ** Notes   :   The Rayleigh spectrum is usually expressed in terms of period, but wave-ice
      !!                fracture is implemented using frequency for all spectra so here it has been
      !!                converted from P(T)dT -> P(f)df. See Bretschneider (1959, Eq. 3.35) for P(T)
      !!                and then the above follows using T = 1/f. For now this is only needed by the
      !!                Yang et al. (2024) 'semi-empirical' fracture scheme (jpfrac_y24b) so this
      !!                saves creating yet another, separate set of coordinate arrays for wave period
      !!                and spectral class widths in terms of period just for this one calculation.
      !!
      !!                Note that in Yang et al. (2024), their Eq. (19) for P(T) is incorrect (suspect
      !!                they used Bretschneider's Eq. 3.34 by mistake, which expresses the Rayleigh
      !!                spectrum in standard form, i.e., in terms of tau = T/T_m).
      !!
      !! ** Callers :   wav_frac_y24b  --> [wav_spec_rayl]   (nn_frac_scheme == jpfrac_y24b)
      !!
      !! ** References
      !!    ----------
      !!    Bretschneider, C. L. (1959).
      !!              Wave variability and wave spectra for wind-generated gravity waves.
      !!              Technical Memorandum No. 118, Beach Erosion Board, U.S. Army Corps of Engineers, Washington, DC, USA.
      !!    Yang, C.-Y., Liu, J., & Chen, D. (2024).
      !!              Understanding the influence of ocean waves on Arctic sea ice simulation:
      !!              a modelling study with an atmosphere-ocean-wave-sea ice coupled model.
      !!              The Cryosphere, 18(3), 1215-1239.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp)                      , INTENT(in)    ::   pwmp_l    ! local wave mean period (s)
      REAL(wp), DIMENSION(nn_nwfreq), INTENT(inout) ::   prayl_l   ! local Rayleigh spectrum P(f) (Hz-1)
      !
      !!-------------------------------------------------------------------

      IF( pwmp_l >= minwmp ) THEN
         prayl_l(:) = 2.7_wp * EXP( -.675_wp / (pwmp_l * wfreq(:))**4 ) / (pwmp_l**4 * wfreq(:)**5)

         ! Normalise (equation is normalised in theory, but discretisation leads to errors):
         IF( SUM(prayl_l(:) * wdfreq(:) ) > 0._wp ) THEN
            prayl_l(:) = prayl_l(:) / SUM( prayl_l(:) * wdfreq(:) )
         ENDIF
      ELSE
         prayl_l(:) = 0._wp
      ENDIF

   END SUBROUTINE wav_spec_rayl


   SUBROUTINE wav_calc_stmer
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE wav_calc_stmer ***
      !!
      !! ** Purpose :   Calculate meridional distances across T cells
      !!
      !! ** Method  :   This calculation needs the sine and cosine of the angle between
      !!                meridians and the j-direction on the T-grid. There is already a
      !!                routine called 'angle' in module geo2ocean.F90 which does this
      !!                and related calculations, but it is not guaranteed to be called
      !!                in all configurations. So the first step of this routine is to
      !!                copy the calculations for the variables gsint and gcost from
      !!                the 'angle' routine. These are then used to compute meridional
      !!                distances across T cells, saving values into the module variable
      !!                stmer which is required by the wave attenuation scheme
      !!                (subroutine ice_wav_attn).
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp) ::   zxnpt, zynpt, znnpt   ! x,y components and norm of the vector: T point to North Pole
      REAL(wp) ::   zxvvt, zyvvt, znvvt   ! x,y components and norm of the vector: between V points below and above a T point
      REAL(wp) ::   zsint, zcost          ! sine and cosine of grid angle at T points
      REAL(wp) ::   zcoststar             ! cosine of T grid angle at which opposite corners (F points) are aligned along a meridian
      INTEGER  ::   ierr                  ! allocate status return value
      INTEGER  ::   ji, jj                ! dummy loop indices
      !
      !!-------------------------------------------------------------------

      ALLOCATE( stmer(jpi,jpj), STAT=ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'wav_calc_stmer: unable to allocate stmer' )

      stmer(:,:) = 0._wp

      DO_2D(0, 1, 0, 1)
         !
         ! Calculate cosine/sine of angle between the northward and grid j-direction at
         ! T points; see subroutine 'angle' of geo2ocean.F90, where this is taken from:
         !
         IF( MOD( ABS( glamv(ji,jj) - glamv(ji,jj-1) ), 360._wp) < 1.e-8_wp ) THEN
            zsint = 0._wp
            zcost = 1._wp
         ELSE
            zxnpt = 0._wp - 2._wp * COS( rad * glamt(ji,jj) ) * TAN( rpi / 4._wp - rad * gphit(ji,jj) / 2._wp )
            zynpt = 0._wp - 2._wp * SIN( rad * glamt(ji,jj) ) * TAN( rpi / 4._wp - rad * gphit(ji,jj) / 2._wp )
            !
            znnpt = zxnpt * zxnpt + zynpt * zynpt
            !
            zxvvt =  2._wp * COS( rad * glamv(ji,jj)   ) * TAN( rpi / 4._wp - rad * gphiv(ji,jj)   / 2._wp )   &
               &  -  2._wp * COS( rad * glamv(ji,jj-1) ) * TAN( rpi / 4._wp - rad * gphiv(ji,jj-1) / 2._wp )
            !
            zyvvt =  2._wp * SIN( rad * glamv(ji,jj)   ) * TAN( rpi / 4._wp - rad * gphiv(ji,jj)   / 2._wp )   &
               &  -  2._wp * SIN( rad * glamv(ji,jj-1) ) * TAN( rpi / 4._wp - rad * gphiv(ji,jj-1) / 2._wp )
            !
            znvvt = MAX( SQRT( znnpt * ( zxvvt * zxvvt + zyvvt * zyvvt ) ), 1.e-14_wp )
            !
            zsint = ( zxnpt * zyvvt - zynpt * zxvvt ) / znvvt
            zcost = ( zxnpt * zxvvt + zynpt * zyvvt ) / znvvt
            !
         ENDIF
         !
         ! Cosine of angle defining aspect ratio of grid cell: this is the value of zcost when
         ! the T grid cell opposite corners (F points) pass through the same meridian. This
         ! determines which component, e1t (i-direction) or e2t (j-direction) is used to work
         ! out the meridional distance along the T grid cell:
         !
         zcoststar = e2t(ji,jj) / SQRT( e1t(ji,jj)**2 + e2t(ji,jj)**2 )
         !
         ! Note: for ORCA configuration zcost > 0 always, but the below
         ! ----- accounts for all possible grid cell orientations
         !
         IF( ABS(zcost) >= zcoststar ) THEN
            stmer(ji,jj) = ABS( e2t(ji,jj) / zcost )
         ELSE
            stmer(ji,jj) = ABS( e1t(ji,jj) / zsint )
         ENDIF
         !
      END_2D

     CALL lbc_lnk( 'wav_calc_stmer', stmer, 'T', 1._wp )

   END SUBROUTINE wav_calc_stmer


   SUBROUTINE ice_wav_init
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE ice_wav_init  ***
      !!
      !! ** Purpose :   Initialise ice wave impacts module.
      !!
      !! ** Method  :   Namelist read.
      !!                Check flags suitably set and stop if not.
      !!                Calculate some module constants.
      !!
      !! ** Callers :   ice_init --> [ice_wav_init]
      !!
      !! ** Note    :   Must be called after ice_fsd_init so that FSD-related
      !!                namelist flags are read and set. Note that some flag
      !!                checking cannot be done until after namelist group
      !!                namsbc_wave is read, which occurs later. Such checks
      !!                are handled in subroutine sbc_wave_init directly.
      !!
      !!-------------------------------------------------------------------
      !
      INTEGER  ::   ji, jj, jf, ji_glo, jj_glo         ! Dummy loop indices
      INTEGER  ::   ierr                               ! Local integer output status for allocate
      INTEGER  ::   ios, ioptio                        ! Local integer output status for namelist read
      !
      !!
      NAMELIST/namwav/ ln_ice_wav     , ln_ice_wav_spec, rn_ice_wav_ecri, nn_frac_scheme,   &
         &             ln_z16_const   , rn_z16_cb      , rn_z16_k       , rn_z16_a      ,   &
         &             rn_z16_b       , rn_z16_hc      , nn_z16_nf      ,                   &
         &             rn_y24a_cw     , rn_y24a_alpha  ,                                    &
         &             nn_ht15_nx1d   , rn_ht15_dx1d   , nn_ht15_rmin   , ln_ht15_rand  ,   &
         &             rn_ht15_tol    , nn_ht15_maxiter, ln_ht15_warn   ,                   &
         &             ln_ice_wav_attn, rn_attn_lam_tol, rn_attn_c0     , rn_attn_ch    ,   &
         &             rn_attn_ch     , rn_attn_ct     , rn_attn_ch2    , rn_attn_ct2   ,   &
         &             rn_attn_cht    , rn_attn_tun
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice, namwav)
      READ_NML_CFG(numnam_ice, namwav)
      IF(lwm) WRITE(numoni, namwav)
      !
      IF(lwp) THEN   ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_wav_init: parameters for wave-ice interactions'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namwav:'
         WRITE(numout,*) '      Wave-ice interactions active or not                        ln_ice_wav = ', ln_ice_wav
         WRITE(numout,*) '         Activate wave-in-ice attenuation scheme or not     ln_ice_wav_attn = ', ln_ice_wav_attn
         WRITE(numout,*) '            Longitude tolerance for meridians (deg E)       rn_attn_lam_tol = ', rn_attn_lam_tol
         WRITE(numout,*) '            Attenuation coefficient, factors in quadratic expression for ln[a(T,h)]:'
         WRITE(numout,*) '               Constant term                                     rn_attn_c0 = ', rn_attn_c0
         WRITE(numout,*) '               Factor of ice thickness (h)                       rn_attn_ch = ', rn_attn_ch
         WRITE(numout,*) '               Factor of wave period (T)                         rn_attn_ct = ', rn_attn_ct
         WRITE(numout,*) '               Factor of h^2                                    rn_attn_ch2 = ', rn_attn_ch2
         WRITE(numout,*) '               Factor of T^2                                    rn_attn_ct2 = ', rn_attn_ct2
         WRITE(numout,*) '               Factor of h*T                                    rn_attn_cht = ', rn_attn_cht
         WRITE(numout,*) '            Tuning factor on a(T,h)                             rn_attn_tun = ', rn_attn_tun
         WRITE(numout,*) '         Read full wave energy spectrum or not              ln_ice_wav_spec = ', ln_ice_wav_spec
         WRITE(numout,*) '         Critical strain at which ice breaks due to waves   rn_ice_wav_ecri = ', rn_ice_wav_ecri
         WRITE(numout,*) '         Wave fracture scheme selection                     nn_frac_scheme  = ', nn_frac_scheme
         WRITE(numout,'(A,I0,A)') '            Zhang et al. (2016) scheme parameters  (nn_frac_scheme = ', jpfrac_z16, '):'
         WRITE(numout,*) '               Use constant participation factor               ln_z16_const = ', ln_z16_const
         WRITE(numout,*) '                  If T, the value is                        rn_z16_cb = ', rn_z16_cb
         WRITE(numout,*) '                  If F, calculate it with:'
         WRITE(numout,*) '                     parameter k                               rn_z16_k = ', rn_z16_k
         WRITE(numout,*) '                     parameter a                               rn_z16_a = ', rn_z16_a
         WRITE(numout,*) '                     parameter b                               rn_z16_b = ', rn_z16_b
         WRITE(numout,*) '                     Cutoff ice thickness                     rn_z16_hc = ', rn_z16_hc
         WRITE(numout,*) '                     Num. surrounding cells for fetch param.  nn_z16_nf = ', nn_z16_nf
         WRITE(numout,'(A,I0,A)') '            Yang et al. (2024) A scheme parameters (nn_frac_scheme = ', jpfrac_y24a, '):'
         WRITE(numout,*) '               Probability function parameter c_w           rn_y24a_cw      = ', rn_y24a_cw
         WRITE(numout,*) '               Probability function parameter alpha         rn_y24a_alpha   = ', rn_y24a_alpha
         WRITE(numout,'(A,I0,A)') '            Horvant & Tziperman (2015) scheme parameters (nn_frac_scheme = ', jpfrac_ht15, '):'
         WRITE(numout,*) '               Size of 1D subdomain for SSH                    nn_ht15_nx1d = ', nn_ht15_nx1d
         WRITE(numout,*) '               Increment of 1D subdomain for (m)               rn_ht15_dx1d = ', rn_ht15_dx1d
         WRITE(numout,*) '               Radius of smallest floes affected by waves (dx1d) nn_ht15_rmin = ', nn_ht15_rmin
         WRITE(numout,*) '               Use random phases or not                        ln_ht15_rand = ', ln_ht15_rand
         WRITE(numout,*) '                  Tolerance for convergence (error fraction)    rn_ht15_tol = ', rn_ht15_tol
         WRITE(numout,*) '                  Max. iterations for convergence           nn_ht15_maxiter = ', nn_ht15_maxiter
         WRITE(numout,*) '                  Warn whenever convergence not reached        ln_ht15_warn = ', ln_ht15_warn
      ENDIF

      IF( ln_ice_wav ) THEN

         ! Checks on flags that do not require knowning SBC wave module flags
         ! (those are handled in subroutine sbc_wave_init).
         !
         ! Wave-ice interactions module requires both wave inputs and FSD
         ! (exception: fracture scheme is Z16, in which case only FSD is needed)
         IF( (nn_frac_scheme /= jpfrac_z16) .AND. .NOT. ln_wave ) THEN
            CALL ctl_stop('ice_wav_init: need to activate SBC wave module (ln_wave=T) for your choice of nn_frac_scheme')
         ENDIF
         IF( .NOT. ln_fsd  ) THEN
            CALL ctl_stop('ice_wav_init: wave-ice interactions require floe size distribution (ln_fsd=T)')
         ENDIF

         ! Warn if both attenuation scheme and reading of full wave spectrum selected
         ! (possible to do so, but usually spectrum comes from a coupled wave model so
         ! should not need to attenuate waves under ice as that is done in the wave model)
         !
         ! (check that spectrum is actually read in is done in sbc_wave_init)
         !
         IF( ln_ice_wav_attn ) THEN
            IF ( ln_ice_wav_spec ) THEN
               CALL ctl_warn('ice_wav_init: ln_ice_wav_attn=T but also using spectrum (ln_ice_wav_spec=T): intentional?')
            ENDIF
            IF ( nn_frac_scheme == jpfrac_z16 ) THEN
               ln_ice_wav_attn = .FALSE.
               CALL ctl_warn('ice_wav_init: ask for attenuation scheme but using Z16 fracture scheme',   &
                  &          '==>> we force ln_ice_wav_attn = F anyway: attenuation scheme is NOT used')
            ENDIF
         ENDIF

         ! Inform SBC wave module which fields need to be read for wave-ice interactions
         !
         ! See sbcwave: l_ini_* determines whether field is allocated/initialised
         !              l_get_* determines whether it is read (from file or coupled wave model)
         !
         ! Never set these flags to F, as they may be needed for other purposes that have already
         ! set them to T. If l_get_VAR=T then the corresponding l_ini_VAR need not be set as it
         ! will be set to T in sbc_wave_init anyway if l_get_VAR=T.
         !
         SELECT CASE( nn_frac_scheme )
            CASE( jpfrac_y24a )
               !
               l_get_hsig = .TRUE.
               l_get_wper = .TRUE.
               !
            CASE( jpfrac_y24b )
               !
               IF( ln_ice_wav_spec ) THEN   ;   l_get_wspec = .TRUE.
               ELSE                         ;   l_get_hsig  = .TRUE.   ;   l_get_wpf   = .TRUE.
               ENDIF
               !
               l_ini_wspec = .TRUE.   ! needed either way (though redundant for ln_ice_wav_spec=T case)
               !
               ! If attenuation scheme active, wmp/wper calculated from attenuated spectrum
               ! If not, need to read it directly:
               IF( .NOT. ln_ice_wav_attn )      l_get_wper  = .TRUE.
               !
            CASE( jpfrac_ht15 )
               !
               IF( ln_ice_wav_spec ) THEN   ;   l_get_wspec = .TRUE.
               ELSE                         ;   l_get_hsig  = .TRUE.   ;   l_get_wpf   = .TRUE.
               ENDIF
               l_ini_wspec = .TRUE.   ! needed either way (though redundant for ln_ice_wav_spec=T case)
               !
            CASE( jpfrac_z16  )
               ! => no wave fields needed; do not update SBC flags
            CASE DEFAULT
               ! Use this opportunity to stop as invalid value given for nn_frac_scheme:
               CALL ctl_stop('STOP', 'ice_wav_init: invalid option for nn_frac_scheme')
         END SELECT

         ! Wave-dependent growth: need wave peak frequency and sig. wave height:
         IF( nn_frac_scheme /= jpfrac_z16 ) THEN
            l_get_hsig  = .TRUE.   ;   l_get_wpf  = .TRUE.
         ENDIF

         ! Wave attenuation scheme: ensure required arrays are at least allocated, and rely on
         ! above to determine which fields are actually read in:
         IF( ln_ice_wav_attn ) THEN
            l_ini_wspec = .TRUE.   ;   l_ini_hsig = .TRUE.
            l_ini_wpf   = .TRUE.   ;   l_ini_wper = .TRUE.
         ENDIF

         ! Allocate and define module constants
         !
         ! Z16 scheme: dimension of local + surrounding domain for sea ice conc./fetch calculation:
         IF( (nn_frac_scheme == jpfrac_z16) .AND. ( nn_z16_nf < 0 ) ) THEN
            CALL ctl_stop('STOP', 'ice_wav_init: nn_z16_nf must be >= 0')
         ENDIF

         ! If fracture scheme uses constant redistributor function, calculate it now
         ! (see Zhang et al. 2015; JGR:O for theory):
         IF( (nn_frac_scheme == jpfrac_z16) .OR. (nn_frac_scheme == jpfrac_y24a) ) THEN
            !
            ALLOCATE( Bfrac_uni(nn_nfsd,nn_nfsd), STAT=ierr )
            !
            IF(ierr /= 0) CALL ctl_stop('ice_wav_init: could not allocate array: Bfrac_uni')
            !
            Bfrac_uni(:,:) = 0._wp
            !
            ! Note: we do not calculate this using formula given in Zhang et al. (2015, Eq. 14)
            ! with cutoff values explicitly because we just assume full range of fractured floe
            ! sizes are possible, which is what Z15, Z16, and Y24A schemes use anyway. Each
            ! Bfrac_uni(jf,:) represents integration of uniform (i.e., constant) B over floe size
            ! category ranges of the second index, so just needs to be scaled in proportion to the
            ! category widths and then normalised at the end.
            !
            ! In general, B needs to be defined on the diagonal -- in principle, floes in category
            ! jf *can* fracture into the same category jf as it represents a range of floe sizes.
            ! In the uniform redistributor case this is obviously valid. The factor of 0.5 accounts
            ! for the fact that not all sub-category floe sizes can result from a given initial
            ! sub-category size in the same range (the 0.5 is intuitive but can also be shown
            ! analytically in this case).
            !
            ! Particularly important: Bfrac(1,1) must always equal 1 (uniform and general case), so
            ! that in the wave fracture equation the smallest category loss term, which can be non-
            ! zero if Q(1) /= 0, always cancels with the corresponding gain term.
            !
            DO jf = 1, nn_nfsd
               Bfrac_uni(jf,1:jf-1) =         floe_ds(1:jf-1)
               Bfrac_uni(jf,jf    ) = .5_wp * floe_ds(jf    )
               !
               ! Normalise so that integral of Bfrac_uni(s,s')*ds' = 1:
               IF( SUM(Bfrac_uni(jf,:)) > 0._wp )   &
                  &   Bfrac_uni(jf,:) = Bfrac_uni(jf,:) / SUM(Bfrac_uni(jf,:))
               !
            ENDDO
            !
         ENDIF

         ! Sub-gridscale domain (1D axis in direction of wave propagation) for
         ! computation of wave fracture distribution in subroutine wav_frac_ht15:
         IF( nn_frac_scheme == jpfrac_ht15 ) THEN
            !
            ALLOCATE( x1d(nn_ht15_nx1d), STAT=ierr )
            !
            IF (ierr /= 0) CALL ctl_stop('ice_wav_init: could not allocate array: x1d')
            !
            x1d(1) = 0._wp   ! x1d = 0., dx, 2*dx, ...
            !
            DO ji = 2, nn_ht15_nx1d
               x1d(ji) = x1d(ji-1) + rn_ht15_dx1d
            ENDDO
         ENDIF

         ! If using wave attenuation scheme, allocate/prepare the global coordinate arrays.
         IF( ln_ice_wav_attn ) THEN
            !
            ALLOCATE( glamt_glo(jpiglo,jpjglo) , gphit_glo(jpiglo,jpjglo) , STAT=ierr )
            !
            IF( ierr /= 0 ) CALL ctl_stop('ice_wav_init: could not allocate global domain coordinates')

            glamt_glo(:,:) = 0._wp
            gphit_glo(:,:) = 0._wp

            ! Calculate global arrays of longitude/latitude:
            DO_2D(0, 0, 0, 0)
               ji_glo = mig(ji,nn_hls)   ! local --> global index
               jj_glo = mjg(jj,nn_hls)   ! local --> global index
               glamt_glo(ji_glo,jj_glo) = glamt(ji,jj)
               gphit_glo(ji_glo,jj_glo) = gphit(ji,jj)
            END_2D
            !
            CALL wav_merge_glo( glamt_glo )   ! merge global arrays to/from all processors
            CALL wav_merge_glo( gphit_glo )   !
            !
            CALL wav_calc_stmer               ! calculate meridional distances across T cells (variable stmer)
            !
         ENDIF

         ! Constant flag for subroutine ice_wav_frac: conditions under which it needs to
         ! compute local wave spectrum (T), otherwise it is already available in wspec (F):
         !
         ! Key point is that attenuation scheme (ln_ice_wav_attn) calculates attenuated wave spectrum
         ! So if NOT reading spectrum from file/model, need to calculate spectrum only if that scheme
         ! is also NOT activated, and this is all only needed for certain fracture schemes.
         !
         l_frac_calc_spec =       (.NOT. ln_ice_wav_spec)           &
            &               .AND. (.NOT. ln_ice_wav_attn)           &
            &               .AND. (nn_frac_scheme == jpfrac_y24b .OR. nn_frac_scheme == jpfrac_ht15)

         ! Similar for routine ice_wav_attn (attenuation scheme): it will need to compute the wave
         ! spectrum everywhere (regardless of ice presence) only if the wave spectrum is NOT read in:
         l_attn_calc_spec = ln_ice_wav_attn .AND. .NOT. ln_ice_wav_spec

         IF(lwp) THEN   ! this was mainly for testing -- not useful to user so could be removed?
            WRITE(numout,*) ''
            WRITE(numout,*) '   Namelist options ==> wave fracture    scheme will calculate spectrum = ', l_frac_calc_spec
            WRITE(numout,*) '                    ==> wave attenuation scheme will calculate spectrum = ', l_attn_calc_spec
         ENDIF

      ELSE
         ! ln_ice_wav = F; force some parameters to be F to avoid issues:
         l_frac_calc_spec = .FALSE.
         l_attn_calc_spec = .FALSE.
         ln_ice_wav_spec  = .FALSE.
         ln_ice_wav_attn  = .FALSE.
      ENDIF

   END SUBROUTINE ice_wav_init


#else
   !!----------------------------------------------------------------------
   !!   Default option          Empty module          NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icewav
