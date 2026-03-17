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
   USE sbc_oce, ONLY :   wndm, ln_wave, ln_wave_spec, nn_nwfreq   ! SBC module
   USE sbcwave, ONLY :   hsw, wpf, wmp, wfreq, wdfreq, wspec      ! SBC: wave variables
   USE ice               ! sea-ice: variables
   USE icefsd , ONLY :   a_ifsd, nf_newice, floe_rl, floe_rc, floe_ru, floe_dr   ! floe size distribution parameters/variables
   USE icefsd , ONLY :   rDt_ice_fsd, fsd_cleanup, ice_fsd_dia                   ! floe size distribution functions/routines

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
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)   ::   wknum       ! angular wave numbers corresponding to wfreq array (m-1)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:)   ::   wlam        ! wavelengths corresponding to wfreq array (m)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   Bfrac_uni   ! uniform fracture redistributor (Y24a scheme)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   stmer       ! meridional distance across T cells (m)

   INTEGER, PARAMETER ::   jpfrac_z16  = 1   ! option for nn_frac_scheme -> Zhang et al. (2016) scheme
   INTEGER, PARAMETER ::   jpfrac_y24a = 2   ! option for nn_frac_scheme -> Yang et al. (2024) scheme A
   INTEGER, PARAMETER ::   jpfrac_y24b = 3   ! option for nn_frac_scheme -> Yang et al. (2024) scheme B
   INTEGER, PARAMETER ::   jpfrac_ht15 = 4   ! option for nn_frac_scheme -> Horvat and Tziperman (2015) scheme

   LOGICAL ::   l_attn_calc_spec   ! whether spectrum needs to be calculated in subroutine ice_wav_attn
   LOGICAL ::   l_frac_calc_spec   ! whether spectrum needs to be calculated in subroutine ice_wav_frac

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
   REAL(wp)        ::   rn_y24a_cw        !: Parameter in Y24A fracture scheme
   REAL(wp)        ::   rn_y24a_alpha     !: Parameter in Y24A fracture scheme
   INTEGER         ::   nn_ht15_nx1d      !: Size of 1D subdomain for wave fracture calculation (HT15 only)
   REAL(wp)        ::   rn_ht15_dx1d      !: Increment of 1D subdomain for wave fracture calculation (m; HT15 only)
   INTEGER         ::   nn_ht15_rmin      !: Radius of smallest floes affected by wave fracture in units of rn_ht15_dx1d
   LOGICAL         ::   ln_ht15_rand      !: Use random phases for SSH in HT15 wave fracture calculation
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
      !!                   2*r_max = sqrt[ (2*C2*Lp^2) / (Wa*g*rhoi*pi^3) ]
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
      !!                   r_max = [1/(2*pi^2*fp^2)] * sqrt[ (C2*g) / (pi*rhoi*Hs) ]
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
      REAL(wp)            ::   zrmax           ! maximum new ice floe size (m)
      INTEGER             ::   jf              ! dummy loop index
      !
      !!-------------------------------------------------------------------

      IF( (nn_frac_scheme == jpfrac_z16) .OR. (phsw < epsi06) .OR. (pwpf < epsi06) ) THEN
         !
         kcat = nf_newice   ! no waves present => set to default new floe size category
         !
      ELSE
         !
         ! --- Calculate maximum floe size from wave conditions, r_max:
         zrmax = SQRT( zC2 * grav / (rpi * rhoi * phsw) ) / (2._wp * rpi**2 * pwpf**2)
         !
         ! --- Find FSD category that r_max belongs to:
         !
         IF( zrmax < floe_rl(1) ) THEN
            ! Smaller than smallest 'resolved' floe size => put in smallest cat. anyway:
            kcat = 1
            !
         ELSE
            ! Find FSD category that r_max belongs to. Note that if r_max exceeds upper
            ! limit of largest floe size category, it goes into that category anyway:
            !
            DO jf = nn_nfsd, 1, -1
               IF( zrmax > floe_rl(jf) ) THEN
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
      !!                     wave spectrum from significant wave height/peak frequency inputs in all
      !!                     grid cells.
      !!
      !!                2.   Copy these terms into global-domain, global (i.e., module scope) arrays. This
      !!                     is necessary because for any given target grid cell the source grid cell may
      !!                     be in a different computational subdomain, so all values must be available to
      !!                     all subdomains before step 3. A call to lib_mpp routine mppsync (a wrapper for
      !!                     MPP_BARRIER) ensures this is possible. Similar global arrays are required for
      !!                     the grid cell longitude and latitudes, but these are compute in advance, in
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
      !!                                           --> ice_wav_calc
      !! ** Invokes :               [ice_wav_attn] --> wav_spec_bret()
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
                  zdmean = zdmean + 2._wp * floe_rc(jf) * a_i(ji,jj,jl) * a_ifsd(ji,jj,jf,jl)
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

         ELSEIF( tmask(ji,jj,1) < 1._wp ) THEN
            zattxp(:) = zf_land    ! set land flag
         ELSE
            zattxp(:) = zf_noice   ! set no-ice flag
         ENDIF

         ! Calculate wave spectrum (*ALL* grid cells), if needed:
         IF( l_attn_calc_spec ) wspec(ji,jj,:) = wav_spec_bret( hsw(ji,jj), wpf(ji,jj) )

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
               ! =========================================================== !
               ! 4  Calculate attenuated wave properties at target grid cell !
               ! =========================================================== !
               !
               CALL ice_wav_calc( wspec(ji,jj,:), hsw(ji,jj), wpf(ji,jj), wmp(ji,jj) )
               !
            ENDIF
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

      ! Moments of the wave spectrum:
      zm0 = SUM(            pwspec(:) * wdfreq(:) )
      zm1 = SUM( wfreq(:) * pwspec(:) * wdfreq(:) )

      ! Significant wave height:
      phsw = 4._wp * SQRT( zm0 )

      ! Wave peak frequency:
      imax = MAXLOC(pwspec(:), DIM=1)
      pwpf = wfreq(imax)

      ! Wave mean period:
      IF( zm1 > 0._wp ) THEN
         pwmp = zm0 / zm1
      ELSE
         pwmp = 0._wp
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
      !!                2.   Calculate the source terms Q(r) and B(s,r) in the equation for the
      !!                     tendency of the floe size distribution, f(r), due to wave fracture:
      !!
      !!                        df(r)/dt = -Q(r)*f(r) + int [ B(s,r)*Q(s)*f(s) ds ]
      !!
      !!                     The first term on the left-hand side represents loss of floes of size r
      !!                     due to them fracturing into smaller floes, and the second term represents
      !!                     gain of floes of size r due to fracturing of floes of size s.
      !!
      !!                     Q(r) is the probability that floes of size r are fractured by waves, per
      !!                     unit time, and B(s,r) is a redistribution function quantifying how floes
      !!                     of size s are transferred to size r; specifically, B(s,r)dr is the fraction
      !!                     of the original area of floes in the size interval [s, s+ds] transferred
      !!                     into the size interval [r,r+dr]. Hence, B(s,r) is normalised for each s
      !!                     such that: int[ B(s,r)dr ] = 1.
      !!
      !!                     This general formulation from Zhang et al. (2015) allows different wave
      !!                     fracture schemes to be used via different choices of Q(r) and B(s,r).
      !!                     These are done by the local subroutines wav_frac_*
      !!
      !!                     In the discretised implementation, f(r) is replaced by L(r) which
      !!                     corresponds to model variables a_ifsd(ji,jj,jf,jl) / floe_dr(jf), while
      !!                     f(s)ds in the integral is replaced directly with a_ifsd(ji,jj,:,jl).
      !!                     Hence there is an additional factor of floe_dr(jf) that is multiplied
      !!                     through and combined with the B(s,r) for which the local variable zBfrac
      !!                     then corresponds to B(s,r)dr (strictly speaking, B integrated over the
      !!                     floe size category interval, representing the total area fraction of
      !!                     fractured floes transferred into that category).
      !!
      !!                3.   Evolve the FSD using adaptive time stepping (Horvat and Tziperman, 2017).
      !!                     Additional time-step restrictions apply when evolving f(r), so a smaller step
      !!                     is (possibly) required, calculated in function rDt_ice_fsd() in module icefsd.
      !!
      !! ** Callers :   ice_stp -> [ice_wav_frac]
      !! ** Calls   :              [ice_wav_frac] -> wav_frac_{z16,y24a,y24b,ht15}
      !! ** Invokes :              [ice_wav_frac] -> wav_spec_bret()   (ln_ice_wav_spec = F AND ln_ice_wav_attn = F)
      !!                           [ice_wav_frac] -> rDt_ice_fsd()
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
      REAL(wp), DIMENSION(nn_nfsd,nn_nfsd)    ::   zBfrac               ! fracture redistribution function B(s,r)dr
      REAL(wp), DIMENSION(nn_nfsd)            ::   zQfrac               ! fracture probability function (s-1)
      REAL(wp), DIMENSION(nn_nfsd)            ::   za_ifsd_tend         ! tendency of FSD due to wave fracture
      REAL(wp)                                ::   zh_i                 ! mean ice thickness
      REAL(wp)                                ::   zfsd_res             ! correction term for area conservation
      REAL(wp)                                ::   zdt_sub              ! adaptive time step (s)
      REAL(wp)                                ::   ztelapsed            ! to track time elapsed during adaptive time stepping (s)
      INTEGER                                 ::   isubt                ! to track number of adaptive time steps
      INTEGER                                 ::   ji, jj, jl, jf       ! dummy loop indices
      !
      REAL(wp), PARAMETER                     ::   zat_i_min = .01_wp   ! minimum concentration for fracture to occur
      INTEGER , PARAMETER                     ::   isubt_max = 100      ! maximum number of adaptive time steps before warning
      !
      !!-------------------------------------------------------------------

      ! Control:
      IF( ln_timing )   CALL timing_start('ice_wav_frac')

      ! Calculate wavelength and angular wave number (2 * pi / wavelength) assuming dispersion
      ! relation of surface deep water gravity waves (these arrays are constants, so calculate once.
      ! Cannot do this in ice_wav_init as that is called before sbc_wave_init..)
      !
      ! ToDo: move this to sbc_wave_init?
      !
      IF( kt == nit000 ) THEN   ! at first time-step
         IF( nn_frac_scheme /= jpfrac_z16 ) THEN   ! frequencies undefined in this case!
            ALLOCATE( wlam(nn_nwfreq), wknum(nn_nwfreq) )
            wlam(:)  = grav / (2._wp * rpi * wfreq(:)**2)
            wknum(:) = 2._wp * rpi / wlam(:)
         ENDIF
      ENDIF

      za_ifsdb(A2D(0),:,:) = a_ifsd(A2D(0),:,:)   ! save a_ifsd before fracture for tendency diagnostics

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
            IF( l_frac_calc_spec ) wspec(ji,jj,:) = wav_spec_bret( hsw(ji,jj), wpf(ji,jj) )

            ! (2) Calculate source terms for the wave fracture equation
            !     This depends on the fracture scheme selected
            !
            ! Terms in common among > 1 schemes:
            zh_i = vt_i(ji,jj) / at_i(ji,jj)   ! mean ice thickness
            !
            SELECT CASE( nn_frac_scheme )
               CASE( jpfrac_z16  )
                  !
                  CALL wav_frac_z16( wndm(ji,jj)      , zh_i     , a_i(ji-1:ji+1,jj-1:jj+1,:),   &
                     &               a_ifsd(ji,jj,:,:), zQfrac(:), zBfrac(:,:) )
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
                  IF( MAXVAL( wspec(ji,jj,:) ) > epsi06 ) THEN
                     CALL wav_frac_ht15( wspec(ji,jj,:), zh_i, zQfrac(:), zBfrac(:,:) )
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
            IF( MAXVAL(zQfrac(:)) > epsi10 ) THEN
               !--------------------------------------!
               ! Begin sub-loop: thickness categories !
               !--------------------------------------!
               DO jl = 1, jpl
                  !
                  CALL fsd_cleanup( a_ifsd(ji,jj,:,jl) )   ! necessary?
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
                     ztelapsed = 0._wp
                     isubt     = 0
                     !
                     DO WHILE (ztelapsed < rDt_ice)
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
                        ! Compute adaptive timestep to increment FSD in each floe size
                        ! category, and make sure we do not overshoot actual time step:
                        zdt_sub = rDt_ice_fsd( a_ifsd(ji,jj,:,jl), za_ifsd_tend(:) )
                        zdt_sub = MIN(zdt_sub, rDt_ice - ztelapsed)
                        !
                        ! Update FSD and time elapsed:
                        a_ifsd(ji,jj,:,jl) = a_ifsd(ji,jj,:,jl) + zdt_sub * za_ifsd_tend(:)
                        ztelapsed          = ztelapsed + zdt_sub
                        isubt              = isubt + 1
                        !
                        IF( isubt == isubt_max ) THEN
                           CALL ctl_warn('ice_wav_frac not converging: ',               &
                              &          'reached maximum number of adaptive time steps')
                        ENDIF
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
                     CALL fsd_cleanup( a_ifsd(ji,jj,:,jl) )
                     !
                  ENDIF ! category jl can fracture
               ENDDO ! -- sub-loop (ice thickness categories)
            ENDIF ! ----- MAXVAL(zQfrac) > 0
         ENDIF ! -------- at_i(ji,jj) > zat_i_min
      END_2D ! ---------- main loop

      ! Write FSD tendency diagnostics due to wave fractue:
      CALL ice_fsd_dia( 'wav', za_ifsdb, a_ifsd(A2D(0),:,:), a_i(A2D(0),:), a_i(A2D(0),:) )

      ! Control:
      IF( ln_timing )   CALL timing_stop('ice_wav_frac')

   END SUBROUTINE ice_wav_frac


   SUBROUTINE wav_frac_z16( puatm, ph_i, pa_i, pa_ifsd, pQfrac, pBfrac )
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE wav_frac_z16 ***
      !!
      !! ** Purpose :   Calculate the probability and redistribution functions for the
      !!                equation evolving the FSD due to wave fracture for the scheme of
      !!                Zhang et al. (2016).
      !!
      !! ** Method  :   Fracture scheme that does *not* use wave forcing data (file or coupled model).
      !!                The probability of a floe of size r undergoing fracture, Q(r), is given by:
      !!
      !!                   Q(r) = MAX[ 0, 1 - int( f(r')dr' )/cb ]
      !!
      !!                where f(r') is the floe size distribution (integrated over thickness) and cb
      !!                is called a 'participation factor' representing the area fraction of ice that
      !!                could participate in fracture. In Zhang et al. (2015) this is set to a constant
      !!                and this behaviour can be achieved by setting namelist parameter ln_z16_const=T
      !!                and rn_z16_cb as required. In Zhang et al. (2016), this scheme was upgraded to
      !!                work in a GCM setting and cb is instead calculated from local wind speed U and
      !!                ice properties, as follows, which is the behaviour when ln_z16_const=F:
      !!
      !!                   c_b = [kU/MAX(hi,hc)] * EXP[ -a(1 - f0) - b(1 - ra/rmax) ] * dt
      !!
      !!                where k, a, and b are dimensionless constants (set in namelist), hi is mean ice
      !!                thickness, hc is a cutoff thickness (namelist), f0 is open water fraction,
      !!                ra is mean floe size, rmax is the largest resolved floe size, and dt is the model
      !!                timestep. The open water fraction is calculated as an average over the grid cell
      !!                and its eight surrounding neighbours.
      !!
      !!                The redistribution function is determined by assuming any floe of size s
      !!                that is undergoing wave fracture is equally likely to fracture into any
      !!                other floe size r < s. This 'unifom redistributor' is calculed in
      !!                subroutine ice_wav_init as it is a constant.
      !!
      !! ** Inputs  :   puatm                   :   local wind speed (m/s)
      !!                ph_i                    :   local mean ice thickness (m)
      !!                pa_i(3,3,jpl)           :   sea ice concentration at local and 8 neighbouring cells
      !!                pa_ifsd(nn_nfsd,jpl)    :   local areal floe size-thickness distribution
      !!
      !! ** Outputs :   pQfrac(nn_nfsd)         :   fracture probability function (s-1)
      !!                pBfrac(nn_nfsd,nn_nfsd) :   fracture redistribution function, B(s,r)dr
      !!                                            (note: first  index corresponds to original floe size s,
      !!                                                   second index corresponds to fractured floe size r)
      !!
      !! ** Callers :   ice_wav_frac --> [wav_frac_z16]
      !!
      !! ** Notes   :   A key distinction of this scheme (beyond it not using wave data) is that it
      !!                technically includes any wind-driven fracture process. In the MIZ, this is
      !!                dominantly wave mediated, but in the pack ice it is mediated by deformation
      !!                and internal stresses. So, with this scheme, there may be a need to adjust
      !!                the brittle fracture (routine ice_fsd_bri) parameters to compensate.
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
      REAL(wp), DIMENSION(3,3,jpl)        , INTENT(in)    ::   pa_i      ! ice concentration, local and 8 surrounding cells
      REAL(wp), DIMENSION(nn_nfsd,jpl)    , INTENT(in)    ::   pa_ifsd   ! local floe size-thickness distribution
      REAL(wp), DIMENSION(nn_nfsd)        , INTENT(inout) ::   pQfrac    ! wave fracture probability function (s-1)
      REAL(wp), DIMENSION(nn_nfsd,nn_nfsd), INTENT(inout) ::   pBfrac    ! wave fracture redistribution function, B(s,r)dr
      !
      REAL(wp), DIMENSION(nn_nfsd) ::   zfsd         ! floe size distribution, integrated over thickness
      REAL(wp)                     ::   zr_a         ! mean floe size
      REAL(wp)                     ::   z1minusf_0   ! 1 minus fetch parameter
      REAL(wp)                     ::   zc_b         ! participation factor
      INTEGER                      ::   jf           ! dummy loop index
      !
      !!-------------------------------------------------------------------

      ! --- Calculate floe size distribution integrated over thickness
      !     and mean floe size (only needed if ln_z16_const = F)
      zfsd(:) = 0._wp
      zr_a    = 0._wp   ! initialise

      ! Note: pa_i(2,2,:) is the current grid cell (other indices are surrounding 8 cells)
      ! Hard-coding this for now; maybe in the future make number of surrounding cells
      ! an option, then calculate this, size of input array, and averaging factor (9 below)
      ! including checks against nn_hls (for now not needed as nn_hls >= 1):
      DO jf = 1, nn_nfsd
         zfsd(jf) = SUM( pa_ifsd(jf,:) * pa_i(2,2,:) )
         zr_a     = zr_a + floe_rc(jf) * zfsd(jf)
      ENDDO

      ! --- Calculate participation factor, zc_b --- !
      !
      IF( ln_z16_const ) THEN
         ! Use constant participation factor (as in Zhang et al. 2015, Eq. 15):
         zc_b = rn_z16_cb
         !
      ELSE
         ! Use variable participation factor (as in Zhang et al. 2016, Eq. 5)
         !
         ! Calculate 1 - f_0, where f_0 open water fraction in local + 8 surrounding grid cells
         ! (see comment above: hard-coding the factor of 9 for now):
         z1minusf_0 = SUM( pa_i(:,:,:) ) / 9._wp
         !
         zc_b = EXP( -rn_z16_a * z1minusf_0 - rn_z16_b * (1._wp - zr_a/floe_rc(nn_nfsd)) )
         zc_b = zc_b * rn_z16_k * puatm * rDt_ice / MAX( ph_i, rn_z16_hc )
      ENDIF

      ! --- Calculate source terms --- !
      DO jf = 1, nn_nfsd
         pQfrac(jf) = MAX( 0._wp, 1._wp - SUM(zfsd(jf:nn_nfsd)) / zc_b )
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
      !!                   Q(r) = c_w * exp[ -alpha * (1 - r/rmax)]
      !!
      !!                where c_w is a constant (sets maximum probability of the largest floes
      !!                fracturing), alpha is a constant (determines strength of probability
      !!                reduction at smaller floe sizes), and rmax is the largest resolved floe size.
      !!                This formula differs slightly from Eq. (13) in Yang et al. (2024) which
      !!                is presumed to have a typesetting mistake in the exponent.
      !!
      !!                The redistribution function is determined by assuming any floe of size s
      !!                that is undergoing wave fracture is equally likely to fracture into any
      !!                other floe size r < s. This 'uniform redistributor' is calculated in
      !!                subroutine ice_wav_init as it is a constant.
      !!
      !! ** Inputs  :   phsw                    :   local significant wave height (m)
      !!                pwmp                    :   local wave mean period (s)
      !!                ph_i                    :   local (grid cell) mean sea ice thickness (m)
      !!
      !! ** Outputs :   pQfrac(nn_nfsd)         :   fracture probability function (s-1)
      !!                pBfrac(nn_nfsd,nn_nfsd) :   fracture redistribution function, B(s,r)dr
      !!                                            (note: first  index corresponds to original floe size s,
      !!                                                   second index corresponds to fractured floe size r)
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
      REAL(wp), DIMENSION(nn_nfsd,nn_nfsd), INTENT(inout) ::   pBfrac   ! wave fracture redistribution function, B(s,r)dr
      !
      INTEGER             ::   jf                ! dummy loop indices
      REAL(wp)            ::   zstrain           ! ice strain due to waves
      !
      !!-------------------------------------------------------------------

      pQfrac(:) = 0._wp   ! reset/initialise, and return value if critical strain not exceeded

      zstrain = 2.5_wp * rpi**4 * ph_i * phsw / (grav**2 * pwmp**4)

      IF( zstrain >= rn_ice_wav_ecri ) THEN
         DO jf = 1, nn_nfsd
            pQfrac(jf) = rn_y24a_cw * EXP( -rn_y24a_alpha * (1._wp - floe_rc(jf) / floe_rc(nn_nfsd)) )
         ENDDO
      ENDIF

      pBfrac(:,:) = Bfrac_uni(:,:)   ! uniform redistribution

   END SUBROUTINE wav_frac_y24a


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
      !!                This scheme also makes assumptions relating wavelength L to floe size r:
      !!
      !!                   1) Wavelength L can only break a floe of diameter 2s if L < 2s
      !!                   2) Fractured floe diameter 2r is half the triggering wavelength, so r = L/4
      !!
      !!                From this follows:
      !!
      !!                   Q(s) = int[ q(L') ] dL'   where the integral is from L' = 0 to 2s
      !!
      !!                            int[ q(L') * Theta(2s - L') * delta(4r - L') ] dL'
      !!                   B(s,r) = --------------------------------------------------
      !!                                    int[ q(L') * Theta(2s - L') ] dL'
      !!
      !!                where Theta(x) = {1 if x >= 0; 0 otherwise}, delta(x) = {1 if x = 0; 0 otherwise},
      !!                and the integral is over all wavelengths L. In the code we just check the conditions
      !!                are met explicitly for each combination of s, r, and L rather than using the last
      !!                formal equation.
      !!
      !! ** Inputs  :   pwmp                    :   local wave mean period (s)
      !!                pWspec(nn_nwfreq)       :   local wave spectrum (spectral energy density; m2.Hz-1)
      !!                ph_i                    :   local (grid cell) mean sea ice thickness (m)
      !!
      !! ** Outputs :   pQfrac(nn_nfsd)         :   fracture probability function (s-1)
      !!                pBfrac(nn_nfsd,nn_nfsd) :   fracture redistribution function, B(s,r)dr
      !!                                            (note: first  index corresponds to original floe size s,
      !!                                                   second index corresponds to fractured floe size r)
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
      REAL(wp), DIMENSION(nn_nfsd,nn_nfsd), INTENT(inout) ::   pBfrac   ! wave fracture redistribution function, B(s,r)dr
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
      zprayl  = wav_spec_rayl( pwmp )                            ! local Rayleigh spectrum (Hz-1)

      ! Probability of waves of each frequency class resulting in fracture:
      zprob(:) = 0._wp
      WHERE( zstrain >= rn_ice_wav_ecri ) zprob = zprayl * wdfreq   ! = q(f)df = q(L)dL

      pQfrac(:)   = 0._wp
      pBfrac(:,:) = 0._wp

      DO jf1 = 1, nn_nfsd
         !
         ! Calculate Q(s) <==> pQfrac(jf1):
         DO jw = 1, nn_nwfreq
            ! Each spectral class can only break ice if wavelength L < initial floe diameter
            ! Note: dL is already (implicitly) multiplied into zprob earlier:
            IF( wlam(jw) < 2._wp * floe_rc(jf1) ) THEN
               pQfrac(jf1) = pQfrac(jf1) + zprob(jw)
            ENDIF
         ENDDO
         !
         ! Calculate B(s,r)dr <==> pBfrac(jf1,jf2):
         DO jf2 = 1, nn_nfsd
            ! Each spectral class can only break ice if wavelength L < initial floe diameter
            ! And it can only break ice into floes of diameter L/2
            ! For the latter, we check fracture floe size is within the floe size category limits:
            DO jw = 1, nn_nwfreq
               IF(          (wlam(jw) < 2._wp * floe_rc(jf1) )   &
                  &   .AND. (wlam(jw) / 4._wp >= floe_rl(jf2))   &
                  &   .AND. (wlam(jw) / 4._wp <  floe_ru(jf2))   ) THEN
                  !
                  pBfrac(jf1,jf2) = pBfrac(jf1,jf2) + zprob(jw)
               ENDIF
            ENDDO
            ! Recall pBfrac includes dr factor; pBfrac(jf1,jf2) <==> B(s,r)*dr:
            pBfrac(jf1,jf2) = pBfrac(jf1,jf2) * floe_dr(jf2)
            !
         ENDDO
         !
         ! Normalise:
         IF( SUM(pBfrac(jf1,:)) > 0._wp ) pBfrac(jf1,:) = pBfrac(jf1,:) / SUM(pBfrac(jf1,:))
         !
      ENDDO

   END SUBROUTINE wav_frac_y24b


   SUBROUTINE wav_frac_ht15( pWspec, ph_i, pQfrac, pBfrac )
      !!-------------------------------------------------------------------
      !!                 *** ROUTINE wav_frac_ht15 ***
      !!
      !! ** Purpose :   Calculate the probability and redistribution functions for the
      !!                equation evolving the FSD due to wave fracture for the scheme
      !!                of Horvat and Tziperman (2015) and Roach et al. (2018)
      !!
      !! ** Method  :   This scheme calculates a distribution of fractured ice lengths from
      !!                the local sea surface height (SSH) field, n(x), defined along a 1D
      !!                sub-gridscale domain and computed from the local wave spectrum.
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
      !!                distribution of these lengths (as radii) is then binned into the FSD category
      !!                bins, and the result is called the 'fracture distribution', W(r), such that
      !!                W(r)dr is the number of fracture radii in the interval [r, r+dr]. From this,
      !!                the probability function in the wave fracture equation is:
      !!
      !!                   Q(r) = 1/(D/2) * int[ r'W(r') dr' ]
      !!
      !!                where the integral limits are from the smallest floe size up to r, and the
      !!                redistribution function:
      !!
      !!                   B(s,r) = rW(r) / int[ r'W(r') dr' ]   for r < s
      !!                          = 0                            otherwise
      !!
      !!                where the integral limits are from the smallest floe size up to s.
      !!                Note that int[ B(s,r)dr ] = 1 as required (see subroutine ice_wav_frac).
      !!
      !!                W(r) satisfies int[ r'W(r') dr' ] = D/2, where the integral is over all floe
      !!                sizes, since the sum of all fracture lengths (twice their radii) must equal
      !!                the domain size, D. So, Q(r) is the fraction of all fracture lengths smaller
      !!                than r. This normalisation factor cancels in the expression for B so there
      !!                are no explicit factors of D.
      !!
      !!                Note also there is an implicit factor of (cg/D) in the expression for Q(r),
      !!                where cg is the wave group velocity, representing the fraction of the domain
      !!                reached by waves. This term comes from Horvat and Tziperman (2015) where the
      !!                model is applied to a single grid box; in this context, the wave field is a
      !!                local quantity for the grid cell assumed to affect all the ice in the grid
      !!                cell, so this factor is set to 1 (per second).
      !!
      !! ** Inputs  :   pWspec(nn_nwfreq)       :   local wave spectrum (spectral energy density; m2.Hz-1)
      !!                ph_i                    :   local (grid cell) mean sea ice thickness (m)
      !!
      !! ** Outputs :   pQfrac(nn_nfsd)         :   fracture probability function (s-1)
      !!                pBfrac(nn_nfsd,nn_nfsd) :   fracture redistribution function, B(s,r)dr
      !!                                            (note: first  index corresponds to original floe size s,
      !!                                                   second index corresponds to fractured floe size r)
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
      REAL(wp), DIMENSION(nn_nfsd,nn_nfsd), INTENT(inout) ::   pBfrac   ! wave fracture redistribution function, B(s,r)dr
      !
      INTEGER                              ::   jx, jy, jf              ! dummy loop indices
      INTEGER                              ::   ixlo, ixhi              ! indices of x1d to locate extrema
      INTEGER                              ::   ixfrac                  ! number of fracture points along x1d
      LOGICAL , DIMENSION(nn_ht15_nx1d)    ::   llmin, llmax, llext     ! sea surface height is a min / is a max / is an extrema
      REAL(wp), DIMENSION(nn_nwfreq)       ::   zphi                    ! phase of wave spectrum components (rad)
      REAL(wp), DIMENSION(nn_ht15_nx1d)    ::   zssh                    ! sea surface height along x1d (m)
      REAL(wp)                             ::   zdx, zdxlo, zdxhi       ! distances between x1d points in finite difference computation (m)
      REAL(wp)                             ::   zstrain                 ! strain experienced by sea ice due to wave field
      REAL(wp), DIMENSION(nn_ht15_nx1d)    ::   zxfrac                  ! distances to points along x1d at which ice fractures
      REAL(wp)                             ::   zfrac_rad               ! floe radius of a piece of fractured ice (m)
      REAL(wp), DIMENSION(nn_nfsd)         ::   zWfrac                  ! fracture distribution (multiplied by dr; dimensionless)
      !
      !!-------------------------------------------------------------------

      ! Control:
      IF( ln_timing )   CALL timing_start('wav_frac_ht15')

      ! Initialisation:
      llmin(:)    = .FALSE.
      llmax(:)    = .FALSE.
      ixfrac      = 1
      zxfrac (:)  = 0._wp
      zWfrac (:)  = 0._wp
      pQfrac(:)   = 0._wp
      pBfrac(:,:) = 0._wp

      ! Spectral phases [constant for now; possibly add (optional!) random phase later]:
      zphi(:) = rpi

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
         !
         DO jx = 2, ixfrac - 1
            !
            zfrac_rad = .5_wp * (zxfrac(jx) - zxfrac(jx-1))   ! factor of 0.5 ==> radius of fractured ice
            !
            ! Populate appropriate floe size category in fracture histogram, zWfrac(:)
            ! Just add 1 for now to get relative proportions in each category; scale whole thing afterwards
            ! Note that zWfrac(:) corresponds to W(r)dr in equation, i.e., there is an implicit factor
            ! of floe_dr(:) which is hence also present in zBfrac(:,:) calculated later.
            !
            DO jf = 1, nn_nfsd - 1
               IF( zfrac_rad < floe_ru(jf) ) THEN
                  zWfrac(jf) = zWfrac(jf) + 1._wp
                  EXIT
               ENDIF
            ENDDO
            !
            ! Separate check for largest fractures (even if it exceeds upper bound of largest
            ! floe size category, it goes into that category anyway; note similar for very small
            ! fractures accounted for in above loop anyway):
            IF( zfrac_rad >= floe_rl(nn_nfsd) ) zWfrac(nn_nfsd) = zWfrac(nn_nfsd) + 1._wp
            !
         ENDDO

         ! Scale fracture histogram with floe size of each category
         ! (noting W only appears multiplied by r in equations for Q and B):
         DO jf = 1, nn_nfsd
            zWfrac(jf) = floe_rc(jf) * zWfrac(jf)
         ENDDO

         ! Calculate the probability (pQfrac) and redistribution (pBfrac) functions
         ! from the fracture distribution [zWfrac, which corresponds to rW(r)dr]:
         !
         DO jf = 2, nn_nfsd
            ! B(s,r)dr = rW(r)dr / int[ r'W(r')dr' ] for r < s and r' < s
            IF( SUM(zWfrac(1:jf-1)) > 0._wp ) THEN
               pBfrac(jf,1:jf-1) = zWfrac(1:jf-1) / SUM(zWfrac(1:jf-1))
            ENDIF
         ENDDO
         !
         DO jf = 1, nn_nfsd
            !
            ! Q(r) = 1/(D/2) * int[ r'W(r')dr' ] for r' < r:
            pQfrac(jf) = 2._wp * SUM(zWfrac(1:jf-1)) / (x1d(nn_ht15_nx1d) - x1d(1))
            !
            ! Divide B(s,r)dr calculated above by the integral/sum over r
            ! This normalises it so that int[ B(s,r)dr ] = 1
            ! (it should be anyway, but in case of discretisation errors this ensures it is so):
            !
            IF( SUM(pBfrac(jf,:)) > 0._wp ) pBfrac(jf,:) = pBfrac(jf,:) / SUM(pBfrac(jf,:))
            !
         ENDDO

      ENDIF

      ! Control:
      IF( ln_timing )   CALL timing_stop('wav_frac_ht15')

   END SUBROUTINE wav_frac_ht15


   FUNCTION wav_spec_bret( phsw_l, pwpf_l )
      !!-------------------------------------------------------------------
      !!                *** FUNCTION wav_spec_bret ***
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
      !! ** Inputs  :   phsw_l :   local significant wave height (m)
      !!                pwpf_l :   local peak frequency (Hz)
      !!
      !! ** Outputs :   local wave energy spectrum (Bretschneider; m2.Hz-1)
      !!
      !! ** Notes   :   The approach of using the Bretschneider formula to estimate the wave spectrum
      !!                for wave-ice interactions follows Horvat and Tziperman (2015) and Roach et al. (2018)
      !!                although the functional form here differs. The above matches Williams et al. (2013, Eq. 21)
      !!                and Bateson et al. (2020, Eq. 5), just converting from their expressions for SB(w) in terms
      !!                of angular frequency, w = 2*pi*f, using SB(w)dw = SB(w)(dw/df)df = 2*pi*SB(w)df,
      !!                hence SB(f) = 2*pi*SB(w). This form can be traced back to Bretschneider (1959).
      !!
      !! ** Invokers:   ice_wav_frac  --> [wav_spec_bret()]   (ln_ice_wav_spec=F AND ln_ice_wav_attn=F)
      !!                wav_attn_spec --> [wav_spec_bret()]   (ln_ice_wav_spec=F AND ln_ice_wav_attn=T)
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
      REAL(wp), INTENT(in)           ::   phsw_l          ! local significant wave height (m)
      REAL(wp), INTENT(in)           ::   pwpf_l          ! local peak frequency (Hz)
      !
      REAL(wp), DIMENSION(nn_nwfreq) ::   wav_spec_bret   ! local wave energy spectrum (m2.Hz-1)
      !
      !!-------------------------------------------------------------------

      wav_spec_bret(:) = .3125_wp * phsw_l**2 * pwpf_l**4 * EXP( -1.25_wp * ( pwpf_l / wfreq(:) )**4 ) / wfreq(:)**5

   END FUNCTION wav_spec_bret


   FUNCTION wav_spec_rayl( pwmp_l )
      !!-------------------------------------------------------------------
      !!                *** FUNCTION wav_spec_rayl ***
      !!-------------------------------------------------------------------
      !!
      !! ** Purpose :   Calculate local Rayleigh spectrum P(f)
      !!
      !! ** Method  :   P(f) = 2.7 * exp[ -0.675 / (T_m * f)^4 ] / (T_m^4 * f^5)
      !!
      !!                where f is frequency and T_m is mean wave period. The equation
      !!                is normalised so that that the integral of P over all f is 1.
      !!
      !! ** Inputs  :   pwmp_l : local wave mean period
      !!
      !! ** Outputs :   local Rayleigh spectrum of frequencies, P(f) (Hz-1)
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
      !! ** Invokers:   wav_frac_y24b  --> [wav_spec_rayl()]   (nn_frac_scheme == jpfrac_y24b)
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
      REAL(wp), INTENT(in)           ::   pwmp_l             ! local wave mean period (s)
      REAL(wp), DIMENSION(nn_nwfreq) ::   wav_spec_rayl      ! local Rayleigh spectrum P(f) (Hz-1)
      !
      !!-------------------------------------------------------------------

      wav_spec_rayl(:) = 2.7_wp * EXP( -.675_wp / (pwmp_l * wfreq(:))**4 ) / (pwmp_l**4 * wfreq(:)**5)

      ! Normalise (equation is normalised in theory, but discretisation leads to errors):
      IF( SUM(wav_spec_rayl(:) * wdfreq(:) ) > 0._wp ) THEN
         wav_spec_rayl(:) = wav_spec_rayl(:) / SUM( wav_spec_rayl(:) * wdfreq(:) )
      ENDIF

   END FUNCTION wav_spec_rayl


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
      INTEGER  ::   ji, jj, jf1, jf2, ji_glo, jj_glo   ! Dummy loop indices
      INTEGER  ::   ierr                               ! Local integer output status for allocate
      INTEGER  ::   ios, ioptio                        ! Local integer output status for namelist read
      !
      REAL(wp) ::   zc1, zc2                           ! Terms for uniform fracture redistributor (Bfrac_uni)
      !
      !!
      NAMELIST/namwav/ ln_ice_wav     , ln_ice_wav_spec, rn_ice_wav_ecri, nn_frac_scheme,   &
         &             ln_z16_const   , rn_z16_cb      , rn_z16_k       , rn_z16_a      ,   &
         &             rn_z16_b       , rn_z16_hc      ,                                    &
         &             rn_y24a_cw     , rn_y24a_alpha  ,                                    &
         &             nn_ht15_nx1d   , rn_ht15_dx1d   , nn_ht15_rmin   , ln_ht15_rand  ,   &
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
         WRITE(numout,'(A,I0,A)') '            Yang et al. (2024) A scheme parameters (nn_frac_scheme = ', jpfrac_y24a, '):'
         WRITE(numout,*) '               Probability function parameter c_w           rn_y24a_cw      = ', rn_y24a_cw
         WRITE(numout,*) '               Probability function parameter alpha         rn_y24a_alpha   = ', rn_y24a_alpha
         WRITE(numout,'(A,I0,A)') '            Horvant & Tziperman (2015) scheme parameters (nn_frac_scheme = ', jpfrac_ht15, '):'
         WRITE(numout,*) '               Size of 1D subdomain for SSH                    nn_ht15_nx1d = ', nn_ht15_nx1d
         WRITE(numout,*) '               Increment of 1D subdomain for (m)               rn_ht15_dx1d = ', rn_ht15_dx1d
         WRITE(numout,*) '               Smallest floe radius affected by waves (dx1d)   nn_ht15_rmin = ', nn_ht15_rmin
         WRITE(numout,*) '               Use random phases or not                        ln_ht15_rand = ', ln_ht15_rand
      ENDIF

      IF( ln_ice_wav ) THEN

         ! Checks on flags that do not require knowning SBC wave module flags
         ! (those are handled in subroutine sbc_wave_init).
         !
         ! Wave-ice interactions module requires both wave inputs and FSD
         ! (exception: fracture scheme is Z16, in which case only FSD is needed)
         IF( (nn_frac_scheme /= jpfrac_z16) .AND. .NOT. ln_wave )   &
            &                CALL ctl_stop('ice_wav_init: ln_ice_wav=T but SBC wave module inactive (ln_wave=F)')
         !
         IF( .NOT. ln_fsd  ) CALL ctl_stop('ice_wav_init: ln_ice_wav=T but FSD inactive (ln_fsd=F)')

         ! Warn if both attenuation scheme and reading of full wave spectrum selected
         ! (possible to do so, but usually spectrum comes from a coupled wave model so
         ! should not need to attenuate waves under ice as that is done in the wave model)
         !
         ! (check that spectrum is actually read in is done in sbc_wave_init)
         !
         IF( ln_ice_wav_attn .AND. ln_ice_wav_spec )   &
            &   CALL ctl_warn('ice_wav_init: ln_ice_wav_attn=T but also using spectrum (ln_ice_wav_spec=T): intentional?')

         ! Allocate and define module constants
         !
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
            zc1 = floe_rc(1) / floe_rc(nn_nfsd)   ! = rmin/rmax
            zc2 = 1._wp - zc1
            !
            DO jf1 = 1, nn_nfsd
               DO jf2 = 1, nn_nfsd
                  IF(          (zc1*floe_rc(jf1) <= floe_rc(jf2))          &
                     &   .AND. (floe_rc(jf2) <= zc2*floe_rc(jf1)) ) THEN
                     !
                     Bfrac_uni(jf1,jf2) = floe_dr(jf2) / ((zc2 - zc1)*floe_rc(jf1))
                  ENDIF
               ENDDO
               ! Normalise so that integral of Bfrac_uni(r1,r2)*dr2 = 1:
               IF( SUM(Bfrac_uni(jf1,:)) > 0._wp )   &
                  &   Bfrac_uni(jf1,:) = Bfrac_uni(jf1,:) / SUM(Bfrac_uni(jf1,:))
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
         l_attn_calc_spec = .NOT. ln_ice_wav_spec

         IF(lwp) THEN
            WRITE(numout,*) ''
            WRITE(numout,*) '   Namelist options ==> wave fracture    scheme will calculate spectrum = ', l_frac_calc_spec
            WRITE(numout,*) '                    ==> wave attenuation scheme will calculate spectrum = ', l_attn_calc_spec
         ENDIF

      ENDIF   ! ln_ice_wav

   END SUBROUTINE ice_wav_init


#else
   !!----------------------------------------------------------------------
   !!   Default option          Empty module          NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icewav
