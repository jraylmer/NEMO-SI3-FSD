MODULE icedyn_rdgrft
   !!======================================================================
   !!                       ***  MODULE icedyn_rdgrft ***
   !!    sea-ice : Mechanical impact on ice thickness distribution
   !!======================================================================
   !! History :       !  2006-02  (M. Vancoppenolle) Original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_dyn_rdgrft       : ridging/rafting of sea ice
   !!   ice_dyn_rdgrft_init  : initialization of ridging/rafting of sea ice
   !!   ice_strength         : ice strength calculation
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE par_icedyn     ! SI3 dynamics parameters
   USE phycst         ! physical constants (ocean directory)
   USE sbc_oce , ONLY : sss_m, sst_m   ! surface boundary condition: ocean fields
   USE ice            ! sea-ice: variables
   USE icevar  , ONLY : ice_var_roundoff
   USE icectl         ! sea-ice: control prints
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_rdgrft        ! called by icestp
   PUBLIC   ice_dyn_rdgrft_init   ! called by icedyn
   PUBLIC   ice_strength          ! called by icedyn_rhg_evp

   INTEGER ::              nice_str   ! choice of the type of strength
   !                                        ! associated indices:
   INTEGER, PARAMETER ::   np_strh79     = 1   ! Hibler 79
   INTEGER, PARAMETER ::   np_strr75     = 2   ! Rothrock 75
   INTEGER, PARAMETER ::   np_strcst     = 3   ! Constant value

   ! Variables shared among ridging subroutines
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   closing_net     ! net rate at which area is removed    (1/s)
      !                                                                 ! (ridging ice area - area of new ridges) / dt
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   opening         ! rate of opening due to divergence/shear
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   closing_gross   ! rate at which area removed, not counting area of new ridges
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   apartf          ! participation function; fraction of ridging/closing associated w/ category n
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hrmin           ! minimum ridge thickness
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hrmax           ! maximum ridge thickness
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hrexp           ! e-folding ridge thickness
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hraft           ! thickness of rafted ice
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   hi_hrdg         ! thickness of ridging ice / mean ridge thickness
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   aridge          ! participating ice ridging
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   araft           ! participating ice rafting
   !
   ! For ridging diagnostics
   LOGICAL                                         ::   ll_diag_rdg     ! activate ridging diagnostics or not
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   airdg1          ! ridging ice area loss
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   airft1          ! rafting ice area loss
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   airdg2          ! new ridged ice area gain
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   airft2          ! new rafted ice area gain
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   virdg           ! new ridged ice volume gain
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   virft           ! new rafted ice volume gain
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_airdg1     ! ridging ice area loss rate diagnostic
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_airft1     ! rafting ice area loss rate diagnostic
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_airdg2     ! new ridged ice area gain rate diagnostic
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_airft2     ! new rafted ice area gain rate diagnostic
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_virdg      ! new ridged ice volume gain rate diagnostic
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_virft      ! new rafted ice volume gain rate diagnostic
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_closing    ! closing rate diagnostic
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:)     ::   diag_opening    ! opening rate diagnostic
   !
   REAL(wp), PARAMETER ::   hrdg_hi_min = 1.1_wp    ! min ridge thickness multiplier: min(hrdg/hi)
   REAL(wp), PARAMETER ::   hi_hrft     = 0.5_wp    ! rafting multiplier: (hi/hraft)
   !
   ! ** namelist (namdyn_rdgrft) **
   LOGICAL  ::   ln_str_smooth    ! ice strength spatial smoothing
   LOGICAL  ::   ln_distf_lin     ! redistribution of ridged ice: linear (Hibler 1980)
   LOGICAL  ::   ln_distf_exp     ! redistribution of ridged ice: exponential (Lipscomb et al 2017)
   REAL(wp) ::   rn_murdg         !    gives e-folding scale of ridged ice (m^.5)
   REAL(wp) ::   rn_csrdg         ! fraction of shearing energy contributing to ridging
   LOGICAL  ::   ln_partf_lin     ! participation function linear (Thorndike et al. (1975))
   REAL(wp) ::   rn_gstar         !    fractional area of young ice contributing to ridging
   LOGICAL  ::   ln_partf_exp     ! participation function exponential (Lipscomb et al. (2007))
   REAL(wp) ::   rn_astar         !    equivalent of G* for an exponential participation function
   LOGICAL  ::   ln_ridging       ! ridging of ice or not
   REAL(wp) ::   rn_hstar         !    thickness that determines the maximal thickness of ridged ice
   REAL(wp) ::   rn_porordg       !    initial porosity of ridges (0.3 regular value)
   REAL(wp) ::   rn_fsnwrdg       !    fractional snow loss to the ocean during ridging
   REAL(wp) ::   rn_fpndrdg       !    fractional pond loss to the ocean during ridging
   LOGICAL  ::   ln_rafting       ! rafting of ice or not
   REAL(wp) ::   rn_hraft         !    threshold thickness (m) for rafting / ridging
   REAL(wp) ::   rn_craft         !    coefficient for smoothness of the hyperbolic tangent in rafting
   REAL(wp) ::   rn_fsnwrft       !    fractional snow loss to the ocean during rafting
   REAL(wp) ::   rn_fpndrft       !    fractional pond loss to the ocean during rafting
   !
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION ice_dyn_rdgrft_alloc()
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_dyn_rdgrft_alloc ***
      !!-------------------------------------------------------------------
      ALLOCATE( closing_net(A2D(0)) , opening(A2D(0))   , closing_gross(A2D(0)),                     &
         &      apartf(A2D(0),0:jpl), hrmin (A2D(0),jpl), hraft(A2D(0),jpl)    , aridge(A2D(0),jpl), &
         &      hrmax (A2D(0),jpl)  , hrexp (A2D(0),jpl), hi_hrdg(A2D(0),jpl)  , araft(A2D(0),jpl) , &
         &      airdg1(A2D(0))      , airft1(A2D(0))    , airdg2(A2D(0))       , airft2(A2D(0))    , &
         &      virdg (A2D(0))      , virft (A2D(0))    ,                                            &
         &      STAT=ice_dyn_rdgrft_alloc )
      !
      CALL mpp_sum ( 'icedyn_rdgrft', ice_dyn_rdgrft_alloc )
      IF( ice_dyn_rdgrft_alloc /= 0 )   CALL ctl_stop( 'STOP',  'ice_dyn_rdgrft_alloc: failed to allocate arrays'  )
      !
   END FUNCTION ice_dyn_rdgrft_alloc


   SUBROUTINE ice_dyn_rdgrft( kt )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_dyn_rdgrft ***
      !!
      !! ** Purpose :   computes the mechanical redistribution of ice thickness
      !!
      !! ** Method  :   Steps :
      !!       0) Identify grid cells with ice
      !!       1) Calculate closing rate, divergence and opening
      !!       2) Identify grid cells with ridging
      !!       3) Start ridging iterations
      !!          - prep = ridged and rafted ice + closing_gross
      !!          - shift = move ice from one category to another
      !!
      !! ** Details
      !!    step1: The net rate of closing is due to convergence and shear, based on Flato and Hibler (1995).
      !!           The energy dissipation rate is equal to the net closing rate times the ice strength.
      !!
      !!    step3: The gross closing rate is equal to the first two terms (open
      !!           water closing and thin ice ridging) without the third term
      !!           (thick, newly ridged ice).
      !!
      !! References :   Flato, G. M., and W. D. Hibler III, 1995, JGR, 100, 18,611-18,626.
      !!                Hibler, W. D. III, 1980, MWR, 108, 1943-1973, 1980.
      !!                Rothrock, D. A., 1975: JGR, 80, 4514-4519.
      !!                Thorndike et al., 1975, JGR, 80, 4501-4513.
      !!                Bitz et al., JGR, 2001
      !!                Amundrud and Melling, JGR 2005
      !!                Babko et al., JGR 2002
      !!
      !!     This routine is based on CICE code and authors William H. Lipscomb,
      !!     and Elizabeth C. Hunke, LANL are gratefully acknowledged
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! number of iteration
      !!
      INTEGER  ::   ji, jj, jk, jl             ! dummy loop index
      INTEGER  ::   iter, iterate_ridging      ! local integer
      REAL(wp) ::   zfac                       ! local scalar
      LOGICAL,  DIMENSION(A2D(0)) ::   ll_ridge       ! To determine whether ridging is required or not within local domain
      REAL(wp), DIMENSION(A2D(0)) ::   zdivu          ! divergence
      !!
      INTEGER, PARAMETER ::   jp_itermax = 20
      !!-------------------------------------------------------------------
      ! controls
      IF( ln_timing    )   CALL timing_start('icedyn_rdgrft')                                                             ! timing
      IF( ln_icediachk )   CALL ice_cons_hsm(0, 'icedyn_rdgrft', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icediachk )   CALL ice_cons2D  (0, 'icedyn_rdgrft',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'ice_dyn_rdgrft: ice ridging and rafting'
         IF(lwp) WRITE(numout,*)'~~~~~~~~~~~~~~'
         ! diagnostics
         IF(  iom_use('lead_open') .OR. iom_use('lead_close') .OR. iom_use('ardg_loss') .OR. iom_use('arft_loss')  .OR. &
            & iom_use('ardg_gain') .OR. iom_use('arft_gain' ) .OR. iom_use('vrdg_gain') .OR. iom_use('vrft_gain') ) THEN
            !
            ll_diag_rdg = .TRUE.
            !
            ALLOCATE( diag_airdg1(A2D(0)) , diag_airft1(A2D(0)) , diag_airdg2 (A2D(0)) , diag_airft2 (A2D(0)) , &
               &      diag_virdg (A2D(0)) , diag_virft (A2D(0)) , diag_closing(A2D(0)) , diag_opening(A2D(0)) )
            !
         ELSE
            ll_diag_rdg = .FALSE.
         ENDIF
         !
      ENDIF

      ! Initialise ridging diagnostics if required
      IF( ll_diag_rdg ) THEN    
         diag_closing(:,:) = 0._wp ; diag_opening(:,:) = 0._wp
         diag_airdg1 (:,:) = 0._wp ; diag_airft1 (:,:) = 0._wp
         diag_airdg2 (:,:) = 0._wp ; diag_airft2 (:,:) = 0._wp
         diag_virdg  (:,:) = 0._wp ; diag_virft  (:,:) = 0._wp
      ENDIF

      !--------------------------------
      ! 0) Identify grid cells with ice
      !--------------------------------
      at_i(A2D(0)) = SUM( a_i(A2D(0),:), dim=3 )

      ! Initialise logical array for ice points where ridging occurs
      ll_ridge(:,:) = .FALSE.

      ! Initialise logical array for ice points
      DO_2D( 0, 0, 0, 0 )
         IF ( at_i(ji,jj) > epsi10 ) THEN
            l_ice_present(ji,jj) = .TRUE.
         ELSE
            l_ice_present(ji,jj) = .FALSE.
         ENDIF
      END_2D
      !--------------------------------------------------------
      ! 1) Dynamical inputs (closing rate, divergence, opening)
      !--------------------------------------------------------
      IF( ANY( l_ice_present(A2D(0)) ) ) THEN

         DO_2D( 0, 0, 0, 0 )
            IF( l_ice_present(ji,jj) ) THEN
               ! closing_net = rate at which open water area is removed + ice area removed by ridging
               !                                                        - ice area added in new ridges
               IF( ln_rhg_EVP .OR. ln_rhg_VP )   closing_net(ji,jj) = rn_csrdg * 0.5_wp * ( delta_i(ji,jj) - ABS( divu_i(ji,jj) ) ) &
                  &                                                                       - MIN( divu_i(ji,jj), 0._wp )
               IF( ln_rhg_EAP )                  closing_net(ji,jj) = rdg_conv(ji,jj)
               !
               IF( divu_i(ji,jj) < 0._wp )       closing_net(ji,jj) = MAX( closing_net(ji,jj), -divu_i(ji,jj) )   ! make sure the closing rate is large enough
               !                                                                                ! to give asum = 1.0 after ridging
               ! Opening rate (non-negative) that will give asum = 1.0 after ridging.
               opening(ji,jj) = closing_net(ji,jj) + divu_i(ji,jj)
            ENDIF
         END_2D
         !
         !------------------------------------
         ! 2) Identify grid cells with ridging
         !------------------------------------
         CALL rdgrft_prep( a_i, v_i, ato_i, l_ice_present(A2D(0)),  &                                            ! <<== in
            &              apartf, aridge, araft, hi_hrdg, hraft, hrmin, hrmax, hrexp, closing_gross, opening, & ! ==>> out
            &              closing_net )                                                                         ! <<== in

         DO_2D( 0, 0, 0, 0 )
            IF( l_ice_present(ji,jj) )  THEN
               IF( SUM( apartf(ji,jj,1:jpl) ) > 0._wp .AND. closing_gross(ji,jj) > 0._wp ) THEN
                  ll_ridge(ji,jj)=.TRUE. ! Set to true only if ridging is required
               ENDIF
            ENDIF
         END_2D

      ENDIF

      !-----------------
      ! 3) Start ridging
      !-----------------
      IF( ANY( ll_ridge(:,:) ) ) THEN

         iter            = 1
         iterate_ridging = 1
         !                                                        !----------------------!
         DO WHILE( iterate_ridging > 0 .AND. iter <= jp_itermax ) !  ridging iterations  !
            !                                                     !----------------------!
            ! Calculate participation function (apartf)
            !       and transfer      function
            !       and closing_gross (+correction on opening)

            CALL rdgrft_prep( a_i, v_i, ato_i, ll_ridge, &                                                          ! <<== in
               &              apartf, aridge, araft, hi_hrdg, hraft, hrmin, hrmax, hrexp, closing_gross, opening, & ! ==>> out
               &              closing_net )                                                                         ! <<== in

            ! Redistribute area, volume, and energy between categories
            CALL rdgrft_shift( apartf, aridge, araft, hi_hrdg, hraft, hrmin, hrmax, hrexp, closing_gross, opening, & ! <<== in
               &               sst_m, sss_m, ll_ridge )                                                              ! <<== in


            ! Do we keep on iterating?
            !-------------------------
            ! Check whether a_i + ato_i = 0
            ! If not, because the closing and opening rates were reduced above, ridge again with new rates

            iterate_ridging = 0
            DO_2D( 0, 0, 0, 0 )
               !
               IF( ll_ridge(ji,jj) ) THEN
                  !
                  ! --- Ridging diagnostics --- !
                  IF( ll_diag_rdg ) THEN
                     diag_opening(ji,jj) = diag_opening(ji,jj) + opening(ji,jj)
                     diag_closing(ji,jj) = diag_closing(ji,jj) + closing_net(ji,jj)
                  ENDIF
                  !
                  ! --- Check convergence --- !
                  zfac = 1._wp - ( ato_i(ji,jj) + SUM( a_i(ji,jj,:) ) )
                  IF( ABS( zfac ) < epsi10 .OR. iter == jp_itermax ) THEN
                     closing_net(ji,jj) = 0._wp
                     opening    (ji,jj) = 0._wp
                     ato_i      (ji,jj) = MAX( 0._wp, 1._wp - SUM( a_i(ji,jj,:) ) )
                  ELSE
                     iterate_ridging  = iterate_ridging + 1
                     zdivu      (ji,jj) = zfac * r1_Dt_ice
                     closing_net(ji,jj) = MAX( 0._wp, -zdivu(ji,jj) )
                     opening    (ji,jj) = MAX( 0._wp,  zdivu(ji,jj) )
                  ENDIF
                  !
               ENDIF
               !
            END_2D
            !
            iter = iter + 1
            !
         END DO
         !clem: check this stop. I think sometimes it does not converge. Do we want to stop the simulation?
!!$         IF( iter  >  jp_itermax )    CALL ctl_stop( 'STOP',  'icedyn_rdgrft: non-converging ridging scheme'  )
         !
      ENDIF  ! high-level ridging check

      ! --- Ridging diagnostics --- !
      IF( ll_diag_rdg ) THEN
         CALL iom_put( 'lead_open' , diag_opening(:,:) )  ! Lead area opening rate
         CALL iom_put( 'lead_close', diag_closing(:,:) )  ! Lead area closing rate
         CALL iom_put( 'ardg_loss' , diag_airdg1 (:,:) )  ! Ridging ice area loss rate
         CALL iom_put( 'arft_loss' , diag_airft1 (:,:) )  ! Rafting ice area loss rate
         CALL iom_put( 'ardg_gain' , diag_airdg2 (:,:) )  ! New ridged ice area gain rate
         CALL iom_put( 'arft_gain' , diag_airft2 (:,:) )  ! New rafted ice area gain rate
         CALL iom_put( 'vrdg_gain' , diag_virdg  (:,:) )  ! New ridged ice volume gain rate
         CALL iom_put( 'vrft_gain' , diag_virft  (:,:) )  ! New rafted ice volume gain rate
      ENDIF

      ! clem: those fields must be updated on the halos: ato_i, a_i, v_i, v_s, sv_i, oa_i, a_ip, v_ip, v_il, e_i, e_s, szv_i

      ! clem: I think we can comment this line but I am not sure it does not change results
!!$      CALL ice_var_agg( 1 )

      ! controls
      IF( sn_cfctl%l_prtctl )   CALL ice_prt3D('icedyn_rdgrft')                                                           ! prints
      IF( ln_icectl    )   CALL ice_prt     (kt, iiceprt, jiceprt,-1, ' - ice dyn rdgrft - ')                             ! prints
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icedyn_rdgrft', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icediachk )   CALL ice_cons2D  (1, 'icedyn_rdgrft',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation
      IF( ln_timing    )   CALL timing_stop ('icedyn_rdgrft')                                                             ! timing
      !
   END SUBROUTINE ice_dyn_rdgrft


   SUBROUTINE rdgrft_prep( pa_i, pv_i, pato_i, lice_pres,  &                                                               ! <<== in
      &                    papartf, paridge, paraft, phi_hrdg, phraft, phrmin, phrmax, phrexp, pclosing_gross, popening, & ! ==>> out
      &                    pclosing_net )                                                                                  ! <<== in
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE rdgrft_prep ***
      !!
      !! ** Purpose :   preparation for ridging calculations
      !!
      !! ** Method  :   Compute the thickness distribution of the ice and open water
      !!                participating in ridging and of the resulting ridges.
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpl) , INTENT(in   ) ::   pa_i, pv_i
      REAL(wp), DIMENSION(jpi,jpj)     , INTENT(in   ) ::   pato_i
      LOGICAL , DIMENSION(A2D(0))      , INTENT(in   ) ::   lice_pres   ! ice presence
      REAL(wp), DIMENSION(A2D(0),0:jpl), INTENT(  out) ::   papartf
      REAL(wp), DIMENSION(A2D(0), jpl) , INTENT(  out) ::   paridge, paraft, phi_hrdg, phraft, phrmin, phrmax, phrexp
      REAL(wp), DIMENSION(A2D(0))      , INTENT(inout), OPTIONAL ::   pclosing_gross, popening
      REAL(wp), DIMENSION(A2D(0))      , INTENT(in   ), OPTIONAL ::   pclosing_net
      !!
      INTEGER  ::   ji, jj, jl                  ! dummy loop indices
      REAL(wp) ::   z1_gstar, z1_astar, zhmean, zfac   ! local scalar
      REAL(wp), DIMENSION(A2D(0))        ::   zasum, z1_asum           ! sum of a_i+ato_i and inverse
      REAL(wp), DIMENSION(A2D(0))        ::   zaksum                   ! normalisation factor
      REAL(wp), DIMENSION(A2D(0),jpl)    ::   zhi                      ! ice thickness
      REAL(wp), DIMENSION(A2D(0),-1:jpl) ::   zGsum                    ! zGsum(n) = sum of areas in categories 0 to n
     !--------------------------------------------------------------------

      z1_gstar = 1._wp / rn_gstar
      z1_astar = 1._wp / rn_astar

      ! Ice thickness needed for rafting
      ! In single precision there were floating point invalids due a sqrt of zhi which happens to have negative values
      ! To solve that an extra check about the value of pv_i was added.
      ! Although adding this condition is safe, the double definition (one for single other for double) has been kept to preserve the results of the sette test.
#if defined key_single

      DO jl = 1, jpl
         DO_2D( 0, 0, 0, 0 )
          IF( lice_pres(ji,jj) ) THEN
             IF (pa_i(ji,jj,jl) > epsi10 .AND. pv_i(ji,jj,jl) > epsi10 ) THEN
                zhi(ji,jj,jl) = pv_i(ji,jj,jl) / pa_i(ji,jj,jl)
#else
      DO jl = 1, jpl
         DO_2D( 0, 0, 0, 0 )
            IF( lice_pres(ji,jj) ) THEN
               IF (pa_i(ji,jj,jl) > epsi10 ) THEN
                  zhi(ji,jj,jl) = pv_i(ji,jj,jl) / pa_i(ji,jj,jl)

#endif
               ELSE
                  zhi(ji,jj,jl) = 0._wp
               ENDIF
            ENDIF
         END_2D
      END DO

      ! 1) Participation function (apartf): a(h) = b(h).g(h)
      !-----------------------------------------------------------------
      ! Compute the participation function = total area lost due to ridging/closing
      ! This is analogous to
      !   a(h) = b(h)g(h) as defined in Thorndike et al. (1975).
      !   assuming b(h) = (2/Gstar) * (1 - G(h)/Gstar).
      !
      ! apartf = integrating b(h)g(h) between the category boundaries
      ! apartf is always >= 0 and SUM(apartf(0:jpl))=1
      !-----------------------------------------------------------------
      !
      ! Compute total area of ice plus open water.
      ! This is in general not equal to one because of divergence during transport
      DO_2D( 0, 0, 0, 0 )
         IF( lice_pres(ji,jj) ) THEN
            zasum(ji,jj) = pato_i(ji,jj) + SUM( pa_i(ji,jj,:))
            IF (zasum(ji,jj) > epsi10 ) THEN
               z1_asum(ji,jj) = 1._wp / zasum(ji,jj)
            ELSE
               z1_asum(ji,jj) = 0._wp
            ENDIF
         ENDIF
      END_2D
      !
      ! Compute cumulative thickness distribution function
      ! Compute the cumulative thickness distribution function zGsum,
      ! where zGsum(n) is the fractional area in categories 0 to n.
      ! initial value (in h = 0) = open water area
      DO_2D( 0, 0, 0, 0 )
         IF( lice_pres(ji,jj) ) THEN
            zGsum(ji,jj,-1) = 0._wp
            zGsum(ji,jj,0 ) = pato_i(ji,jj) * z1_asum(ji,jj)
         ENDIF
      END_2D
      DO jl = 1, jpl
         DO_2D( 0, 0, 0, 0 )
            IF( lice_pres(ji,jj) ) THEN
               zGsum(ji,jj,jl) = ( pato_i(ji,jj) + SUM( pa_i(ji,jj,1:jl)) ) * z1_asum(ji,jj)  ! sum(1:jl) is correct (and not jpl)
            ENDIF
         END_2D
      END DO
      !
      IF( ln_partf_lin ) THEN          !--- Linear formulation (Thorndike et al., 1975)
         DO jl = 0, jpl
            DO_2D( 0, 0, 0, 0 )
               IF( lice_pres(ji,jj) ) THEN
                  IF    ( zGsum(ji,jj,jl)   < rn_gstar ) THEN
                     papartf(ji,jj,jl) = z1_gstar * ( zGsum(ji,jj,jl)   - zGsum(ji,jj,jl-1) ) * &
                        &                 ( 2._wp - ( zGsum(ji,jj,jl-1) + zGsum(ji,jj,jl)   ) * z1_gstar )
                  ELSEIF( zGsum(ji,jj,jl-1) < rn_gstar ) THEN
                     papartf(ji,jj,jl) = z1_gstar * ( rn_gstar          - zGsum(ji,jj,jl-1) ) *  &
                        &                 ( 2._wp - ( zGsum(ji,jj,jl-1) + rn_gstar          ) * z1_gstar )
                  ELSE
                     papartf(ji,jj,jl) = 0._wp
                  ENDIF
               ENDIF
            END_2D
         END DO
         !
      ELSEIF( ln_partf_exp ) THEN      !--- Exponential, more stable formulation (Lipscomb et al, 2007)
         !
         zfac = 1._wp / ( 1._wp - EXP(-z1_astar) )
         DO_2D( 0, 0, 0, 0 )
            IF( lice_pres(ji,jj) ) THEN
               zGsum(ji,jj,-1) = EXP( -zGsum(ji,jj,-1) * z1_astar ) * zfac
            ENDIF
         END_2D
         DO jl = 0, jpl
            DO_2D( 0, 0, 0, 0 )
               IF( lice_pres(ji,jj) ) THEN
                  zGsum  (ji,jj,jl) = EXP( -zGsum(ji,jj,jl) * z1_astar ) * zfac
                  papartf(ji,jj,jl) = zGsum(ji,jj,jl-1) - zGsum(ji,jj,jl)
               ENDIF
            END_2D
         END DO
         !
      ENDIF
      !                                !--- Ridging and rafting participation concentrations
      IF( ln_rafting .AND. ln_ridging ) THEN             !- ridging & rafting
         DO jl = 1, jpl
            DO_2D( 0, 0, 0, 0 )
               IF( lice_pres(ji,jj) ) THEN
                  paridge(ji,jj,jl) = ( 1._wp + TANH ( rn_craft * ( zhi(ji,jj,jl) - rn_hraft ) ) ) * 0.5_wp * papartf(ji,jj,jl)
                  paraft (ji,jj,jl) = papartf(ji,jj,jl) - paridge(ji,jj,jl)
               ENDIF
            END_2D
         END DO
      ELSEIF( ln_ridging .AND. .NOT. ln_rafting ) THEN   !- ridging alone
         DO jl = 1, jpl
            DO_2D( 0, 0, 0, 0 )
               IF( lice_pres(ji,jj) ) THEN
                  paridge(ji,jj,jl) = papartf(ji,jj,jl)
                  paraft (ji,jj,jl) = 0._wp
               ENDIF
            END_2D
         END DO
      ELSEIF( ln_rafting .AND. .NOT. ln_ridging ) THEN   !- rafting alone
         DO jl = 1, jpl
            DO_2D( 0, 0, 0, 0 )
               IF( lice_pres(ji,jj) ) THEN
                  paridge(ji,jj,jl) = 0._wp
                  paraft (ji,jj,jl) = papartf(ji,jj,jl)
               ENDIF
            END_2D
         END DO
      ELSE                                               !- no ridging & no rafting
         DO jl = 1, jpl
            DO_2D( 0, 0, 0, 0 )
               IF( lice_pres(ji,jj) ) THEN
                  paridge(ji,jj,jl) = 0._wp
                  paraft (ji,jj,jl) = 0._wp
               ENDIF
            END_2D
         END DO
      ENDIF

      ! 2) Transfer function
      !-----------------------------------------------------------------
      ! If assuming ridged ice is uniformly distributed between hrmin and  
      ! hrmax (ln_distf_lin):
      !
      ! Compute max and min ridged ice thickness for each ridging category.
      !
      ! This parameterization is a modified version of Hibler (1980).
      ! The mean ridging thickness, zhmean, is proportional to hi^(0.5)
      !  and for very thick ridging ice must be >= hrdg_hi_min*hi
      !
      ! The minimum ridging thickness, hrmin, is equal to 2*hi
      !  (i.e., rafting) and for very thick ridging ice is
      !  constrained by hrmin <= (zhmean + hi)/2.
      !
      ! The maximum ridging thickness, hrmax, is determined by zhmean and hrmin.
      !
      ! These modifications have the effect of reducing the ice strength
      ! (relative to the Hibler formulation) when very thick ice is ridging.
      !
      !-----------------------------------------------------------------
      ! If assuming ridged ice ITD is a negative exponential 
      ! (ln_distf_exp) and following CICE implementation: 
      ! 
      !  g(h) ~ exp[-(h-hrmin)/hrexp], h >= hrmin 
      ! 
      ! where hrmin is the minimum thickness of ridging ice and 
      ! hrexp is the e-folding thickness.
      ! 
      ! Here, assume as above that hrmin = min(2*hi, hi+maxraft).
      ! That is, the minimum ridge thickness results from rafting,
      !  unless the ice is thicker than maxraft.
      !
      ! Also, assume that hrexp = mu_rdg*sqrt(hi).
      ! The parameter mu_rdg is tuned to give e-folding scales mostly
      !  in the range 2-4 m as observed by upward-looking sonar.
      !
      ! Values of mu_rdg in the right column give ice strengths
      !  roughly equal to values of Hstar in the left column
      !  (within ~10 kN/m for typical ITDs):
      !
      !   Hstar     mu_rdg
      !
      !     25        3.0
      !     50        4.0
      !     75        5.0
      !    100        6.0
      !
      ! zaksum = net area removed/ total area participating
      ! where total area participating = area of ice that ridges
      !         net area removed = total area participating - area of new ridges
      !-----------------------------------------------------------------
      zfac = 1._wp / hi_hrft
      DO_2D( 0, 0, 0, 0 )
         IF( lice_pres(ji,jj) ) THEN
            zaksum(ji,jj) = papartf(ji,jj,0)
         ENDIF
      END_2D
      !
      DO jl = 1, jpl
         DO_2D( 0, 0, 0, 0 )
            IF( lice_pres(ji,jj) ) THEN
               IF ( papartf(ji,jj,jl) > 0._wp ) THEN
                  zhmean           = MAX( SQRT( rn_hstar * zhi(ji,jj,jl) ), zhi(ji,jj,jl) * hrdg_hi_min )
                  phrmin(ji,jj,jl) = MIN( 2._wp * zhi(ji,jj,jl), 0.5_wp * ( zhmean + zhi(ji,jj,jl) ) )
                  phraft(ji,jj,jl) = zhi(ji,jj,jl) * zfac
                  !
                  IF( ln_distf_lin ) THEN
                     phrmax  (ji,jj,jl) = 2._wp * zhmean - phrmin(ji,jj,jl)
                     phi_hrdg(ji,jj,jl) = zhi(ji,jj,jl) / MAX( zhmean, epsi20 )
                  ELSEIF( ln_distf_exp ) THEN
                     phrexp  (ji,jj,jl) = rn_murdg * SQRT( zhi(ji,jj,jl) )
                     phi_hrdg(ji,jj,jl) = zhi(ji,jj,jl) / MAX( epsi20, phrmin(ji,jj,jl) + phrexp(ji,jj,jl) )
                  ENDIF
                  !
                  ! Normalization factor : zaksum, ensures mass conservation
                  zaksum(ji,jj) = zaksum(ji,jj) + paridge(ji,jj,jl) * ( 1._wp - phi_hrdg(ji,jj,jl) )    &
                     &                          + paraft (ji,jj,jl) * ( 1._wp - hi_hrft )
               ELSE
                  phrmin  (ji,jj,jl) = 0._wp
                  phrmax  (ji,jj,jl) = 0._wp
                  phrexp  (ji,jj,jl) = 0._wp
                  phraft  (ji,jj,jl) = 0._wp
                  phi_hrdg(ji,jj,jl) = 1._wp
               ENDIF
            ENDIF
         END_2D
      END DO

      !
      IF( PRESENT( pclosing_net ) ) THEN
         !
         ! 3) closing_gross
         !-----------------
         ! Based on the ITD of ridging and ridged ice, convert the net closing rate to a gross closing rate.
         ! NOTE: 0 < aksum <= 1
         DO_2D( 0, 0, 0, 0 )
            IF( lice_pres(ji,jj) ) THEN
               IF (zaksum(ji,jj) > epsi10 ) THEN
                  pclosing_gross(ji,jj) = pclosing_net(ji,jj) / zaksum(ji,jj)
               ELSE
                  pclosing_gross(ji,jj) = 0._wp
               ENDIF
            ENDIF
         END_2D

         ! correction to closing rate if excessive ice removal
         !----------------------------------------------------
         ! Reduce the closing rate if more than 100% of any ice category would be removed
         ! Reduce the opening rate in proportion
         DO jl = 1, jpl
            DO_2D( 0, 0, 0, 0 )
               IF( lice_pres(ji,jj) ) THEN
                  zfac = papartf(ji,jj,jl) * pclosing_gross(ji,jj) * rDt_ice
                  IF( zfac > pa_i(ji,jj,jl) .AND. papartf(ji,jj,jl) /= 0._wp ) THEN
                     pclosing_gross(ji,jj) = pa_i(ji,jj,jl) / papartf(ji,jj,jl) * r1_Dt_ice
                  ENDIF
               ENDIF
            END_2D
         END DO
         
         ! 4) correction to opening if excessive open water removal
         !---------------------------------------------------------
         ! Reduce the closing rate if more than 100% of the open water would be removed
         ! Reduce the opening rate in proportion
         DO_2D( 0, 0, 0, 0 )
            IF( lice_pres(ji,jj) ) THEN
               zfac = pato_i(ji,jj) + ( popening(ji,jj) - papartf(ji,jj,0) * pclosing_gross(ji,jj) ) * rDt_ice
               IF( zfac < 0._wp ) THEN           ! would lead to negative ato_i
                  popening(ji,jj) = papartf(ji,jj,0) * pclosing_gross(ji,jj) - pato_i(ji,jj) * r1_Dt_ice
               ELSEIF( zfac > zasum(ji,jj) ) THEN   ! would lead to ato_i > asum
                  popening(ji,jj) = papartf(ji,jj,0) * pclosing_gross(ji,jj) + ( zasum(ji,jj) - pato_i(ji,jj) ) * r1_Dt_ice
               ENDIF
            ENDIF
         END_2D
         !

      ENDIF
      !
   END SUBROUTINE rdgrft_prep


   SUBROUTINE rdgrft_shift( papartf, paridge, paraft, phi_hrdg, phraft, phrmin, phrmax, phrexp, pclosing_gross, popening, & ! <<== in
      &                     psst, psss, lice_pres )                                                                         ! <<== in
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE rdgrft_shift ***
      !!
      !! ** Purpose :   shift ridging ice among thickness categories of ice thickness
      !!
      !! ** Method  :   Remove area, volume, and energy from each ridging category
      !!                and add to thicker ice categories.
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(A2D(0),0:jpl)     , INTENT(in) ::   papartf
      REAL(wp), DIMENSION(A2D(0),jpl)       , INTENT(in) ::   paridge, paraft, phi_hrdg, phraft, phrmin, phrmax, phrexp
      REAL(wp), DIMENSION(A2D(0))           , INTENT(in) ::   pclosing_gross, popening
      REAL(wp), DIMENSION(jpi,jpj)          , INTENT(in) ::   psst, psss
      LOGICAL , DIMENSION(A2D(0))           , INTENT(in) ::   lice_pres  ! ice presence
     !
      INTEGER  ::   ji, jj, jl, jl1, jl2, jk       ! dummy loop indices
      REAL(wp) ::   hL, hR, farea              ! left and right limits of integration and new area going to jl2
      REAL(wp) ::   expL, expR                 ! exponentials involving hL, hR
      REAL(wp) ::   vsw                        ! vol of water trapped into ridges
      REAL(wp) ::   afrdg, afrft               ! fraction of category area ridged/rafted
      REAL(wp), DIMENSION(A2D(0)) ::   oirdg, aprdg, vsrdg, vprdg, vlrdg  ! area etc of new ridges
      REAL(wp), DIMENSION(A2D(0)) ::   oirft, aprft, vsrft, vprft, vlrft  ! area etc of rafted ice
      REAL(wp), DIMENSION(A2D(0),nlay_i) ::   sirdg, sirft
      !
      REAL(wp) ::   ersw             ! enthalpy of water trapped into ridges
      REAL(wp) ::   zswitch, fvol    ! new ridge volume going to jl2
      REAL(wp) ::   z1_ai            ! 1 / a
      REAL(wp), DIMENSION(A2D(0)) ::   zvti             ! sum(v_i)
      !
      REAL(wp), DIMENSION(A2D(0),nlay_s) ::   esrft     ! snow energy of rafting ice
      REAL(wp), DIMENSION(A2D(0),nlay_i) ::   eirft     ! ice  energy of rafting ice
      REAL(wp), DIMENSION(A2D(0),nlay_s) ::   esrdg     ! enth*volume of new ridges
      REAL(wp), DIMENSION(A2D(0),nlay_i) ::   eirdg     ! enth*volume of new ridges
      !
      INTEGER , DIMENSION(A2D(0)) ::   itest_rdg, itest_rft   ! test for conservation
      LOGICAL , DIMENSION(A2D(0)) ::   ll_shift         ! logical for doing calculation or not
      !!-------------------------------------------------------------------
      !
      ll_shift(:,:) = .false.

      zvti(:,:) = SUM( v_i(A2D(0),:), dim=3 )   ! total ice volume
      !
      ! 1) Change in open water area due to closing and opening
      !--------------------------------------------------------
      DO_2D( 0, 0, 0, 0 )
         IF( lice_pres(ji,jj) ) THEN
            ato_i(ji,jj) = MAX( 0._wp, ato_i(ji,jj) + ( popening(ji,jj) - papartf(ji,jj,0) * pclosing_gross(ji,jj) ) * rDt_ice )
         ENDIF
      END_2D

      DO jl1 = 1, jpl

         DO_2D( 0, 0, 0, 0 )
            !
            IF( lice_pres(ji,jj) ) THEN
               !
               ! set logical to true when ridging
               IF( papartf(ji,jj,jl1) > 0._wp .AND. pclosing_gross(ji,jj) > 0._wp ) THEN   ;   ll_shift(ji,jj) = .TRUE.
               ELSE                                                                        ;   ll_shift(ji,jj) = .FALSE.
               ENDIF

               IF( ll_shift(ji,jj) ) THEN   ! only if ice is ridging

                  IF( a_i(ji,jj,jl1) > epsi10 ) THEN   ;   z1_ai = 1._wp / a_i(ji,jj,jl1)
                  ELSE                                 ;   z1_ai = 0._wp
                  ENDIF

                  ! area of ridging / rafting ice (airdg1) and of new ridge (airdg2)
                  airdg1(ji,jj) = paridge(ji,jj,jl1) * pclosing_gross(ji,jj) * rDt_ice
                  airft1(ji,jj) = paraft (ji,jj,jl1) * pclosing_gross(ji,jj) * rDt_ice

                  airdg2(ji,jj) = airdg1(ji,jj) * phi_hrdg(ji,jj,jl1)
                  airft2(ji,jj) = airft1(ji,jj) *  hi_hrft

                  ! ridging /rafting fractions
                  afrdg = airdg1(ji,jj) * z1_ai
                  afrft = airft1(ji,jj) * z1_ai

                  ! volume and enthalpy (J/m2, >0) of seawater trapped into ridges
                  IF    ( zvti(ji,jj) <= 10. ) THEN ; vsw = v_i(ji,jj,jl1) * afrdg * rn_porordg                                           ! v <= 10m then porosity = rn_porordg
                  ELSEIF( zvti(ji,jj) >= 20. ) THEN ; vsw = 0._wp                                                                         ! v >= 20m then porosity = 0
                  ELSE                           ; vsw = v_i(ji,jj,jl1) * afrdg * rn_porordg * MAX( 0._wp, 2._wp - 0.1_wp * zvti(ji,jj) ) ! v > 10m and v < 20m then porosity = linear transition to 0
                  ENDIF
                  ersw = -rhoi * vsw * rcp * psst(ji,jj)   ! clem: if sst>0, then ersw <0 (is that possible?)

                  ! volume etc of ridging / rafting ice and new ridges (vi, vs, sm, oi, es, ei)
                  virdg(ji,jj) = v_i (ji,jj,jl1)    * afrdg + vsw
                  vsrdg(ji,jj) = v_s (ji,jj,jl1)    * afrdg
                  oirdg(ji,jj) = oa_i(ji,jj,jl1)    * afrdg * phi_hrdg(ji,jj,jl1)

                  virft(ji,jj) = v_i (ji,jj,jl1)    * afrft
                  vsrft(ji,jj) = v_s (ji,jj,jl1)    * afrft
                  oirft(ji,jj) = oa_i(ji,jj,jl1)    * afrft * hi_hrft

                  IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                     aprdg(ji,jj) = a_ip(ji,jj,jl1) * afrdg * phi_hrdg(ji,jj,jl1)
                     vprdg(ji,jj) = v_ip(ji,jj,jl1) * afrdg
                     aprft(ji,jj) = a_ip(ji,jj,jl1) * afrft *  hi_hrft
                     vprft(ji,jj) = v_ip(ji,jj,jl1) * afrft
                     IF ( ln_pnd_lids ) THEN
                        vlrdg(ji,jj) = v_il(ji,jj,jl1) * afrdg
                        vlrft(ji,jj) = v_il(ji,jj,jl1) * afrft
                     ENDIF
                  ENDIF

                  esrdg(ji,jj,:) = e_s(ji,jj,:,jl1) * afrdg
                  esrft(ji,jj,:) = e_s(ji,jj,:,jl1) * afrft
                  eirdg(ji,jj,:) = e_i(ji,jj,:,jl1) * afrdg + ersw * r1_nlay_i
                  eirft(ji,jj,:) = e_i(ji,jj,:,jl1) * afrft

                  IF( nn_icesal == 4 ) THEN
                     sirdg(ji,jj,:) = szv_i(ji,jj,:,jl1) * afrdg + vsw * psss(ji,jj) * r1_nlay_i
                     sirft(ji,jj,:) = szv_i(ji,jj,:,jl1) * afrft
                  ELSE
                     sirdg(ji,jj,1) = sv_i(ji, jj, jl1) * afrdg + vsw * psss(ji,jj)
                     sirft(ji,jj,1) = sv_i(ji, jj, jl1) * afrft
                  ENDIF

                  ! Ice-ocean exchanges associated with ice porosity
                  wfx_dyn(ji,jj) = wfx_dyn(ji,jj) - vsw * rhoi * r1_Dt_ice   ! increase in ice volume due to seawater frozen in voids
                  sfx_dyn(ji,jj) = sfx_dyn(ji,jj) - vsw * psss(ji,jj) * rhoi * r1_Dt_ice
                  hfx_dyn(ji,jj) = hfx_dyn(ji,jj) + ersw * r1_Dt_ice          ! > 0 [W.m-2]

                  ! Put the snow and pond lost by ridging into the ocean
                  !  Note that esrdg > 0; the ocean must cool to melt snow. If the ocean temp = Tf already, new ice must grow.
                  wfx_snw_dyn(ji,jj) = wfx_snw_dyn(ji,jj) + ( rhos * vsrdg(ji,jj) * ( 1._wp - rn_fsnwrdg )   &   ! fresh water source for ocean
                     &                                      + rhos * vsrft(ji,jj) * ( 1._wp - rn_fsnwrft ) ) * r1_Dt_ice
                  DO jk = 1, nlay_s
                     hfx_dyn(ji,jj) = hfx_dyn(ji,jj) + ( - esrdg(ji,jj,jk) * ( 1._wp - rn_fsnwrdg )   &          ! heat sink for ocean (<0, W.m-2)
                        &                                - esrft(ji,jj,jk) * ( 1._wp - rn_fsnwrft ) ) * r1_Dt_ice
                  END DO

                  IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                     wfx_pnd(ji,jj)    = wfx_pnd(ji,jj)   + ( rhow * vprdg(ji,jj) * ( 1._wp - rn_fpndrdg )   &   ! fresh water source for ocean
                        &                                   + rhow * vprft(ji,jj) * ( 1._wp - rn_fpndrft ) ) * r1_Dt_ice
                     IF ( ln_pnd_lids ) THEN
                        wfx_pnd(ji,jj) = wfx_pnd(ji,jj)   + ( rhow * vlrdg(ji,jj) * ( 1._wp - rn_fpndrdg )   &   ! fresh water source for ocean
                           &                                + rhow * vlrft(ji,jj) * ( 1._wp - rn_fpndrft ) ) * r1_Dt_ice
                     ENDIF
                  ENDIF

                  ! virtual salt flux to keep salinity constant
                  IF( nn_icesal == 1 .OR. nn_icesal == 3 )  THEN
                     sirdg(ji,jj,1) = sirdg(ji,jj,1) - ( psss(ji,jj) - s_i(ji,jj,jl1) ) * vsw                      ! ridge salinity = s_i
                     sfx_bri(ji,jj) = sfx_bri(ji,jj) + ( psss(ji,jj) - s_i(ji,jj,jl1) ) * vsw * rhoi * r1_Dt_ice   ! put back sss_m into the ocean
                     !                                                                                             ! and get  s_i  from the ocean
                  ENDIF

                  ! Remove area, volume of new ridge to each category jl1
                  !------------------------------------------------------
                  a_i (ji,jj,jl1) = a_i (ji,jj,jl1) - airdg1(ji,jj) - airft1(ji,jj)
                  v_i (ji,jj,jl1) = v_i (ji,jj,jl1)     * ( 1._wp - afrdg - afrft )
                  v_s (ji,jj,jl1) = v_s (ji,jj,jl1)     * ( 1._wp - afrdg - afrft )
                  oa_i(ji,jj,jl1) = oa_i(ji,jj,jl1)     * ( 1._wp - afrdg - afrft )
                  IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                     a_ip(ji,jj,jl1)  = a_ip(ji,jj,jl1) * ( 1._wp - afrdg - afrft )
                     v_ip(ji,jj,jl1)  = v_ip(ji,jj,jl1) * ( 1._wp - afrdg - afrft )
                     IF ( ln_pnd_lids ) v_il(ji,jj,jl1) = v_il(ji,jj,jl1) * ( 1._wp - afrdg - afrft )
                  ENDIF
                  !
                  e_s(ji,jj,:,jl1) = e_s(ji,jj,:,jl1)   * ( 1._wp - afrdg - afrft )
                  e_i(ji,jj,:,jl1) = e_i(ji,jj,:,jl1)   * ( 1._wp - afrdg - afrft )
                  !
                  IF( nn_icesal == 4 ) THEN   ;   szv_i(ji,jj,:,jl1) = szv_i(ji,jj,:,jl1) * ( 1._wp - afrdg - afrft )
                  ELSE                        ;   sv_i (ji,jj,  jl1) = sv_i (ji,jj,  jl1) * ( 1._wp - afrdg - afrft )
                  ENDIF

                  ! --- Ridging diagnostics ---!
                  IF( ll_diag_rdg ) THEN
                     diag_airdg1 (ji,jj) = diag_airdg1 (ji,jj) + airdg1 (ji,jj) * r1_Dt_ice
                     diag_airft1 (ji,jj) = diag_airft1 (ji,jj) + airft1 (ji,jj) * r1_Dt_ice
                     diag_airdg2 (ji,jj) = diag_airdg2 (ji,jj) + airdg2 (ji,jj) * r1_Dt_ice
                     diag_airft2 (ji,jj) = diag_airft2 (ji,jj) + airft2 (ji,jj) * r1_Dt_ice
                     diag_virdg  (ji,jj) = diag_virdg  (ji,jj) + virdg  (ji,jj) * r1_Dt_ice
                     diag_virft  (ji,jj) = diag_virft  (ji,jj) + virft  (ji,jj) * r1_Dt_ice
                  ENDIF

               ENDIF
               !
            ENDIF ! lice_pres
            !
         END_2D
         
 
         ! 2) compute categories in which ice is added (jl2)
         !--------------------------------------------------
         itest_rdg(:,:) = 0
         itest_rft(:,:) = 0
         DO jl2  = 1, jpl
            !
            DO_2D( 0, 0, 0, 0 )

               IF( lice_pres(ji,jj) .AND. ll_shift(ji,jj) ) THEN

                  ! Compute the fraction of ridged ice area and volume going to thickness category jl2
                  IF( ln_distf_lin ) THEN ! Hibler (1980) linear formulation
                     !
                     IF( phrmin(ji,jj,jl1) <= hi_max(jl2) .AND. phrmax(ji,jj,jl1) > hi_max(jl2-1) ) THEN
                        hL = MAX( phrmin(ji,jj,jl1), hi_max(jl2-1) )
                        hR = MIN( phrmax(ji,jj,jl1), hi_max(jl2)   )
                        farea = ( hR      - hL      ) / ( phrmax(ji,jj,jl1)                  - phrmin(ji,jj,jl1)                  )
                        fvol  = ( hR * hR - hL * hL ) / ( phrmax(ji,jj,jl1) * phrmax(ji,jj,jl1) - phrmin(ji,jj,jl1) * phrmin(ji,jj,jl1) )
                        !
                        itest_rdg(ji,jj) = 1   ! test for conservation
                     ELSE
                        farea = 0._wp
                        fvol  = 0._wp
                     ENDIF
                     !
                  ELSEIF( ln_distf_exp ) THEN ! Lipscomb et al. (2007) exponential formulation                      
                     !
                     IF( jl2 < jpl ) THEN
                        !
                        IF( phrmin(ji,jj,jl1) <= hi_max(jl2) ) THEN
                           hL    = MAX( phrmin(ji,jj,jl1), hi_max(jl2-1) )
                           hR    = hi_max(jl2)
                           expL  = EXP( -( hL - phrmin(ji,jj,jl1) ) / MAX( epsi20, phrexp(ji,jj,jl1) ) )
                           expR  = EXP( -( hR - phrmin(ji,jj,jl1) ) / MAX( epsi20, phrexp(ji,jj,jl1) ) )
                           farea = expL - expR
                           fvol  = ( ( hL + phrexp(ji,jj,jl1) ) * expL  &
                              &    - ( hR + phrexp(ji,jj,jl1) ) * expR ) / MAX( epsi20, phrmin(ji,jj,jl1) + phrexp(ji,jj,jl1) )
                        ELSE
                           farea = 0._wp
                           fvol  = 0._wp
                        END IF
                        !                 
                     ELSE             ! jl2 = jpl
                        !
                        hL    = MAX( phrmin(ji,jj,jl1), hi_max(jl2-1) )
                        expL  = EXP(-( hL - phrmin(ji,jj,jl1) ) / MAX( epsi20, phrexp(ji,jj,jl1) ) )
                        farea = expL
                        fvol  = ( hL + phrexp(ji,jj,jl1) ) * expL / MAX( epsi20, phrmin(ji,jj,jl1) + phrexp(ji,jj,jl1) )
                        !
                     END IF            ! jl2 < jpl
                     ! 
                     itest_rdg(ji,jj) = 1   ! test for conservation => clem: I am not sure about that
                     !
                  END IF             ! ridge redistribution
                     
                  ! Compute the fraction of rafted ice area and volume going to thickness category jl2
                  IF( phraft(ji,jj,jl1) <= hi_max(jl2) .AND. phraft(ji,jj,jl1) >  hi_max(jl2-1) ) THEN
                     zswitch = 1._wp
                     !
                     itest_rft(ji,jj) = 1   ! test for conservation
                  ELSE
                     zswitch = 0._wp
                  ENDIF
                  !
                  ! Patch to ensure perfect conservation if ice thickness goes mad
                  ! Sometimes thickness is larger than hi_max(jpl) because of advection scheme (for very small areas)
                  ! Then ice volume is removed from one category but the ridging/rafting scheme
                  ! does not know where to move it, leading to a conservation issue.
                  IF( itest_rdg(ji,jj) == 0 .AND. jl2 == jpl ) THEN   ;   farea = 1._wp   ;   fvol = 1._wp   ;   ENDIF
                  IF( itest_rft(ji,jj) == 0 .AND. jl2 == jpl ) zswitch = 1._wp
                  !
                  ! Add area, volume of new ridge to category jl2
                  !----------------------------------------------
                  a_i (ji,jj,jl2) = a_i (ji,jj,jl2) + ( airdg2(ji,jj)              * farea + airft2(ji,jj)              * zswitch )
                  oa_i(ji,jj,jl2) = oa_i(ji,jj,jl2) + ( oirdg (ji,jj)              * farea + oirft (ji,jj)              * zswitch )
                  v_i (ji,jj,jl2) = v_i (ji,jj,jl2) + ( virdg (ji,jj)              * fvol  + virft (ji,jj)              * zswitch )
                  v_s (ji,jj,jl2) = v_s (ji,jj,jl2) + ( vsrdg (ji,jj) * rn_fsnwrdg * fvol  + vsrft (ji,jj) * rn_fsnwrft * zswitch )
                  IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                     v_ip (ji,jj,jl2) = v_ip(ji,jj,jl2)    + ( vprdg(ji,jj) * rn_fpndrdg * fvol  + vprft(ji,jj) * rn_fpndrft * zswitch )
                     a_ip (ji,jj,jl2) = a_ip(ji,jj,jl2)    + ( aprdg(ji,jj) * rn_fpndrdg * farea + aprft(ji,jj) * rn_fpndrft * zswitch )
                     IF ( ln_pnd_lids ) THEN
                        v_il (ji,jj,jl2) = v_il(ji,jj,jl2) + ( vlrdg(ji,jj) * rn_fpndrdg * fvol  + vlrft(ji,jj) * rn_fpndrft * zswitch )
                     ENDIF
                  ENDIF
                  e_s(ji,jj,:,jl2) = e_s(ji,jj,:,jl2) + ( esrdg(ji,jj,:) * rn_fsnwrdg * fvol + esrft(ji,jj,:) * rn_fsnwrft * zswitch )
                  e_i(ji,jj,:,jl2) = e_i(ji,jj,:,jl2) + ( eirdg(ji,jj,:)              * fvol + eirft(ji,jj,:)              * zswitch )
                  IF( nn_icesal == 4 ) THEN
                     szv_i(ji,jj,:,jl2) = szv_i(ji,jj,:,jl2) + ( sirdg(ji,jj,:) * fvol + sirft(ji,jj,:) * zswitch )
                  ELSE
                     sv_i (ji,jj,  jl2) = sv_i (ji, jj, jl2) + ( sirdg(ji,jj,1) * fvol + sirft(ji,jj,1) * zswitch )
                  ENDIF
               ENDIF

            END_2D
            !
         END DO ! jl2
         !
      END DO ! jl1
      !
      ! roundoff errors
      !----------------
      ! In case ridging/rafting lead to very small negative values (sometimes it happens)
      CALL ice_var_roundoff( a_i, v_i, v_s, sv_i, oa_i, a_ip, v_ip, v_il, e_s, e_i, szv_i, lice_pres )

      !
   END SUBROUTINE rdgrft_shift


   SUBROUTINE ice_strength
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE ice_strength ***
      !!
      !! ** Purpose :   computes ice strength used in dynamics routines of ice thickness
      !!
      !! ** Method  :   Compute the strength of the ice pack, defined as the energy (J m-2) 
      !!              dissipated per unit area removed from the ice pack under compression,
      !!              and assumed proportional to the change in potential energy caused
      !!              by ridging. Note that ice strength using Hibler's formulation must be
      !!              smoothed.
      !!----------------------------------------------------------------------
      INTEGER             ::   ji, jj, jl  ! dummy loop indices
      REAL(wp)            ::   z1_3        ! local scalars
      REAL(wp), DIMENSION(A2D(0))      ::   zworka         ! temporary array used here
      REAL(wp), DIMENSION(jpi,jpj,jpl) ::   za_i_cap       ! local capped ice concentration
      !!
      LOGICAL             ::   ln_str_R75
      REAL(wp)            ::   zhi, zcp
      REAL(wp)            ::   h2rdg                     ! mean value of h^2 for new ridge
      REAL(wp), DIMENSION(A2D(0)) ::   zaksum              ! normalisation factor
      !!----------------------------------------------------------------------
      ! at_i needed for strength
      at_i(:,:) = SUM( a_i, dim=3 )
      !
      SELECT CASE( nice_str )          !--- Set which ice strength is chosen

      CASE ( np_strr75 )           !== Rothrock(1975)'s method ==!

         ! this should be defined once for all at the 1st time step
         zcp = 0.5_wp * grav * (rho0-rhoi) * rhoi * r1_rho0   ! proport const for PE
         !
         strength(:,:) = 0._wp
         !
         ! Initialise local capped a_i to zero
         ! Note that if 0 < a_i < epsi10, can end up with zhi=0 but apartf>0 rdgrft_prep; za_i_cap avoids this here
         DO jl = 1, jpl
           DO_2D( 0, 0, 0, 0 )
              za_i_cap(ji,jj,jl) = 0.0
           END_2D
         END DO
         !
         l_ice_present(A2D(0)) = .FALSE.
         ! Identify grid cells with ice
         DO_2D( 0, 0, 0, 0 )
            IF ( at_i(ji,jj) > epsi10 ) THEN
               l_ice_present(ji,jj) = .TRUE.
            ENDIF
         END_2D

         IF( ANY( l_ice_present(A2D(0)) ) ) THEN

            ! Cap a_i to avoid zhi in rdgrft_prep going below minimum
            za_i_cap(A2D(0),:) = a_i(A2D(0),:)
            DO jl = 1, jpl
               DO_2D( 0, 0, 0, 0 )
                  IF( l_ice_present(ji,jj) ) THEN
                     IF ( a_i(ji,jj,jl) > epsi10 .AND. ( v_i(ji,jj,jl) / a_i(ji,jj,jl) ) .LT. rn_himin ) THEN
                        za_i_cap(ji,jj,jl) = a_i(ji,jj,jl) * ( v_i(ji,jj,jl) / a_i(ji,jj,jl) ) / rn_himin
                     ENDIF
                  ENDIF
               END_2D
            END DO

            CALL rdgrft_prep( za_i_cap, v_i, ato_i, l_ice_present(A2D(0)), &                 ! <<== in
               &              apartf, aridge, araft, hi_hrdg, hraft, hrmin, hrmax, hrexp )   ! ==>> out
            !
            zaksum(:,:) = apartf(:,:,0) !clem: aksum should be defined in the header => local to module
            DO jl = 1, jpl
               DO_2D( 0, 0, 0, 0 )
                  IF( l_ice_present(ji,jj) ) THEN
                     IF ( apartf(ji,jj,jl) > 0._wp ) THEN
                        zaksum(ji,jj) = zaksum(ji,jj) + aridge(ji,jj,jl) * ( 1._wp - hi_hrdg(ji,jj,jl) )    &
                           &                          + araft (ji,jj,jl) * ( 1._wp - hi_hrft )
                     ENDIF
                  ENDIF
               END_2D
            END DO
            !
            z1_3 = 1._wp / 3._wp
            DO jl = 1, jpl
               DO_2D( 0, 0, 0, 0 )
                  IF( l_ice_present(ji,jj) ) THEN
                     !
                     IF( apartf(ji,jj,jl) > 0._wp ) THEN
                        !
                        IF( ln_distf_lin ) THEN       ! Uniform redistribution of ridged ice
                           h2rdg = z1_3 * ( hrmax(ji,jj,jl) * hrmax(ji,jj,jl) +     & ! (a**3-b**3)/(a-b) = a*a+ab+b*b
                              &             hrmin(ji,jj,jl) * hrmin(ji,jj,jl) +     &
                              &             hrmax(ji,jj,jl) * hrmin(ji,jj,jl) )
                           !
                        ELSEIF( ln_distf_exp ) THEN   ! Exponential redistribution of ridged ice
                           h2rdg =          hrmin(ji,jj,jl) * hrmin(ji,jj,jl)   &
                              &   + 2._wp * hrmin(ji,jj,jl) * hrexp(ji,jj,jl)   &
                              &   + 2._wp * hrexp(ji,jj,jl) * hrexp(ji,jj,jl)
                        END IF
                        !
                        IF( a_i(ji,jj,jl) > epsi10 ) THEN   ;   zhi = v_i(ji,jj,jl) / a_i(ji,jj,jl)
                        ELSE                                ;   zhi = 0._wp
                        ENDIF

                        ! Make sure ice thickness is not below the minimum
                        ! Do not adjust concentration as don't want strength routine to be able to do this
                        IF( a_i(ji,jj,jl) > epsi10 .AND. zhi < rn_himin ) THEN
                           zhi = rn_himin
                        ENDIF


!!$                     zstrength(ii) = zstrength(ii) -         apartf(ii,jl) * zhi * zhi                  ! PE loss from deforming ice
!!$                     zstrength(ii) = zstrength(ii) + 2._wp * araft (ii,jl) * zhi * zhi                  ! PE gain from rafting ice
!!$                     zstrength(ii) = zstrength(ii) +         aridge(ii,jl) * hi_hrdg(ii,jl) * z1_3 *  & ! PE gain from ridging ice
!!$                        &                                   ( hrmax(ii,jl) * hrmax(ii,jl) +           & ! (a**3-b**3)/(a-b) = a*a+ab+b*b
!!$                        &                                     hrmin(ii,jl) * hrmin(ii,jl) +           &
!!$                        &                                     hrmax(ii,jl) * hrmin(ii,jl) )
                        strength(ji,jj) = strength(ji,jj) - apartf(ji,jj,jl) * zhi * zhi                  ! PE loss
                        strength(ji,jj) = strength(ji,jj) + 2._wp * araft(ji,jj,jl) * zhi * zhi           ! PE gain (rafting)
                        strength(ji,jj) = strength(ji,jj) + aridge(ji,jj,jl) * h2rdg *  hi_hrdg(ji,jj,jl)    ! PE gain (ridging)

                     ENDIF
                     !
                  ENDIF
               END_2D
            END DO
            !
            DO_2D( 0, 0, 0, 0 )
               IF( l_ice_present(ji,jj) ) THEN
                  strength(ji,jj) = rn_pe_rdg * zcp * strength(ji,jj) / zaksum(ji,jj)
               ENDIF
            END_2D
            !
         ENDIF

         CALL lbc_lnk( 'icedyn_rdgrft', strength, 'T', 1.0_wp ) ! this call could be removed if calculations were done on the full domain
         !                                                      ! but we decided it is more efficient this way
         !
      CASE ( np_strh79 )           !== Hibler(1979)'s method ==!
         !
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            IF( at_i(ji,jj) > epsi10 ) THEN
               strength(ji,jj) = rn_pstar * SUM( v_i(ji,jj,:) ) * EXP( -rn_crhg * ( 1._wp - at_i(ji,jj) ) )
            ELSE
               strength(ji,jj) = 0._wp
            ENDIF
         END_2D
         !
      CASE ( np_strcst )           !== Constant strength ==!
         !
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            IF( at_i(ji,jj) > epsi10 ) THEN
               strength(ji,jj) = rn_str
            ELSE
               strength(ji,jj) = 0._wp
            ENDIF
         END_2D
         !
      END SELECT
      !
      IF( ln_str_smooth ) THEN         !--- Spatial smoothing
         DO_2D( 0, 0, 0, 0 )
            IF( at_i(ji,jj) > epsi10 ) THEN
               zworka(ji,jj) = ( 4._wp * strength(ji,jj)              &
                  &                    + ( ( strength(ji-1,jj) * tmask(ji-1,jj,1) + strength(ji+1,jj) * tmask(ji+1,jj,1) ) &
                  &                      + ( strength(ji,jj-1) * tmask(ji,jj-1,1) + strength(ji,jj+1) * tmask(ji,jj+1,1) ) ) &
                  &            ) / ( 4._wp + tmask(ji-1,jj,1) + tmask(ji+1,jj,1) + tmask(ji,jj-1,1) + tmask(ji,jj+1,1) )
            ELSE
               zworka(ji,jj) = 0._wp
            ENDIF
         END_2D

         DO_2D( 0, 0, 0, 0 )
            strength(ji,jj) = zworka(ji,jj)
         END_2D
         CALL lbc_lnk( 'icedyn_rdgrft', strength, 'T', 1.0_wp )
         !
      ENDIF
      !
   END SUBROUTINE ice_strength


   SUBROUTINE ice_dyn_rdgrft_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_rdgrft_init ***
      !!
      !! ** Purpose :   Physical constants and parameters linked
      !!                to the mechanical ice redistribution
      !!
      !! ** Method  :   Read the namdyn_rdgrft namelist
      !!                and check the parameters values
      !!                called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namdyn_rdgrft
      !!-------------------------------------------------------------------
      INTEGER :: ios, ioptio                ! Local integer output status for namelist read
      !!
      NAMELIST/namdyn_rdgrft/ ln_str_H79, rn_pstar, rn_crhg, ln_str_R75, rn_pe_rdg, ln_str_CST, rn_str, ln_str_smooth, &
         &                    ln_distf_lin, ln_distf_exp, rn_murdg, rn_csrdg,            &
         &                    ln_partf_lin, rn_gstar, ln_partf_exp, rn_astar,            &
         &                    ln_ridging, rn_hstar, rn_porordg, rn_fsnwrdg, rn_fpndrdg,  &
         &                    ln_rafting, rn_hraft, rn_craft  , rn_fsnwrft, rn_fpndrft
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namdyn_rdgrft)
      READ_NML_CFG(numnam_ice,namdyn_rdgrft)
      IF(lwm) WRITE ( numoni, namdyn_rdgrft )
      !
      IF (lwp) THEN                          ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_rdgrft_init: ice parameters for ridging/rafting '
         WRITE(numout,*) '~~~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_rdgrft:'
         WRITE(numout,*) '      ice strength parameterization Hibler (1979)              ln_str_H79   = ', ln_str_H79
         WRITE(numout,*) '            1st bulk-rheology parameter                        rn_pstar     = ', rn_pstar
         WRITE(numout,*) '            2nd bulk-rhelogy parameter                         rn_crhg      = ', rn_crhg
         WRITE(numout,*) '      ice strength parameterization Rothrock (1975)            ln_str_R75   = ', ln_str_R75
         WRITE(numout,*) '            coef accounting for frictional dissipation         rn_pe_rdg    = ', rn_pe_rdg         
         WRITE(numout,*) '      ice strength parameterization Constant                   ln_str_CST   = ', ln_str_CST
         WRITE(numout,*) '            ice strength value                                 rn_str       = ', rn_str
         WRITE(numout,*) '      spatial smoothing of the strength                        ln_str_smooth= ', ln_str_smooth
         WRITE(numout,*) '      redistribution of ridged ice: linear (Hibler 1980)       ln_distf_lin = ', ln_distf_lin
         WRITE(numout,*) '      redistribution of ridged ice: exponential(Lipscomb 2017) ln_distf_exp = ', ln_distf_exp
         WRITE(numout,*) '            e-folding scale of ridged ice                      rn_murdg     = ', rn_murdg
         WRITE(numout,*) '      Fraction of shear energy contributing to ridging         rn_csrdg     = ', rn_csrdg
         WRITE(numout,*) '      linear ridging participation function                    ln_partf_lin = ', ln_partf_lin
         WRITE(numout,*) '            Fraction of ice coverage contributing to ridging   rn_gstar     = ', rn_gstar
         WRITE(numout,*) '      Exponential ridging participation function               ln_partf_exp = ', ln_partf_exp
         WRITE(numout,*) '            Equivalent to G* for an exponential function       rn_astar     = ', rn_astar
         WRITE(numout,*) '      Ridging of ice sheets or not                             ln_ridging   = ', ln_ridging
         WRITE(numout,*) '            max ridged ice thickness                           rn_hstar     = ', rn_hstar
         WRITE(numout,*) '            Initial porosity of ridges                         rn_porordg   = ', rn_porordg
         WRITE(numout,*) '            Fraction of snow volume conserved during ridging   rn_fsnwrdg   = ', rn_fsnwrdg
         WRITE(numout,*) '            Fraction of pond volume conserved during ridging   rn_fpndrdg   = ', rn_fpndrdg
         WRITE(numout,*) '      Rafting of ice sheets or not                             ln_rafting   = ', ln_rafting
         WRITE(numout,*) '            Parmeter thickness (threshold between ridge-raft)  rn_hraft     = ', rn_hraft
         WRITE(numout,*) '            Rafting hyperbolic tangent coefficient             rn_craft     = ', rn_craft
         WRITE(numout,*) '            Fraction of snow volume conserved during rafting   rn_fsnwrft   = ', rn_fsnwrft
         WRITE(numout,*) '            Fraction of pond volume conserved during rafting   rn_fpndrft   = ', rn_fpndrft
      ENDIF
      !
      ioptio = 0
      IF( ln_str_H79    ) THEN   ;   ioptio = ioptio + 1   ;   nice_str = np_strh79       ;   ENDIF
      IF( ln_str_R75    ) THEN   ;   ioptio = ioptio + 1   ;   nice_str = np_strr75       ;   ENDIF
      IF( ln_str_CST    ) THEN   ;   ioptio = ioptio + 1   ;   nice_str = np_strcst       ;   ENDIF
      IF( ioptio /= 1 )   &
         &   CALL ctl_stop( 'ice_dyn_rdgrft_init: one and only one ice strength option has to be defined ' )
      !
      ioptio = 0
      IF( ln_distf_lin ) THEN   ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ln_distf_exp ) THEN   ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ioptio /= 1 )   &
         &   CALL ctl_stop( 'ice_dyn_rdgrft_init: choose one and only one redistribution function (ln_distf_lin or ln_distf_exp)' )
      !
      ioptio = 0
      IF( ln_partf_lin ) THEN   ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ln_partf_exp ) THEN   ;   ioptio = ioptio + 1   ;   ENDIF
      IF( ioptio /= 1 )   &
         &   CALL ctl_stop( 'ice_dyn_rdgrft_init: choose one and only one participation function (ln_partf_lin or ln_partf_exp)' )
      !
      IF( .NOT. ln_icethd ) THEN
         rn_porordg = 0._wp
         rn_fsnwrdg = 1._wp ; rn_fsnwrft = 1._wp
         rn_fpndrdg = 1._wp ; rn_fpndrft = 1._wp
         IF( lwp ) THEN
            WRITE(numout,*) '      ==> only ice dynamics is activated, thus some parameters must be changed'
            WRITE(numout,*) '            rn_porordg   = ', rn_porordg
            WRITE(numout,*) '            rn_fsnwrdg   = ', rn_fsnwrdg
            WRITE(numout,*) '            rn_fpndrdg   = ', rn_fpndrdg
            WRITE(numout,*) '            rn_fsnwrft   = ', rn_fsnwrft
            WRITE(numout,*) '            rn_fpndrft   = ', rn_fpndrft
         ENDIF
      ENDIF
      !                              ! allocate arrays
      IF( ice_dyn_rdgrft_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'ice_dyn_rdgrft_init: unable to allocate arrays' )
      !
  END SUBROUTINE ice_dyn_rdgrft_init

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icedyn_rdgrft
