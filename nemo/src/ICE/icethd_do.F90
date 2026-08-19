MODULE icethd_do
   !!======================================================================
   !!                       ***  MODULE icethd_do   ***
   !!   sea-ice: sea ice growth in the leads (open water)  
   !!======================================================================
   !! History :       !  2005-12  (M. Vancoppenolle) Original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_thd_do        : ice growth in open water (=lateral accretion of ice)
   !!   ice_thd_do_init   : initialization
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE par_oce
   USE dom_oce , ONLY : umask, vmask, smask0
   USE phycst         ! physical constants
   USE ice            ! sea-ice: variables
   USE sbc_oce , ONLY : sss_m
   USE sbcwave , ONLY : hsw, wpf
   USE sbc_ice , ONLY : utau_ice, vtau_ice
   USE icectl         ! sea-ice: conservation
   USE icevar  , ONLY : ice_var_vremap
   USE icethd_sal     ! sea-ice: salinity profiles
   USE icefsd  , ONLY : ice_fsd_part_newice, ice_fsd_add_newice, ice_fsd_thd, ice_fsd_weld, ice_fsd_dia
   USE icefsd  , ONLY : a_ifsd, nf_newice
   USE icewav  , ONLY : ice_wav_newice
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd_do        ! called by ice_thd
   PUBLIC   ice_thd_frazil    ! called by ice_thd
   PUBLIC   ice_thd_do_init   ! called by ice_stp
   !
   !                             !!** namelist (namthd_do) **
   REAL(wp) ::   rn_hinew         !  thickness for new ice formation (m)
   LOGICAL  ::   ln_frazil        !  use of frazil ice collection as function of wind (T) or not (F)
   REAL(wp) ::   rn_maxfraz       !  maximum portion of frazil ice collecting at the ice bottom
   REAL(wp) ::   rn_vfraz         !  threshold drift speed for collection of bottom frazil ice
   REAL(wp) ::   rn_Cfraz         !  squeezing coefficient for collection of bottom frazil ice

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd_do
      !!-------------------------------------------------------------------
      !!               ***   ROUTINE ice_thd_do  ***
      !!  
      !! ** Purpose : Computation of the evolution of the ice thickness and 
      !!              concentration as a function of the heat balance in the leads
      !!       
      !! ** Method  : Ice is formed in the open water when ocean looses heat
      !!              (heat budget of open water is negative) following
      !!
      !!       (dA/dt)acc = F[ (1-A)/(1-a) ] * [ Bl / (Li*h0) ]
      !!          where - h0 is the thickness of ice created in the lead
      !!                - a is a minimum fraction for leads
      !!                - F is a monotonic non-increasing function defined as:
      !!                  F(X)=( 1 - X**exld )**(1.0/exld)
      !!                - exld is the exponent closure rate (=2 default val.)
      !! 
      !! ** Action : - Adjustment of snow and ice thicknesses and heat
      !!                content in brine pockets
      !!             - Updating ice internal temperature
      !!             - Computation of variation of ice volume and mass
      !!             - Computation of a_i after lateral accretion and 
      !!               update h_s, h_i      
      !!------------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      !
      REAL(wp) ::   ztmelts
      REAL(wp) ::   zdE
      REAL(wp) ::   zQm          ! enthalpy exchanged with the ocean (J/m2, >0 towards ocean)
      REAL(wp) ::   zEi          ! sea ice specific enthalpy (J/kg)
      REAL(wp) ::   zEw          ! seawater specific enthalpy (J/kg)
      REAL(wp) ::   zfmdt        ! mass flux x time step (kg/m2, >0 towards ocean)
      !
      INTEGER  ::   jcat        ! indexes of categories where new ice grows
      INTEGER  ::   npti
      !
      REAL(wp) ::   zv_newfra
      REAL(wp) ::   zv_newice   ! volume of accreted ice
      REAL(wp) ::   za_newice   ! fractional area of accreted ice
      REAL(wp) ::   ze_newice   ! heat content of accreted ice
      REAL(wp) ::   zo_newice   ! age of accreted ice
      REAL(wp) ::   zdv_res     ! residual volume in case of excessive heat budget
      REAL(wp) ::   zda_res     ! residual area in case of excessive heat budget
      REAL(wp) ::   zv_frazb    ! accretion of frazil ice at the ice bottom
      !
      REAL(wp), DIMENSION(jpl) ::   zv_b    ! old volume of ice in category jl
      REAL(wp), DIMENSION(jpl) ::   za_b    ! old area of ice in category jl
      !
      REAL(wp), DIMENSION(A2D(0)) ::   zs_newice     ! salinity of accreted ice
      !
      REAL(wp), DIMENSION(0:nlay_i+1) ::   zh_i_old, ze_i_old, zs_i_old
      !
      ! For floe size distribution:
      REAL(wp), DIMENSION(jpl) ::   zda_latgro        ! fsd: change in ice area due to lateral growth in cat. jl
      REAL(wp), DIMENSION(jpl) ::   zv_latgro_cat     ! fsd: lateral growth volume in cat. jl
      REAL(wp)                 ::   zv_latgro         ! fsd: lateral growth volume, total
      REAL(wp)                 ::   zv_newice_total   ! fsd: new ice + lat. growth, used when updating e_i and szv_i
      INTEGER                  ::   jcat_fsd          ! fsd: new ice floe size category
      !
      REAL(wp), DIMENSION(A2D(0),nn_nfsd,jpl) ::   za_ifsdb_latgro   ! a_ifsd before lateral growth, for diagnostics
      REAL(wp), DIMENSION(A2D(0),nn_nfsd,jpl) ::   za_ifsda_latgro   ! "      after  lateral growth  "
      REAL(wp), DIMENSION(A2D(0),nn_nfsd,jpl) ::   za_ifsda_newice   ! "      after  new ice growth  "
      REAL(wp), DIMENSION(A2D(0),jpl)         ::   za_ib_latgro      ! a_i before lateral growth, for diagnostics
      REAL(wp), DIMENSION(A2D(0),jpl)         ::   za_ia_latgro      ! "   after  lateral growth  "
      REAL(wp), DIMENSION(A2D(0),jpl)         ::   za_ia_newice      ! "   after  new ice growth  "
      !
      !!-----------------------------------------------------------------------!
      !
      IF( ln_timing    )   CALL timing_start('icethd_do')
      IF( ln_icediachk )   CALL ice_cons_hsm( 0, 'icethd_do', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft )
      IF( ln_icediachk )   CALL ice_cons2D  ( 0, 'icethd_do',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft )

      !------------------------------------------------------------------------------!
      ! 1) Compute thickness, salinity, enthalpy, age, area and volume of new ice
      !------------------------------------------------------------------------------!
      ! it occurs if cooling
      at_i(A2D(0)) = SUM( a_i(A2D(0),:), dim=3 )

      ! Identify grid points where new ice forms
      npti = 0
      DO_2D( 0, 0, 0, 0 )
         IF ( qlead(ji,jj)  <  0._wp ) THEN
            npti = npti + 1
         ENDIF
      END_2D

      IF( ln_fsd ) THEN
         ! Prepare before/after arrays for diagnostics computed externally in ice_fsd_dia:
         ! (set all equal to initial values as some grid cells may not change):
         za_ifsdb_latgro(A2D(0),:,:) = a_ifsd(A2D(0),:,:)
         za_ifsda_latgro(A2D(0),:,:) = a_ifsd(A2D(0),:,:)
         za_ifsda_newice(A2D(0),:,:) = a_ifsd(A2D(0),:,:)
         za_ib_latgro(A2D(0),:) = a_i(A2D(0),:)
         za_ia_latgro(A2D(0),:) = a_i(A2D(0),:)
         za_ia_newice(A2D(0),:) = a_i(A2D(0),:)
      ENDIF

      IF ( npti > 0 ) THEN

         ! Convert units for ice internal energy and salt content
         DO jl = 1, jpl
            DO jk = 1, nlay_i
               DO_2D(0, 0, 0, 0)
                  IF (qlead(ji,jj) < 0._wp) THEN
                     IF (v_i(ji,jj,jl) > 0._wp) THEN
                        e_i(ji,jj,jk,jl) = e_i(ji,jj,jk,jl) / v_i(ji,jj,jl) * REAL( nlay_i )
                        szv_i(ji,jj,jk,jl) = szv_i(ji,jj,jk,jl) / v_i(ji,jj,jl) * REAL( nlay_i )
                     ELSE
                        e_i(ji,jj,jk,jl) = 0._wp
                        szv_i(ji,jj,jk,jl) = 0._wp
                     ENDIF   
                  ENDIF
               END_2D
            END DO
         END DO

         ! --- Salinity of new ice --- ! 
         SELECT CASE ( nn_icesal )
         CASE ( 1 )                    ! Sice = constant 
            zs_newice(:,:) = rn_icesal
         CASE ( 2 , 4 )                ! Sice = F(z,t) [Griewank and Notz 2013 ; Rees Jones and Worster 2014]
             DO_2D(0, 0, 0, 0)
                IF (qlead(ji,jj) < 0._wp) THEN
                   zs_newice(ji,jj) = rn_sinew * sss_m(ji,jj)
                ENDIF
             END_2D
         CASE ( 3 )                    ! Sice = F(z) [multiyear ice]
            zs_newice(:,:) =   2.3
         END SELECT
         !
         !                       ! ==================== !
         !                       ! Start main loop here !
         !                       ! ==================== !
         DO_2D(0,0,0,0)
            IF( qlead(ji,jj) < 0._wp ) THEN ! qlead is the heat budget in the first ocean level. Only grow ice when it is negative
               ! Keep old ice areas and volume in memory
               DO jl = 1, jpl
                  zv_b(jl) = v_i(ji,jj,jl) 
                  za_b(jl) = a_i(ji,jj,jl)
               ENDDO
            
               ! --- Heat content of new ice --- !
               ! We assume that new ice is formed at the seawater freezing point
               ztmelts   = - rTmlt * zs_newice(ji,jj)                  ! Melting point (C)
               ze_newice =   rhoi * (  rcpi  * ( ztmelts - ( t_bo(ji,jj) - rt0 ) )                               &
                  &                  + rLfus * MAX( 0._wp, 1._wp - ztmelts / MIN( t_bo(ji,jj) - rt0, -epsi10 ) ) & ! clem: max to deal with different eq. freezing in ice and ocean (but useless for now)
                  &                  - rcp   * ztmelts )
            
               ! --- Age of new ice --- !
               zo_newice = 0._wp

               ! --- Volume of new ice --- !
               zEi           = - ze_newice * r1_rhoi                  ! specific enthalpy of forming ice [J/kg]

               zEw           = rcp * ( t_bo(ji,jj) - rt0 )            ! specific enthalpy of seawater at t_bo [J/kg]
                                                                      ! clem: we suppose we are already at the freezing point (condition qlead<0 is satisfyied) 
                                                                   
               zdE           = zEi - zEw                              ! specific enthalpy difference [J/kg] (<0)
                                              
               zfmdt         = - qlead(ji,jj) / zdE                   ! Fm.dt [kg/m2] (<0) 
                                                                      ! clem: we use qlead instead of zqld (icethd) because we suppose we are at the freezing point   
               zv_newice     = - zfmdt * r1_rhoi

               zQm           = zfmdt * zEw                            ! heat to the ocean >0 associated with mass flux  

               ! Contribution to heat flux to the ocean [W.m-2], >0  
               hfx_thd(ji,jj) = hfx_thd(ji,jj) + zfmdt * zEw * r1_Dt_ice
               ! Total heat flux used in this process [W.m-2]  
               hfx_opw(ji,jj) = hfx_opw(ji,jj) - zfmdt * zdE * r1_Dt_ice
               ! mass flux
               wfx_opw(ji,jj) = wfx_opw(ji,jj) - zv_newice * rhoi * r1_Dt_ice
               ! salt flux
               sfx_opw(ji,jj) = sfx_opw(ji,jj) - zv_newice * rhoi * zs_newice(ji,jj) * r1_Dt_ice
         
               IF( ln_fsd ) THEN
                  ! --- floe size distribution --- !
                  !
                  ! Partition new ice growth (zv_newice) into open water new ice
                  ! growth and lateral growth at floe edges. The latter is
                  ! assigned to zv_latgro, zv_newice is updated accordingly, then
                  ! the latter is treated as usual regardless of ln_fsd:
                  !
                  CALL ice_fsd_part_newice( za_b(:), zv_b(:), a_ifsd(ji,jj,:,:), zv_newice, zv_latgro, zda_latgro )
                  !
               ELSE
                  zv_latgro     = 0._wp
                  zda_latgro(:) = 0._wp   ! area changes due to lateral growth
               ENDIF

               ! Lateral growth volume per category is calculated during the loop
               ! below where they are added to v_i in place, but will need to
               ! save them anyway (to this array) for later update of e_i:
               zv_latgro_cat(:) = 0._wp

               ! A fraction fraz_frac of frazil ice is accreted at the ice bottom
               IF( at_i(ji,jj) > 0._wp ) THEN
                  zv_frazb  =           fraz_frac(ji,jj)   * zv_newice
                  zv_newice = ( 1._wp - fraz_frac(ji,jj) ) * zv_newice
               ELSE
                  zv_frazb  = 0._wp
               ENDIF
               ! --- Area of new ice --- !
               za_newice = zv_newice / ht_i_new(ji,jj)

               ! --- Redistribute new ice area and volume into ice categories --- !

               ! --- lateral ice growth --- !
               ! If lateral ice growth gives an ice concentration > amax, then
               ! we keep the excessive volume in memory and attribute it later to bottom accretion
               IF ( za_newice > MAX( 0._wp, rn_amax_2d(ji,jj) - at_i(ji,jj) - SUM(zda_latgro(:)) ) ) THEN ! max is for roundoff error
                  zda_res   = za_newice - MAX( 0._wp, rn_amax_2d(ji,jj) - at_i(ji,jj) - SUM(zda_latgro(:)) )
                  zdv_res   = zda_res * ht_i_new(ji,jj) 
                  za_newice = MAX( 0._wp, za_newice - zda_res )
                  zv_newice = MAX( 0._wp, zv_newice - zdv_res )
               ELSE
                  zda_res = 0._wp
                  zdv_res = 0._wp
               ENDIF

               ! find which category to fill
               at_i(ji,jj) = 0._wp
               DO jl = 1, jpl
                  IF( ht_i_new(ji,jj) > hi_max(jl-1) .AND. ht_i_new(ji,jj) <= hi_max(jl) ) THEN
                     a_i(ji,jj,jl) = a_i(ji,jj,jl) + za_newice
                     v_i(ji,jj,jl) = v_i(ji,jj,jl) + zv_newice
                     jcat = jl
                  ENDIF

                  ! --- floe size distribution --- !
                  !
                  ! Lateral growth of existing ice in all thickness categories.
                  ! FSD is updated with new ice in ice_fsd_add_newice, called later.
                  ! Note if ln_fsd = .false. then zda_latgro(:) = 0.
                  !
                  IF( zda_latgro(jl) > 0._wp ) THEN
                     !
                     a_i(ji,jj,jl) = a_i(ji,jj,jl) + zda_latgro(jl)
                     !
                     IF( a_i(ji,jj,jl) > 0._wp ) THEN
                        !
                        ! Lateral growth volume for this thickness cat. (save for updating e_i later):
                        ! NOTE: use zv_b/za_b, not v_i/a_i: latter already updated above with
                        ! new ice for one of the categories!
                        !
                        zv_latgro_cat(jl) = zda_latgro(jl) * zv_b(jl) / za_b(jl)
                        v_i(ji,jj,jl) = v_i(ji,jj,jl) + zv_latgro_cat(jl)
                     ENDIF
                     !
                     ! Update FSD due to lateral growth:
                     CALL ice_fsd_thd( a_ifsd(ji,jj,:,jl), zv_latgro / rDt_ice )
                     !
                  ENDIF
                  !
                  ! Save FSD and a_i after lateral growth for diagnostics:
                  IF( ln_fsd ) THEN
                     za_ifsda_latgro(ji,jj,:,jl) = a_ifsd(ji,jj,:,jl)
                     !
                     ! CAREFUL: a_i has now evolved due to new ice AND lateral growth. So, get what it
                     ! would be only due to lateral growth from za_b (start of time step) and zda_latgro:
                     za_ia_latgro(ji,jj,jl) = za_b(jl) + zda_latgro(jl)
                  ENDIF
                  ! ------------------------------ !

                  at_i(ji,jj) = at_i(ji,jj) + a_i(ji,jj,jl)
               END DO

               ! --- floe size distribution --- !
               !
               ! For new ice cat (jcat), this needs to be AFTER lateral growth of FSD
               ! (i.e., subroutine ice_thd_evolve). It also requires a_i BEFORE new ice
               ! growth, but AFTER lateral growth, which are both done above but we can
               ! recover the correct value using za_b (a_i at beginning of ice_thd_do)
               ! and zda_latgro (FSD lateral area growth per thickness category):
               !
               IF( ln_fsd ) THEN
                  !
                  ! Floe size category that new ice is added to, jcat_fsd, can be modified by
                  ! ocean waves if ln_ice_wav=T, else it is set in FSD module (nf_newice):
                  !
                  jcat_fsd = nf_newice
                  IF( ln_ice_wav ) CALL ice_wav_newice( hsw(ji,jj), wpf(ji,jj), jcat_fsd )
                  !
                  CALL ice_fsd_add_newice( a_ifsd(ji,jj,:,jcat)         , za_newice,   &
                     &                     za_b(jcat) + zda_latgro(jcat), jcat_fsd     )
                  !
                  ! Save FSD and a_i after adding new ice for diagnostics
                  !
                  ! Note diagnostics are done in the order that FSD is updated (latgro then newice, and later welding)
                  ! za_ifsda_newice and za_ia_newice really mean 'after new ice AND lateral growth'
                  ! Both variables at their current state now include both processes:
                  za_ifsda_newice(ji,jj,:,:) = a_ifsd(ji,jj,:,:)
                  za_ia_newice   (ji,jj,  :) = a_i   (ji,jj,  :)
               ENDIF

               ! Heat content
               !
               ! With floe size distribution, we have added new ice area to all categories
               ! --> update enthalpy (e_i) and salinity content (szv_i) in each category
               !     using zv_latgro_cat calculated above.
               !
               ! Without floe size distribution, we only add new ice area to category jcat
               ! --> update in category jcat only; other jl in loop below therefore does
               !     nothing, as zv_latgro_cat will be 0, recovering original implementation
               !     prior to adding FSD.
               !
               DO jl = 1, jpl
                  ! Total new ice volume added laterally (from FSD) and from new ice (if jl == jcat):
                  zv_newice_total = zv_latgro_cat(jl)
                  IF( jl == jcat ) zv_newice_total = zv_newice_total + zv_newice
                  !
                  IF( zv_newice_total > 0._wp ) THEN
                     IF( za_b(jl) > 0._wp ) THEN
                        e_i(ji,jj,:,jl) = ( ze_newice     * zv_newice_total + e_i(ji,jj,:,jl) * zv_b(jl) ) / MAX( v_i(ji,jj,jl), epsi20 )
                        szv_i(ji,jj,:,jl) = ( zs_newice(ji,jj) * zv_newice_total + szv_i(ji,jj,:,jl) * zv_b(jl) ) / MAX( v_i(ji,jj,jl), epsi20 )
                     ELSE
                        e_i(ji,jj,:,jl) = ze_newice
                        szv_i(ji,jj,:,jl) = zs_newice(ji,jj)
                     ENDIF
                  ENDIF
               ENDDO

               ! --- bottom ice growth + ice enthalpy remapping + FSD floe welding --- !
               DO jl = 1, jpl
               
                  ! for remapping
                  zh_i_old(0:nlay_i+1) = 0._wp
                  ze_i_old(0:nlay_i+1) = 0._wp
                  zs_i_old(0:nlay_i+1) = 0._wp
                  DO jk = 1, nlay_i
                     zh_i_old(jk) =                      v_i(ji,jj,jl) * r1_nlay_i
                     ze_i_old(jk) = e_i(ji,jj,jk,jl) * v_i(ji,jj,jl) * r1_nlay_i
                     zs_i_old(jk) = szv_i(ji,jj,jk,jl) * v_i(ji,jj,jl) * r1_nlay_i
                  END DO

                  ! new volumes including lateral/bottom accretion + residual
                  IF( at_i(ji,jj) >= epsi20 ) THEN
                     zv_newfra     = ( zdv_res + zv_frazb ) * a_i(ji,jj,jl) / MAX( at_i(ji,jj) , epsi20 )
                  ELSE                  
                     zv_newfra     = 0._wp
                     a_i(ji,jj,jl) = 0._wp
                  ENDIF
                  v_i(ji,jj,jl) = v_i(ji,jj,jl) + zv_newfra
                  ! for remapping
                  zh_i_old(nlay_i+1) = zv_newfra
                  ze_i_old(nlay_i+1) = ze_newice     * zv_newfra
                  zs_i_old(nlay_i+1) = zs_newice(ji,jj) * zv_newfra
           
                  ! --- Update bulk salinity --- !
                  sv_i(ji,jj,jl) = sv_i(ji,jj,jl) + zs_newice(ji,jj) * ( v_i(ji,jj,jl) - zv_b(jl) )
              
                  ! --- Ice enthalpy and salt remapping --- !
                                         CALL ice_var_vremap( zh_i_old, ze_i_old, e_i(ji,jj,:,jl) ) 
                  IF( nn_icesal == 4 )   CALL ice_var_vremap( zh_i_old, zs_i_old, szv_i(ji,jj,:,jl) ) 

                  ! --- Floe welding (only changes FSD) --- !
                  IF( ln_fsd ) CALL ice_fsd_weld( a_ifsd(ji,jj,:,jl), a_i(ji,jj,jl) )
                  !
               END DO
            ENDIF ! qlead < 0   
         END_2D
         !                       ! ================== !
         !                       ! End main loop here !
         !                       ! ================== !
         !
         ! Change units for e_i/szv_i
         DO jl = 1, jpl
            DO jk = 1, nlay_i
               DO_2D(0, 0, 0, 0)
                  IF (qlead(ji,jj) < 0._wp) THEN
                     e_i(ji,jj,jk,jl) = e_i(ji,jj,jk,jl) * v_i(ji,jj,jl) * r1_nlay_i 
                     szv_i(ji,jj,jk,jl) = szv_i(ji,jj,jk,jl) * v_i(ji,jj,jl) * r1_nlay_i 
                  ENDIF
               END_2D
            END DO
         END DO
         !
      ENDIF ! npti > 0
      !
      ! Floe size distribution: tendency diagnostics (lateral growth 'lag'; new ice/open water 'opw'; welding 'wel')
      IF( ln_fsd ) THEN
         !               !-------------------------------------------------------------------------------!
         !               ! Process |   FSTD before   |    FSTD after     |  a_i before  |   a_i after    !
         !               !---------+-----------------+-------------------+--------------+----------------!
         CALL ice_fsd_dia(  'lag'  , za_ifsdb_latgro ,  za_ifsda_latgro  , za_ib_latgro , za_ia_latgro   )
         CALL ice_fsd_dia(  'opw'  , za_ifsda_latgro ,  za_ifsda_newice  , za_ia_latgro , za_ia_newice   )
         CALL ice_fsd_dia(  'wel'  , za_ifsda_newice , a_ifsd(A2D(0),:,:), a_i(A2D(0),:), a_i (A2D(0),:) )
      ENDIF
      !
      ! the following fields need to be updated on the halos (done in icethd): a_i, v_i, sv_i, e_i 
      !
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icethd_do', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      IF( ln_icediachk )   CALL ice_cons2D  (1, 'icethd_do',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)
      IF( ln_timing    )   CALL timing_stop ('icethd_do')
      !
   END SUBROUTINE ice_thd_do


   SUBROUTINE ice_thd_frazil
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_frazil ***
      !!
      !! ** Purpose :   frazil ice collection thickness and fraction
      !!
      !! ** Inputs  :   u_ice, v_ice, utau_ice, vtau_ice
      !! ** Ouputs  :   ht_i_new, fraz_frac
      !!-----------------------------------------------------------------------
      INTEGER  ::   ji, jj             ! dummy loop indices
      INTEGER  ::   iter
      REAL(wp) ::   zvfrx, zvgx, ztaux, zf, ztenagm, zvfry, zvgy, ztauy, zvrel2, zfp, ztwogp
      REAL(wp), PARAMETER ::   zcai    = 1.4e-3_wp                       ! ice-air drag (clem: should be dependent on coupling/forcing used)
      REAL(wp), PARAMETER ::   zhicrit = 0.04_wp                         ! frazil ice thickness
      REAL(wp), PARAMETER ::   zsqcd   = 1.0_wp / SQRT( 1.3_wp * zcai )  ! 1/SQRT(airdensity*drag)
      REAL(wp), PARAMETER ::   zgamafr = 0.03_wp
      !!-----------------------------------------------------------------------
      !
      !---------------------------------------------------------!
      ! Collection thickness of ice formed in leads and polynyas
      !---------------------------------------------------------!    
      ! ht_i_new is the thickness of new ice formed in open water
      ! ht_i_new can be either prescribed (ln_frazil=F) or computed (ln_frazil=T)
      ! Frazil ice forms in open water, is transported by wind, accumulates at the edge of the consolidated ice edge
      ! where it forms aggregates of a specific thickness called collection thickness.
      !
      fraz_frac(:,:) = 0._wp
      !
      ! Default new ice thickness
      WHERE( qlead(:,:) < 0._wp ) ! cooling
         ht_i_new(:,:) = rn_hinew
      ELSEWHERE
         ht_i_new(:,:) = 0._wp
      END WHERE

      IF( ln_frazil ) THEN
         ztwogp  = 2._wp * rho0 / ( grav * 0.3_wp * ( rho0 - rhoi ) )  ! reduced grav
         !
         DO_2D( 0, 0, 0, 0 )
            IF ( qlead(ji,jj) < 0._wp ) THEN ! cooling
               ! -- Wind stress -- !
               ztaux = utau_ice(ji,jj) * smask0(ji,jj)
               ztauy = vtau_ice(ji,jj) * smask0(ji,jj)
               ! Square root of wind stress
               ztenagm = SQRT( SQRT( ztaux * ztaux + ztauy * ztauy ) )

               ! -- Frazil ice velocity -- !
               IF( ztenagm >= epsi10 ) THEN
                  zvfrx = zgamafr * zsqcd * ztaux / MAX( ztenagm, epsi10 )
                  zvfry = zgamafr * zsqcd * ztauy / MAX( ztenagm, epsi10 )
               ELSE
                  zvfrx = 0._wp
                  zvfry = 0._wp
               ENDIF
               ! -- Pack ice velocity -- !
               zvgx = ( u_ice(ji-1,jj  ) * umask(ji-1,jj  ,1)  + u_ice(ji,jj) * umask(ji,jj,1) ) * 0.5_wp
               zvgy = ( v_ice(ji  ,jj-1) * vmask(ji  ,jj-1,1)  + v_ice(ji,jj) * vmask(ji,jj,1) ) * 0.5_wp

               ! -- Relative frazil/pack ice velocity & fraction of frazil ice-- !
               IF( at_i(ji,jj) >= epsi10 ) THEN
                  zvrel2 = MAX( (zvfrx - zvgx)*(zvfrx - zvgx) + (zvfry - zvgy)*(zvfry - zvgy), 0.15_wp*0.15_wp )
                  fraz_frac(ji,jj) = ( TANH( rn_Cfraz * ( SQRT(zvrel2) - rn_vfraz ) ) + 1._wp ) * 0.5_wp * rn_maxfraz
               ELSE
                  zvrel2 = 0._wp
                  fraz_frac(ji,jj) = 0._wp
               ENDIF
               
               ! -- new ice thickness (iterative loop) -- !
               ht_i_new(ji,jj) = zhicrit +   ( zhicrit + 0.1_wp )    &
                  &                      / ( ( zhicrit + 0.1_wp ) * ( zhicrit + 0.1_wp ) -  zhicrit * zhicrit ) * ztwogp * zvrel2
               iter = 1
               DO WHILE ( iter < 20 ) 
                  zf  = ( ht_i_new(ji,jj) - zhicrit ) * ( ht_i_new(ji,jj) * ht_i_new(ji,jj) - zhicrit * zhicrit ) -   &
                     &    ht_i_new(ji,jj) * zhicrit * ztwogp * zvrel2
                  zfp = ( ht_i_new(ji,jj) - zhicrit ) * ( 3.0_wp * ht_i_new(ji,jj) + zhicrit ) - zhicrit * ztwogp * zvrel2

                  ht_i_new(ji,jj) = ht_i_new(ji,jj) - zf / MAX( zfp, epsi20 )
                  iter = iter + 1
               END DO
               !
               ! bound ht_i_new (though I don't see why it should be necessary)
               ht_i_new(ji,jj) = MAX( 0.01_wp, MIN( ht_i_new(ji,jj), hi_max(jpl) ) )
               !
            ELSE
               ht_i_new(ji,jj) = 0._wp
            ENDIF
            !
         END_2D
         ! 
      ENDIF
   END SUBROUTINE ice_thd_frazil

   SUBROUTINE ice_thd_do_init
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_do_init *** 
      !!                 
      !! ** Purpose :   Physical constants and parameters associated with
      !!                ice growth in the leads
      !!
      !! ** Method  :   Read the namthd_do namelist and check the parameters
      !!                called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namthd_do
      !!-------------------------------------------------------------------
      INTEGER  ::   ios   ! Local integer 
      !!
      NAMELIST/namthd_do/ rn_hinew, ln_frazil, rn_maxfraz, rn_vfraz, rn_Cfraz
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namthd_do)
      READ_NML_CFG(numnam_ice,namthd_do)
      IF(lwm) WRITE( numoni, namthd_do )
      !
      IF(lwp) THEN                          ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_thd_do_init: Ice growth in open water'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namthd_do:'
         WRITE(numout,*) '      ice thickness for lateral accretion                       rn_hinew   = ', rn_hinew
         WRITE(numout,*) '      Frazil ice thickness as a function of wind or not         ln_frazil  = ', ln_frazil
         WRITE(numout,*) '      Maximum proportion of frazil ice collecting at bottom     rn_maxfraz = ', rn_maxfraz
         WRITE(numout,*) '      Threshold relative drift speed for collection of frazil   rn_vfraz   = ', rn_vfraz
         WRITE(numout,*) '      Squeezing coefficient for collection of frazil            rn_Cfraz   = ', rn_Cfraz
      ENDIF
      !
      IF ( rn_hinew < rn_himin )   CALL ctl_stop( 'ice_thd_do_init : rn_hinew should be >= rn_himin' )
      !
   END SUBROUTINE ice_thd_do_init
   
#else
   !!----------------------------------------------------------------------
   !!   Default option                                NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icethd_do
