MODULE iceistate
   !!======================================================================
   !!                     ***  MODULE  iceistate  ***
   !!   sea-ice : Initialization of ice variables
   !!======================================================================
   !! History :  2.0  !  2004-01  (C. Ethe, G. Madec) Original code
   !!            3.0  !  2007     (M. Vancoppenolle)  Rewrite for ice cats
   !!            4.0  !  2018     (many people)       SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_istate       :  initialization of diagnostics ice variables
   !!   ice_istate_init  :  initialization of ice state and namelist read
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE phycst         ! physical constant
   USE oce     , ONLY : ssh
   USE sbc_oce , ONLY : sst_m, sss_m, ln_ice_embd 
   USE sbc_ice , ONLY : tn_ice, snwice_mass, snwice_mass_b
   USE eosbn2         ! equation of state
# if defined key_qco
   USE domqco         ! Quasi-Eulerian coord.
# endif
   USE ice            ! sea-ice: variables
   USE icevar  , ONLY : ice_var_salprof, ice_var_itd
   USE icefsd  , ONLY : ice_fsd_istate      ! floe size distribution
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran    ! Fortran routines library
   USE fldread        ! read input fields
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)

# if defined key_agrif
   USE agrif_oce
   USE agrif_ice_interp 
# endif   

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_istate        ! called by icestp.F90
   PUBLIC   ice_istate_init   ! called by icestp.F90
   !
   !                                     !! ** namelist (namini) **
   REAL(wp) ::   rn_thres_sst
   REAL(wp) ::   rn_hti_ini_n, rn_hts_ini_n, rn_ati_ini_n, rn_smi_ini_n, rn_tmi_ini_n, rn_tsu_ini_n, rn_tms_ini_n
   REAL(wp) ::   rn_hti_ini_s, rn_hts_ini_s, rn_ati_ini_s, rn_smi_ini_s, rn_tmi_ini_s, rn_tsu_ini_s, rn_tms_ini_s
   REAL(wp) ::   rn_apd_ini_n, rn_hpd_ini_n, rn_hld_ini_n
   REAL(wp) ::   rn_apd_ini_s, rn_hpd_ini_s, rn_hld_ini_s
   !
   !                              ! if nn_iceini_file = 1
   INTEGER , PARAMETER ::   jpfldi = 10          ! maximum number of files to read
   INTEGER , PARAMETER ::   jp_hti = 1           ! index of ice thickness    (m)
   INTEGER , PARAMETER ::   jp_hts = 2           ! index of snw thickness    (m)
   INTEGER , PARAMETER ::   jp_ati = 3           ! index of ice fraction     (-)
   INTEGER , PARAMETER ::   jp_smi = 4           ! index of ice salinity     (g/kg)
   INTEGER , PARAMETER ::   jp_tmi = 5           ! index of ice temperature  (K)
   INTEGER , PARAMETER ::   jp_tsu = 6           ! index of ice surface temp (K)
   INTEGER , PARAMETER ::   jp_tms = 7           ! index of snw temperature  (K)
   INTEGER , PARAMETER ::   jp_apd = 8           ! index of pnd fraction     (-)
   INTEGER , PARAMETER ::   jp_hpd = 9           ! index of pnd depth        (m)
   INTEGER , PARAMETER ::   jp_hld = 10          ! index of pnd lid depth    (m)
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   si  ! structure of input fields (file informations, fields read)
   !
#if defined key_agrif
   REAL(wp), PUBLIC ::   rsshadj   !: initial mean ssh adjustment due to initial ice+snow mass
#endif
   !
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_istate( kt, Kbb, Kmm, Kaa )
      !!-------------------------------------------------------------------
      !!                    ***  ROUTINE ice_istate  ***
      !!
      !! ** Purpose :   defined the sea-ice initial state
      !!
      !! ** Method  :   This routine will put some ice where ocean
      !!                is at the freezing point, then fill in ice 
      !!                state variables using prescribed initial 
      !!                values in the namelist            
      !!
      !! ** Steps   :   1) Set initial surface and basal temperatures
      !!                2) Recompute or read sea ice state variables
      !!                3) Fill in space-dependent arrays for state variables
      !!                4) snow-ice mass computation
      !!
      !! ** Notes   : o_i, t_su, t_s, t_i, sz_i must be filled everywhere, even
      !!              where there is no ice
      !!--------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt            ! time step 
      INTEGER, INTENT(in) :: Kbb, Kmm, Kaa ! ocean time level indices
      !
      INTEGER  ::   ji, jj, jk, jl         ! dummy loop indices
      REAL(wp) ::   ztmelts, zsshadj, area
      INTEGER , DIMENSION(4)       ::   itest
      REAL(wp), DIMENSION(jpi,jpj) ::   ztfrz      ! freezing temperature (degC)
      REAL(wp), DIMENSION(jpi,jpj) ::   zswitch    ! ice indicator
      REAL(wp), DIMENSION(jpi,jpj) ::   zmsk       ! ice indicator
      REAL(wp), DIMENSION(jpi,jpj) ::   zht_i_ini, zat_i_ini, ztm_s_ini            !data from namelist or nc file
      REAL(wp), DIMENSION(jpi,jpj) ::   zt_su_ini, zht_s_ini, zsm_i_ini, ztm_i_ini !data from namelist or nc file
      REAL(wp), DIMENSION(jpi,jpj) ::   zapnd_ini, zhpnd_ini, zhlid_ini            !data from namelist or nc file
      !--------------------------------------------------------------------

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'ice_istate: sea-ice initialization '
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'

      !---------------------------
      ! 1) 1st init. of the fields
      !---------------------------
      !
      ! basal temperature (considered at freezing point)   [Kelvin]
      CALL eos_fzp( sss_m(:,:), t_bo(:,:), kbnd=0 )
      t_bo(:,:) = ( t_bo(:,:) + rt0 ) * smask0(:,:) 
      !
      ! == reduced arrays == !
      DO jl = 1, jpl
         DO_2D( 0, 0, 0, 0 )
            !
            cnd_ice(ji,jj,jl) = 0._wp                  ! conductivity at the ice top
            !
            tn_ice(ji,jj,jl) = t_i (ji,jj,1,jl)        ! temp for coupled runs
            t1_ice(ji,jj,jl) = t_i (ji,jj,1,jl)        ! temp for coupled runs
            !
            a_ip_eff(ji,jj,jl) = 0._wp   ! melt pond effective fraction
         END_2D
         !
      ENDDO
      
      ! == full arrays == !
      DO jl = 1, jpl
         !
         DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, nlay_i )
            ! heat
            e_i(ji,jj,jk,jl) = 0._wp
            t_i(ji,jj,jk,jl) = rt0 * tmask(ji,jj,1) ! ice temp
            ! salt
            szv_i(ji,jj,jk,jl) = 0._wp
            sz_i (ji,jj,jk,jl) = rn_simin * tmask(ji,jj,1)
         END_3D
         
         DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, nlay_s )
            ! heat
            e_s(ji,jj,jk,jl) = 0._wp
            t_s(ji,jj,jk,jl) = rt0 * tmask(ji,jj,1) ! snw temp
         END_3D
         !
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            ! general fields
            a_i (ji,jj,jl) = 0._wp
            v_i (ji,jj,jl) = 0._wp
            v_s (ji,jj,jl) = 0._wp
            sv_i(ji,jj,jl) = 0._wp
            oa_i(ji,jj,jl) = 0._wp
            h_i (ji,jj,jl) = 0._wp
            h_s (ji,jj,jl) = 0._wp
            s_i (ji,jj,jl) = 0._wp
            o_i (ji,jj,jl) = 0._wp
            t_su(ji,jj,jl) = rt0 * tmask(ji,jj,1)
            !
            ! melt ponds
            a_ip(ji,jj,jl) = 0._wp
            v_ip(ji,jj,jl) = 0._wp
            v_il(ji,jj,jl) = 0._wp
            h_ip(ji,jj,jl) = 0._wp
            h_il(ji,jj,jl) = 0._wp
         END_2D
         !
      ENDDO
      !
      ! ice velocities
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         u_ice(ji,jj) = 0._wp
         v_ice(ji,jj) = 0._wp
      END_2D
      !
      !------------------------------------------------------------------------
      ! 2) overwrite some of the fields with namelist parameters or netcdf file
      !------------------------------------------------------------------------
      IF( ln_iceini ) THEN
         !
#if defined key_agrif
         IF ( ( Agrif_Root() ).OR.(.NOT.ln_init_chfrpar ) ) THEN
#endif
            !                             !---------------!
            IF( nn_iceini_file == 1 )THEN ! Read a file   !
               !                          !---------------!
               WHERE( ff_t(:,:) >= 0._wp )   ;   zmsk(:,:) = 1._wp
               ELSEWHERE                     ;   zmsk(:,:) = 0._wp
               END WHERE
               !
               CALL fld_read( kt, 1, si ) ! input fields provided at the current time-step
               !
               ! -- mandatory fields -- !
               zht_i_ini(:,:) = si(jp_hti)%fnow(:,:,1) * tmask(:,:,1)
               zht_s_ini(:,:) = si(jp_hts)%fnow(:,:,1) * tmask(:,:,1)
               zat_i_ini(:,:) = si(jp_ati)%fnow(:,:,1) * tmask(:,:,1)

               ! -- optional fields -- !
               !    if fields do not exist then set them to the values present in the namelist (except for temperatures)
               !
               ! ice salinity
               IF( TRIM(si(jp_smi)%clrootname) == 'NOT USED' ) &
                  &     si(jp_smi)%fnow(:,:,1) = ( rn_smi_ini_n * zmsk + rn_smi_ini_s * (1._wp - zmsk) ) * tmask(:,:,1)
               !
               ! temperatures
               IF    ( TRIM(si(jp_tmi)%clrootname) == 'NOT USED' .AND. TRIM(si(jp_tsu)%clrootname) == 'NOT USED' .AND. &
                  &    TRIM(si(jp_tms)%clrootname) == 'NOT USED' ) THEN
                  si(jp_tmi)%fnow(:,:,1) = ( rn_tmi_ini_n * zmsk + rn_tmi_ini_s * (1._wp - zmsk) ) * tmask(:,:,1)
                  si(jp_tsu)%fnow(:,:,1) = ( rn_tsu_ini_n * zmsk + rn_tsu_ini_s * (1._wp - zmsk) ) * tmask(:,:,1)
                  si(jp_tms)%fnow(:,:,1) = ( rn_tms_ini_n * zmsk + rn_tms_ini_s * (1._wp - zmsk) ) * tmask(:,:,1)
               ENDIF
               IF( TRIM(si(jp_tmi)%clrootname) == 'NOT USED' .AND. TRIM(si(jp_tms)%clrootname) /= 'NOT USED' ) & ! if T_s is read and not T_i, set T_i = (T_s + T_freeze)/2
                  &     si(jp_tmi)%fnow(:,:,1) = 0.5_wp * ( si(jp_tms)%fnow(:,:,1) + 271.15 )
               IF( TRIM(si(jp_tmi)%clrootname) == 'NOT USED' .AND. TRIM(si(jp_tsu)%clrootname) /= 'NOT USED' ) & ! if T_su is read and not T_i, set T_i = (T_su + T_freeze)/2
                  &     si(jp_tmi)%fnow(:,:,1) = 0.5_wp * ( si(jp_tsu)%fnow(:,:,1) + 271.15 )
               IF( TRIM(si(jp_tsu)%clrootname) == 'NOT USED' .AND. TRIM(si(jp_tms)%clrootname) /= 'NOT USED' ) & ! if T_s is read and not T_su, set T_su = T_s
                  &     si(jp_tsu)%fnow(:,:,1) = si(jp_tms)%fnow(:,:,1)
               IF( TRIM(si(jp_tsu)%clrootname) == 'NOT USED' .AND. TRIM(si(jp_tmi)%clrootname) /= 'NOT USED' ) & ! if T_i is read and not T_su, set T_su = T_i
                  &     si(jp_tsu)%fnow(:,:,1) = si(jp_tmi)%fnow(:,:,1)
               IF( TRIM(si(jp_tms)%clrootname) == 'NOT USED' .AND. TRIM(si(jp_tsu)%clrootname) /= 'NOT USED' ) & ! if T_su is read and not T_s, set T_s = T_su
                  &     si(jp_tms)%fnow(:,:,1) = si(jp_tsu)%fnow(:,:,1)
               IF( TRIM(si(jp_tms)%clrootname) == 'NOT USED' .AND. TRIM(si(jp_tmi)%clrootname) /= 'NOT USED' ) & ! if T_i is read and not T_s, set T_s = T_i
                  &     si(jp_tms)%fnow(:,:,1) = si(jp_tmi)%fnow(:,:,1)
               !
               ! pond concentration
               IF( TRIM(si(jp_apd)%clrootname) == 'NOT USED' ) &
                  &     si(jp_apd)%fnow(:,:,1) = ( rn_apd_ini_n * zmsk + rn_apd_ini_s * (1._wp - zmsk) ) * tmask(:,:,1) & ! rn_apd = pond fraction => rn_apnd * a_i = pond conc.
                  &                              * si(jp_ati)%fnow(:,:,1) 
               !
               ! pond depth
               IF( TRIM(si(jp_hpd)%clrootname) == 'NOT USED' ) &
                  &     si(jp_hpd)%fnow(:,:,1) = ( rn_hpd_ini_n * zmsk + rn_hpd_ini_s * (1._wp - zmsk) ) * tmask(:,:,1)
               !
               ! pond lid depth
               IF( TRIM(si(jp_hld)%clrootname) == 'NOT USED' ) &
                  &     si(jp_hld)%fnow(:,:,1) = ( rn_hld_ini_n * zmsk + rn_hld_ini_s * (1._wp - zmsk) ) * tmask(:,:,1)
               !
               zsm_i_ini(:,:) = si(jp_smi)%fnow(:,:,1) * tmask(:,:,1)
               ztm_i_ini(:,:) = si(jp_tmi)%fnow(:,:,1) * tmask(:,:,1)
               zt_su_ini(:,:) = si(jp_tsu)%fnow(:,:,1) * tmask(:,:,1)
               ztm_s_ini(:,:) = si(jp_tms)%fnow(:,:,1) * tmask(:,:,1)
               zapnd_ini(:,:) = si(jp_apd)%fnow(:,:,1) * tmask(:,:,1)
               zhpnd_ini(:,:) = si(jp_hpd)%fnow(:,:,1) * tmask(:,:,1)
               zhlid_ini(:,:) = si(jp_hld)%fnow(:,:,1) * tmask(:,:,1)
               !
               !                          !---------------!
            ELSE                          ! Read namelist !
               !                          !---------------!
               ! no ice if (sst - Tfreez) >= thresold
               CALL eos_fzp( sss_m(:,:), ztfrz(:,:), kbnd=nn_hls )
               WHERE( ( sst_m(:,:) - ztfrz(:,:) ) * tmask(:,:,1) >= rn_thres_sst )   ;   zmsk(:,:) = 0._wp 
               ELSEWHERE                                                             ;   zmsk(:,:) = tmask(:,:,1)
               END WHERE
               !
               ! assign initial thickness, concentration, snow depth and salinity to an hemisphere-dependent array
               WHERE( ff_t(:,:) >= 0._wp )
                  zht_i_ini(:,:) = rn_hti_ini_n * zmsk(:,:)
                  zht_s_ini(:,:) = rn_hts_ini_n * zmsk(:,:)
                  zat_i_ini(:,:) = rn_ati_ini_n * zmsk(:,:)
                  zsm_i_ini(:,:) = rn_smi_ini_n * zmsk(:,:)
                  ztm_i_ini(:,:) = rn_tmi_ini_n * zmsk(:,:)
                  zt_su_ini(:,:) = rn_tsu_ini_n * zmsk(:,:)
                  ztm_s_ini(:,:) = rn_tms_ini_n * zmsk(:,:)
                  zapnd_ini(:,:) = rn_apd_ini_n * zmsk(:,:) * zat_i_ini(:,:) ! rn_apd = pond fraction => rn_apd * a_i = pond conc. 
                  zhpnd_ini(:,:) = rn_hpd_ini_n * zmsk(:,:)
                  zhlid_ini(:,:) = rn_hld_ini_n * zmsk(:,:)
               ELSEWHERE
                  zht_i_ini(:,:) = rn_hti_ini_s * zmsk(:,:)
                  zht_s_ini(:,:) = rn_hts_ini_s * zmsk(:,:)
                  zat_i_ini(:,:) = rn_ati_ini_s * zmsk(:,:)
                  zsm_i_ini(:,:) = rn_smi_ini_s * zmsk(:,:)
                  ztm_i_ini(:,:) = rn_tmi_ini_s * zmsk(:,:)
                  zt_su_ini(:,:) = rn_tsu_ini_s * zmsk(:,:)
                  ztm_s_ini(:,:) = rn_tms_ini_s * zmsk(:,:)
                  zapnd_ini(:,:) = rn_apd_ini_s * zmsk(:,:) * zat_i_ini(:,:) ! rn_apd = pond fraction => rn_apd * a_i = pond conc.
                  zhpnd_ini(:,:) = rn_hpd_ini_s * zmsk(:,:)
                  zhlid_ini(:,:) = rn_hld_ini_s * zmsk(:,:)
               END WHERE
               !
            ENDIF

            ! make sure ponds = 0 if no ponds scheme
            IF ( .NOT.ln_pnd ) THEN
               zapnd_ini(:,:) = 0._wp
               zhpnd_ini(:,:) = 0._wp
               zhlid_ini(:,:) = 0._wp
            ENDIF
            
            IF ( .NOT.ln_pnd_lids ) THEN
               zhlid_ini(:,:) = 0._wp
            ENDIF
            
            !----------------!
            ! 3) fill fields !
            !----------------!
            ! distribute 1-cat into jpl-cat: (jpi*jpj) -> (jpi*jpj,jpl)
            DO jj = 1, jpj
               CALL ice_var_itd( zht_i_ini(:,jj)     , zht_s_ini(:,jj)     , zat_i_ini(:,jj)   ,                     & ! <<= in
                  &                h_i    (:,jj,:)   ,   h_s    (:,jj,:)   ,   a_i    (:,jj,:) ,                     & ! =>> out
                  &              ztm_i_ini(:,jj)     , ztm_s_ini(:,jj)     , zt_su_ini(:,jj)   , zsm_i_ini(:,jj)   , & ! <<= in
                  &              zapnd_ini(:,jj)     , zhpnd_ini(:,jj)     , zhlid_ini(:,jj)   ,                     & ! <<= in
                  &                t_i    (:,jj,1,:) ,   t_s    (:,jj,1,:) ,   t_su   (:,jj,:) ,       s_i(:,jj,:) , & ! =>> out
                  &                a_ip   (:,jj,:)   ,   h_ip   (:,jj,:)   ,   h_il   (:,jj,:) )                       ! =>> out
            ENDDO
            DO jl = 1, jpl
               DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 2, nlay_s )
                  t_s(ji,jj,jk,jl) = t_s(ji,jj,1,jl)
               END_3D
            ENDDO
            DO jl = 1, jpl
               DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 2, nlay_i )
                  t_i(ji,jj,jk,jl) = t_i(ji,jj,1,jl)
               END_3D
            ENDDO

            ! switch for the following
            WHERE( SUM(a_i(:,:,:),dim=3) > 0._wp )   ;   zswitch(:,:) = tmask(:,:,1) 
            ELSEWHERE                                ;   zswitch(:,:) = 0._wp
            END WHERE
            
            ! calculate extensive and intensive variables
            DO jl = 1, jpl
               DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                  s_i(ji,jj,jl) = MIN( MAX( rn_simin , s_i(ji,jj,jl) ) , rn_sinew * sss_m(ji,jj) )
               END_2D
            END DO
            CALL ice_var_salprof ! for sz_i

            DO jl = 1, jpl
               DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                  v_i (ji,jj,jl) = h_i(ji,jj,jl) * a_i(ji,jj,jl)
                  v_s (ji,jj,jl) = h_s(ji,jj,jl) * a_i(ji,jj,jl)
                  sv_i(ji,jj,jl) = s_i(ji,jj,jl) * v_i(ji,jj,jl)
               END_2D
            END DO
            !
            DO jl = 1, jpl
               DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, nlay_s )
                  e_s(ji,jj,jk,jl) = zswitch(ji,jj) * v_s(ji,jj,jl) * r1_nlay_s * &
                     &               rhos * ( rcpi * ( rt0 - t_s(ji,jj,jk,jl) ) + rLfus )
               END_3D
            END DO
            !
            DO jl = 1, jpl
               DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, nlay_i )
                  ! salt
                  szv_i(ji,jj,jk,jl) = sz_i(ji,jj,jk,jl) * v_i(ji,jj,jl) * r1_nlay_i
                  ! heat
                  ztmelts          = - rTmlt * sz_i(ji,jj,jk,jl) + rt0 ! melting temperature in K
                  e_i(ji,jj,jk,jl) = zswitch(ji,jj) * v_i(ji,jj,jl) * r1_nlay_i * &
                     &               rhoi * (  rcpi  * ( ztmelts - t_i(ji,jj,jk,jl) ) + &
                     &                         rLfus * ( 1._wp - (ztmelts-rt0) / MIN( (t_i(ji,jj,jk,jl)-rt0), -epsi20 ) ) &
                     &                       - rcp   * ( ztmelts - rt0 ) )
               END_3D
            END DO
            !
#if defined key_agrif
         ELSE
            CALL  agrif_istate_ice
         ENDIF
#endif
         ! Melt ponds
         WHERE( a_i(A2D(0),:) > epsi10 )   ;   a_ip_eff(:,:,:) = a_ip(A2D(0),:) / a_i(A2D(0),:)
         ELSEWHERE                         ;   a_ip_eff(:,:,:) = 0._wp
         END WHERE
         v_ip(:,:,:) = h_ip(:,:,:) * a_ip(:,:,:)
         v_il(:,:,:) = h_il(:,:,:) * a_ip(:,:,:)
         
         ! specific temperatures for coupled runs
         DO_2D( 0, 0, 0, 0 )
            tn_ice(ji,jj,:) = t_su(ji,jj,:)
            t1_ice(ji,jj,:) = t_i (ji,jj,1,:)
         END_2D
         !
         ! ice concentration should not exceed amax
         at_i(:,:) = SUM( a_i, dim=3 )
         DO jl = 1, jpl
            WHERE( at_i(:,:) > rn_amax_2d(:,:) )   a_i(:,:,jl) = a_i(:,:,jl) * rn_amax_2d(:,:) / at_i(:,:)
         END DO
        at_i(:,:) = SUM( a_i, dim=3 )
         !
      ENDIF ! ln_iceini
      !
      ! === floe size distribution === !
      !
      ! Multiple options to initialise FSD handled by ice_fsd_istate (in icefsd.F90)
      ! (it checks for the case ln_iceini = F and forces initial FSD = 0 if so):
      IF( ln_fsd )   CALL ice_fsd_istate
      !
      !----------------------------------------------------------
      ! 4) Adjust ssh and vertical scale factors to snow-ice mass
      !----------------------------------------------------------
      snwice_mass  (:,:) = tmask(:,:,1) * SUM( rhos * v_s + rhoi * v_i + rhow * ( v_ip + v_il ), dim=3  )   ! snow+ice mass
      snwice_mass_b(:,:) = snwice_mass(:,:)
      !
      IF( ln_ice_embd ) THEN            ! embedded sea-ice: deplete the initial ssh below sea-ice area
         !                              ! ----------------
         ssh(:,:,Kmm) = ssh(:,:,Kmm) - snwice_mass(:,:) * r1_rho0
         ssh(:,:,Kbb) = ssh(:,:,Kbb) - snwice_mass(:,:) * r1_rho0
         !
      ELSE                              ! levitating sea-ice: deplete the initial ssh over the whole domain
         !                              ! ------------------
         area    = glob_2Dsum( 'iceistate', e1e2t(:,:) * ssmask(:,:) )
         zsshadj = glob_2Dsum( 'iceistate', snwice_mass(:,:) * r1_rho0 * e1e2t(:,:) ) / area
#if defined key_agrif
         ! Override ssh adjustment in nested domains by the root-domain ssh adjustment;
         ! store the adjustment value in a global module variable to make it retrievable in nested domains
         IF( .NOT.Agrif_Root() ) THEN
            zsshadj = Agrif_Parent(rsshadj)
         ELSE
            rsshadj = zsshadj
         ENDIF
#endif
         IF(lwp) WRITE(numout,'(A23,F10.6,A20)') ' sea level adjusted by ', -zsshadj, ' m to compensate for'
         IF(lwp) WRITE(numout,*) ' the initial snow+ice mass'
         !
         WHERE( ssmask(:,:) == 1._wp )
            ssh(:,:,Kmm) = ssh(:,:,Kmm) - zsshadj
            ssh(:,:,Kbb) = ssh(:,:,Kbb) - zsshadj
         ENDWHERE
         !
      ENDIF
      !
      IF( .NOT.lk_linssh ) THEN
#if defined key_qco
         CALL dom_qco_zgr( Kbb, Kmm )        ! upadte of r3=ssh/h0 ratios
#elif defined key_linssh
         !                                   ! Fix in time : key_linssh case, set through domzgr_substitute.h90
#endif
      ENDIF
      !
      !!clem: output of initial state should be written here but it is impossible because
      !!      the ocean and ice are in the same file
      !!      CALL dia_wri_state( 'output.init' )
      !
   END SUBROUTINE ice_istate


   SUBROUTINE ice_istate_init
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_istate_init  ***
      !!        
      !! ** Purpose :   Definition of initial state of the ice 
      !!
      !! ** Method  :   Read the namini namelist and check the parameter 
      !!              values called at the first timestep (nit000)
      !!
      !! ** input   :  Namelist namini
      !!
      !!-----------------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer output status for namelist read
      INTEGER ::   ifpr, ierror
      !
      CHARACTER(len=256) ::  cn_dir          ! Root directory for location of ice files
      TYPE(FLD_N)                    ::   sn_hti, sn_hts, sn_ati, sn_smi, sn_tmi, sn_tsu, sn_tms, sn_apd, sn_hpd, sn_hld
      TYPE(FLD_N), DIMENSION(jpfldi) ::   slf_i                 ! array of namelist informations on the fields to read
      !
      NAMELIST/namini/ ln_iceini, nn_iceini_file, rn_thres_sst, &
         &             rn_hti_ini_n, rn_hti_ini_s, rn_hts_ini_n, rn_hts_ini_s, &
         &             rn_ati_ini_n, rn_ati_ini_s, rn_smi_ini_n, rn_smi_ini_s, &
         &             rn_tmi_ini_n, rn_tmi_ini_s, rn_tsu_ini_n, rn_tsu_ini_s, rn_tms_ini_n, rn_tms_ini_s, &
         &             rn_apd_ini_n, rn_apd_ini_s, rn_hpd_ini_n, rn_hpd_ini_s, rn_hld_ini_n, rn_hld_ini_s, &
         &             sn_hti, sn_hts, sn_ati, sn_tsu, sn_tmi, sn_smi, sn_tms, sn_apd, sn_hpd, sn_hld, cn_dir
      !!-----------------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namini)
      READ_NML_CFG(numnam_ice,namini)
      IF(lwm) WRITE ( numoni, namini )
      !
      slf_i(jp_hti) = sn_hti  ;  slf_i(jp_hts) = sn_hts
      slf_i(jp_ati) = sn_ati  ;  slf_i(jp_smi) = sn_smi
      slf_i(jp_tmi) = sn_tmi  ;  slf_i(jp_tsu) = sn_tsu   ;   slf_i(jp_tms) = sn_tms
      slf_i(jp_apd) = sn_apd  ;  slf_i(jp_hpd) = sn_hpd   ;   slf_i(jp_hld) = sn_hld
      !
      IF(lwp) THEN                          ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_istate_init: ice parameters inititialisation '
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namini:'
         WRITE(numout,*) '      ice initialization (T) or not (F)                ln_iceini      = ', ln_iceini
         WRITE(numout,*) '      ice initialization from a netcdf file            nn_iceini_file = ', nn_iceini_file
         WRITE(numout,*) '      max ocean temp. above Tfreeze with initial ice   rn_thres_sst   = ', rn_thres_sst
         IF( ln_iceini .AND. nn_iceini_file == 0 ) THEN
            WRITE(numout,*) '      initial snw thickness in the north-south         rn_hts_ini     = ', rn_hts_ini_n,rn_hts_ini_s 
            WRITE(numout,*) '      initial ice thickness in the north-south         rn_hti_ini     = ', rn_hti_ini_n,rn_hti_ini_s
            WRITE(numout,*) '      initial ice concentr  in the north-south         rn_ati_ini     = ', rn_ati_ini_n,rn_ati_ini_s
            WRITE(numout,*) '      initial ice salinity  in the north-south         rn_smi_ini     = ', rn_smi_ini_n,rn_smi_ini_s
            WRITE(numout,*) '      initial surf temperat in the north-south         rn_tsu_ini     = ', rn_tsu_ini_n,rn_tsu_ini_s
            WRITE(numout,*) '      initial ice temperat  in the north-south         rn_tmi_ini     = ', rn_tmi_ini_n,rn_tmi_ini_s
            WRITE(numout,*) '      initial snw temperat  in the north-south         rn_tms_ini     = ', rn_tms_ini_n,rn_tms_ini_s
            WRITE(numout,*) '      initial pnd fraction  in the north-south         rn_apd_ini     = ', rn_apd_ini_n,rn_apd_ini_s
            WRITE(numout,*) '      initial pnd depth     in the north-south         rn_hpd_ini     = ', rn_hpd_ini_n,rn_hpd_ini_s
            WRITE(numout,*) '      initial pnd lid depth in the north-south         rn_hld_ini     = ', rn_hld_ini_n,rn_hld_ini_s
         ENDIF
      ENDIF
      !
      IF( nn_iceini_file == 1 ) THEN                      ! Ice initialization using input file
         !
         ! set si structure
         ALLOCATE( si(jpfldi), STAT=ierror )
         IF( ierror > 0 ) THEN
            CALL ctl_stop( 'ice_istate_ini in iceistate: unable to allocate si structure' )   ;   RETURN
         ENDIF
         !
         DO ifpr = 1, jpfldi
            ALLOCATE( si(ifpr)%fnow(jpi,jpj,1) )
            IF( slf_i(ifpr)%ln_tint )  ALLOCATE( si(ifpr)%fdta(jpi,jpj,1,2) )
         END DO
         !
         ! fill si with slf_i and control print
         CALL fld_fill( si, slf_i, cn_dir, 'ice_istate_ini', 'initialization of sea ice fields', 'numnam_ice' )
         !
      ENDIF
      !
      IF( .NOT.ln_pnd ) THEN
         rn_apd_ini_n = 0. ; rn_apd_ini_s = 0.
         rn_hpd_ini_n = 0. ; rn_hpd_ini_s = 0.
         rn_hld_ini_n = 0. ; rn_hld_ini_s = 0.
         CALL ctl_warn( 'rn_apd_ini & rn_hpd_ini = 0 & rn_hld_ini = 0 when no ponds' )
      ENDIF
      !
      IF( .NOT.ln_pnd_lids ) THEN
         rn_hld_ini_n = 0. ; rn_hld_ini_s = 0.
      ENDIF
      !
   END SUBROUTINE ice_istate_init

#else
   !!----------------------------------------------------------------------
   !!   Default option :         Empty module         NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE iceistate
