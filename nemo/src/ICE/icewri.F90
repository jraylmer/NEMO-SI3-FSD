MODULE icewri
   !!======================================================================
   !!                     ***  MODULE  icewri  ***
   !!   sea-ice : output ice variables
   !!======================================================================
   !! History :  4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_wri       : write of the diagnostics variables in ouput file
   !!   ice_wri_state : write for initial state or/and abandon
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE phycst         ! physical constant
   USE sbc_oce , ONLY : sprecip, qns, qsr, sst_m, sss_m
   USE sbc_ice        ! Surface boundary condition: ice fields
   USE ice            ! sea-ice: variables
   USE icevar  , ONLY : ice_var_brine
   !
   USE ioipsl         !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lib_fortran
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC ice_wri        ! called by ice_stp
   PUBLIC ice_wri_state  ! called by dia_wri_state

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_wri( kt )
      !!-------------------------------------------------------------------
      !!  This routine ouputs some (most?) of the sea ice fields
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! time-step
      !
      INTEGER  ::   ji, jj, jk, jl  ! dummy loop indices
      REAL(wp) ::   z2da, z2db, zrho1, zrho2
      REAL(wp) ::   zmiss           ! missing value retrieved from xios
      REAL(wp), DIMENSION(A2D(0))               ::   z2d                            ! 2D workspace
      REAL(wp), DIMENSION(A2D(0))               ::   zmsk00, zmsk05, zmsk15, zmsksn ! O%, 5% and 15% concentration mask and snow mask
      REAL(wp), DIMENSION(A2D(0),jpl)           ::   zmsk00c, zmsksnc               ! categories masks
      REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE ::   zmsk00l                        ! layers masks
      REAL(wp), DIMENSION(:,:)    , ALLOCATABLE ::   zfast, zalb, zmskalb           ! 2D workspace
      !
      ! Global ice diagnostics (SIMIP)
      REAL(wp), DIMENSION(A2D(0),6) ::   ztmp
      COMPLEX(dp), DIMENSION(6)     ::   ysum
      REAL(wp), DIMENSION(6)        ::   zdiag
      !!-------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('icewri')

      ! get missing value from xml
      CALL iom_miss_val( 'icetemp', zmiss )

      ! brine volume
      IF( iom_use('icebrv') .OR. iom_use('icebrv_cat') )   CALL ice_var_brine

      ! tresholds for outputs
      zmsk00 (:,:)   = MERGE( 1._wp, 0._wp, at_i(A2D(0))  >= epsi06  )
      zmsk05 (:,:)   = MERGE( 1._wp, 0._wp, at_i(A2D(0))  >= 0.05_wp )
      zmsk15 (:,:)   = MERGE( 1._wp, 0._wp, at_i(A2D(0))  >= 0.15_wp )
      zmsksn (:,:)   = MERGE( 1._wp, 0._wp, vt_s(A2D(0))  >= epsi06  )
      zmsk00c(:,:,:) = MERGE( 1._wp, 0._wp, a_i(A2D(0),:) >= epsi06  )
      zmsksnc(:,:,:) = MERGE( 1._wp, 0._wp, v_s(A2D(0),:) >= epsi06  )

      IF( iom_use('icetemp_lay' ) .OR. iom_use('icesalt_lay' ) ) THEN
         ALLOCATE( zmsk00l(A2D(0),nlay_i,jpl) )
         DO jl = 1, jpl
            DO_3D( 0, 0, 0, 0, 1, nlay_i )
               zmsk00l(ji,jj,jk,jl) = zmsk00c(ji,jj,jl)
            END_3D
         ENDDO
      ENDIF
   
      !-----------------
      ! Standard outputs
      !-----------------
      zrho1 = ( rho0 - rhoi ) * r1_rho0 ; zrho2 = rhos * r1_rho0
      ! masks
      CALL iom_put( 'icemask'  , zmsk00 )   ! ice mask 0%
      CALL iom_put( 'icemask05', zmsk05 )   ! ice mask 5%
      CALL iom_put( 'icemask15', zmsk15 )   ! ice mask 15%
      CALL iom_put( 'icepres'  , zmsk00 )   ! Ice presence (1 or 0)
      !
      ! general fields
      IF( iom_use('icemass' ) )   CALL iom_put( 'icemass', vt_i(A2D(0)) * rhoi * zmsk00 )                                     ! Ice mass per cell area
      IF( iom_use('snwmass' ) )   CALL iom_put( 'snwmass', vt_s(A2D(0)) * rhos * zmsksn )                                     ! Snow mass per cell area
      IF( iom_use('iceconc' ) )   CALL iom_put( 'iceconc', at_i(A2D(0))        * zmsk00 )                                     ! ice concentration
      IF( iom_use('icevolu' ) )   CALL iom_put( 'icevolu', vt_i(A2D(0))        * zmsk00 )                                     ! ice volume = mean ice thickness over the cell
      IF( iom_use('icethic' ) )   CALL iom_put( 'icethic', hm_i(:,:)           * zmsk00 )                                     ! ice thickness
      IF( iom_use('snwthic' ) )   CALL iom_put( 'snwthic', hm_s(:,:)           * zmsk00 )                                     ! snw thickness
      IF( iom_use('icebrv'  ) )   CALL iom_put( 'icebrv' , vm_ibr(:,:)* 100.   * zmsk00 )                                     ! brine volume
      IF( iom_use('iceage'  ) )   CALL iom_put( 'iceage' , om_i(:,:) / rday    * zmsk15 + zmiss * ( 1._wp - zmsk15 ) )        ! ice age
      IF( iom_use('icehnew' ) )   CALL iom_put( 'icehnew', ht_i_new(:,:)                )                                     ! new ice thickness formed in the leads
      IF( iom_use('snwvolu' ) )   CALL iom_put( 'snwvolu', vt_s(A2D(0))        * zmsksn )                                     ! snow volume
      IF( iom_use('icefrb'  ) ) THEN                                                                                          ! Ice freeboard
         z2d(:,:) = zrho1 * hm_i(:,:) - zrho2 * hm_s(:,:)
         WHERE( z2d < 0._wp )   z2d = 0._wp
                                  CALL iom_put( 'icefrb' , z2d * zmsk00 )
      ENDIF
      ! melt ponds
      IF( iom_use('iceapnd' ) )   CALL iom_put( 'iceapnd' , at_ip(A2D(0))     * zmsk00 )                                      ! melt pond total fraction
      IF( iom_use('iceapeff' ) )  CALL iom_put( 'iceapeff', at_ip_eff(A2D(0)) * zmsk00 )                                      ! melt pond total fraction (effective)
      IF( iom_use('icehpnd' ) )   CALL iom_put( 'icehpnd' , hm_ip(:,:)        * zmsk00 )                                      ! melt pond depth
      IF( iom_use('icevpnd' ) )   CALL iom_put( 'icevpnd' , vt_ip(A2D(0))     * zmsk00 )                                      ! melt pond total volume per unit area
      IF( iom_use('icehlid' ) )   CALL iom_put( 'icehlid' , hm_il(:,:)        * zmsk00 )                                      ! melt pond lid depth
      IF( iom_use('icevlid' ) )   CALL iom_put( 'icevlid' , vt_il(A2D(0))     * zmsk00 )                                      ! melt pond lid total volume per unit area
      ! salt
      IF( iom_use('icesalt' ) )   CALL iom_put( 'icesalt', sm_i(:,:)                 * zmsk00 + zmiss * ( 1._wp - zmsk00 ) )  ! mean ice salinity
      IF( iom_use('icesalm' ) )   CALL iom_put( 'icesalm', st_i(:,:) * rhoi * 1.0e-3 * zmsk00 )                               ! Mass of salt in sea ice per cell area
      ! heat
      IF( iom_use('icetemp' ) )   CALL iom_put( 'icetemp', ( tm_i (:,:) - rt0 ) * zmsk00 + zmiss * ( 1._wp - zmsk00 ) )       ! ice mean temperature
      IF( iom_use('snwtemp' ) )   CALL iom_put( 'snwtemp', ( tm_s (:,:) - rt0 ) * zmsksn + zmiss * ( 1._wp - zmsksn ) )       ! snw mean temperature
      IF( iom_use('icettop' ) )   CALL iom_put( 'icettop', ( tm_su(:,:) - rt0 ) * zmsk00 + zmiss * ( 1._wp - zmsk00 ) )       ! temperature at the ice surface
      IF( iom_use('icetbot' ) )   CALL iom_put( 'icetbot', ( t_bo (:,:) - rt0 ) * zmsk00 + zmiss * ( 1._wp - zmsk00 ) )       ! temperature at the ice bottom
      IF( iom_use('icetsni' ) )   CALL iom_put( 'icetsni', ( tm_si(:,:) - rt0 ) * zmsk00 + zmiss * ( 1._wp - zmsk00 ) )       ! temperature at the snow-ice interface
      IF( iom_use('icehc'   ) )   CALL iom_put( 'icehc'  ,  -et_i (:,:)         * zmsk00 )                                    ! ice heat content
      IF( iom_use('snwhc'   ) )   CALL iom_put( 'snwhc'  ,  -et_s (:,:)         * zmsksn )                                    ! snow heat content
      ! momentum
      IF( iom_use('uice'    ) )   CALL iom_put( 'uice'   ,   u_ice(:,:) )                                                     ! ice velocity u
      IF( iom_use('vice'    ) )   CALL iom_put( 'vice'   ,   v_ice(:,:) )                                                     ! ice velocity v
      !
      IF( iom_use('icevel') .OR. iom_use('fasticepres') ) THEN                                                                ! module of ice velocity & fast ice
         DO_2D( 0, 0, 0, 0 )
            z2da  = u_ice(ji,jj) + u_ice(ji-1,jj)
            z2db  = v_ice(ji,jj) + v_ice(ji,jj-1)
            z2d(ji,jj) = 0.5_wp * SQRT( z2da * z2da + z2db * z2db )
         END_2D
         CALL iom_put( 'icevel', z2d )

         ALLOCATE( zfast(A2D(0)) )
         zfast(:,:) = MERGE( 1._wp, 0._wp,  z2d(:,:) < 5.e-04_wp .AND. zmsk15(:,:) == 1._wp ) ! record presence of fast ice
         CALL iom_put( 'fasticepres', zfast )
         DEALLOCATE( zfast )
      ENDIF
      !
      IF( iom_use('icealb') .OR. iom_use('albedo') ) THEN                                                                     ! ice albedo and surface albedo
         ALLOCATE( zalb(A2D(0)), zmskalb(A2D(0)) )
         ! ice albedo
         WHERE( at_i_b(:,:) < 1.e-03 )
            zmskalb(:,:) = 0._wp
            zalb   (:,:) = ralb_oce
         ELSEWHERE
            zmskalb(:,:) = 1._wp
            zalb   (:,:) = SUM( alb_ice(:,:,:) * a_i_b(:,:,:), dim=3 ) / at_i_b(:,:)
         END WHERE
         CALL iom_put( 'icealb' , zalb * zmskalb + zmiss * ( 1._wp - zmskalb ) )
         ! ice+ocean albedo
         zalb(:,:) = SUM( alb_ice(:,:,:) * a_i_b(:,:,:), dim=3 ) + ralb_oce * ( 1._wp - at_i_b(:,:) )
         CALL iom_put( 'albedo' , zalb )
         DEALLOCATE( zalb, zmskalb )
      ENDIF
      !
      ! --- category-dependent fields --- !
      IF( iom_use('icemask_cat' ) )   CALL iom_put( 'icemask_cat' ,                        zmsk00c                               ) ! ice mask 0%
      IF( iom_use('iceconc_cat' ) )   CALL iom_put( 'iceconc_cat' , a_i(A2D(0),:)        * zmsk00c                               ) ! area for categories
      IF( iom_use('icethic_cat' ) )   CALL iom_put( 'icethic_cat' , h_i(A2D(0),:)        * zmsk00c + zmiss * ( 1._wp - zmsk00c ) ) ! thickness for categories
      IF( iom_use('snwthic_cat' ) )   CALL iom_put( 'snwthic_cat' , h_s(A2D(0),:)        * zmsksnc + zmiss * ( 1._wp - zmsksnc ) ) ! snow depth for categories
      IF( iom_use('icesalt_cat' ) )   CALL iom_put( 'icesalt_cat' , s_i(A2D(0),:)        * zmsk00c + zmiss * ( 1._wp - zmsk00c ) ) ! salinity for categories
      IF( iom_use('iceage_cat'  ) )   CALL iom_put( 'iceage_cat'  , o_i(A2D(0),:) / rday * zmsk00c + zmiss * ( 1._wp - zmsk00c ) ) ! ice age
      IF( iom_use('icebrv_cat'  ) )   CALL iom_put( 'icebrv_cat'  , v_ibr(:,:,:) * 100.  * zmsk00c + zmiss * ( 1._wp - zmsk00c ) ) ! brine volume
      IF( iom_use('iceapnd_cat' ) )   CALL iom_put( 'iceapnd_cat' , a_ip(A2D(0),:)       * zmsk00c                               ) ! melt pond frac for categories
      IF( iom_use('icevpnd_cat' ) )   CALL iom_put( 'icevpnd_cat' , v_ip(A2D(0),:)       * zmsk00c                               ) ! melt pond volume for categories
      IF( iom_use('icehpnd_cat' ) )   CALL iom_put( 'icehpnd_cat' , h_ip(A2D(0),:)       * zmsk00c + zmiss * ( 1._wp - zmsk00c ) ) ! melt pond thickness for categories
      IF( iom_use('icevlid_cat' ) )   CALL iom_put( 'icevlid_cat' , v_il(A2D(0),:)       * zmsk00c + zmiss * ( 1._wp - zmsk00c ) ) ! melt pond lid volume for categories
      IF( iom_use('icehlid_cat' ) )   CALL iom_put( 'icehlid_cat' , h_il(A2D(0),:)       * zmsk00c + zmiss * ( 1._wp - zmsk00c ) ) ! melt pond lid thickness for categories
      IF( iom_use('iceafpnd_cat') )   CALL iom_put( 'iceafpnd_cat', a_ip_frac(:,:,:)     * zmsk00c                               ) ! melt pond frac per ice area for categories
      IF( iom_use('iceaepnd_cat') )   CALL iom_put( 'iceaepnd_cat', a_ip_eff(A2D(0),:)   * zmsk00c                               ) ! melt pond effective frac for categories
      IF( iom_use('icealb_cat'  ) )   CALL iom_put( 'icealb_cat'  , alb_ice(:,:,:)       * zmsk00c + zmiss * ( 1._wp - zmsk00c ) ) ! ice albedo for categories
      IF( iom_use('icettop_cat' ) )   CALL iom_put( 'icettop_cat' , (t_su(A2D(0),:)-rt0) * zmsk00c + zmiss * ( 1._wp - zmsk00c ) ) ! surface temperature
      IF( iom_use('icetemp_cat' ) )   CALL iom_put( 'icetemp_cat' , (SUM( t_i(A2D(0),:,:), dim=3 ) * r1_nlay_i - rt0) * zmsk00c  &
         &                                                                                         + zmiss * ( 1._wp - zmsk00c ) ) ! ice temperature
      IF( iom_use('snwtemp_cat' ) )   CALL iom_put( 'snwtemp_cat' , (SUM( t_s(A2D(0),:,:), dim=3 ) * r1_nlay_s - rt0) * zmsksnc  &
         &                                                                                         + zmiss * ( 1._wp - zmsksnc ) ) ! snow temperature

      ! --- layer-dependent fields --- !
      IF( iom_use('icetemp_lay' ) .OR. iom_use('icesalt_lay' ) ) THEN
         CALL iom_put( 'icetemp_lay' , (t_i(A2D(0),:,:)-rt0) * zmsk00l + zmiss * ( 1._wp - zmsk00l ) ) ! ice temperature
         CALL iom_put( 'icesalt_lay' , sz_i(A2D(0),:,:)      * zmsk00l + zmiss * ( 1._wp - zmsk00l ) ) ! ice salinity
         DEALLOCATE( zmsk00l )
      ENDIF

      !------------------
      ! Add-ons for SIMIP
      !------------------
      ! trends
      IF( iom_use('dmithd') )   CALL iom_put( 'dmithd', - wfx_bog(:,:) - wfx_bom(:,:) - wfx_sum(:,:) - wfx_sni(:,:) &
         &                                              - wfx_opw(:,:) - wfx_lam(:,:) - wfx_res(:,:)                        ) ! Sea-ice mass change from thermodynamics
      IF( iom_use('dmidyn') )   CALL iom_put( 'dmidyn', - wfx_dyn + rhoi * diag_trp_vi                                      ) ! Sea-ice mass change from dynamics(kg/m2/s)
      IF( iom_use('dmiopw') )   CALL iom_put( 'dmiopw', - wfx_opw                                                           ) ! Sea-ice mass change through growth in open water
      IF( iom_use('dmibog') )   CALL iom_put( 'dmibog', - wfx_bog                                                           ) ! Sea-ice mass change through basal growth
      IF( iom_use('dmisni') )   CALL iom_put( 'dmisni', - wfx_sni                                                           ) ! Sea-ice mass change through snow-to-ice conversion
      IF( iom_use('dmisum') )   CALL iom_put( 'dmisum', - wfx_sum                                                           ) ! Sea-ice mass change through surface melting
      IF( iom_use('dmibom') )   CALL iom_put( 'dmibom', - wfx_bom                                                           ) ! Sea-ice mass change through bottom melting
      IF( iom_use('dmilam') )   CALL iom_put( 'dmilam', - wfx_lam                                                           ) ! Sea-ice mass change through lateral melting
      IF( iom_use('dmtsub') )   CALL iom_put( 'dmtsub', - wfx_sub                                                           ) ! Sea-ice mass change through evaporation and sublimation
      IF( iom_use('dmssub') )   CALL iom_put( 'dmssub', - wfx_snw_sub                                                       ) ! Snow mass change through sublimation
      IF( iom_use('dmisub') )   CALL iom_put( 'dmisub', - wfx_ice_sub                                                       ) ! Sea-ice mass change through sublimation
      IF( iom_use('dmsspr') )   CALL iom_put( 'dmsspr', - wfx_spr                                                           ) ! Snow mass change through snow fall
      IF( iom_use('dmsssi') )   CALL iom_put( 'dmsssi',   wfx_sni*rhos*r1_rhoi                                              ) ! Snow mass change through snow-to-ice conversion
      IF( iom_use('dmsmel') )   CALL iom_put( 'dmsmel', - wfx_snw_sum                                                       ) ! Snow mass change through melt
      IF( iom_use('dmsdyn') )   CALL iom_put( 'dmsdyn', - wfx_snw_dyn + rhos * diag_trp_vs                                  ) ! Snow mass change through dynamics(kg/m2/s)

      ! Global ice diagnostics
      IF(  iom_use('NH_icearea') .OR. iom_use('NH_icevolu') .OR. iom_use('NH_iceextt') .OR. &
         & iom_use('SH_icearea') .OR. iom_use('SH_icevolu') .OR. iom_use('SH_iceextt') ) THEN
         !
         WHERE( ff_t(A2D(0)) > 0._wp )   ;   z2d(:,:) = 1._wp
         ELSEWHERE                       ;   z2d(:,:) = 0.
         END WHERE
         !
         ztmp(:,:,1) = at_i(A2D(0)) *           z2d   * e1e2t(A2D(0)) * 1.e-12
         ztmp(:,:,2) = vt_i(A2D(0)) *           z2d   * e1e2t(A2D(0)) * 1.e-12
         ztmp(:,:,3) =                          z2d   * e1e2t(A2D(0)) * 1.e-12 * zmsk15
         ztmp(:,:,4) = at_i(A2D(0)) * ( 1._wp - z2d ) * e1e2t(A2D(0)) * 1.e-12 
         ztmp(:,:,5) = vt_i(A2D(0)) * ( 1._wp - z2d ) * e1e2t(A2D(0)) * 1.e-12
         ztmp(:,:,6) =                ( 1._wp - z2d ) * e1e2t(A2D(0)) * 1.e-12 * zmsk15
         !
         ysum(:) = local_2Dsum( ztmp )
         CALL mpp_sum( 'icewri', ysum, cdelay = 'icediags' )
         zdiag(:) = REAL(ysum, wp)
         !
         CALL iom_put( 'NH_icearea' , zdiag(1) )
         CALL iom_put( 'NH_icevolu' , zdiag(2) )
         CALL iom_put( 'NH_iceextt' , zdiag(3) )
         CALL iom_put( 'SH_icearea' , zdiag(4) )
         CALL iom_put( 'SH_icevolu' , zdiag(5) )
         CALL iom_put( 'SH_iceextt' , zdiag(6) )
         !
      ENDIF
      !
!!CR      !     !  Create an output files (output.lim.abort.nc) if S < 0 or u > 20 m/s
!!CR      !     IF( kindic < 0 )   CALL ice_wri_state( 'output.abort' )
!!CR      !     not yet implemented
!!gm  idem for the ocean...  Ask Seb how to get rid of ioipsl....
      !
      IF( ln_timing )  CALL timing_stop('icewri')
      !
   END SUBROUTINE ice_wri


   SUBROUTINE ice_wri_state( kid )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE ice_wri_state  ***
      !!
      !! ** Purpose :   create a NetCDF file named cdfile_name which contains
      !!      the instantaneous ice state and forcing fields for ice model
      !!        Used to find errors in the initial state or save the last
      !!      ocean state in case of abnormal end of a simulation
      !!
      !! History :   4.0  !  2013-06  (C. Rousset)
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kid
      !!----------------------------------------------------------------------
      !
      !! The file is open in dia_wri_state (ocean routine)

      CALL iom_rstput( 0, 0, kid, 'sithic'  , hm_i     )   ! Ice thickness
      CALL iom_rstput( 0, 0, kid, 'siconc'  , at_i     )   ! Ice concentration
      CALL iom_rstput( 0, 0, kid, 'sitemp'  , tm_i     )   ! Ice temperature
      CALL iom_rstput( 0, 0, kid, 'sittop'  , t_su     )   ! Surface Ice temperature
      CALL iom_rstput( 0, 0, kid, 'sitbot'  , t_bo     )   ! Bottom Ice temperature
      CALL iom_rstput( 0, 0, kid, 'sivelu'  , u_ice    )   ! i-Ice speed
      CALL iom_rstput( 0, 0, kid, 'sivelv'  , v_ice    )   ! j-Ice speed
      CALL iom_rstput( 0, 0, kid, 'utau_ice', utau_ice )   ! i-Wind stress over ice
      CALL iom_rstput( 0, 0, kid, 'vtau_ice', vtau_ice )   ! i-Wind stress over ice
      CALL iom_rstput( 0, 0, kid, 'snowpre' , sprecip  )   ! Snow precipitation
      CALL iom_rstput( 0, 0, kid, 'sisali'  , sm_i     )   ! Ice salinity
      CALL iom_rstput( 0, 0, kid, 'sivolu'  , vt_i     )   ! Ice volume
      CALL iom_rstput( 0, 0, kid, 'siapnd'  , at_ip    )   ! Melt pond fraction
      CALL iom_rstput( 0, 0, kid, 'sivpnd'  , vt_ip    )   ! Melt pond volume
      CALL iom_rstput( 0, 0, kid, 'sithicat', h_i      )   ! Ice thickness
      CALL iom_rstput( 0, 0, kid, 'siconcat', a_i      )   ! Ice concentration
      CALL iom_rstput( 0, 0, kid, 'sisalcat', s_i      )   ! Ice salinity
      CALL iom_rstput( 0, 0, kid, 'snthicat', h_s      )   ! Snw thickness
      CALL iom_rstput( 0, 0, kid, 'qsr'     , qsr      )   ! Solar flx over ocean
      CALL iom_rstput( 0, 0, kid, 'qns'     , qns      )   ! NonSolar flx over ocean
      CALL iom_rstput( 0, 0, kid, 'sst'     , sst_m    )   ! sst
      CALL iom_rstput( 0, 0, kid, 'sss'     , sss_m    )   ! sss

    END SUBROUTINE ice_wri_state

#else
   !!----------------------------------------------------------------------
   !!   Default option :         Empty module         NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icewri
