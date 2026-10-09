MODULE icevar
   !!======================================================================
   !!                       ***  MODULE icevar ***
   !!   sea-ice:  series of functions to transform or compute ice variables
   !!======================================================================
   !! History :   -   !  2006-01  (M. Vancoppenolle) Original code
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!
   !!                 There are three sets of variables
   !!                 VGLO : global variables of the model
   !!                        - v_i (jpi,jpj,jpl)
   !!                        - v_s (jpi,jpj,jpl)
   !!                        - a_i (jpi,jpj,jpl)
   !!                        - t_s (jpi,jpj,jpl)
   !!                        - e_i (jpi,jpj,nlay_i,jpl)
   !!                        - e_s (jpi,jpj,nlay_s,jpl)
   !!                        - sv_i(jpi,jpj,jpl)
   !!                        - oa_i(jpi,jpj,jpl)
   !!                 VEQV : equivalent variables sometimes used in the model
   !!                        - h_i(jpi,jpj,jpl)
   !!                        - h_s(jpi,jpj,jpl)
   !!                        - t_i(jpi,jpj,nlay_i,jpl)
   !!                        ...
   !!                 VAGG : aggregate variables, averaged/summed over all
   !!                        thickness categories
   !!                        - vt_i(jpi,jpj)
   !!                        - vt_s(jpi,jpj)
   !!                        - at_i(jpi,jpj)
   !!                        - st_i(jpi,jpj)
   !!                        - et_s(jpi,jpj)  total snow heat content
   !!                        - et_i(jpi,jpj)  total ice thermal content
   !!                        - sm_i(jpi,jpj)  mean ice salinity
   !!                        - tm_i(jpi,jpj)  mean ice temperature
   !!                        - tm_s(jpi,jpj)  mean snw temperature
   !!----------------------------------------------------------------------
   !!   ice_var_agg       : integrate variables over layers and categories
   !!   ice_var_glo2eqv   : transform from VGLO to VEQV
   !!   ice_var_eqv2glo   : transform from VEQV to VGLO
   !!   ice_var_salprof   : salinity profile in the ice
   !!   ice_var_salprof1d : salinity profile in the ice 1D
   !!   ice_var_zapsmall  : remove very small area and volume
   !!   ice_var_zapneg    : remove negative ice fields
   !!   ice_var_roundoff  : remove negative values arising from roundoff erros
   !!   ice_var_brine     : brine volume
   !!   ice_var_enthalpy  : compute ice and snow enthalpies from temperature
   !!   ice_var_sshdyn    : compute equivalent ssh in lead
   !!   ice_var_itd       : convert N-cat to M-cat
   !!   ice_var_snwfra    : fraction of ice covered by snow
   !!   ice_var_snwblow   : distribute snow fall between ice and ocean
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE phycst         ! physical constants (ocean directory)
   USE sbc_oce , ONLY : sss_m, sst_m, ln_ice_embd, nn_fsbc
   USE ice            ! sea-ice: variables
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_var_agg
   PUBLIC   ice_var_glo2eqv
!!$   PUBLIC   ice_var_eqv2glo
   PUBLIC   ice_var_salprof
!   PUBLIC   ice_var_salprof1d
   PUBLIC   ice_var_zapsmall
   PUBLIC   ice_var_zapneg
   PUBLIC   ice_var_roundoff
   PUBLIC   ice_var_brine
   PUBLIC   ice_var_enthalpy
   PUBLIC   ice_var_vremap
   PUBLIC   snw_var_vremap
   PUBLIC   ice_var_sshdyn
   PUBLIC   ice_var_itd
   PUBLIC   ice_var_snwfra
   PUBLIC   ice_var_snwblow

   INTERFACE ice_var_itd
      MODULE PROCEDURE ice_var_itd_1c1c, ice_var_itd_Nc1c, ice_var_itd_1cMc, ice_var_itd_NcMc
   END INTERFACE

   !! * Substitutions
#  include "do_loop_substitute.h90"

   INTERFACE ice_var_snwfra
      MODULE PROCEDURE ice_var_snwfra_2d, ice_var_snwfra_3d
   END INTERFACE

   INTERFACE ice_var_snwblow
      MODULE PROCEDURE ice_var_snwblow_2d
   END INTERFACE

   !!----------------------------------------------------------------------
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_var_agg( kn )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_var_agg  ***
      !!
      !! ** Purpose :   aggregates ice-thickness-category variables to
      !!              all-ice variables, i.e. it turns VGLO into VAGG
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kn     ! =1 state variables only
      !                                 ! >1 state variables + others
      !
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      REAL(wp) ::   z1_vt_i, z1_vt_s, z1_at_i
      !!-------------------------------------------------------------------
      !
      ! full    arrays: vt_i, vt_s, at_i, vt_ip, vt_il, at_ip
      ! reduced arrays: the rest
      !
      ! --- integrated values
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         vt_i(ji,jj)  = SUM( v_i (ji,jj,:) )
         vt_s(ji,jj)  = SUM( v_s (ji,jj,:) )
         at_i(ji,jj)  = SUM( a_i (ji,jj,:) )
         !
         at_ip    (ji,jj) = SUM( a_ip    (ji,jj,:) ) ! melt ponds
         vt_ip    (ji,jj) = SUM( v_ip    (ji,jj,:) )
         vt_il    (ji,jj) = SUM( v_il    (ji,jj,:) )
         !
         ato_i(ji,jj) = 1._wp - at_i(ji,jj)  ! open water fraction
      END_2D
      !
      DO_2D( 0, 0, 0, 0 )
         et_s(ji,jj)  = SUM( SUM( e_s (ji,jj,:,:), dim=2 ) )
         et_i(ji,jj)  = SUM( SUM( e_i (ji,jj,:,:), dim=2 ) )
         !
         at_ip_eff(ji,jj) = SUM( a_ip_eff(ji,jj,:) * a_i(ji,jj,:) )
         !
         !!GS: tm_su always needed by ABL over sea-ice
         IF( at_i(ji,jj) <= epsi20 ) THEN
            tm_su  (ji,jj) = rt0
         ELSE
            tm_su  (ji,jj) = SUM( t_su(ji,jj,:) * a_i(ji,jj,:) ) / at_i(ji,jj)
         ENDIF
      END_2D
      !
      IF( nn_icesal == 4 ) THEN
         DO_2D( 0, 0, 0, 0 )
            st_i(ji,jj) = SUM( SUM( szv_i (ji,jj,:,:), dim=2 ) )
         END_2D
      ELSE
         DO_2D( 0, 0, 0, 0 )
            st_i(ji,jj) = SUM( sv_i(ji,jj,:) )
         END_2D
      ENDIF
      
      ! The following fields are calculated for diagnostics and outputs only
      ! ==> Do not use them for other purposes
      IF( kn > 1 ) THEN
         !
         DO_2D( 0, 0, 0, 0 )
            IF( at_i(ji,jj) > epsi20 ) THEN   ;   z1_at_i = 1._wp / at_i(ji,jj)
            ELSE                              ;   z1_at_i = 0._wp
            ENDIF
            IF( vt_i(ji,jj) > epsi20 ) THEN   ;   z1_vt_i = 1._wp / vt_i(ji,jj)
            ELSE                              ;   z1_vt_i = 0._wp
            ENDIF

            ! mean ice/snow thickness
            hm_i(ji,jj) = vt_i(ji,jj) * z1_at_i
            hm_s(ji,jj) = vt_s(ji,jj) * z1_at_i
            !
            ! mean temperature (K), salinity and age
            tm_si(ji,jj) = SUM( t_si(ji,jj,:) * a_i(ji,jj,:)  ) * z1_at_i
            om_i (ji,jj) = SUM( oa_i(ji,jj,:)                 ) * z1_at_i
            sm_i (ji,jj) =      st_i(ji,jj)                     * z1_vt_i
         END_2D
            !
         tm_i(:,:) = 0._wp
         tm_s(:,:) = 0._wp
         DO jl = 1, jpl
            DO_3D( 0, 0, 0, 0, 1, nlay_i )
               IF( vt_i(ji,jj) > epsi20 ) THEN
                  tm_i(ji,jj) = tm_i(ji,jj) + r1_nlay_i * t_i (ji,jj,jk,jl) * v_i(ji,jj,jl) / vt_i(ji,jj)
               ELSE
                  tm_i(ji,jj) = rt0
               ENDIF
            END_3D
         END DO
         DO jl = 1, jpl
            DO_3D( 0, 0, 0, 0, 1, nlay_s )
               IF( vt_s(ji,jj) > epsi20 ) THEN
                  tm_s(ji,jj) = tm_s(ji,jj) + r1_nlay_s * t_s (ji,jj,jk,jl) * v_s(ji,jj,jl) / vt_s(ji,jj)
               ELSE
                  tm_s(ji,jj) = rt0
               ENDIF
            END_3D
         END DO
            !            
!!$         DO_2D( 0, 0, 0, 0 )
!!$            IF( at_i(ji,jj) > epsi20 ) THEN   ;   z1_at_i = 1._wp / at_i(ji,jj)
!!$            ELSE                              ;   z1_at_i = 0._wp
!!$            ENDIF
!!$            IF( vt_i(ji,jj) > epsi20 ) THEN   ;   z1_vt_i = 1._wp / vt_i(ji,jj)
!!$            ELSE                              ;   z1_vt_i = 0._wp
!!$            ENDIF
!!$            IF( vt_s(ji,jj) > epsi20 ) THEN   ;   z1_vt_s = 1._wp / vt_s(ji,jj)
!!$            ELSE                              ;   z1_vt_s = 0._wp
!!$            ENDIF
!!$
!!$            ! mean ice/snow thickness
!!$            hm_i(ji,jj) = vt_i(ji,jj) * z1_at_i
!!$            hm_s(ji,jj) = vt_s(ji,jj) * z1_at_i
!!$            !
!!$            ! mean temperature (K), salinity and age
!!$            tm_si(ji,jj) = SUM( t_si(ji,jj,:) * a_i(ji,jj,:)  ) * z1_at_i
!!$            om_i (ji,jj) = SUM( oa_i(ji,jj,:)                 ) * z1_at_i
!!$            sm_i (ji,jj) =      st_i(ji,jj)                     * z1_vt_i
!!$            !
!!$            tm_i(ji,jj) = 0._wp
!!$            tm_s(ji,jj) = 0._wp
!!$            DO jl = 1, jpl
!!$               DO jk = 1, nlay_i
!!$                  tm_i(ji,jj) = tm_i(ji,jj) + r1_nlay_i * t_i (ji,jj,jk,jl) * v_i(ji,jj,jl) * z1_vt_i
!!$               END DO
!!$               DO jk = 1, nlay_s
!!$                  tm_s(ji,jj) = tm_s(ji,jj) + r1_nlay_s * t_s (ji,jj,jk,jl) * v_s(ji,jj,jl) * z1_vt_s
!!$               END DO
!!$            END DO
!!$            !            
!!$         END_2D
         ! put rt0 where there is no ice
         WHERE( at_i(A2D(0)) <= epsi20 )
            tm_si(:,:) = rt0
            tm_i (:,:) = rt0
            tm_s (:,:) = rt0
         END WHERE
         !
         ! mean melt pond depth
         WHERE( at_ip(A2D(0)) > epsi20 )
            hm_ip(:,:) = vt_ip(A2D(0)) / at_ip(A2D(0))
            hm_il(:,:) = vt_il(A2D(0)) / at_ip(A2D(0))
         ELSEWHERE
            hm_ip(:,:) = 0._wp
            hm_il(:,:) = 0._wp
         END WHERE
         !
      ENDIF
      !
   END SUBROUTINE ice_var_agg


   SUBROUTINE ice_var_glo2eqv( kn )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_var_glo2eqv ***
      !!
      !! ** Purpose :   computes equivalent variables as function of
      !!              global variables, i.e. it turns VGLO into VEQV
      !!-------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kn     ! =1 everything including ponds (necessary for init)
      !                                 ! =2            excluding ponds if ln_pnd=F
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      REAL(wp) ::   ze_i, zs_i                      ! local scalars
      REAL(wp) ::   ze_s, ztmelts, zbbb, zccc       !   -      -
      REAL(wp) ::   zhmax, z1_hmax                  !   -      -
      REAL(wp) ::   zlay_i, zlay_s                  !   -      -
      REAL(wp) ::   z1_hl, z1_a_i, z1_a_ip
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   za_s_fra
      !!-------------------------------------------------------------------

!!gm Question 2:  It is possible to define existence of sea-ice in a common way between
!!                ice area and ice volume ?
!!                the idea is to be able to define one for all at the begining of this routine
!!                a criteria for icy area (i.e. a_i > epsi20 and v_i > epsi20 )

      !---------------------------------------------------------------
      ! Ice thickness, snow thickness, ice salinity, ice age and ponds
      !---------------------------------------------------------------
      !
      ! bound h_i by hi_max (i.e. 99 m) with associated update of ice area
      ! clem: if a>1 then do something
      zhmax   =          hi_max(jpl)
      z1_hmax =  1._wp / hi_max(jpl)
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         IF( v_i(ji,jj,jpl) > ( zhmax*a_i(ji,jj,jpl) ) )   a_i(ji,jj,jpl) = MIN( 1._wp, v_i(ji,jj,jpl) * z1_hmax )
      END_2D
      
      DO jl = 1, jpl
         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
            !                                            !--- inverse of the ice area
            IF( a_i(ji,jj,jl) > epsi20 ) THEN   ;   z1_a_i = 1._wp / a_i(ji,jj,jl)
            ELSE                                ;   z1_a_i = 0._wp
            ENDIF
            !                                            !--- ice thickness
            h_i(ji,jj,jl) = v_i (ji,jj,jl) * z1_a_i
            !                                            !--- snow thickness
            h_s(ji,jj,jl) = v_s (ji,jj,jl) * z1_a_i
            !                                            !--- ice age
            o_i(ji,jj,jl) = oa_i(ji,jj,jl) * z1_a_i
            !
         END_2D
      ENDDO
      !
      IF( kn == 1 .OR. ln_pnd ) THEN
         !
         z1_hl = 1._wp / ( rn_pnd_hl_max - rn_pnd_hl_min )
         DO jl = 1, jpl
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               IF( a_ip(ji,jj,jl) > epsi20 ) THEN   ;   z1_a_ip = 1._wp / a_ip(ji,jj,jl)
               ELSE                                 ;   z1_a_ip = 0._wp
               ENDIF
               !                                         !--- pond and lid thickness
               h_ip(ji,jj,jl) = v_ip(ji,jj,jl) * z1_a_ip
               h_il(ji,jj,jl) = v_il(ji,jj,jl) * z1_a_ip
            END_2D
            !                                            !--- melt pond effective area (used for albedo)
            DO_2D( 0, 0, 0, 0 )
               IF( a_i(ji,jj,jl) > epsi20 ) THEN   ;   a_ip_frac(ji,jj,jl) = a_ip(ji,jj,jl) / a_i(ji,jj,jl)
               ELSE                                ;   a_ip_frac(ji,jj,jl) = 0._wp
               ENDIF
               IF    ( h_il(ji,jj,jl) <= rn_pnd_hl_min ) THEN   ;   a_ip_eff(ji,jj,jl) = a_ip_frac(ji,jj,jl)       ! lid is very thin.  Expose all the pond
               ELSEIF( h_il(ji,jj,jl) >= rn_pnd_hl_max ) THEN   ;   a_ip_eff(ji,jj,jl) = 0._wp                     ! lid is very thick. Cover all the pond up with ice and snow
               ELSE                                             ;   a_ip_eff(ji,jj,jl) = a_ip_frac(ji,jj,jl) * &   ! lid is in between. Expose part of the pond
                  &                                                                    ( rn_pnd_hl_max - h_il(ji,jj,jl) ) * z1_hl 

               ENDIF
               IF( h_ip(ji,jj,jl) < h_s(ji,jj,jl) )   a_ip_eff(ji,jj,jl) = 0._wp
               !
            END_2D
         ENDDO
         !
         ALLOCATE( za_s_fra(A2D(nn_hls),jpl) )   ! must have the same shape as h_s
         CALL ice_var_snwfra( h_s(:,:,:), za_s_fra(:,:,:) )               ! calculate ice fraction covered by snow
         a_ip_eff(A2D(0),:) = MIN( a_ip_eff(A2D(0),:), 1._wp - za_s_fra(A2D(0),:) )   ! make sure (a_ip_eff + a_s_fra) <= 1
         !
         DEALLOCATE( za_s_fra )
      ENDIF
      !
      !-------------------
      ! Ice salinity         (with a min value rn_simin and a max value rn_sinew*sss)
      !-------------------
      IF( nn_icesal == 1 .OR. nn_icesal == 3 ) THEN
         !
         CALL ice_var_salprof   ! salinity profile
         !
      ELSEIF( nn_icesal == 2 ) THEN
         !
         DO jl = 1, jpl
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               IF( v_i(ji,jj,jl) > epsi20 ) THEN
!!clem test                  s_i(ji,jj,jl) = MAX( rn_simin , MIN( rn_sinew * sss_m(ji,jj), sv_i(ji,jj,jl) / v_i(ji,jj,jl) ) )
                  s_i(ji,jj,jl) = sv_i(ji,jj,jl) / v_i(ji,jj,jl)
               ELSE
                  s_i(ji,jj,jl) = MIN( rn_simin, rn_sinew * sss_m(ji,jj) )
               ENDIF
            END_2D
         ENDDO
         CALL ice_var_salprof   ! salinity profile
         !
      ELSEIF( nn_icesal == 4 ) THEN
         !
         s_i(:,:,:) = 0._wp
         zlay_i = REAL( nlay_i , wp )    ! number of layers
         DO jl = 1, jpl
            DO_3D( nn_hls, nn_hls, nn_hls, nn_hls, 1, nlay_i )
               IF ( v_i(ji,jj,jl) > epsi20 ) THEN     !--- icy area
                  zs_i = szv_i(ji,jj,jk,jl) / ( v_i(ji,jj,jl) * r1_nlay_i )
!!clem test                  sz_i(ji,jj,jk,jl) = MAX( rn_simin , MIN( rn_sinew * sss_m(ji,jj), zs_i ) )
                  sz_i(ji,jj,jk,jl) = zs_i
               ELSE                                   !--- no ice
                  sz_i(ji,jj,jk,jl) = MIN( rn_simin, rn_sinew * sss_m(ji,jj) )
               ENDIF
               s_i(ji,jj,jl) = s_i(ji,jj,jl) + sz_i(ji,jj,jk,jl) * r1_nlay_i
            END_3D
         END DO
         !         
      ENDIF
         
      !-------------------
      ! Ice temperature   [K]   (with a minimum value (rt0 - 100.))
      !-------------------
      zlay_i = REAL( nlay_i , wp )    ! number of layers
      DO jl = 1, jpl
         DO jk = 1, nlay_i
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               IF ( v_i(ji,jj,jl) > epsi20 ) THEN     !--- icy area
                  !
                  ze_i             =   e_i (ji,jj,jk,jl) / v_i(ji,jj,jl) * zlay_i             ! Energy of melting e(S,T) [J.m-3]
                  ztmelts          = - sz_i(ji,jj,jk,jl) * rTmlt                              ! Ice layer melt temperature [C]
                  ! Conversion q(S,T) -> T (second order equation)
                  zbbb             = ( rcp - rcpi ) * ztmelts + ze_i * r1_rhoi - rLfus
                  zccc             = SQRT( MAX( zbbb * zbbb - 4._wp * rcpi * rLfus * ztmelts , 0._wp) )
                  t_i(ji,jj,jk,jl) = MAX( -100._wp , MIN( -( zbbb + zccc ) * 0.5_wp * r1_rcpi , ztmelts ) ) + rt0   ! [K] with bounds: -100 < t_i < ztmelts

                  !
               ELSE                                   !--- no ice
                  t_i(ji,jj,jk,jl) = rt0
               ENDIF
            END_2D
         END DO
      END DO

      !--------------------
      ! Snow temperature   [K]   (with a minimum value (rt0 - 100.))
      !--------------------
      zlay_s = REAL( nlay_s , wp )
      DO jl = 1, jpl
         DO jk = 1, nlay_s
            DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
               IF ( v_s(ji,jj,jl) > epsi20 ) THEN     !--- icy area
                  t_s(ji,jj,jk,jl) = rt0 + MAX( -100._wp ,  &
                     &               MIN( r1_rcpi*( -r1_rhos*( e_s(ji,jj,jk,jl) / v_s(ji,jj,jl) * zlay_s ) + rLfus ) , 0._wp ) )
               ELSE                                   !--- no ice
                  t_s(ji,jj,jk,jl) = rt0
               ENDIF
            END_2D
         END DO
      END DO
      !
      ! integrated values
      vt_i (:,:) = SUM( v_i , dim=3 )
      vt_ip(:,:) = SUM( v_ip, dim=3 )
      vt_s (:,:) = SUM( v_s , dim=3 )
      at_i (:,:) = SUM( a_i , dim=3 )
      at_ip(:,:) = SUM( a_ip, dim=3 )
      !
   END SUBROUTINE ice_var_glo2eqv


!!$   SUBROUTINE ice_var_eqv2glo
!!$      !!-------------------------------------------------------------------
!!$      !!                ***  ROUTINE ice_var_eqv2glo ***
!!$      !!
!!$      !! ** Purpose :   computes global variables as function of
!!$      !!              equivalent variables,  i.e. it turns VEQV into VGLO
!!$      !!-------------------------------------------------------------------
!!$      INTEGER  ::   ji, jj, jl   ! dummy loop indices
!!$      !!-------------------------------------------------------------------
!!$      !
!!$      DO jl = 1, jpl
!!$         DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
!!$            v_i (ji,jj,jl) = h_i (ji,jj,jl) * a_i (ji,jj,jl)
!!$            v_s (ji,jj,jl) = h_s (ji,jj,jl) * a_i (ji,jj,jl)
!!$            sv_i(ji,jj,jl) = s_i (ji,jj,jl) * v_i (ji,jj,jl)
!!$            v_ip(ji,jj,jl) = h_ip(ji,jj,jl) * a_ip(ji,jj,jl)
!!$            v_il(ji,jj,jl) = h_il(ji,jj,jl) * a_ip(ji,jj,jl)
!!$         END_2D
!!$      ENDDO
!!$      !
!!$   END SUBROUTINE ice_var_eqv2glo


   SUBROUTINE ice_var_salprof
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_var_salprof ***
      !!
      !! ** Purpose :   computes salinity profile in function of bulk salinity
      !!
      !! ** Method  : If bulk salinity greater than zsi1,
      !!              the profile is assumed to be constant (S_inf)
      !!              If bulk salinity lower than zsi0,
      !!              the profile is linear with 0 at the surface (S_zero)
      !!              If it is between zsi0 and zsi1, it is a
      !!              alpha-weighted linear combination of s_inf and s_zero
      !!
      !! ** References : Vancoppenolle et al., 2007
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop index
      REAL(wp) ::   z1_dS, ztmp1, ztmp2, zalpha
      REAL(wp), PARAMETER ::   zsi0 = 3.5_wp
      REAL(wp), PARAMETER ::   zsi1 = 4.5_wp
      !!-------------------------------------------------------------------

!!gm Question: Remove the option 3 ?  How many years since it last use ?

      SELECT CASE ( nn_icesal )
      !
      !               !---------------------------------------!
      CASE( 1 )       !  constant salinity in time and space  !
         !            !---------------------------------------!
         sz_i(:,:,:,:) = rn_icesal
         s_i (:,:,:)   = rn_icesal
         !
         !            !---------------------------------------------!
      CASE( 2 , 4 )   !  time varying salinity with linear profile  !
         !            !---------------------------------------------!
         z1_dS = 1._wp / ( zsi1 - zsi0 )
         !
         DO jl = 1, jpl
            DO jk = 1, nlay_i
               DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                  !
                  IF    ( s_i(ji,jj,jl) >= zsi1 ) THEN    ;    zalpha = 0._wp
                  ELSEIF( s_i(ji,jj,jl) <= zsi0 ) THEN    ;    zalpha = 1._wp
                  ELSE                                    ;    zalpha = ( zsi1 - s_i(ji,jj,jl) ) * z1_dS
                  ENDIF
                  !
                  IF( s_i(ji,jj,jl) >=                             ( 0.5_wp * rn_sinew * sss_m(ji,jj) ) )   zalpha = 0._wp ! force constant profile when SSS too low (Baltic Sea)
                  IF( s_i(ji,jj,jl) <= ( REAL( nlay_i , wp ) * MIN( rn_simin, rn_sinew * sss_m(ji,jj) ) ) ) zalpha = 0._wp ! force constant profile when ice surface salinity too small
                  !                                                                                                        !       it depends on nlay_i which is not ideal
                  !
                  ! linear profile with 0 surface value
                  sz_i(ji,jj,jk,jl) =             zalpha   * s_i(ji,jj,jl) * 2._wp * ( REAL(jk,wp) - 0.5_wp ) * r1_nlay_i &
                     &                + ( 1._wp - zalpha ) * s_i(ji,jj,jl)
                  !  => mean(sz_i(jk)) = s_i
               END_2D
            ENDDO
         ENDDO
         !
         !            !-------------------------------------------!
      CASE( 3 )       ! constant salinity with a fix profile      ! (Schwarzacher (1959) multiyear salinity profile
         !            !-------------------------------------------!                                   (mean = 2.30)
         !
         s_i(:,:,:) = 2.30_wp
!!gm Remark: if we keep the case 3, then compute an store one for all time-step
!!           a array  S_prof(1:nlay_i) containing the calculation and just do:
!         DO jk = 1, nlay_i
!            sz_i(:,:,jk,:) = S_prof(jk)
!         END DO
!!gm end
         DO jl = 1, jpl
            DO jk = 1, nlay_i
               ztmp1 = ( REAL(jk,wp) - 0.5_wp ) * r1_nlay_i
               ztmp2 = 1.6_wp * (  1._wp - COS( rpi * ztmp1**(0.407_wp/(0.573_wp+ztmp1)) ) )
               DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
                  sz_i(ji,jj,jk,jl) =  ztmp2
               END_2D
            END DO
         END DO
          !
      END SELECT
      !
   END SUBROUTINE ice_var_salprof


   SUBROUTINE ice_var_zapsmall
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_zapsmall ***
      !!
      !! ** Purpose :   Remove too small sea ice areas and correct fluxes
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jl, jk   ! dummy loop indices
      REAL(wp) ::   zsmall
      !!-------------------------------------------------------------------
      !
      WHERE( a_i(A2D(0),:) > epsi10 )   ;   h_i(A2D(0),:) = v_i(A2D(0),:) / a_i(A2D(0),:)
      ELSEWHERE                         ;   h_i(A2D(0),:) = 0._wp
      END WHERE
      !
      !-----------------------------------------------------------------
      ! Zap ice volume, add salt to ocean
      !-----------------------------------------------------------------
      IF( nn_icesal == 4 ) THEN
         DO jl = 1, jpl
            DO_3D( 0, 0, 0, 0, 1, nlay_i )
               !
               zsmall = MIN( a_i(ji,jj,jl), v_i(ji,jj,jl),  h_i(ji,jj,jl) )
               !
               IF( zsmall < epsi10 ) THEN
                  ! update exchanges with ocean
                  sfx_res(ji,jj)  = sfx_res(ji,jj) + szv_i(ji,jj,jk,jl) * rhoi * r1_Dt_ice
                  szv_i(ji,jj,jk,jl) = 0._wp
                  sz_i (ji,jj,jk,jl) = MIN( rn_simin, rn_sinew * sss_m(ji,jj) )
               ENDIF
            END_3D
         ENDDO
      ELSE
         DO jl = 1, jpl
            DO_2D( 0, 0, 0, 0 )
               !
               zsmall = MIN( a_i(ji,jj,jl), v_i(ji,jj,jl),  h_i(ji,jj,jl) )
               !
               IF( zsmall < epsi10 ) THEN
                  ! update exchanges with ocean
                  sfx_res(ji,jj)  = sfx_res(ji,jj) + sv_i(ji,jj,jl)   * rhoi * r1_Dt_ice
                  sv_i (ji,jj,jl) = 0._wp
               ENDIF
            END_2D
         END DO
      ENDIF

      !-----------------------------------------------------------------
      ! Zap ice energy and use ocean heat to melt ice
      !-----------------------------------------------------------------
      DO jl = 1, jpl
         DO_3D( 0, 0, 0, 0, 1, nlay_i )
            !
            zsmall = MIN( a_i(ji,jj,jl), v_i(ji,jj,jl),  h_i(ji,jj,jl) )
            !
            IF( zsmall < epsi10 ) THEN
               ! update exchanges with ocean
               hfx_res(ji,jj)   = hfx_res(ji,jj) - e_i(ji,jj,jk,jl) * r1_Dt_ice ! W.m-2 <0
               e_i(ji,jj,jk,jl) = 0._wp
               t_i(ji,jj,jk,jl) = rt0
            ENDIF
         END_3D
      ENDDO
      
      DO jl = 1, jpl
         DO_3D( 0, 0, 0, 0, 1, nlay_s )
            !
            zsmall = MIN( a_i(ji,jj,jl), v_i(ji,jj,jl),  h_i(ji,jj,jl) )
            !
            IF( zsmall < epsi10 ) THEN
               ! update exchanges with ocean
               hfx_res(ji,jj)   = hfx_res(ji,jj) - e_s(ji,jj,jk,jl) * r1_Dt_ice ! W.m-2 <0
               e_s(ji,jj,jk,jl) = 0._wp
               t_s(ji,jj,jk,jl) = rt0
            ENDIF
         END_3D
      ENDDO
      !
      !-----------------------------------------------------------------
      ! zap ice and snow volume, add water to ocean
      !-----------------------------------------------------------------
      DO jl = 1, jpl
         DO_2D( 0, 0, 0, 0 )
            !
            zsmall = MIN( a_i(ji,jj,jl), v_i(ji,jj,jl),  h_i(ji,jj,jl) )
            !
            IF( zsmall < epsi10 ) THEN
               ! update exchanges with ocean
               wfx_res(ji,jj)  = wfx_res(ji,jj) + v_i (ji,jj,jl)   * rhoi * r1_Dt_ice
               wfx_res(ji,jj)  = wfx_res(ji,jj) + v_s (ji,jj,jl)   * rhos * r1_Dt_ice
               wfx_res(ji,jj)  = wfx_res(ji,jj) + ( v_ip(ji,jj,jl)+v_il(ji,jj,jl) ) * rhow * r1_Dt_ice
               !
               a_i  (ji,jj,jl) = 0._wp
               v_i  (ji,jj,jl) = 0._wp
               v_s  (ji,jj,jl) = 0._wp
               t_su (ji,jj,jl) = sst_m(ji,jj) + rt0
               oa_i (ji,jj,jl) = 0._wp
               !
               h_i (ji,jj,jl) = 0._wp
               h_s (ji,jj,jl) = 0._wp
               !
               a_ip (ji,jj,jl) = 0._wp
               v_ip (ji,jj,jl) = 0._wp
               v_il (ji,jj,jl) = 0._wp
               h_ip (ji,jj,jl) = 0._wp
               h_il (ji,jj,jl) = 0._wp
            ENDIF
         END_2D
      END DO
      
      ! to be sure that at_i is the sum of a_i(jl)
      at_i (A2D(0)) = SUM( a_i (A2D(0),:), dim=3 )
      vt_i (A2D(0)) = SUM( v_i (A2D(0),:), dim=3 )
!!clem add?
!      vt_s (:,:) = SUM( v_s (:,:,:), dim=3 )
!      st_i (:,:) = SUM( sv_i(:,:,:), dim=3 )
!      et_s(:,:)  = SUM( SUM( e_s (:,:,:,:), dim=4 ), dim=3 )
!      et_i(:,:)  = SUM( SUM( e_i (:,:,:,:), dim=4 ), dim=3 )
!!clem

      ! open water = 1 if at_i=0
      WHERE( at_i(A2D(0)) == 0._wp )   ato_i(A2D(0)) = 1._wp
      !
   END SUBROUTINE ice_var_zapsmall


   SUBROUTINE ice_var_zapneg( ihls, pdt, pv_i, pv_s, psv_i, poa_i, pa_i, pa_ip, pv_ip, pv_il, pe_s, pe_i, pszv_i )
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_zapneg ***
      !!
      !! ** Purpose :   Remove negative sea ice fields and correct fluxes
      !!-------------------------------------------------------------------
      INTEGER                     , INTENT(in   ) ::   ihls       ! loop index
      REAL(wp)                    , INTENT(in   ) ::   pdt        ! tracer time-step
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pv_i       ! ice volume
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pv_s       ! snw volume
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   psv_i      ! salt content
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   poa_i      ! age content
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pa_i       ! ice concentration
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pa_ip      ! melt pond fraction
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pv_ip      ! melt pond volume
      REAL(wp), DIMENSION(:,:,:)  , INTENT(inout) ::   pv_il      ! melt pond lid volume
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pe_s       ! snw heat content
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pe_i       ! ice heat content
      REAL(wp), DIMENSION(:,:,:,:), INTENT(inout) ::   pszv_i     ! ice salt content
      !
      INTEGER  ::   ji, jj, jl, jk   ! dummy loop indices
      REAL(wp) ::   z1_dt
      REAL(wp), DIMENSION(jpi,jpj) ::   zwfx_res, zhfx_res, zsfx_res ! needed since loop is not (0,0,0,0)
      !!-------------------------------------------------------------------
      !
      DO_2D( ihls, ihls, ihls, ihls )
         zwfx_res(ji,jj) = 0._wp
         zhfx_res(ji,jj) = 0._wp
         zsfx_res(ji,jj) = 0._wp
      END_2D
      
      z1_dt = 1._wp / pdt
      !
      ! make sure a_i=0 where v_i<=0
      WHERE( pv_i(:,:,:) <= 0._wp )   pa_i(:,:,:) = 0._wp
      
      !--------------------------------------
      ! zap ice salt and send it to the ocean
      !--------------------------------------
      IF( nn_icesal == 4 ) THEN
         DO jl = 1, jpl
            DO_3D( ihls, ihls, ihls, ihls, 1, nlay_i )
               IF( pszv_i(ji,jj,jk,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp .OR. pv_i(ji,jj,jl) <= 0._wp ) THEN
                  zsfx_res(ji,jj)     = zsfx_res(ji,jj) + pszv_i(ji,jj,jk,jl) * rhoi * z1_dt
                  pszv_i(ji,jj,jk,jl) = 0._wp
               ENDIF
            END_3D
         ENDDO
      ELSE
         DO jl = 1, jpl
            DO_2D( ihls, ihls, ihls, ihls )
               IF( psv_i(ji,jj,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp .OR. pv_i(ji,jj,jl) <= 0._wp ) THEN
                  zsfx_res(ji,jj)    = zsfx_res(ji,jj) + psv_i(ji,jj,jl) * rhoi * z1_dt
                  psv_i   (ji,jj,jl) = 0._wp
               ENDIF
            END_2D
         END DO
      ENDIF
      !
      !----------------------------------------
      ! zap ice energy and send it to the ocean
      !----------------------------------------
      DO jl = 1, jpl
         DO_3D( ihls, ihls, ihls, ihls, 1, nlay_i )
            IF( pe_i(ji,jj,jk,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp .OR. pv_i(ji,jj,jl) <= 0._wp ) THEN
               zhfx_res(ji,jj)   = zhfx_res(ji,jj) - pe_i(ji,jj,jk,jl) * z1_dt ! W.m-2 >0
               pe_i(ji,jj,jk,jl) = 0._wp
            ENDIF
         END_3D
         !
         DO_3D( ihls, ihls, ihls, ihls, 1, nlay_s )
            IF( pe_s(ji,jj,jk,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp .OR. pv_s(ji,jj,jl) <= 0._wp ) THEN
               zhfx_res(ji,jj)   = zhfx_res(ji,jj) - pe_s(ji,jj,jk,jl) * z1_dt ! W.m-2 <0
               pe_s(ji,jj,jk,jl) = 0._wp
            ENDIF
         END_3D
      ENDDO
      !
      !--------------------------------------------
      ! zap ice and snow volume, add water to ocean
      !--------------------------------------------
      DO jl = 1, jpl
         DO_2D( ihls, ihls, ihls, ihls )
            IF( pv_i(ji,jj,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp ) THEN
               zwfx_res(ji,jj)    = zwfx_res(ji,jj) + pv_i (ji,jj,jl) * rhoi * z1_dt
               pv_i    (ji,jj,jl) = 0._wp
            ENDIF
            IF( pv_s(ji,jj,jl) < 0._wp .OR. pa_i(ji,jj,jl) <= 0._wp ) THEN
               zwfx_res(ji,jj)    = zwfx_res(ji,jj) + pv_s (ji,jj,jl) * rhos * z1_dt
               pv_s    (ji,jj,jl) = 0._wp
            ENDIF
            IF( pv_ip(ji,jj,jl) < 0._wp .OR. pv_il(ji,jj,jl) < 0._wp .OR. pa_ip(ji,jj,jl) <= 0._wp ) THEN
               zwfx_res(ji,jj)    = zwfx_res(ji,jj) + pv_il(ji,jj,jl) * rhow * z1_dt
               pv_il   (ji,jj,jl) = 0._wp
            ENDIF
            IF( pv_ip(ji,jj,jl) < 0._wp .OR. pa_ip(ji,jj,jl) <= 0._wp ) THEN
               zwfx_res(ji,jj)    = zwfx_res(ji,jj) + pv_ip(ji,jj,jl) * rhow * z1_dt
               pv_ip   (ji,jj,jl) = 0._wp
            ENDIF
         END_2D
      END DO
      !
      ! record residual fluxes
      DO_2D( 0, 0, 0, 0 )
         wfx_res(ji,jj) = wfx_res(ji,jj) + zwfx_res(ji,jj)
         hfx_res(ji,jj) = hfx_res(ji,jj) + zhfx_res(ji,jj)
         sfx_res(ji,jj) = sfx_res(ji,jj) + zsfx_res(ji,jj)
      END_2D
      !
      WHERE( poa_i (:,:,:) < 0._wp )   poa_i (:,:,:) = 0._wp
      WHERE( pa_i  (:,:,:) < 0._wp )   pa_i  (:,:,:) = 0._wp
      WHERE( pa_ip (:,:,:) < 0._wp )   pa_ip (:,:,:) = 0._wp
      !
   END SUBROUTINE ice_var_zapneg


   SUBROUTINE ice_var_roundoff( pa_i, pv_i, pv_s, psv_i, poa_i, pa_ip, pv_ip, pv_il, pe_s, pe_i, pszv_i, lice_pres )
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_roundoff ***
      !!
      !! ** Purpose :   Remove negative sea ice values arising from roundoff errors
      !!-------------------------------------------------------------------
      INTEGER ::   ji, jj
      !
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pa_i       ! ice concentration
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pv_i       ! ice volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pv_s       ! ice volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   psv_i      ! salt content
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   poa_i      ! age content
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pa_ip      ! melt pond fraction
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pv_ip      ! melt pond volume
      REAL(wp), DIMENSION(jpi,jpj,jpl)  , INTENT(inout) ::   pv_il      ! melt pond lid volume
      REAL(wp), DIMENSION(jpi,jpj,nlay_s,jpl), INTENT(inout) ::   pe_s       ! snw heat content
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pe_i       ! ice heat content
      REAL(wp), DIMENSION(jpi,jpj,nlay_i,jpl), INTENT(inout) ::   pszv_i     ! ice salt content
      LOGICAL , DIMENSION(A2D(0))       , INTENT(in)    ::   lice_pres  ! ice presence
      !!-------------------------------------------------------------------
      DO_2D( 0, 0, 0, 0 )
         IF( lice_pres(ji,jj) ) THEN
            pa_i (ji,jj,:)   = MAX(pa_i (ji,jj,:),   0._wp)
            pv_i (ji,jj,:)   = MAX(pv_i (ji,jj,:),   0._wp)
            pv_s (ji,jj,:)   = MAX(pv_s (ji,jj,:),   0._wp)
            poa_i(ji,jj,:)   = MAX(poa_i(ji,jj,:),   0._wp)
            pe_i (ji,jj,:,:) = MAX(pe_i (ji,jj,:,:), 0._wp)
            pe_s (ji,jj,:,:) = MAX(pe_s (ji,jj,:,:), 0._wp)
         ENDIF
      END_2D
      IF( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
         DO_2D( 0, 0, 0, 0 )
            IF( lice_pres(ji,jj) ) THEN
               pa_ip(ji,jj,:) = MAX(pa_ip(ji,jj,:), 0._wp)
               pv_ip(ji,jj,:) = MAX(pv_ip(ji,jj,:), 0._wp)
            ENDIF
         END_2D
         IF( ln_pnd_lids ) THEN
            DO_2D( 0, 0, 0, 0 )
               IF( lice_pres(ji,jj) ) THEN
                  pv_il(ji,jj,:) = MAX(pv_il(ji,jj,:), 0._wp)
               ENDIF
            END_2D
         ENDIF
      ENDIF
      IF( nn_icesal == 4 ) THEN
         DO_2D( 0, 0, 0, 0 )
            IF( lice_pres(ji,jj) ) THEN
               pszv_i(ji,jj,:,:) = MAX(pszv_i(ji,jj,:,:), 0._wp)
            ENDIF
         END_2D
      ELSE
         DO_2D( 0, 0, 0, 0 )
            IF( lice_pres(ji,jj) ) THEN
               psv_i(ji,jj,:) = MAX(psv_i(ji,jj,:), 0._wp)
            ENDIF
         END_2D
      ENDIF
      !
   END SUBROUTINE ice_var_roundoff


   SUBROUTINE ice_var_brine
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_var_brine ***
      !!
      !! ** Purpose :   computes brine volume fraction (%)
      !!                         and salinity of the brine in sea ice
      !!
      !! ** Method  : e = - 0.054 * S (ppt) / T (C)
      !!
      !! References : Vancoppenolle et al., JGR, 2007
      !!-------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk, jl   ! dummy loop indices
      REAL(wp) ::   zt1, zt2, zt3, zs_br
      !!-------------------------------------------------------------------
      !
      v_ibr(:,:,:) = 0._wp
      DO jl = 1, jpl
         DO_3D( 0, 0, 0, 0, 1, nlay_i )
            ! brine salinity
            zt1 = t_i(ji,jj,jk,jl) - rt0
            zt2 = zt1 * zt1
            zt3 = zt2 * zt1
            IF    ( nn_liquidus == 1 ) THEN ; zs_br = - zt1 / rTmlt                                      ! --- Linear liquidus
            ELSEIF( nn_liquidus == 2 ) THEN ; zs_br = -18.7_wp * zt1 - 0.519_wp * zt2 - 0.00535_wp * zt3 ! --- 3rd order liquidus, VC19
            ELSEIF( nn_liquidus == 3 ) THEN ; zs_br = -17.6_wp * zt1 - 0.389_wp * zt2 - 0.00362_wp * zt3 ! --- Weast 71 liquidus in RJW14
            ENDIF
            ! brine volume fraction
            IF( zt1 < - epsi10 )   v_ibr(ji,jj,jl) = v_ibr(ji,jj,jl) + r1_nlay_i * sz_i(ji,jj,jk,jl) / zs_br
         END_3D
      ENDDO
      !
      ! mean brine volume fraction
      DO_2D( 0, 0, 0, 0 )
         IF( vt_i(ji,jj) > epsi20 ) THEN
            vm_ibr(ji,jj) = SUM( v_ibr(ji,jj,:) * v_i(ji,jj,:) ) / vt_i(ji,jj)
         ELSE
            vm_ibr(ji,jj) = 0._wp
         ENDIF
      END_2D
      !
   END SUBROUTINE ice_var_brine

   
   SUBROUTINE ice_var_enthalpy( jl_cat )
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_enthalpy ***
      !!
      !! ** Purpose :   Computes sea ice energy of melting q_i (J.m-3) from temperature
      !!
      !! ** Method  :   Formula (Bitz and Lipscomb, 1999)
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in)  ::  jl_cat
      INTEGER              ::  ji, jj, jk          ! dummy loop indices
      REAL(wp)             ::  ztmelts             ! local scalar
      !!-------------------------------------------------------------------
      !
      DO jk = 1, nlay_i             ! Sea ice energy of melting
         DO_2D( 0, 0, 0, 0 )
            IF( l_ice_present(ji,jj) ) THEN
               ztmelts       = - rTmlt  * sz_i(ji,jj,jk,jl_cat)
               t_i(ji,jj,jk,jl_cat) = MIN( t_i(ji,jj,jk,jl_cat), ztmelts + rt0 ) ! Force t_i to be lower than melting point => likely conservation issue
               !                                                                 !   (sometimes zdf scheme produces abnormally high temperatures)
               e_i(ji,jj,jk,jl_cat) = rhoi * ( rcpi  * ( ztmelts - ( t_i(ji,jj,jk,jl_cat) - rt0 ) )                       &
                  &                  + rLfus * MAX( 0._wp, 1._wp - ztmelts / MIN( t_i(ji,jj,jk,jl_cat) - rt0, -epsi10 ) ) & ! clem: max to deal with different eq. freezing in ice and ocean (but useless for now)
                  &                  - rcp   * ztmelts )
          ENDIF
         END_2D
      END DO
      !
      DO jk = 1, nlay_s             ! Snow energy of melting
         DO_2D( 0, 0, 0, 0 )
            IF( l_ice_present(ji,jj) ) THEN
               e_s(ji,jj,jk,jl_cat) = rhos * ( rcpi * ( rt0 - t_s(ji,jj,jk,jl_cat) ) + rLfus )
            ENDIF
         END_2D
      END DO
      !
   END SUBROUTINE ice_var_enthalpy

   SUBROUTINE ice_var_vremap( ph_old, pts_old, pts_i )
      !!-------------------------------------------------------------------
      !!               ***   ROUTINE ice_var_vremap  ***
      !!
      !! ** Purpose :
      !!           This routine computes new vertical grids in the ice, 
      !!           and consistently redistributes temperatures and salinities 
      !!           Redistribution is made so as to ensure energy/salt conservation
      !!
      !!
      !! ** Method  : linear conservative remapping
      !!           
      !! ** Steps : 1) cumulative integrals of old enthalpies/salinities/thicknesses
      !!            2) linear remapping on the new layers
      !!
      !! ------------ cum0(0)                        ------------- cum1(0)
      !!                                    NEW      -------------
      !! ------------ cum0(1)               ==>      -------------
      !!     ...                                     -------------
      !! ------------                                -------------
      !! ------------ cum0(nlay_i+2)                 ------------- cum1(nlay_i)
      !!
      !!
      !! References : Bitz & Lipscomb, JGR 99; Vancoppenolle et al., GRL, 2005
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(0:nlay_i+1), INTENT(in)    ::   ph_old, pts_old  ! old tickness (m), enthlapy (J.m-2) or salt (m.g/kg)
      REAL(wp), DIMENSION(1:nlay_i)  , INTENT(inout) ::   pts_i            ! new enthlapies (J.m-3, remapped) or salt (g/kg)
      !
      INTEGER  ::   jk0, jk1   !  old/new layer indices
      !
      REAL(wp), DIMENSION(0:nlay_i+2) ::   zts_cum0, zh_cum0   ! old cumulative enthlapies/salinities and layers interfaces
      REAL(wp), DIMENSION(0:nlay_i)   ::   zts_cum1, zh_cum1   ! new cumulative enthlapies/salinities and layers interfaces
      REAL(wp)                        ::   zhnew               ! new layers thicknesses
      !!-------------------------------------------------------------------

      !-------------------------------------------------------------------------------
      !  1) Cumulative integral of old enthalpy/salt * thickness and layers interfaces
      !-------------------------------------------------------------------------------
      zts_cum0(0) = 0._wp 
      zh_cum0 (0) = 0._wp
      DO jk0 = 1, nlay_i+2
         zts_cum0(jk0) = zts_cum0(jk0-1) + pts_old(jk0-1)
         zh_cum0 (jk0) = zh_cum0 (jk0-1) + ph_old(jk0-1)
      END DO

      !------------------------------------
      !  2) Interpolation on the new layers
      !------------------------------------
      ! new layer thickesses
      zhnew = SUM( ph_old(0:nlay_i+1) ) * r1_nlay_i  

      ! new layers interfaces
      zh_cum1(0) = 0._wp
      DO jk1 = 1, nlay_i
         zh_cum1(jk1) = zh_cum1(jk1-1) + zhnew
      END DO

      zts_cum1(0:nlay_i) = 0._wp 
      ! new cumulative q*h => linear interpolation
      DO jk0 = 1, nlay_i+2
         DO jk1 = 1, nlay_i-1
            IF( zh_cum1(jk1) <= zh_cum0(jk0) .AND. zh_cum1(jk1) > zh_cum0(jk0-1) )   THEN
               zts_cum1(jk1) = ( zts_cum0(jk0-1) * ( zh_cum0(jk0) - zh_cum1(jk1  ) ) +  &
                  &              zts_cum0(jk0  ) * ( zh_cum1(jk1) - zh_cum0(jk0-1) ) )  &
                  &            / ( zh_cum0(jk0) - zh_cum0(jk0-1) )
            ENDIF
         END DO
      END DO
      ! to ensure that total heat/salt content is strictly conserved, set:
      zts_cum1(nlay_i) = zts_cum0(nlay_i+2) 

      ! new enthalpies/salinities
      DO jk1 = 1, nlay_i
         pts_i(jk1) = MAX( 0._wp, zts_cum1(jk1) - zts_cum1(jk1-1) ) / MAX( zhnew, epsi10 ) ! max for roundoff error
      END DO
      
   END SUBROUTINE ice_var_vremap

   SUBROUTINE snw_var_vremap( ph_old, pe_old, pe_s )
      !!-------------------------------------------------------------------
      !!               ***   ROUTINE snw_var_vremap  ***
      !! ** Purpose :
      !!           This routine computes new vertical grids in the snow,
      !!           and consistently redistributes temperatures.
      !!           Redistribution is made so as to ensure to energy conservation
      !!
      !!
      !! ** Method  : linear conservative remapping
      !!
      !! ** Steps : 1) cumulative integrals of old enthalpies/thicknesses
      !!            2) linear remapping on the new layers
      !!
      !! ------------ cum0(0)                        ------------- cum1(0)
      !!                                    NEW      -------------
      !! ------------ cum0(1)               ==>      -------------
      !!     ...                                     -------------
      !! ------------                                -------------
      !! ------------ cum0(nlay_s+1)                 ------------- cum1(nlay_s)
      !!
      !!
      !! References : Bitz & Lipscomb, JGR 99; Vancoppenolle et al., GRL, 2005
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(0:nlay_s), INTENT(in)    ::   ph_old      ! old thicknesses (m)
      REAL(wp), DIMENSION(0:nlay_s), INTENT(in)    ::   pe_old      ! old enthlapies  (J.m-3)
      REAL(wp), DIMENSION(1:nlay_s), INTENT(inout) ::   pe_s        ! new enthlapies  (J.m-3, remapped)
      !
      INTEGER  :: jk0, jk1   !  old/new layer indices
      !
      REAL(wp), DIMENSION(0:nlay_s+1) ::   zeh_cum0, zh_cum0   ! old cumulative enthlapies and layers interfaces
      REAL(wp), DIMENSION(0:nlay_s)   ::   zeh_cum1, zh_cum1   ! new cumulative enthlapies and layers interfaces
      REAL(wp)                        ::   zhnew               ! new layers thicknesses
      !!-------------------------------------------------------------------

      !--------------------------------------------------------------------------
      !  1) Cumulative integral of old enthalpy * thickness and layers interfaces
      !--------------------------------------------------------------------------
      zeh_cum0(0) = 0._wp
      zh_cum0 (0) = 0._wp
      DO jk0 = 1, nlay_s+1
         zeh_cum0(jk0) = zeh_cum0(jk0-1) + pe_old(jk0-1) * ph_old(jk0-1)
         zh_cum0 (jk0) = zh_cum0 (jk0-1) + ph_old(jk0-1)
      END DO

      !------------------------------------
      !  2) Interpolation on the new layers
      !------------------------------------
      ! new layer thickesses
      zhnew = SUM( ph_old(0:nlay_s) ) * r1_nlay_s

      ! new layers interfaces
      zh_cum1(0) = 0._wp
      DO jk1 = 1, nlay_s
         zh_cum1(jk1) = zh_cum1(jk1-1) + zhnew
      END DO

      zeh_cum1(0:nlay_s) = 0._wp
      ! new cumulative q*h => linear interpolation
      DO jk0 = 1, nlay_s+1
         DO jk1 = 1, nlay_s-1
            IF( zh_cum1(jk1) <= zh_cum0(jk0) .AND. zh_cum1(jk1) > zh_cum0(jk0-1) )   THEN
               zeh_cum1(jk1) = ( zeh_cum0(jk0-1) * ( zh_cum0(jk0) - zh_cum1(jk1  ) ) +  &
                  &              zeh_cum0(jk0  ) * ( zh_cum1(jk1) - zh_cum0(jk0-1) ) )  &
                  &            / ( zh_cum0(jk0) - zh_cum0(jk0-1) )
            ENDIF
         END DO
      END DO
      ! to ensure that total heat content is strictly conserved, set:
      zeh_cum1(nlay_s) = zeh_cum0(nlay_s+1)

      ! new enthalpies
      DO jk1 = 1, nlay_s
         pe_s(jk1) = MAX( 0._wp, zeh_cum1(jk1) - zeh_cum1(jk1-1) ) / MAX( zhnew, epsi20 ) ! max for roundoff error
      END DO

   END SUBROUTINE snw_var_vremap

   FUNCTION ice_var_sshdyn(pssh, psnwice_mass, psnwice_mass_b)
      !!---------------------------------------------------------------------
      !!                   ***  ROUTINE ice_var_sshdyn  ***
      !!
      !! ** Purpose :  compute the equivalent ssh in lead when sea ice is embedded
      !!
      !! ** Method  :  ssh_lead = ssh + (Mice + Msnow) / rho0
      !!
      !! ** Reference : Jean-Michel Campin, John Marshall, David Ferreira,
      !!                Sea ice-ocean coupling using a rescaled vertical coordinate z*,
      !!                Ocean Modelling, Volume 24, Issues 1-2, 2008
      !!----------------------------------------------------------------------
      !
      ! input
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: pssh            !: ssh [m]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: psnwice_mass    !: mass of snow and ice at current  ice time step [Kg/m2]
      REAL(wp), DIMENSION(jpi,jpj), INTENT(in) :: psnwice_mass_b  !: mass of snow and ice at previous ice time step [Kg/m2]
      !
      ! output
      REAL(wp), DIMENSION(jpi,jpj) :: ice_var_sshdyn  ! equivalent ssh in lead [m]
      !
      ! temporary
      REAL(wp) :: zintn, zintb                     ! time interpolation weights []
      !
      ! compute ice load used to define the equivalent ssh in lead
      IF( ln_ice_embd ) THEN
         !
         ! average interpolation coeff as used in dynspg = (1/nn_fsbc)   * {SUM[n/nn_fsbc], n=0,nn_fsbc-1}
         !                                               = (1/nn_fsbc)^2 * {SUM[n]        , n=0,nn_fsbc-1}
         zintn = REAL( nn_fsbc - 1 ) / REAL( nn_fsbc ) * 0.5_wp
         !
         ! average interpolation coeff as used in dynspg = (1/nn_fsbc)   *    {SUM[1-n/nn_fsbc], n=0,nn_fsbc-1}
         !                                               = (1/nn_fsbc)^2 * (nn_fsbc^2 - {SUM[n], n=0,nn_fsbc-1})
         zintb = REAL( nn_fsbc + 1 ) / REAL( nn_fsbc ) * 0.5_wp
         !
         ! compute equivalent ssh in lead
         ice_var_sshdyn(:,:) = pssh(:,:) + ( zintn * psnwice_mass(:,:) + zintb * psnwice_mass_b(:,:) ) * r1_rho0
         !
      ELSE
         ! compute equivalent ssh in lead
         ice_var_sshdyn(:,:) = pssh(:,:)
      ENDIF
      !
   END FUNCTION ice_var_sshdyn


   !!-------------------------------------------------------------------
   !!                ***  INTERFACE ice_var_itd   ***
   !!
   !! ** Purpose :  converting N-cat ice to jpl ice categories
   !!-------------------------------------------------------------------
   SUBROUTINE ice_var_itd_1c1c( phti, phts, pati ,                             ph_i, ph_s, pa_i, &
      &                         ptmi, ptms, ptmsu, psmi, patip, phtip, phtil,  pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il )
      !!-------------------------------------------------------------------
      !! ** Purpose :  converting 1-cat ice to 1 ice category
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(:), INTENT(in)    ::   phti, phts, pati    ! input  ice/snow variables
      REAL(wp), DIMENSION(:), INTENT(inout) ::   ph_i, ph_s, pa_i    ! output ice/snow variables
      REAL(wp), DIMENSION(:), INTENT(in)    ::   ptmi, ptms, ptmsu, psmi, patip, phtip, phtil    ! input  ice/snow temp & sal & ponds
      REAL(wp), DIMENSION(:), INTENT(inout) ::   pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il    ! output ice/snow temp & sal & ponds
      !!-------------------------------------------------------------------
      ! == thickness and concentration == !
      ph_i(:) = phti(:)
      ph_s(:) = phts(:)
      pa_i(:) = pati(:)
      !
      ! == temperature and salinity and ponds == !
      pt_i (:) = ptmi (:)
      pt_s (:) = ptms (:)
      pt_su(:) = ptmsu(:)
      ps_i (:) = psmi (:)
      pa_ip(:) = patip(:)
      ph_ip(:) = phtip(:)
      ph_il(:) = phtil(:)

   END SUBROUTINE ice_var_itd_1c1c

   SUBROUTINE ice_var_itd_Nc1c( phti, phts, pati ,                             ph_i, ph_s, pa_i, &
      &                         ptmi, ptms, ptmsu, psmi, patip, phtip, phtil,  pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il )
      !!-------------------------------------------------------------------
      !! ** Purpose :  converting N-cat ice to 1 ice category
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)    ::   phti, phts, pati    ! input  ice/snow variables
      REAL(wp), DIMENSION(:)  , INTENT(inout) ::   ph_i, ph_s, pa_i    ! output ice/snow variables
      REAL(wp), DIMENSION(:,:), INTENT(in)    ::   ptmi, ptms, ptmsu, psmi, patip, phtip, phtil    ! input  ice/snow temp & sal & ponds
      REAL(wp), DIMENSION(:)  , INTENT(inout) ::   pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il    ! output ice/snow temp & sal & ponds
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:) ::   z1_ai, z1_vi, z1_vs
      !
      INTEGER ::   idim
      !!-------------------------------------------------------------------
      !
      idim = SIZE( phti, 1 )
      !
      ! == thickness and concentration == !
      ALLOCATE( z1_ai(idim), z1_vi(idim), z1_vs(idim) )
      !
      pa_i(:) = SUM( pati(:,:), dim=2 )

      WHERE( ( pa_i(:) ) /= 0._wp )   ;   z1_ai(:) = 1._wp / pa_i(:)
      ELSEWHERE                       ;   z1_ai(:) = 0._wp
      END WHERE

      ph_i(:) = SUM( phti(:,:) * pati(:,:), dim=2 ) * z1_ai(:)
      ph_s(:) = SUM( phts(:,:) * pati(:,:), dim=2 ) * z1_ai(:)
      !
      ! == temperature and salinity == !
      WHERE( ( pa_i(:) * ph_i(:) ) /= 0._wp )   ;   z1_vi(:) = 1._wp / ( pa_i(:) * ph_i(:) )
      ELSEWHERE                                 ;   z1_vi(:) = 0._wp
      END WHERE
      WHERE( ( pa_i(:) * ph_s(:) ) /= 0._wp )   ;   z1_vs(:) = 1._wp / ( pa_i(:) * ph_s(:) )
      ELSEWHERE                                 ;   z1_vs(:) = 0._wp
      END WHERE
      pt_i (:) = SUM( ptmi (:,:) * pati(:,:) * phti(:,:), dim=2 ) * z1_vi(:)
      pt_s (:) = SUM( ptms (:,:) * pati(:,:) * phts(:,:), dim=2 ) * z1_vs(:)
      pt_su(:) = SUM( ptmsu(:,:) * pati(:,:)            , dim=2 ) * z1_ai(:)
      ps_i (:) = SUM( psmi (:,:) * pati(:,:) * phti(:,:), dim=2 ) * z1_vi(:)

      ! == ponds == !
      pa_ip(:) = SUM( patip(:,:), dim=2 )
      WHERE( pa_ip(:) /= 0._wp )
         ph_ip(:) = SUM( phtip(:,:) * patip(:,:), dim=2 ) / pa_ip(:)
         ph_il(:) = SUM( phtil(:,:) * patip(:,:), dim=2 ) / pa_ip(:)
      ELSEWHERE
         ph_ip(:) = 0._wp
         ph_il(:) = 0._wp
      END WHERE
      !
      DEALLOCATE( z1_ai, z1_vi, z1_vs )
      !
   END SUBROUTINE ice_var_itd_Nc1c

   SUBROUTINE ice_var_itd_1cMc( phti, phts, pati ,                             ph_i, ph_s, pa_i, &
      &                         ptmi, ptms, ptmsu, psmi, patip, phtip, phtil,  pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il )
      !!-------------------------------------------------------------------
      !!
      !! ** Purpose :  converting 1-cat ice to jpl ice categories
      !!
      !!
      !! ** Method:   ice thickness distribution follows a gamma function from Abraham et al. (2015)
      !!              it has the property of conserving total concentration and volume
      !!
      !!
      !! ** Arguments : phti: 1-cat ice thickness
      !!                phts: 1-cat snow depth
      !!                pati: 1-cat ice concentration
      !!
      !! ** Output    : jpl-cat
      !!
      !!  Abraham, C., Steiner, N., Monahan, A. and Michel, C., 2015.
      !!               Effects of subgrid‐scale snow thickness variability on radiative transfer in sea ice.
      !!               Journal of Geophysical Research: Oceans, 120(8), pp.5597-5614
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(:),   INTENT(in)    ::   phti, phts, pati    ! input  ice/snow variables
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   ph_i, ph_s, pa_i    ! output ice/snow variables
      REAL(wp), DIMENSION(:)  , INTENT(in)    ::   ptmi, ptms, ptmsu, psmi, patip, phtip, phtil    ! input  ice/snow temp & sal & ponds
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il    ! output ice/snow temp & sal & ponds
      !
      REAL(wp), ALLOCATABLE, DIMENSION(:) ::   zfra, z1_hti
      INTEGER  ::   ji, jk, jl
      INTEGER  ::   idim
      REAL(wp) ::   zv, zdh
      !!-------------------------------------------------------------------
      !
      idim = SIZE( phti , 1 )
      !
      ph_i(1:idim,1:jpl) = 0._wp
      ph_s(1:idim,1:jpl) = 0._wp
      pa_i(1:idim,1:jpl) = 0._wp
      !
      ALLOCATE( z1_hti(idim) )
      WHERE( phti(:) /= 0._wp )   ;   z1_hti(:) = 1._wp / phti(:)
      ELSEWHERE                   ;   z1_hti(:) = 0._wp
      END WHERE
      !
      ! == thickness and concentration == !
      ! for categories 1:jpl-1, integrate the gamma function from hi_max(jl-1) to hi_max(jl)
      DO jl = 1, jpl-1
         DO ji = 1, idim
            !
            IF( phti(ji) > 0._wp ) THEN
               ! concentration : integrate ((4A/H^2)xexp(-2x/H))dx from x=hi_max(jl-1) to hi_max(jl)
               pa_i(ji,jl) = pati(ji) * z1_hti(ji) * (  ( phti(ji) + 2.*hi_max(jl-1) ) * EXP( -2.*hi_max(jl-1)*z1_hti(ji) ) &
                  &                                   - ( phti(ji) + 2.*hi_max(jl  ) ) * EXP( -2.*hi_max(jl  )*z1_hti(ji) ) )
               !
               ! volume : integrate ((4A/H^2)x^2exp(-2x/H))dx from x=hi_max(jl-1) to hi_max(jl)
               zv = pati(ji) * z1_hti(ji) * (  ( phti(ji)*phti(ji) + 2.*phti(ji)*hi_max(jl-1) + 2.*hi_max(jl-1)*hi_max(jl-1) ) &
                  &                            * EXP( -2.*hi_max(jl-1)*z1_hti(ji) ) &
                  &                          - ( phti(ji)*phti(ji) + 2.*phti(ji)*hi_max(jl) + 2.*hi_max(jl)*hi_max(jl) ) &
                  &                            * EXP(-2.*hi_max(jl)*z1_hti(ji)) )
               ! thickness
               IF( pa_i(ji,jl) > epsi06 ) THEN
                  ph_i(ji,jl) = zv / pa_i(ji,jl)
               ELSE
                  ph_i(ji,jl) = 0.
                  pa_i(ji,jl) = 0.
               ENDIF
            ENDIF
            !
         ENDDO
      ENDDO
      !
      ! for the last category (jpl), integrate the gamma function from hi_max(jpl-1) to infinity
      DO ji = 1, idim
         !
         IF( phti(ji) > 0._wp ) THEN
            ! concentration : integrate ((4A/H^2)xexp(-2x/H))dx from x=hi_max(jpl-1) to infinity
            pa_i(ji,jpl) = pati(ji) * z1_hti(ji) * ( phti(ji) + 2.*hi_max(jpl-1) ) * EXP( -2.*hi_max(jpl-1)*z1_hti(ji) )

            ! volume : integrate ((4A/H^2)x^2exp(-2x/H))dx from x=hi_max(jpl-1) to infinity
            zv = pati(ji) * z1_hti(ji) * ( phti(ji)*phti(ji) + 2.*phti(ji)*hi_max(jpl-1) + 2.*hi_max(jpl-1)*hi_max(jpl-1) ) &
               &                         * EXP( -2.*hi_max(jpl-1)*z1_hti(ji) )
            ! thickness
            IF( pa_i(ji,jpl) > epsi06 ) THEN
               ph_i(ji,jpl) = zv / pa_i(ji,jpl)
            else
               ph_i(ji,jpl) = 0.
               pa_i(ji,jpl) = 0.
            ENDIF
         ENDIF
         !
      ENDDO
      !
      ! Add Snow in each category where pa_i is not 0
      DO jl = 1, jpl
         DO ji = 1, idim
            IF( pa_i(ji,jl) > 0._wp ) THEN
               ph_s(ji,jl) = ph_i(ji,jl) * phts(ji) * z1_hti(ji)
               ! In case snow load is in excess that would lead to transformation from snow to ice
               ! Then, transfer the snow excess into the ice (different from icethd_dh)
               zdh = MAX( 0._wp, ( rhos * ph_s(ji,jl) + ( rhoi - rho0 ) * ph_i(ji,jl) ) * r1_rho0 )
               ! recompute h_i, h_s avoiding out of bounds values
               ph_i(ji,jl) = MIN( hi_max(jl), ph_i(ji,jl) + zdh )
               ph_s(ji,jl) = MAX( 0._wp, ph_s(ji,jl) - zdh * rhoi * r1_rhos )
            ENDIF
         END DO
      END DO
      !
      DEALLOCATE( z1_hti )
      !
      ! == temperature and salinity == !
      DO jl = 1, jpl
         pt_i (:,jl) = ptmi (:)
         pt_s (:,jl) = ptms (:)
         pt_su(:,jl) = ptmsu(:)
         ps_i (:,jl) = psmi (:)
      END DO
      !
      ! == ponds == !
      ALLOCATE( zfra(idim) )
      ! keep the same pond fraction atip/ati for each category
      WHERE( pati(:) /= 0._wp )   ;   zfra(:) = patip(:) / pati(:)
      ELSEWHERE                   ;   zfra(:) = 0._wp
      END WHERE
      DO jl = 1, jpl
         pa_ip(:,jl) = zfra(:) * pa_i(:,jl)
      END DO
      ! keep the same v_ip/v_i ratio for each category
      WHERE( ( phti(:) * pati(:) ) /= 0._wp )   ;   zfra(:) = ( phtip(:) * patip(:) ) / ( phti(:) * pati(:) )
      ELSEWHERE                                 ;   zfra(:) = 0._wp
      END WHERE
      DO jl = 1, jpl
         WHERE( pa_ip(:,jl) /= 0._wp )   ;   ph_ip(:,jl) = zfra(:) * ( ph_i(:,jl) * pa_i(:,jl) ) / pa_ip(:,jl)
         ELSEWHERE                       ;   ph_ip(:,jl) = 0._wp
         END WHERE
      END DO
      ! keep the same v_il/v_i ratio for each category
      WHERE( ( phti(:) * pati(:) ) /= 0._wp )   ;   zfra(:) = ( phtil(:) * patip(:) ) / ( phti(:) * pati(:) )
      ELSEWHERE                                 ;   zfra(:) = 0._wp
      END WHERE
      DO jl = 1, jpl
         WHERE( pa_ip(:,jl) /= 0._wp )   ;   ph_il(:,jl) = zfra(:) * ( ph_i(:,jl) * pa_i(:,jl) ) / pa_ip(:,jl)
         ELSEWHERE                       ;   ph_il(:,jl) = 0._wp
         END WHERE
      END DO
      DEALLOCATE( zfra )
      !
   END SUBROUTINE ice_var_itd_1cMc

   SUBROUTINE ice_var_itd_NcMc( phti, phts, pati ,                             ph_i, ph_s, pa_i, &
      &                         ptmi, ptms, ptmsu, psmi, patip, phtip, phtil,  pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il )
      !!-------------------------------------------------------------------
      !!
      !! ** Purpose :  converting N-cat ice to jpl ice categories
      !!
      !!                  ice thickness distribution follows a gaussian law
      !!               around the concentration of the most likely ice thickness
      !!                           (similar as iceistate.F90)
      !!
      !! ** Method:   Iterative procedure
      !!
      !!               1) Fill ice cat that correspond to input thicknesses
      !!                  Find the lowest(jlmin) and highest(jlmax) cat that are filled
      !!
      !!               2) Expand the filling to the cat jlmin-1 and jlmax+1
      !!                   by removing 25% ice area from jlmin and jlmax (resp.)
      !!
      !!               3) Expand the filling to the empty cat between jlmin and jlmax
      !!                   by a) removing 25% ice area from the lower cat (ascendant loop jlmin=>jlmax)
      !!                      b) removing 25% ice area from the higher cat (descendant loop jlmax=>jlmin)
      !!
      !! ** Arguments : phti: N-cat ice thickness
      !!                phts: N-cat snow depth
      !!                pati: N-cat ice concentration
      !!
      !! ** Output    : jpl-cat
      !!
      !!  (Example of application: BDY forcings when inputs have N-cat /= jpl)
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(in)    ::   phti, phts, pati    ! input  ice/snow variables
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   ph_i, ph_s, pa_i    ! output ice/snow variables
      REAL(wp), DIMENSION(:,:), INTENT(in)    ::   ptmi, ptms, ptmsu, psmi, patip, phtip, phtil    ! input  ice/snow temp & sal & ponds
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pt_i, pt_s, pt_su, ps_i, pa_ip, ph_ip, ph_il    ! output ice/snow temp & sal & ponds
      !
      INTEGER , ALLOCATABLE, DIMENSION(:,:) ::   jlfil, jlfil2
      INTEGER , ALLOCATABLE, DIMENSION(:)   ::   jlmax, jlmin
      REAL(wp), ALLOCATABLE, DIMENSION(:)   ::   z1_ai, z1_vi, z1_vs, ztmp, zfra
      !
      REAL(wp), PARAMETER ::   ztrans = 0.25_wp
      INTEGER  ::   ji, jl, jl1, jl2
      INTEGER  ::   idim, icat
      !!-------------------------------------------------------------------
      !
      idim = SIZE( phti, 1 )
      icat = SIZE( phti, 2 )
      !
      ! == thickness and concentration == !
      !                                 ! ---------------------- !
      IF( icat == jpl ) THEN            ! input cat = output cat !
         !                              ! ---------------------- !
         ph_i(:,:) = phti(:,:)
         ph_s(:,:) = phts(:,:)
         pa_i(:,:) = pati(:,:)
         !
         ! == temperature and salinity and ponds == !
         pt_i (:,:) = ptmi (:,:)
         pt_s (:,:) = ptms (:,:)
         pt_su(:,:) = ptmsu(:,:)
         ps_i (:,:) = psmi (:,:)
         pa_ip(:,:) = patip(:,:)
         ph_ip(:,:) = phtip(:,:)
         ph_il(:,:) = phtil(:,:)
         !                              ! ---------------------- !
      ELSEIF( icat == 1 ) THEN          ! input cat = 1          !
         !                              ! ---------------------- !
         CALL  ice_var_itd_1cMc( phti(:,1), phts(:,1), pati (:,1), &
            &                    ph_i(:,:), ph_s(:,:), pa_i (:,:), &
            &                    ptmi(:,1), ptms(:,1), ptmsu(:,1), psmi(:,1), patip(:,1), phtip(:,1), phtil(:,1), &
            &                    pt_i(:,:), pt_s(:,:), pt_su(:,:), ps_i(:,:), pa_ip(:,:), ph_ip(:,:), ph_il(:,:)  )
         !                              ! ---------------------- !
      ELSEIF( jpl == 1 ) THEN           ! output cat = 1         !
         !                              ! ---------------------- !
         CALL  ice_var_itd_Nc1c( phti(:,:), phts(:,:), pati (:,:), &
            &                    ph_i(:,1), ph_s(:,1), pa_i (:,1), &
            &                    ptmi(:,:), ptms(:,:), ptmsu(:,:), psmi(:,:), patip(:,:), phtip(:,:), phtil(:,:), &
            &                    pt_i(:,1), pt_s(:,1), pt_su(:,1), ps_i(:,1), pa_ip(:,1), ph_ip(:,1), ph_il(:,1)  )
         !                              ! ----------------------- !
      ELSE                              ! input cat /= output cat !
         !                              ! ----------------------- !

         ALLOCATE( jlfil(idim,jpl), jlfil2(idim,jpl) )       ! allocate arrays
         ALLOCATE( jlmin(idim), jlmax(idim) )

         ! --- initialize output fields to 0 --- !
         ph_i(1:idim,1:jpl) = 0._wp
         ph_s(1:idim,1:jpl) = 0._wp
         pa_i(1:idim,1:jpl) = 0._wp
         !
         ! --- fill the categories --- !
         !     find where cat-input = cat-output and fill cat-output fields
         jlmax(:) = 0
         jlmin(:) = 999
         jlfil(:,:) = 0
         DO jl1 = 1, jpl
            DO jl2 = 1, icat
               DO ji = 1, idim
                  IF( hi_max(jl1-1) <= phti(ji,jl2) .AND. hi_max(jl1) > phti(ji,jl2) ) THEN
                     ! fill the right category
                     ph_i(ji,jl1) = phti(ji,jl2)
                     ph_s(ji,jl1) = phts(ji,jl2)
                     pa_i(ji,jl1) = pati(ji,jl2)
                     ! record categories that are filled
                     jlmax(ji) = MAX( jlmax(ji), jl1 )
                     jlmin(ji) = MIN( jlmin(ji), jl1 )
                     jlfil(ji,jl1) = jl1
                  ENDIF
               END DO
            END DO
         END DO
         !
         ! --- fill the gaps between categories --- !
         !     transfer from categories filled at the previous step to the empty ones in between
         DO ji = 1, idim
            jl1 = jlmin(ji)
            jl2 = jlmax(ji)
            IF( jl1 > 1 ) THEN
               ! fill the lower cat (jl1-1)
               pa_i(ji,jl1-1) = ztrans * pa_i(ji,jl1)
               ph_i(ji,jl1-1) = hi_mean(jl1-1)
               ! remove from cat jl1
               pa_i(ji,jl1  ) = ( 1._wp - ztrans ) * pa_i(ji,jl1)
            ENDIF
            IF( jl2 < jpl ) THEN
               ! fill the upper cat (jl2+1)
               pa_i(ji,jl2+1) = ztrans * pa_i(ji,jl2)
               ph_i(ji,jl2+1) = hi_mean(jl2+1)
               ! remove from cat jl2
               pa_i(ji,jl2  ) = ( 1._wp - ztrans ) * pa_i(ji,jl2)
            ENDIF
         END DO
         !
         jlfil2(:,:) = jlfil(:,:)
         ! fill categories from low to high
         DO jl = 2, jpl-1
            DO ji = 1, idim
               IF( jlfil(ji,jl-1) /= 0 .AND. jlfil(ji,jl) == 0 ) THEN
                  ! fill high
                  pa_i(ji,jl) = ztrans * pa_i(ji,jl-1)
                  ph_i(ji,jl) = hi_mean(jl)
                  jlfil(ji,jl) = jl
                  ! remove low
                  pa_i(ji,jl-1) = ( 1._wp - ztrans ) * pa_i(ji,jl-1)
               ENDIF
            END DO
         END DO
         !
         ! fill categories from high to low
         DO jl = jpl-1, 2, -1
            DO ji = 1, idim
               IF( jlfil2(ji,jl+1) /= 0 .AND. jlfil2(ji,jl) == 0 ) THEN
                  ! fill low
                  pa_i(ji,jl) = pa_i(ji,jl) + ztrans * pa_i(ji,jl+1)
                  ph_i(ji,jl) = hi_mean(jl)
                  jlfil2(ji,jl) = jl
                  ! remove high
                  pa_i(ji,jl+1) = ( 1._wp - ztrans ) * pa_i(ji,jl+1)
               ENDIF
            END DO
         END DO
         !
         DEALLOCATE( jlfil, jlfil2 )      ! deallocate arrays
         DEALLOCATE( jlmin, jlmax )
         !
         ! == temperature and salinity == !
         !
         ALLOCATE( z1_ai(idim), z1_vi(idim), z1_vs(idim), ztmp(idim) )
         !
         WHERE( SUM( pa_i(:,:), dim=2 ) /= 0._wp )               ;   z1_ai(:) = 1._wp / SUM( pa_i(:,:), dim=2 )
         ELSEWHERE                                               ;   z1_ai(:) = 0._wp
         END WHERE
         WHERE( SUM( pa_i(:,:) * ph_i(:,:), dim=2 ) /= 0._wp )   ;   z1_vi(:) = 1._wp / SUM( pa_i(:,:) * ph_i(:,:), dim=2 )
         ELSEWHERE                                               ;   z1_vi(:) = 0._wp
         END WHERE
         WHERE( SUM( pa_i(:,:) * ph_s(:,:), dim=2 ) /= 0._wp )   ;   z1_vs(:) = 1._wp / SUM( pa_i(:,:) * ph_s(:,:), dim=2 )
         ELSEWHERE                                               ;   z1_vs(:) = 0._wp
         END WHERE
         !
         ! fill all the categories with the same value
         ztmp(:) = SUM( ptmi (:,:) * pati(:,:) * phti(:,:), dim=2 ) * z1_vi(:)
         DO jl = 1, jpl
            pt_i (:,jl) = ztmp(:)
         END DO
         ztmp(:) = SUM( ptms (:,:) * pati(:,:) * phts(:,:), dim=2 ) * z1_vs(:)
         DO jl = 1, jpl
            pt_s (:,jl) = ztmp(:)
         END DO
         ztmp(:) = SUM( ptmsu(:,:) * pati(:,:)            , dim=2 ) * z1_ai(:)
         DO jl = 1, jpl
            pt_su(:,jl) = ztmp(:)
         END DO
         ztmp(:) = SUM( psmi (:,:) * pati(:,:) * phti(:,:), dim=2 ) * z1_vi(:)
         DO jl = 1, jpl
            ps_i (:,jl) = ztmp(:)
         END DO
         !
         DEALLOCATE( z1_ai, z1_vi, z1_vs, ztmp )
         !
         ! == ponds == !
         ALLOCATE( zfra(idim) )
         ! keep the same pond fraction atip/ati for each category
         WHERE( SUM( pati(:,:), dim=2 ) /= 0._wp )   ;   zfra(:) = SUM( patip(:,:), dim=2 ) / SUM( pati(:,:), dim=2 )
         ELSEWHERE                                   ;   zfra(:) = 0._wp
         END WHERE
         DO jl = 1, jpl
            pa_ip(:,jl) = zfra(:) * pa_i(:,jl)
         END DO
         ! keep the same v_ip/v_i ratio for each category
         WHERE( SUM( phti(:,:) * pati(:,:), dim=2 ) /= 0._wp )
            zfra(:) = SUM( phtip(:,:) * patip(:,:), dim=2 ) / SUM( phti(:,:) * pati(:,:), dim=2 )
         ELSEWHERE
            zfra(:) = 0._wp
         END WHERE
         DO jl = 1, jpl
            WHERE( pa_ip(:,jl) /= 0._wp )   ;   ph_ip(:,jl) = zfra(:) * ( ph_i(:,jl) * pa_i(:,jl) ) / pa_ip(:,jl)
            ELSEWHERE                       ;   ph_ip(:,jl) = 0._wp
            END WHERE
         END DO
         ! keep the same v_il/v_i ratio for each category
         WHERE( SUM( phti(:,:) * pati(:,:), dim=2 ) /= 0._wp )
            zfra(:) = SUM( phtil(:,:) * patip(:,:), dim=2 ) / SUM( phti(:,:) * pati(:,:), dim=2 )
         ELSEWHERE
            zfra(:) = 0._wp
         END WHERE
         DO jl = 1, jpl
            WHERE( pa_ip(:,jl) /= 0._wp )   ;   ph_il(:,jl) = zfra(:) * ( ph_i(:,jl) * pa_i(:,jl) ) / pa_ip(:,jl)
            ELSEWHERE                       ;   ph_il(:,jl) = 0._wp
            END WHERE
         END DO
         DEALLOCATE( zfra )
         !
      ENDIF
      !
   END SUBROUTINE ice_var_itd_NcMc

   !!-------------------------------------------------------------------
   !! INTERFACE ice_var_snwfra
   !!
   !! ** Purpose :  fraction of ice covered by snow
   !!
   !! ** Method  :  In absence of proper snow model on top of sea ice,
   !!               we argue that snow does not cover the whole ice because
   !!               of wind blowing...
   !!
   !! ** Arguments : ph_s: snow thickness
   !!
   !! ** Output    : pa_s_fra: fraction of ice covered by snow
   !!
   !!-------------------------------------------------------------------
   SUBROUTINE ice_var_snwfra_3d( ph_s, pa_s_fra )
      REAL(wp), DIMENSION(A2D(nn_hls),jpl), INTENT(in   ) ::   ph_s        ! snow thickness
      REAL(wp), DIMENSION(A2D(nn_hls),jpl), INTENT(  out) ::   pa_s_fra    ! ice fraction covered by snow
      INTEGER  ::   ji, jj, jl
      IF    ( nn_snwfra == 0 ) THEN   ! basic 0 or 1 snow cover
         DO jl = 1, jpl
            DO_2D( 0, 0, 0, 0 )
               pa_s_fra(ji,jj,jl) = MERGE( 1._wp, 0._wp, ph_s(ji,jj,jl) > 0._wp )
            END_2D
         END DO
      ELSEIF( nn_snwfra == 1 ) THEN   ! snow cover depends on hsnow (met-office style)
         DO jl = 1, jpl
            DO_2D( 0, 0, 0, 0 )
               pa_s_fra(ji,jj,jl) = 1._wp - EXP( -0.2_wp * rhos * ph_s(ji,jj,jl)  )
            END_2D
         END DO
      ELSEIF( nn_snwfra == 2 ) THEN   ! snow cover depends on hsnow (cice style)
         DO jl = 1, jpl
            DO_2D( 0, 0, 0, 0 )
               pa_s_fra(ji,jj,jl) = ph_s(ji,jj,jl) / ( ph_s(ji,jj,jl) + 0.02_wp )
            END_2D
         END DO
      ENDIF
   END SUBROUTINE ice_var_snwfra_3d

   SUBROUTINE ice_var_snwfra_2d( ph_s, pa_s_fra )
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   ph_s        ! snow thickness
      REAL(wp), DIMENSION(A2D(0)), INTENT(  out) ::   pa_s_fra    ! ice fraction covered by snow
      IF    ( nn_snwfra == 0 ) THEN   ! basic 0 or 1 snow cover
         WHERE( ph_s > 0._wp ) ; pa_s_fra = 1._wp
         ELSEWHERE             ; pa_s_fra = 0._wp
         END WHERE
      ELSEIF( nn_snwfra == 1 ) THEN   ! snow cover depends on hsnow (met-office style)
         pa_s_fra = 1._wp - EXP( -0.2_wp * rhos * ph_s )
      ELSEIF( nn_snwfra == 2 ) THEN   ! snow cover depends on hsnow (cice style)
         pa_s_fra = ph_s / ( ph_s + 0.02_wp )
      ENDIF
   END SUBROUTINE ice_var_snwfra_2d

   !!--------------------------------------------------------------------------
   !! INTERFACE ice_var_snwblow
   !!
   !! ** Purpose :   Compute distribution of precip over the ice
   !!
   !!                Snow accumulation in one thermodynamic time step
   !!                snowfall is partitionned between leads and ice.
   !!                If snow fall was uniform, a fraction (1-at_i) would fall into leads
   !!                but because of the winds, more snow falls on leads than on sea ice
   !!                and a greater fraction (1-at_i)^beta of the total mass of snow
   !!                (beta < 1) falls in leads.
   !!                In reality, beta depends on wind speed,
   !!                and should decrease with increasing wind speed but here, it is
   !!                considered as a constant. an average value is 0.66
   !!--------------------------------------------------------------------------
!!gm  I think it can be usefull to set this as a FUNCTION, not a SUBROUTINE....
   SUBROUTINE ice_var_snwblow_2d( pin, pout )
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) :: pin   ! previous fraction lead ( 1. - a_i_b )
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) :: pout
      pout = ( 1._wp - ( pin )**rn_snwblow )
   END SUBROUTINE ice_var_snwblow_2d

#else
   !!----------------------------------------------------------------------
   !!   Default option         Dummy module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icevar
