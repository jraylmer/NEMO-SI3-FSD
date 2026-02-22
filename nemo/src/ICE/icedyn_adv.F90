MODULE icedyn_adv
   !!======================================================================
   !!                       ***  MODULE icedyn_adv   ***
   !!   sea-ice: advection
   !!======================================================================
   !! History :  4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_dyn_adv   : advection of sea ice variables
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE phycst         ! physical constant
   USE sbc_oce , ONLY : nn_fsbc   ! frequency of sea-ice call
   USE ice            ! sea-ice: variables
   USE icevar         ! sea-ice: operations
   USE icedyn_adv_pra ! sea-ice: advection scheme (Prather)
   USE icedyn_adv_umx ! sea-ice: advection scheme (ultimate-macho)
   USE icectl         ! sea-ice: control prints
   USE icefsd  , ONLY : a_ifsd, ice_fsd_dia   ! sea-ice: floe size distribution
   !
   USE in_out_manager ! I/O manager
   USE iom            ! I/O manager library
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_dyn_adv        ! called by icestp
   PUBLIC   ice_dyn_adv_init   ! called by icedyn

   INTEGER ::              nice_adv   ! choice of the type of advection scheme
   !                                       ! associated indices:
   INTEGER, PARAMETER ::   np_advPRA = 1   ! Prather scheme
   INTEGER, PARAMETER ::   np_advUMx = 2   ! Ultimate-Macho scheme
   !
   !                            !!** ice-advection namelist (namdyn_adv) **
   LOGICAL ::   ln_adv_Pra       ! Prather        advection scheme
   LOGICAL ::   ln_adv_UMx       ! Ultimate-Macho advection scheme
   INTEGER ::   nn_UMx           ! order of the UMx advection scheme   
   !
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL licence     (./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_dyn_adv( kt ) 
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_dyn_adv ***
      !!                    
      !! ** purpose : advection of sea ice
      !!
      !! ** method  : One can choose between 
      !!     a) an Ultimate-Macho scheme (with order defined by nn_UMx) => ln_adv_UMx
      !!     b) and a second order Prather scheme => ln_adv_Pra
      !!
      !! ** action :
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! number of iteration
      !
      ! Floe size distribution: local variables for tendency diagnostics:
      REAL(wp), DIMENSION(A2D(0),nn_nfsd,jpl) ::   za_ifsdb   ! floe size distribution before advection
      REAL(wp), DIMENSION(A2D(0),jpl)         ::   za_ib      ! sea ice concentration before advection
      !
      !!---------------------------------------------------------------------
      !
      ! controls
      IF( ln_timing    )   CALL timing_start('icedyn_adv')                                                             ! timing
      IF( ln_icediachk )   CALL ice_cons_hsm(0, 'icedyn_adv', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_adv: sea-ice advection'
         WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      ! Floe size distribution: save before values for diagnostics:
      IF( ln_fsd ) THEN
         za_ib(A2D(0),:)      = a_i(A2D(0),:)
         za_ifsdb(A2D(0),:,:) = a_ifsd(A2D(0),:,:)
      ENDIF
      !
      !---------------!
      !== Advection ==!
      !---------------!
      SELECT CASE( nice_adv )
      !                                !-----------------------!
      CASE( np_advUMx )                ! ULTIMATE-MACHO scheme !
         !                             !-----------------------!
         CALL ice_dyn_adv_umx( nn_UMx, kt, u_ice, v_ice, h_i, h_s, h_ip, &
            &                          ato_i, v_i, v_s, sv_i, oa_i, a_i, a_ip, v_ip, v_il, e_s, e_i, szv_i )
         !                             !-----------------------!
      CASE( np_advPRA )                ! PRATHER scheme        !
         !                             !-----------------------!
         CALL ice_dyn_adv_pra(         kt, u_ice, v_ice, h_i, h_s, h_ip, &
            &                          ato_i, v_i, v_s, sv_i, oa_i, a_i, a_ip, v_ip, v_il, e_s, e_i, szv_i, &
            &                          a_ifsd )
      END SELECT

      !------------
      ! diagnostics
      !------------
      diag_trp_ei(:,:) =    SUM(SUM( e_i  (A2D(0),:,:) - e_i_b  (A2D(0),:,:), dim=4 ), dim=3 ) * r1_Dt_ice
      diag_trp_es(:,:) =    SUM(SUM( e_s  (A2D(0),:,:) - e_s_b  (A2D(0),:,:), dim=4 ), dim=3 ) * r1_Dt_ice
      IF( nn_icesal == 4 ) THEN
         diag_trp_sv(:,:) = SUM(SUM( szv_i(A2D(0),:,:) - szv_i_b(A2D(0),:,:), dim=4 ), dim=3 ) * r1_Dt_ice
      ELSE
         diag_trp_sv(:,:) = SUM(     sv_i (A2D(0),:)   - sv_i_b (A2D(0),:)           , dim=3 ) * r1_Dt_ice
      ENDIF
      diag_trp_vi(:,:) =    SUM(     v_i  (A2D(0),:)   - v_i_b  (A2D(0),:)           , dim=3 ) * r1_Dt_ice
      diag_trp_vs(:,:) =    SUM(     v_s  (A2D(0),:)   - v_s_b  (A2D(0),:)           , dim=3 ) * r1_Dt_ice
      IF( iom_use('icemtrp') )   CALL iom_put( 'icemtrp' ,  diag_trp_vi * rhoi          )   ! ice mass transport
      IF( iom_use('snwmtrp') )   CALL iom_put( 'snwmtrp' ,  diag_trp_vs * rhos          )   ! snw mass transport
      IF( iom_use('salmtrp') )   CALL iom_put( 'salmtrp' ,  diag_trp_sv * rhoi * 1.e-03 )   ! salt mass transport (kg/m2/s)
      IF( iom_use('dihctrp') )   CALL iom_put( 'dihctrp' , -diag_trp_ei                 )   ! advected ice heat content (W/m2)
      IF( iom_use('dshctrp') )   CALL iom_put( 'dshctrp' , -diag_trp_es                 )   ! advected snw heat content (W/m2)

      ! Floe size distribution diagnostics: tendencies due to advection:
      IF( ln_fsd )   CALL ice_fsd_dia( 'adv', za_ifsdb, a_ifsd(A2D(0),:,:), za_ib, a_i(A2D(0),:) )

      ! controls
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icedyn_adv', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icectl    )   CALL ice_prt     (kt, iiceprt, jiceprt,-1, ' - ice dyn & trp - ')                           ! prints
      IF( sn_cfctl%l_prtctl )   &
         &                 CALL ice_prt3D('icedyn_adv')                                                                ! prints
      IF( ln_timing    )   CALL timing_stop ('icedyn_adv')                                                             ! timing
      !
   END SUBROUTINE ice_dyn_adv


   SUBROUTINE ice_dyn_adv_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_dyn_adv_init  ***
      !!
      !! ** Purpose :   Physical constants and parameters linked to the ice
      !!                dynamics
      !!
      !! ** Method  :   Read the namdyn_adv namelist and check the ice-dynamic
      !!                parameter values called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namdyn_adv
      !!-------------------------------------------------------------------
      INTEGER ::   ios, ioptio   ! Local integer output status for namelist read
      !!
      NAMELIST/namdyn_adv/ ln_adv_Pra, ln_adv_UMx, nn_UMx
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namdyn_adv)
      READ_NML_CFG(numnam_ice,namdyn_adv)
      IF(lwm) WRITE( numoni, namdyn_adv )
      !
      IF(lwp) THEN                     ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_dyn_adv_init: ice parameters for ice dynamics '
         WRITE(numout,*) '~~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_adv:'
         WRITE(numout,*) '      type of advection scheme (Prather)             ln_adv_Pra = ', ln_adv_Pra 
         WRITE(numout,*) '      type of advection scheme (Ulimate-Macho)       ln_adv_UMx = ', ln_adv_UMx 
         WRITE(numout,*) '         order of the Ultimate-Macho scheme          nn_UMx     = ', nn_UMx
      ENDIF
      !
      !                             !== set the choice of ice advection ==!
      ioptio = 0 
      IF( ln_adv_Pra ) THEN   ;   ioptio = ioptio + 1   ;   nice_adv = np_advPRA    ;   ENDIF
      IF( ln_adv_UMx ) THEN   ;   ioptio = ioptio + 1   ;   nice_adv = np_advUMx    ;   ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'ice_dyn_adv_init: choose one and only one ice adv. scheme (ln_adv_Pra or ln_adv_UMx)' )
      !
      IF( ln_adv_Pra )   CALL adv_pra_init  !* read or initialize all required files
      !
   END SUBROUTINE ice_dyn_adv_init

#else
   !!----------------------------------------------------------------------
   !!   Default option         Empty Module           NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icedyn_adv

