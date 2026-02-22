MODULE icethd
   !!======================================================================
   !!                  ***  MODULE icethd   ***
   !!   sea-ice : master routine for thermodynamics
   !!======================================================================
   !! History :  1.0  !  2000-01  (M.A. Morales Maqueda, H. Goosse, T. Fichefet) original code 1D
   !!            4.0  !  2018     (many people)       SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_thd       : thermodynamics of sea ice
   !!   ice_thd_init  : initialisation of sea-ice thermodynamics
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE phycst         ! physical constants
   USE par_oce
   USE ice            ! sea-ice: variables
   USE sbc_oce , ONLY : sss_m, sst_m, frq_m, sprecip
   USE sbc_ice , ONLY : qsr_ice, qns_ice, dqns_ice, evap_ice, qprec_ice, qml_ice, qcn_ice, qtr_ice_top
   USE snwthd_prec    ! sea-ice: snow fall, sublimation & deposition
   USE icethd_zdf     ! sea-ice: vertical heat diffusion
   USE icethd_dh      ! sea-ice: ice growth and melt
   USE snwthd_dh      ! sea-ice: snow melt
   USE icethd_da      ! sea-ice: lateral melting
   USE icethd_sal     ! sea-ice: salinity
   USE icethd_do      ! sea-ice: growth in open water
   USE icethd_pnd     ! sea-ice: melt ponds
   USE iceitd  , ONLY : ice_itd_rem
   USE icefsd  , ONLY : a_ifsd, ice_fsd_dia
   USE icecor         ! sea-ice: corrections
   USE icectl         ! sea-ice: control print
   !
   USE in_out_manager ! I/O manager
   USE iom            , ONLY : iom_miss_val, iom_put       ! I/O manager library
   USE lib_mpp        ! MPP library
   USE lbclnk         ! lateral boundary conditions (or mpp links)
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_thd         ! called by limstp module
   PUBLIC   ice_thd_init    ! called by ice_init

   LOGICAL , ALLOCATABLE, DIMENSION(:,:,:) ::   llmsk
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   zq_top      ! heat for surface ablation                    (J.m-2)
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   zq_bot      ! heat for bottom ablation                     (J.m-2)
   REAL(wp), ALLOCATABLE, DIMENSION(:,:)   ::   zf_tt       ! Heat budget to determine melting or freezing (W.m-2)

   ! sanity checks
   CHARACTER(LEN=50)      ::   clname="cfl_icesalt.ascii"    ! ascii filename
   INTEGER , DIMENSION(3) ::   iloc
   REAL(wp)               ::   zcfl_drain_max, zcfl_flush_max
   INTEGER                ::   numcfl                        ! outfile unit
   
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_thd( kt )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_thd  ***
      !!
      !! ** Purpose : This routine manages ice thermodynamics
      !!
      !! ** Action : - computation of oceanic sensible heat flux at the ice base
      !!                              energy budget in the leads
      !!                              net fluxes on top of ice and of ocean
      !!             - selection of grid cells with ice
      !!                - call ice_thd_zdf  for vertical heat diffusion
      !!                - call ice_thd_dh   for vertical ice growth and melt
      !!                - call ice_thd_pnd  for melt ponds
      !!                - call ice_thd_temp to  retrieve temperature from ice enthalpy
      !!                - call ice_thd_sal  for ice desalination
      !!                - call ice_thd_temp to  retrieve temperature from ice enthalpy
      !!                - call ice_thd_mono for extra lateral ice melt if active virtual thickness distribution
      !!                - call ice_thd_da   for lateral ice melt
      !!             - back to the geographic grid
      !!                - call ice_thd_rem  for remapping thickness distribution
      !!                - call ice_thd_do   for ice growth in leads
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! number of iteration
      !
      REAL(wp), DIMENSION(A2D(0)) ::   zevap_rema   ! remaining of evaporation after snow sublimation (in kg/m2)
      !
      REAL(wp), DIMENSION(A2D(0),nn_nfsd,jpl) ::   za_ifsdb_da   ! FSD before lateral melt for FSD diagnostics
      REAL(wp), DIMENSION(A2D(0),jpl)         ::   za_ib_da      ! a_i before lateral melt for FSD diagnostics
      !
      INTEGER ::   ji, jj, jk, jl   ! dummy loop indices
      !!-------------------------------------------------------------------

      ! controls
      IF( ln_timing    )   CALL timing_start('icethd')                                                             ! timing
      IF( ln_icediachk )   CALL ice_cons_hsm(0, 'icethd', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft) ! conservation
      IF( ln_icediachk )   CALL ice_cons2D  (0, 'icethd',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft) ! conservation

      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'ice_thd: sea-ice thermodynamics'
         WRITE(numout,*) '~~~~~~~'
      ENDIF

      ALLOCATE( zq_top(A2D(0)), zq_bot(A2D(0)), zf_tt(A2D(0)) )
      
      ! convergence tests
      IF( ln_zdf_chkcvg ) THEN
         ALLOCATE( ztice_cvgerr(A2D(0),jpl) , ztice_cvgstp(A2D(0),jpl) )
         ztice_cvgerr = 0._wp ; ztice_cvgstp = 0._wp
      ENDIF
      !
      IF( ln_sal_chk )   CALL ice_thd_salchk( kt, 1 )
      !
      !-------------------------------------------------------------------------------------------!
      ! Thermodynamic computation (only on grid points covered by ice) => loop over ice categories
      !-------------------------------------------------------------------------------------------!
      !
      CALL ice_thd_frazil             !--- frazil ice: collection thickness (ht_i_new) & fraction of frazil (fraz_frac)
      !
      ! Save a_ifsd before lateral melt (ice_thd_da; only thing affecting this variable in jl loop below):
      IF( ln_fsd )   za_ifsdb_da(A2D(0),:,:) = a_ifsd(A2D(0),:,:)
      !
      DO jl = 1, jpl

         l_ice_present(A2D(0)) = .FALSE.
         DO_2D( 0, 0, 0, 0 )
            IF( a_i(ji,jj,jl) > epsi10 )   l_ice_present(ji,jj) = .TRUE. ! select ice covered grid points
         END_2D
        
         IF( ANY( l_ice_present(A2D(0)) ) ) THEN 
            !
                              CALL ice_thd_unit_convert( jl, 1 )            ! --- & Change units of e_i, e_s from J/m2 to J/m3 --- !
            !
            dh_s_tot  (:,:) = 0._wp                              ! --- some init --- !  (important to have them here)
            dh_i_bom  (:,:) = 0._wp ; dh_i_itm(:,:) = 0._wp
            dh_i_sub  (:,:) = 0._wp ; dh_i_bog(:,:) = 0._wp
            dh_snowice(:,:) = 0._wp ; dh_s_itm(:,:) = 0._wp
            zevap_rema(:,:) = 0._wp
            !
            IF( ln_icedH )    CALL snw_thd_prec( jl, zevap_rema )                        ! --- Snow fall and sublimation/deposition --- !
            !                                               ==> out: zevap_rema
            !
                              CALL ice_thd_zdf ( jl )                                    ! --- Ice & Snow temperature diffusion --- !
            !                                               ==> out: qtr_ice_bot, qcn_ice_bot, qcn_ice_top, cnd_ice, qcn_ice
            !                  
                              CALL ice_thd_bdg ( jl )                                    ! --- Ice & Snow Top/Bottom heat budgets --- !
            !                                               ==> out: qml_ice, zq_top, zq_bot, zf_tt
            !                  
            IF( ln_icedH )    CALL snw_thd_dh  ( jl, zq_top )                            ! --- Snow melt --- !
            !                                            <==> inout: zq_top 
            !
            IF( ln_icedH )    CALL ice_thd_dh  ( jl, zq_top, zq_bot, zf_tt, zevap_rema ) ! --- Ice Growing/Melting & snow-ice --- !
            !                                                ==> in: zq_top, zq_bot, zf_tt, zevap_rema
            !                                                                              
                              CALL ice_thd_temp( jl )                                    ! --- Ice Temperature update --- !
            !
            !                  
                              CALL ice_thd_sal ( jl )                                    ! --- Ice Salinity --- !
            !
            !
                              CALL ice_thd_temp( jl )                                    ! --- Ice Temperature update --- !
            !
            IF( ln_fsd )      za_ib_da(A2D(0),jl) = a_i(A2D(0),jl)                       ! --- Save a_i before lateral melt for FSD diagnostics --- !
            !
            IF( ln_icedH .AND. ln_virtual_itd ) &
               &              CALL ice_thd_mono( jl )                                    ! --- Extra lateral melting if virtual_itd --- !
            !
            IF( ln_icedA )    CALL ice_thd_da  ( jl )                                    ! --- Ice Lateral melting --- !
            !
                              CALL ice_thd_unit_convert( jl, 2 )            ! --- Change units of e_i, e_s from J/m3 to J/m2 --- !
            !
         ENDIF ! l_ice_present
         !
         !
      END DO ! jl loop 
      !
      ! Floe size distribution: tendency diagnostics
      !
      ! NOTE: if anything is added after ice_thd_da in loop above (or order changed)
      ! ----- will need to save a_i after lateral melt as well. Currently it is not needed;
      !       we can just use a_i because it is the most recent process.
      !
      IF( ln_fsd )         CALL ice_fsd_dia( 'lam', za_ifsdb_da, a_ifsd(A2D(0),:,:), za_ib_da, a_i(A2D(0),:) )
      !
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'icethd', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      IF( ln_icediachk )   CALL ice_cons2D  (1, 'icethd',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)
      !
      IF ( ln_pnd .AND. ln_icedH ) &
         &                    CALL ice_thd_pnd                      ! --- Melt ponds --- !
      !
      IF( jpl > 1  )          CALL ice_itd_rem( kt )                ! --- Transport ice between thickness categories --- !
      !
      IF( ln_icedO )          CALL ice_thd_do                       ! --- Frazil ice growth in leads --- !
      !
                              CALL ice_cor( kt , 2 )                ! --- Corrections --- !
      !
      oa_i(A2D(0),:) = oa_i(A2D(0),:) + a_i(A2D(0),:) * rDt_ice     ! --- Ice natural aging incrementation
      !
      !                                                             ! --- LBC for the halos --- !
      CALL lbc_lnk( 'icethd', a_i , 'T', 1._wp, v_i , 'T', 1._wp, v_s , 'T', 1._wp, sv_i, 'T', 1._wp, oa_i, 'T', 1._wp, &
         &                    t_su, 'T', 1._wp, a_ip, 'T', 1._wp, v_ip, 'T', 1._wp, v_il, 'T', 1._wp )
      CALL lbc_lnk( 'icethd', e_i , 'T', 1._wp, e_s , 'T', 1._wp, szv_i,'T', 1._wp )
      !
      IF( ln_fsd ) CALL lbc_lnk( 'icethd', a_ifsd, 'T', 1._wp )
      !
      at_i(:,:) = SUM( a_i, dim=3 )
      DO_2D( 0, 0, 0, 0 )                                           ! --- Ice velocity corrections
         IF( at_i(ji,jj) == 0._wp ) THEN   ! if ice has melted
            IF( at_i(ji+1,jj) == 0._wp )   u_ice(ji  ,jj) = 0._wp   ! right side
            IF( at_i(ji-1,jj) == 0._wp )   u_ice(ji-1,jj) = 0._wp   ! left side
            IF( at_i(ji,jj+1) == 0._wp )   v_ice(ji,jj  ) = 0._wp   ! upper side
            IF( at_i(ji,jj-1) == 0._wp )   v_ice(ji,jj-1) = 0._wp   ! bottom side
         ENDIF
      END_2D
      CALL lbc_lnk( 'icethd', u_ice, 'U', -1.0_wp, v_ice, 'V', -1.0_wp )
      !
      ! convergence tests
      IF( ln_zdf_chkcvg ) THEN
         CALL iom_put( 'tice_cvgerr', ztice_cvgerr ) ; DEALLOCATE( ztice_cvgerr )
         CALL iom_put( 'tice_cvgstp', ztice_cvgstp ) ; DEALLOCATE( ztice_cvgstp )
      ENDIF
      !
      ! sanity checks for salt drainage and flushing
      IF( ln_sal_chk )   CALL ice_thd_salchk( kt, 2 )
      
      DEALLOCATE( zq_top, zq_bot, zf_tt )

      ! controls
      IF( ln_icectl )   CALL ice_prt    (kt, iiceprt, jiceprt, 1, ' - ice thermodyn. - ') ! prints
      IF( sn_cfctl%l_prtctl )   &
        &               CALL ice_prt3D  ('icethd')                                        ! prints
      IF( ln_timing )   CALL timing_stop('icethd')                                        ! timing
      !
   END SUBROUTINE ice_thd


   SUBROUTINE ice_thd_bdg( jl_cat )
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_bdg ***
      !!
      !! ** Purpose :   Computes heat 
      !!
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   jl_cat
      INTEGER  ::   ji, jj   ! dummy loop indices
      !!-------------------------------------------------------------------
      IF( .NOT.ln_cndflx .OR. ln_cndemulate ) THEN
         DO_2D( 0, 0, 0, 0 )
            IF( l_ice_present(ji,jj) ) THEN
               IF( t_su(ji,jj,jl_cat) >= rt0 ) THEN
                  qml_ice(ji,jj,jl_cat) =   qns_ice    (ji,jj,jl_cat) + qsr_ice    (ji,jj,jl_cat)  &
                     &                    - qtr_ice_top(ji,jj,jl_cat) - qcn_ice_top(ji,jj,jl_cat)
               ELSE
                  qml_ice(ji,jj,jl_cat) = 0._wp
               ENDIF
               !
            ENDIF
         END_2D
      ENDIF
      !
      DO_2D( 0, 0, 0, 0 )      
         IF( l_ice_present(ji,jj) ) THEN
            zq_top(ji,jj) = MAX( 0._wp, qml_ice(ji,jj,jl_cat) * rDt_ice )
            zf_tt (ji,jj) = qcn_ice_bot(ji,jj,jl_cat) + qsb_ice_bot(ji,jj) + fhld(ji,jj) + qtr_ice_bot(ji,jj,jl_cat) * frq_m(ji,jj)
            zq_bot(ji,jj) = MAX( 0._wp, zf_tt(ji,jj) * rDt_ice )
         ENDIF
      END_2D
      !
   END SUBROUTINE ice_thd_bdg


   SUBROUTINE ice_thd_temp( jl_cat )
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_temp ***
      !!
      !! ** Purpose :   Computes sea ice temperature (Kelvin) from enthalpy
      !!
      !! ** Method  :   Formula (Bitz and Lipscomb, 1999)
      !!-------------------------------------------------------------------
      INTEGER, INTENT(in) ::   jl_cat
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   ztmelts, zbbb, zccc  ! local scalar
      !!-------------------------------------------------------------------
      ! Recover ice temperature
      DO jk = 1, nlay_i
         DO_2D( 0, 0, 0, 0 )
            IF( l_ice_present(ji,jj) ) THEN
               ztmelts = -rTmlt * sz_i(ji,jj,jk,jl_cat)
               ! Conversion q(S,T) -> T (second order equation)
               zbbb = ( rcp - rcpi ) * ztmelts + e_i(ji,jj,jk,jl_cat) * r1_rhoi - rLfus
               zccc = SQRT( MAX( zbbb * zbbb - 4._wp * rcpi * rLfus * ztmelts, 0._wp ) )
               t_i(ji,jj,jk,jl_cat) = rt0 - ( zbbb + zccc ) * 0.5_wp * r1_rcpi
            ELSE
               t_i(ji,jj,jk,jl_cat) = rt0
            ENDIF
         END_2D
      END DO
      !
   END SUBROUTINE ice_thd_temp


   SUBROUTINE ice_thd_mono( jl_cat )
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_mono ***
      !!
      !! ** Purpose :   Lateral melting in case virtual_itd
      !!                          ( dA = A/2h dh )
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   jl_cat
      INTEGER  ::   ji,jj              ! dummy loop indices
      REAL(wp) ::   zhi_bef            ! ice thickness before thermo
      REAL(wp) ::   zdh_mel, zda_mel   ! net melting
      REAL(wp) ::   zvi, zvs           ! ice/snow volumes
      !!-----------------------------------------------------------------------
      !
      DO_2D( 0, 0, 0, 0 )
         !
         IF( l_ice_present(ji,jj) ) THEN
            !
            zdh_mel = MIN( 0._wp, dh_i_itm(ji,jj) + dh_i_sum_3d(ji,jj,jl_cat) + dh_i_bom(ji,jj) &
               &                                          + dh_snowice(ji,jj) + dh_i_sub(ji,jj) )
            !
            IF( zdh_mel < 0._wp .AND. a_i(ji,jj,jl_cat) > 0._wp )  THEN
               zvi     = a_i(ji,jj,jl_cat) * h_i(ji,jj,jl_cat)
               zvs     = a_i(ji,jj,jl_cat) * h_s(ji,jj,jl_cat)
               ! lateral melting = concentration change
               zhi_bef = h_i(ji,jj,jl_cat) - zdh_mel
               zda_mel = MAX( -a_i(ji,jj,jl_cat) , a_i(ji,jj,jl_cat) * zdh_mel / ( 2._wp * MAX( zhi_bef, epsi20 ) ) )
               a_i(ji,jj,jl_cat) = MAX( epsi20, a_i(ji,jj,jl_cat) + zda_mel )
               ! adjust thickness
               h_i(ji,jj,jl_cat) = zvi / a_i(ji,jj,jl_cat)
               h_s(ji,jj,jl_cat) = zvs / a_i(ji,jj,jl_cat)
               ! retrieve total concentration
               at_i(ji,jj) = a_i(ji,jj,jl_cat)
            END IF
            !
         ENDIF
         !
      END_2D
      !
   END SUBROUTINE ice_thd_mono

   SUBROUTINE ice_thd_salchk( kt, kn )
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_salchk ***
      !!
      !! ** Purpose :   checking salt drainage and flushing
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, kn 
      !
      !INTEGER ::   jk   ! dummy loop indices
      !!-----------------------------------------------------------------------
     
      ! sanity checks for salt flushing and drainage
      IF( kn == 1 ) THEN
         !
         ALLOCATE( llmsk(A2D(0),jpl) )
         ALLOCATE( zcfl_flush(A2D(0),jpl) , zcfl_drain(A2D(0),jpl), zsneg_flush(A2D(0),jpl) , zsneg_drain(A2D(0),jpl) )
         !
         zcfl_flush = 0._wp ; zcfl_drain = 0._wp
         zsneg_flush = 0._wp ; zsneg_drain = 0._wp
         !
      ELSEIF( kn == 2 ) THEN
         !
         CALL iom_put( 'sice_flush_neg', zsneg_flush ) ; DEALLOCATE( zsneg_flush )
         CALL iom_put( 'sice_drain_neg', zsneg_drain ) ; DEALLOCATE( zsneg_drain )
         !
         CALL iom_put( 'cfl_flush', zcfl_flush )
         CALL iom_put( 'cfl_drain', zcfl_drain )

         !                    ! calculate maximum values and locations
         llmsk(Nis0:Nie0,Njs0:Nje0,:) = h_i(Nis0:Nie0,Njs0:Nje0,:) > rn_himin        ! define only where h > 0.10m
         CALL mpp_maxloc( 'icethd', zcfl_drain, llmsk, zcfl_drain_max, iloc )
         CALL mpp_maxloc( 'icethd', zcfl_flush, llmsk, zcfl_flush_max, iloc )

         IF( lwp ) THEN       ! write out to file
            WRITE(numcfl,FMT='(2x,i6,3x,a10,4x,f8.4,1x,i4,1x,i4,1x,i4)') kt, 'Max Cdrain', zcfl_drain_max, iloc(1), iloc(2), iloc(3)
            WRITE(numcfl,FMT='(11x,     a10,4x,f8.4,1x,i4,1x,i4,1x,i4)')     'Max Cflush', zcfl_flush_max, iloc(1), iloc(2), iloc(3)
         ENDIF
         DEALLOCATE( zcfl_flush, zcfl_drain )
         DEALLOCATE( llmsk )

         IF( kt == nitend .AND. lwp )   CLOSE( numcfl )
      ENDIF
      
   END SUBROUTINE ice_thd_salchk


   SUBROUTINE ice_thd_unit_convert( kl, kn )
      !!-----------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_unit_convert ***
      !!
      !! ** Purpose :   to handle the conversion of units from J/m2 to J/m3 and reverse
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kl   ! index of the ice category
      INTEGER, INTENT(in) ::   kn   ! 1= from J/m2 to J/m3   ;   2= from J/m3 to J/m2
      !
      INTEGER ::   jk   ! dummy loop indices
      !!-----------------------------------------------------------------------
      !
      SELECT CASE( kn )
      !                    !-------------------------!
      CASE( 1 )            !==  from J/m2 to J/m3  ==!
         !                 !-------------------------!
         ! --- Change units of e_i, e_s from J/m2 to J/m3 --- !
         ! Here we make sure that we don't divide by very small, but physically
         ! meaningless, products of sea ice thicknesses/snow depths and sea ice concentration
         DO jk = 1, nlay_i
            WHERE( (h_i(:,:,kl) * a_i(:,:,kl)) > epsi20 )
               e_i(:,:,jk,kl) = e_i(:,:,jk,kl) / (h_i(:,:,kl) * a_i(:,:,kl)) * nlay_i
            ELSEWHERE
               e_i(:,:,jk,kl) = 0._wp
            ENDWHERE
         END DO
         DO jk = 1, nlay_s
            WHERE( (h_s(:,:,kl) * a_i(:,:,kl)) > epsi20 )
               e_s(:,:,jk,kl) = e_s(:,:,jk,kl) / (h_s(:,:,kl) * a_i(:,:,kl)) * nlay_s
            ELSEWHERE
               e_s(:,:,jk,kl) = 0._wp
            ENDWHERE
         END DO
         !
         !                 !-------------------------!
      CASE( 2 )            !==  from J/m3 to J/m2  ==!
         !                 !-------------------------!
         ! --- Change units of e_i, e_s from J/m3 to J/m2 --- !
         DO jk = 1, nlay_i
            e_i(:,:,jk,kl) = e_i(:,:,jk,kl) * h_i(:,:,kl) * a_i(:,:,kl) * r1_nlay_i
         END DO
         DO jk = 1, nlay_s
            e_s(:,:,jk,kl) = e_s(:,:,jk,kl) * h_s(:,:,kl) * a_i(:,:,kl) * r1_nlay_s
         END DO
         !
         ! Change thickness to volume (replaces routine ice_var_eqv2glo)
         v_i (:,:,kl) = h_i(:,:,kl) * a_i(:,:,kl)
         v_s (:,:,kl) = h_s(:,:,kl) * a_i(:,:,kl)
         sv_i(:,:,kl) = s_i(:,:,kl) * v_i(:,:,kl)
         oa_i(:,:,kl) = o_i(:,:,kl) * a_i(:,:,kl)
         DO jk = 1, nlay_i
            szv_i(:,:,jk,kl) = sz_i(:,:,jk,kl) * v_i(:,:,kl) * r1_nlay_i
         END DO
      END SELECT
      !
   END SUBROUTINE ice_thd_unit_convert


   SUBROUTINE ice_thd_init
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE ice_thd_init ***
      !!
      !! ** Purpose :   Physical constants and parameters associated with
      !!                ice thermodynamics
      !!
      !! ** Method  :   Read the namthd namelist and check the parameters
      !!                called at the first timestep (nit000)
      !!
      !! ** input   :   Namelist namthd
      !!-------------------------------------------------------------------
      INTEGER  ::   ios   ! Local integer output status for namelist read
      !!
      NAMELIST/namthd/ ln_icedH, ln_icedA, ln_icedO, ln_leadhfx
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice,namthd)
      READ_NML_CFG(numnam_ice,namthd)
      IF(lwm) WRITE( numoni, namthd )
      !
      IF(lwp) THEN                          ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_thd_init: Ice Thermodynamics'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namthd:'
         WRITE(numout,*) '      activate ice thick change from top/bot (T) or not (F)                ln_icedH   = ', ln_icedH
         WRITE(numout,*) '      activate lateral melting (T) or not (F)                              ln_icedA   = ', ln_icedA
         WRITE(numout,*) '      activate ice growth in open-water (T) or not (F)                     ln_icedO   = ', ln_icedO
         WRITE(numout,*) '      heat in the leads is used to melt sea-ice before warming the ocean   ln_leadhfx = ', ln_leadhfx
     ENDIF
      !
                       CALL ice_thd_zdf_init   ! set ice heat diffusion parameters
      IF( ln_icedA )   CALL ice_thd_da_init    ! set ice lateral melting parameters
      IF( ln_icedO )   CALL ice_thd_do_init    ! set ice growth in open water parameters
                       CALL ice_thd_sal_init   ! set ice salinity parameters
                       CALL ice_thd_pnd_init   ! set melt ponds parameters
      !
      IF( ln_sal_chk ) THEN
         ! create output ascii file
         CALL ctl_opn( numcfl, clname, 'UNKNOWN', 'FORMATTED', 'SEQUENTIAL', 1, numout, lwp, 1 )
         WRITE(numcfl,*) 'Timestep  Direction   Max C     i    j    k'
         WRITE(numcfl,*) '*******************************************'
      ENDIF
      !
   END SUBROUTINE ice_thd_init

#else
   !!----------------------------------------------------------------------
   !!   Default option         Dummy module          NO  SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icethd
