MODULE icefsd
   !!======================================================================
   !!                       ***  MODULE icefsd ***
   !!   sea-ice : floe size distribution
   !!======================================================================
   !! History :  5.0  !  2024     (J.R. Aylmer)         Original code based
   !!                                                   on CPOM-CICE and
   !!                                                   CICE/Icepack
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3' :                                     SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_fsd_init : namelist read
   !!----------------------------------------------------------------------
   USE par_ice          ! SI3 parameters
   USE ice              ! sea-ice: variables

   USE in_out_manager   ! I/O manager (needed for lwm and lwp logicals)
   USE iom              ! I/O manager library (needed for iom_put)
   USE lib_mpp          ! MPP library (needed for read_nml_substitute.h90)

   IMPLICIT NONE
   PRIVATE

   INTERFACE ice_fsd_cor
      !!----------------------------------------------------------------------
      !!                  ***  INTERFACE ice_fsd_cor ***
      !!----------------------------------------------------------------------
      !! ** Purpose :   Generic interface for applying numerical corrections and
      !!                renormalisation for the floe size distribution (FSD)
      !!----------------------------------------------------------------------
      MODULE PROCEDURE fsd_cor_1d   ! e.g., a_ifsd(ji,jj,:,jl)
      MODULE PROCEDURE fsd_cor_2d   ! e.g., a_ifsd(ji,jj,:,:)
      MODULE PROCEDURE fsd_cor_4d   ! e.g., a_ifsd(:,:,:,:)
   END INTERFACE

   PUBLIC ::   ice_fsd_init               ! routine called by ice_init
   PUBLIC ::   ice_fsd_istate             ! routine called by ice_istate, ice_rst_read
   PUBLIC ::   ice_fsd_wri                ! routine called by ice_stp
   PUBLIC ::   ice_fsd_dia                ! routine called by various routines
   PUBLIC ::   ice_fsd_brit               ! routine called by ice_dyn
   PUBLIC ::   ice_fsd_part_newice        ! routine called by ice_thd_do
   PUBLIC ::   ice_fsd_add_newice         ! routine called by ice_thd_do
   PUBLIC ::   ice_fsd_weld               ! routine called by ice_thd_do
   PUBLIC ::   ice_fsd_thd                ! routine called by ice_thd_d{a,o}
   PUBLIC ::   ice_fsd_tstep              ! routine called by ice_wav_frac
   PUBLIC ::   fsd_eff_size               ! function called by ice_frm
   PUBLIC ::   floe_size_dist             ! function called by ice_frm
   PUBLIC ::   fsd_peri_dens              ! function called by ice_thd_da
   PUBLIC ::   ice_fsd_cor                ! generic interface: small/negative value corrections and re-normalisation

   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: floe_sl      !: FSD floe size, lower bounds of categories (m)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: floe_sc      !: FSD floe size, centre       of categories (m)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: floe_su      !: FSD floe size, upper bounds of categories (m)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:)   :: floe_ds      !: FSD category widths (m)
   REAL(wp),         ALLOCATABLE, DIMENSION(:)   :: floe_al      !: FSD floe areas, floes of size floe_sl (m2)
   REAL(wp),         ALLOCATABLE, DIMENSION(:)   :: floe_ac      !: FSD floe areas, floes of size floe_sc (m2)
   REAL(wp),         ALLOCATABLE, DIMENSION(:)   :: floe_au      !: FSD floe areas, floes of size floe_su (m2)
   REAL(wp),         ALLOCATABLE, DIMENSION(:)   :: floe_dlog_sc !: FSD cat. centre spacing in log(size) space
   INTEGER ,         ALLOCATABLE, DIMENSION(:,:) :: floe_iweld   !: index of FSD cat. two given FSD cats. can weld to
   INTEGER , PUBLIC                              :: nf_newice    !: index of FSD cat. for new ice in absence of waves (m)

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   a_ifsd      !: FSD per ice thickness category
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   a_ifsd_b    !: FSD at "before" time step (see icestp)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:,:) ::   a_ifsd_b0   !: FSD truly at before time step
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:)   ::   a_i_b0      !: a_i truly at before time step

   ! ** namelist (namfsd) **
   REAL(wp), DIMENSION(100) ::      rn_fsd_catbnd   ! User-defined category limits if nn_nfsd_catini = 0 (below)
   INTEGER  ::   nn_fsd_catini      ! FSD category definition option
   REAL(wp) ::   rn_fsd_smin        ! Minimum floe size (m; nn_fsd_catini >= 1)
   REAL(wp) ::   rn_fsd_smax        ! Minimum floe size (m; nn_fsd_catini >= 1)
   REAL(wp) ::   rn_fsd_spc         ! Category spacing non-linearity parameter (nn_fsd_catini = 2,3)
   INTEGER  ::   nn_fsd_ini         ! FSD init. options (1 = all in largest FSD cat; 2 = imposed power law)
   REAL(wp) ::   rn_fsd_ini_alpha   ! Parameter used for power law initial FSD with nn_fsd_ini = 2 only
   REAL(wp) ::   rn_fsd_s_newice    ! Floe size of new ice in absence of wave field [m]
   REAL(wp) ::   rn_fsd_brit_grad   ! Maximum gradient of number density FSD in log-log space
   REAL(wp) ::   rn_fsd_brit_tres   ! Restoring timescale (units: s)
   REAL(wp) ::   rn_fsd_amin_weld   ! Minimum concentration required for floe welding to take effect
   REAL(wp) ::   rn_fsd_c_weld      ! Floe welding coefficient [m-2.s-1]

   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"

CONTAINS

   FUNCTION floe_size_dist( pa_ifsd, pa_i )
      !!-------------------------------------------------------------------
      !!                   *** FUNCTION floe_size_dist ***
      !! ** Purpose :   Compute floe size distribution (FSD) from prognostic variables
      !! ** Method  :   f(s)ds = int[ L(s,h)ds * g(h)dh ] (integrate over thickness h)
      !! ** Input   :   pa_ifsd(nn_nfsd,jpl)    :  modified-areal floe size-thickness distribution L(s,h)ds
      !!                pa_i   (        jpl)    :  ice thickness distribution, g(h)dh
      !! ** Output  :   floe_size_dist(nn_nfsd) :  floe size distribution, f(s)ds
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(nn_nfsd,jpl), INTENT(in) :: pa_ifsd          ! modified floe size-thickness distribution, L(s,h)ds
      REAL(wp), DIMENSION(jpl)        , INTENT(in) :: pa_i             ! ice thickness distribution, g(h)dh
      REAL(wp), DIMENSION(nn_nfsd)                 :: floe_size_dist   ! floe size distribution, f(s)ds
      INTEGER                                      :: jf               ! dummy loop index
      !!-------------------------------------------------------------------
      floe_size_dist(:) = 0._wp
      DO jf = 1, nn_nfsd
         floe_size_dist(jf) = SUM( pa_ifsd(jf,:) * pa_i(:) )
      ENDDO
   END FUNCTION floe_size_dist


   FUNCTION peri_dens_dist( pa_ifsd, pa_i )
      !!-------------------------------------------------------------------
      !!                   *** FUNCTION peri_dens_dist ***
      !! ** Purpose :   Compute perimeter-density floe size distribution from prognostic variables
      !! ** Method  :   p(s)ds = (pi / floeshape * c) * int[ (L(s,h)/s)ds * g(h)dh ] (integrate over thickness h)
      !! ** Input   :   pa_ifsd(nn_nfsd,jpl)    : modified-areal floe size-thickness distribution L(s,h)ds
      !!                pa_i   (        jpl)    : ice thickness distribution, g(h)dh
      !! ** Output  :   peri_dens_dist(nn_nfsd) : perimeter-density floe size distribution, p(s)ds
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(nn_nfsd,jpl), INTENT(in) :: pa_ifsd          ! modified floe size-thickness distribution, L(s,h)ds
      REAL(wp), DIMENSION(jpl)        , INTENT(in) :: pa_i             ! ice thickness distribution, g(h)dh
      REAL(wp)                                     :: zc               ! sea ice concentration, c
      REAL(wp), DIMENSION(nn_nfsd)                 :: peri_dens_dist   ! perimeter density floe size distribution, p(s)ds
      INTEGER                                      :: jf               ! dummy loop index
      !!-------------------------------------------------------------------
      peri_dens_dist(:) = 0._wp
      zc = SUM( pa_i(:) )  ! sea ice concentration
      IF( zc > epsi10 ) THEN
         DO jf = 1, nn_nfsd
            peri_dens_dist(jf) = rpi * SUM( pa_ifsd(jf,:) * pa_i(:) ) / (rn_floeshape * zc * floe_sc(jf))
         ENDDO
      ENDIF
   END FUNCTION peri_dens_dist


   FUNCTION fsd_peri_dens( pfsd )
      !!-------------------------------------------------------------------
      !!                   *** FUNCTION ice_fsd_peri ***
      !! ** Purpose :   Compute perimeter density from floe size distribution
      !! ** Method  :   rho = (pi / floeshape * a) * int[ (f(s)/s)ds ] (integrate over floe size s)
      !!                where 'a' is the area fraction of ice over the thickness range relevant to f(s)ds
      !! ** Input   :   pfsd(nn_nfsd) : floe size distribution f(s)ds OR L(s,h)ds
      !! ** Output  :   fsd_peri_dens : perimeter density (units: m-1)
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(nn_nfsd), INTENT(in)  ::   pfsd            ! floe size distribution, f(s)ds
      REAL(wp)                                  ::   fsd_peri_dens   ! perimeter density, rho (units: m-1)
      REAL(wp)                                  ::   za              ! sea ice area fraction
      INTEGER                                   ::   jf              ! dummy loop index
      !!-------------------------------------------------------------------
      fsd_peri_dens = 0._wp   ! initialise
      za = SUM(pfsd(:))
      IF( za > epsi10 ) THEN
         DO jf = 1, nn_nfsd
            fsd_peri_dens = fsd_peri_dens + pfsd(jf) / floe_sc(jf)
         ENDDO
         ! Input pfsd can either be 'floe size distribution', f(s)ds, the sum (above to give za) of
         ! which equals the sea ice concentration, or it can be the 'modified-areal floe size-thickness
         ! distribution', i.e., the prognostic variable a_ifsd, the sum of which is 1. Thus either way,
         ! the result below is correctly normalised to the relevant sea ice area (za):
         fsd_peri_dens = fsd_peri_dens * rpi / (rn_floeshape * za)
      ENDIF
   END FUNCTION fsd_peri_dens


   FUNCTION fsd_eff_size( pfsd )
      !!-------------------------------------------------------------------
      !!                   *** FUNCTION fsd_eff_size ***
      !! ** Purpose :   Compute effective floe size from floe size distribution
      !! ** Method  :   seff = a / int[ (f(s)/s)ds ] (integrate over floe size, s)
      !!                where 'a' is the area fraction of ice over the thickness range relevant to f(s)ds
      !! ** Input   :   pfsd(nn_nfsd) : floe size distribution f(s)ds OR L(s,h)ds
      !! ** Output  :   fsd_eff_size  : effective floe size, seff (units: m)
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(nn_nfsd), INTENT(in) ::   pfsd           ! floe size distribution, f(s)ds
      REAL(wp)                                 ::   fsd_eff_size   ! effective floe size (units: m)
      REAL(wp)                                 ::   za             ! sea ice area fraction
      INTEGER                                  ::   jf             ! dummy loop index
      !!-------------------------------------------------------------------
      fsd_eff_size = 0._wp   ! initialise
      za = SUM(pfsd(:))
      IF( za > epsi10 ) THEN
         DO jf = 1, nn_nfsd
            fsd_eff_size = fsd_eff_size + pfsd(jf) / floe_sc(jf)
         ENDDO
         ! Input pfsd can either be 'floe size distribution', f(s)ds, the sum (above to give za) of
         ! which equals the sea ice concentration, or it can be the 'modified-areal floe size-thickness
         ! distribution', i.e., the prognostic variable a_ifsd, the sum of which is 1. Thus either way,
         ! the result below is correctly normalised to the relevant sea ice area (za):
         fsd_eff_size = za / fsd_eff_size
      ENDIF
   END FUNCTION fsd_eff_size


   SUBROUTINE fsd_cor_1d( pa_ifsd_jl )
      !!-------------------------------------------------------------------
      !!                   ***  ROUTINE fsd_cor_1d  ***
      !! ** Purpose :   Remove small/negative values and re-normalise mFSTD
      !! ** Input   :   a_ifsd(ji,jj,:,jl) (i.e., at one grid cell and one thickness category)
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(nn_nfsd), INTENT(inout) ::   pa_ifsd_jl   ! mFSTD (one grid cell, one ITD cat.)
      REAL(wp)                                    ::   ztotfrac     ! for normalisation
      INTEGER                                     ::   jf           ! dummy loop index
      !!-------------------------------------------------------------------
      !
      ! Remove negative and/or very small values in each floe size category:
      WHERE( pa_ifsd_jl <= epsi10 )   pa_ifsd_jl = 0._wp
      !
      ztotfrac = SUM(pa_ifsd_jl(:))   ! should = 1 when properly normalised
      IF(ztotfrac > epsi10) THEN
         DO jf = 1, nn_nfsd
            pa_ifsd_jl(jf) = pa_ifsd_jl(jf) / ztotfrac   ! re-normalise
         ENDDO
      ELSE
         pa_ifsd_jl(:) = 0._wp   ! => ice-free grid cell, set to exactly 0
      ENDIF
      !
   END SUBROUTINE fsd_cor_1d


   SUBROUTINE fsd_cor_2d( pa_ifsd )
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE fsd_cor_2d  ***
      !! ** Purpose :   Remove small/negative values and re-normalise mFSTD
      !! ** Input   :   a_ifsd(ji,jj,:,:) (i.e., at one grid cell)
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(nn_nfsd,jpl), INTENT(inout) ::   pa_ifsd   ! mFSTD (one grid cell, all ITD cats.)
      INTEGER                                         ::   jl        ! dummy loop index
      !!-------------------------------------------------------------------
      DO jl = 1, jpl
         CALL fsd_cor_1d( pa_ifsd(:,jl) )
      ENDDO
   END SUBROUTINE fsd_cor_2d


   SUBROUTINE fsd_cor_4d( pa_ifsd )
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE fsd_cor_4d  ***
      !! ** Purpose :   Remove small/negative values and re-normalise mFSTD
      !! ** Input   :   a_ifsd(:,:,:,:) (i.e., full prognostic mFSTD array a_ifsd)
      !!-------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,nn_nfsd,jpl), INTENT(inout) ::   pa_ifsd      ! Full mFSTD (a_ifsd) array
      INTEGER                                                 ::   ji, jj, jl   ! dummy loop indices
      !!-------------------------------------------------------------------
      DO jl = 1, jpl
         DO_2D(0, 0, 0, 0)
            CALL fsd_cor_1d( pa_ifsd(ji,jj,:,jl) )
         END_2D
      ENDDO
   END SUBROUTINE fsd_cor_4d


   SUBROUTINE ice_fsd_tstep( cdcrn, pa_ifsd, ptendency, pt_elapsed, ksubt )
      !!-------------------------------------------------------------------
      !!                    *** ROUTINE ice_fsd_tstep ***
      !!
      !! ** Purpose :   Evolve the mFSTD under a given tendency using adaptive time stepping
      !!
      !! ** Method  :   Calculate time step restrictions for incrementing the current 'modified-
      !!                -areal' floe size thickness distribution (mFSTD, L(s,h); i.e., the
      !!                prognostic variable a_ifsd) under a specified tendency.
      !!
      !!                Since L(s,h) must be bounded by 0 and 1, there are additional time step
      !!                restrictions which may not be satisfied by the global time step rDt_ice:
      !!
      !!                   0 <= L(s,h) + (dL/dt) * dt <= 1
      !!
      !!                Depending on the sign of dL/dt, one end of the inequality implies the
      !!                maximum allowed dt to update the current L with that tendency (which
      !!                differs for each floe size category). Here, dt is calculated as the most
      !!                restrictive across categories, and then the mFSTD is evolved under the
      !!                given tendency over the time interval dt (or, if it is smaller, the
      !!                global time step minus the time already elapsed from previous iterations).
      !!
      !! ** Notes   :   After evolving L, the time elapsed is increased by dt. The tendency then
      !!                needs to be recomputed and the process repeated until the full time step,
      !!                rDt_ice, has elapsed. This routine carries out one adaptive sub-time step
      !!                and manages general stability checks and warnings. It is called from other
      !!                routines within the following general template structure:
      !!
      !!                   ...
      !!                   zt_elapsed = 0._wp ; isubt = 0   ! <-- start adaptive time stepping
      !!                   DO WHILE( zt_elapsed < rDt_ice )
      !!                      ! <-- calculate tendency -->
      !!                      CALL ice_fsd_tstep( 'name_of_routine', a_ifsd, tendency, zt_elapsed, isubt)
      !!                      ! => here, isubt increased by 1
      !!                      !    and zt_elapsed increased by dt
      !!                      ! <-- ad-hoc checks-->
      !!                   ENDDO
      !!                   ...
      !!
      !!                See Horvat and Tziperman (2017; App. A), for further details on time stepping.
      !!
      !! ** Input   :   cdcrn              : name of calling subroutine (for warning prints)
      !!                pa_ifsd(nn_nfsd)   : current value of mFSTD
      !!                ptendency(nn_nfsd) : required tendency of mFSTD (units: s-1)
      !!                pt_elapsed         : total time elapsed from previous iterations (units: s)
      !!                ksubt              : number of sub-time steps used so far
      !!
      !! ** Output  :   pa_ifsd(nn_nfsd)   : updated (incremented by sub-time step * tendency)
      !!                pt_elapsed         : updated (increased   by sub-time step)
      !!                ksubt              : updated (increased   by 1)
      !!
      !! ** References
      !!    ----------
      !!    Horvat, C. & Tziperman, E. (2017).
      !!              The evolution of scaling laws in the sea ice floe size distribution.
      !!              Journal of Geophysical Research: Oceans, 122(9), 7630-7650.
      !!-------------------------------------------------------------------
      !
      CHARACTER(len=*)            , INTENT(in)    ::   cdcrn        ! calling routine name
      REAL(wp), DIMENSION(nn_nfsd), INTENT(inout) ::   pa_ifsd      ! current mFSTD
      REAL(wp), DIMENSION(nn_nfsd), INTENT(in)    ::   ptendency    ! mFSTD tendency (units: s-1)
      REAL(wp)                    , INTENT(inout) ::   pt_elapsed   ! time elapsed from previous iterations (units: s)
      INTEGER                     , INTENT(inout) ::   ksubt        ! number of adaptive time steps used
      !
      CHARACTER(len=3)                            ::   cl_warn      ! for warning print
      REAL(wp), DIMENSION(nn_nfsd)                ::   zdt_restr    ! time step restrictions (units: s)
      REAL(wp)                                    ::   zt_remain    ! remaining time to evolve (units: s)
      REAL(wp)                                    ::   zdt          ! largest allowed time step (units: s)
      INTEGER                                     ::   jf           ! dummy loop index
      !
      INTEGER, PARAMETER :: isubt_warn = 100   ! number of iterations at which warning raised
      !
      !!-------------------------------------------------------------------

      ! Determine time remaining to evolve (global time step minus time already evolved)
      ! (note: should never be here if pt_elapsed > rDt_ice, so zt_remain > 0):
      zt_remain = rDt_ice - pt_elapsed

      ! Calculate maximum possible time step in each floe category and save to zdt_restr
      !
      ! Afterwards we use MINVAL to select maximum allowed time step, but it cannot be larger
      ! than the remaining time. So, we safely use that as the initial/default value:
      zdt_restr(:) = zt_remain
      !
      DO jf = 1, nn_nfsd
         IF( ptendency(jf) >  epsi10 ) zdt_restr(jf) = (1._wp - pa_ifsd(jf)) /     ptendency(jf)
         IF( ptendency(jf) < -epsi10 ) zdt_restr(jf) =          pa_ifsd(jf)  / ABS(ptendency(jf))
      ENDDO

      ! Time step = most restricting value; note we still need to check against time remaining
      ! despite default value above, as ALL(zdt_restr(:) > rDt_ice - pt_elapsed) is possible:
      zdt = MIN(zt_remain, MINVAL(zdt_restr))

      IF( (zdt / zt_remain ) < epsi10 ) THEN   ! stop if adaptive time step is too small
         CALL ctl_stop('STOP', '',                                                             &
            &  '   FSD tendency has become unstable during subroutine: '//TRIM(cdcrn),         &
            &  '   (suggestion: reducing width of floe size categories may overcome the issue)')
      ENDIF

      ! Update FSTD, reduce time increment, and update number of iterations:
      pa_ifsd(:) = pa_ifsd(:) + zdt * ptendency(:)
      pt_elapsed = pt_elapsed + zdt
      ksubt      = ksubt + 1

      IF( ksubt == isubt_warn ) THEN
         WRITE(cl_warn,'(I3)') isubt_warn
         CALL ctl_warn(TRIM(cdcrn)//': reached '//TRIM(cl_warn)//' adaptive sub-time steps')
      ENDIF

   END SUBROUTINE ice_fsd_tstep


   SUBROUTINE ice_fsd_brit
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE ice_fsd_brit  ***
      !!
      !! ** Purpose :   In-plane (brittle) fracture of sea ice
      !!
      !! ** Method  :   Quasi-restoring of the number density floe size distribution (NDFSD)
      !!                towards a theoretical gradient in log(floe size)--log(NDFSD) space whenever
      !!                the actual NDFSD locally exceeds that gradient.
      !!
      !!                The NDFSD is related to the prognostic 'modified-areal' floe size-thickness
      !!                distribution (mFSTD, L(s,h); variable: a_ifsd) for each ice thickness
      !!                category by division of the latter by floe area, which is proportional to
      !!                floe size squared. Thus, the maximum allowed gradient in log-log space of
      !!                the mFSTD is that of the NDFSD plus 2 (see external docs).
      !!
      !!                When restoring is needed, based on the above condition estimated between
      !!                floe size categories j and j-1 using a backward-in-space finite difference,
      !!                area fraction from category j is transferred into category j-1 following a
      !!                (quasi-)exponential decay:
      !!
      !!                   d(ln L)/d(ln s)|_j = -L_j / t_res
      !!
      !!                where t_res is a restoring timescale (namelist rn_fsd_brit_tres) and the
      !!                same tendency (right-hand side) is applied to category j-1 with opposite
      !!                sign. Such tendency contributions are summed to give the net tendency in
      !!                each category and evolved using adaptive time stepping (see ice_fsd_tstep).
      !!
      !!                See docs and Bateson et al. (2022) for theory/motivation.
      !!
      !! ** References
      !!    ----------
      !!    Bateson, A. W., Feltham, D. L., Schroeder, D. L., Wang, Y., Hwang, B., Ridley, J. K., & Aksenov, Y. (2022).
      !!              Sea ice floe size: its impact on pan-Arctic and local ice mass and required model complexity.
      !!              The Cryosphere, 16, 2565-2593.
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(A2D(0),nn_nfsd,jpl) ::   zfstd_b           ! FSTD before brittle fracture (for diagnostics)
      REAL(wp), DIMENSION(nn_nfsd)            ::   ztendency         ! tendency of FSTD
      REAL(wp), DIMENSION(nn_nfsd)            ::   zloss, zgain      ! loss and gain terms to compute tendency
      REAL(wp)                                ::   zlogfsd_grad      ! forward-in-space gradient of FSD in log-log space
      REAL(wp)                                ::   zfsd_res          ! correction term for area conservation
      REAL(wp)                                ::   zt_elapsed        ! time elapsed during adaptive time stepping (units: s)
      INTEGER                                 ::   isubt             ! number of iterations used in adaptive time stepping
      INTEGER                                 ::   ji, jj, jl, jf    ! dummy loop indices
      !
      !!-------------------------------------------------------------------

      zfstd_b(A2D(0),:,:) = a_ifsd(A2D(0),:,:)   ! FSTD before brittle fracture (for diagnostics)

      DO jl = 1, jpl
         DO_2D( 0, 0, 0, 0 )
            IF( (a_i(ji,jj,jl) > epsi10) .AND. (ALL(a_ifsd(ji,jj,:,jl) > epsi10)) ) THEN

               ! Start adaptive time stepping:
               zt_elapsed = 0._wp   ! time elapsed during adaptive time stepping
               isubt      = 0       ! number of sub time steps taken
               !
               DO WHILE (zt_elapsed < rDt_ice)

                  zloss(:) = 0._wp  ! initialise or reset
                  zgain(:) = 0._wp

                  DO jf = 2, nn_nfsd
                     ! Backward-in-(log)-space gradient (denominator pre-computed in fsd_initbounds):
                     zlogfsd_grad = ( LOG(a_ifsd(ji,jj,jf,jl)) - LOG(a_ifsd(ji,jj,jf-1,jl)) )   &
                        &           / floe_dlog_sc(jf)
                     !
                     ! If gradient is too large here, restore toward maximum allowed gradient by
                     ! transferring area fraction to the smaller category as exponential decay
                     !
                     ! Offset of 2 is to transform from area-density FSD (a_ifsd)
                     !                             to number-density FSD (gradient rn_fsd_brit_grad)
                     IF ( zlogfsd_grad > rn_fsd_brit_grad + 2._wp ) THEN
                        zloss(jf)   = zloss(jf)   + a_ifsd(ji,jj,jf,jl)   ! divide by timescale after
                        zgain(jf-1) = zgain(jf-1) + a_ifsd(ji,jj,jf,jl)
                     ENDIF
                  ENDDO

                  ! Net tendency per category:
                  ztendency(:) = (zgain(:) - zloss(:)) / rn_fsd_brit_tres

                  ! Evolve a_ifsd over maximum stable time step and increase zt_elapsed accordingly:
                  CALL ice_fsd_tstep('ice_fsd_brit', a_ifsd(ji,jj,:,jl), ztendency(:), zt_elapsed, isubt)

                  ! Break adaptive time stepping loop if all ice now in smallest category
                  ! (=> all possible fracture occurred):
                  IF( a_ifsd(ji,jj,1,jl) > (1._wp - epsi10)) EXIT

               ENDDO

               ! Brittle fracture physically cannot/should not directly lead to loss of ice area
               ! So, correct for any numerical residual by adding it back to the smallest floe size
               ! category (if some area is lost) or by taking it away from the largest category
               ! that has at least that residual available (if some area has been gained).
               !
               zfsd_res = SUM(a_ifsd(ji,jj,:,jl)) - 1._wp
               !
               IF( zfsd_res <= 0._wp ) THEN   ! area lost
                  a_ifsd(ji,jj,1,jl) = a_ifsd(ji,jj,1,jl) + ABS(zfsd_res)
               ELSE   ! area gained
                  DO jf = nn_nfsd, 1, -1
                     IF( a_ifsd(ji,jj,jf,jl) > zfsd_res) THEN
                        a_ifsd(ji,jj,jf,jl) = a_ifsd(ji,jj,jf,jl) - zfsd_res
                        EXIT
                     ENDIF
                  ENDDO
               ENDIF
               !
               ! Now ensure normalisation and [0-1] bounding:
               CALL ice_fsd_cor( a_ifsd(ji,jj,:,jl) )
               !
            ENDIF
         END_2D
      ENDDO

      ! Write FSD tendency diagnostics due to brittle fracture:
      CALL ice_fsd_dia( 'bfr', zfstd_b, a_ifsd(A2D(0),:,:), a_i(A2D(0),:), a_i(A2D(0),:) )

   END SUBROUTINE ice_fsd_brit


   SUBROUTINE ice_fsd_part_newice( pa_i, pv_i, pa_ifstd, pv_newice, pv_latgro, pda_latgro )
      !!-------------------------------------------------------------------
      !!            ***  ROUTINE ice_fsd_part_newice  ***
      !!
      !! ** Purpose :   Partition total new ice volume into new ice formation
      !!                in open water and lateral growth of existing ice
      !!
      !! ** Method  :   Calculate lead area and lateral surface area of floes
      !!                following Horvat and Tziperman (2015):
      !!
      !!                A_lead = int [ f(s,h) * (4*r_lw/s + 4*(r_lw/s)**2) ] ds dh
      !!                A_lat  = int [ f(s,h) * pi * h / (a_shape * s) ] ds dh
      !!
      !!                where A_lead = lead area (per unit ocean area)
      !!                      A_lat  = total area of vertical edges of floes
      !!                               (per unit ocean area)
      !!                      int    = integral over floe size s and thickness h
      !!                      f(s,h) = g(h)L(s,h) floe size-thickness distribution
      !!                      a_shape= floe shape parameter
      !!
      !!                The lead width r_lw is the annulus surrounding floes for
      !!                which freezing of existing floes occurs, distinguished
      !!                from 'open water area' in which new ice forms away from
      !!                existing ice. Hence, 'lead area' plus 'open water area'
      !!                equals one minus sea ice concentration. See Horvat and
      !!                Tziperman (2015), sect. 2.1 and Fig. 1 for details.
      !!                Following Roach et al. (2018), the smallest floe size
      !!                is used for r_lw.
      !!
      !!                The lateral growth volume (in the 'lead area') is then:
      !!
      !!                v_latgro = v_newice * A_lead / [ 1 + (at_i / A_lat) ]
      !!
      !!                where v_newice = total new ice growth
      !!                      at_i     = sea ice concentration
      !!
      !!                v_newice, initially calculated in ice_thd_do from which
      !!                this routine is called, is then updated by subtracting
      !!                v_latgro. In ice_thd_do, it is then used to add new ice
      !!                in open water (as usual since there is no lateral growth
      !!                by default, i.e., without FSD) and the total new ice
      !!                volume added, now (v_newice + v_latgro), is unchanged.
      !!
      !! ** Input   :   pa_i, pv_i      : local ice concentration and volume (per category)
      !!                pa_ifstd        : local floe size-thickness distribution
      !!                pv_newice       : volume of new ice to grow in total as
      !!                                  calculated in ice_thd_do, before
      !!                                  accounting for FSD (ln_fsd) and/or
      !!                                  frazil ice collection (ln_frazil)
      !!
      !! ** Output  :   pv_newice       : input updated by subtracting pv_latgro
      !!                                  (so it now represents new ice volume
      !!                                  grown in open water).
      !!                pv_latgro       : volume of ice to grow laterally on
      !!                                  existing floes in the lead area
      !!                pda_latgro(jpl) : a_i change due to lateral growth of
      !!                                  existing floes in each thickness cat.
      !!
      !! ** Note    :   no updates to ice concentration, volume, or FSD
      !!                prognostic variables are made in this routine. That is
      !!                done in the routines: ice_thd_do, ice_fsd_thd, and ice_fsd_add_newice.
      !!
      !! ** References
      !!    ----------
      !!    Horvat, C., & Tziperman, E. (2015).
      !!              A prognostic model of the sea-ice floe size and thickness distribution.
      !!              The Cryosphere, 9, 2119-2134.
      !!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018).
      !!              An emergent sea ice floe size distribution in a global coupled ocean-sea ice model
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(jpl)        , INTENT(in)    ::   pa_i         ! local ice concentration (per category)
      REAL(wp), DIMENSION(jpl)        , INTENT(in)    ::   pv_i         ! local ice volume (per category)
      REAL(wp), DIMENSION(nn_nfsd,jpl), INTENT(in)    ::   pa_ifstd     ! local floe size-thickness distribution
      REAL(wp)                        , INTENT(inout) ::   pv_newice    ! local total new ice volume (from ice_thd_do)
      REAL(wp)                        , INTENT(out)   ::   pv_latgro    ! lateral growth volume
      REAL(wp), DIMENSION(jpl)        , INTENT(out)   ::   pda_latgro   ! a_i change due to lateral growth
      !
      INTEGER  ::   jl, jf        ! dummy loop indices
      REAL(wp) ::   za_lead       ! lead area for open water growth (per unit ocean area)
      REAL(wp) ::   za_lat_surf   ! lateral surface area of floes (per unit ocean area)
      REAL(wp) ::   zat_i         ! total ice concentration in grid cell
      REAL(wp) ::   zh_i          ! ice thickness (m)
      REAL(wp) ::   zr_lw         ! width of lead region (m)
      !
      !!-------------------------------------------------------------------

      pv_latgro     = 0._wp   ! initialise
      pda_latgro(:) = 0._wp
      za_lead       = 0._wp
      za_lat_surf   = 0._wp
      zat_i         = SUM(pa_i)
      zr_lw         = floe_sc(1)   ! smallest floe size for width of lead region

      ! --- Calculate za_lead and za_lat_surf (integrate/sum over thickness
      !     and floe size cats.):
      DO jl = 1, jpl

         ! need ice thickness (m) for za_lat_surf:
         IF ( pa_i(jl) > 0._wp ) THEN
            zh_i = pv_i(jl) / pa_i(jl)
         ELSE
            zh_i = 0._wp
         ENDIF

         DO jf = 1, nn_nfsd
            !
            za_lead = za_lead + pa_i(jl) * pa_ifstd(jf,jl)   &
               &                         * 4._wp * (zr_lw / floe_sc(jf) + zr_lw**2 / floe_sc(jf)**2)
            !
            za_lat_surf = za_lat_surf + pa_ifstd(jf,jl) * pa_i(jl)   &
               &                        * rpi * zh_i / (rn_floeshape * floe_sc(jf))
            !
         ENDDO
      ENDDO

      ! --- Lead area cannot exceed open water fraction and must be > 0:
      za_lead = MAX( 0._wp, MIN( za_lead, 1._wp - zat_i ) )

      ! --- Calculate lateral growth volume, zv_latgro, and the change in a_i
      !     due to lateral growth, or leave both as 0 if no lateral growth:
      IF (za_lat_surf > epsi10) THEN

         pv_latgro = pv_newice * za_lead / (1._wp + zat_i / za_lat_surf)

         DO jl = 1, jpl
            DO jf = 1, nn_nfsd

               ! note lateral growth rate = zv_latgro / rDt_ice, but here we
               ! calculate the growth over time step which is just zv_latgro:
               pda_latgro(jl) = pda_latgro(jl) + 4._wp * pa_i(jl) * pa_ifstd(jf,jl)   &
                  &                                    * pv_latgro / floe_sc(jf)

            ENDDO
         ENDDO

         IF ( SUM(pda_latgro) >= za_lead ) THEN
            ! --- Cannot expand ice laterally beyond the lead region
            !     so normalise net lateral area growth to equal lead area:
            pda_latgro(:) = pda_latgro(:) / SUM(pda_latgro)
            pda_latgro(:) = pda_latgro(:) * za_lead
         ENDIF

      ENDIF

      ! --- Update volume of new ice to grow in open water:
      pv_newice = pv_newice - pv_latgro

   END SUBROUTINE ice_fsd_part_newice


   SUBROUTINE ice_fsd_add_newice( pa_ifsd, pa_newice, pa_i_before, kcat )
      !!-------------------------------------------------------------------
      !!                ***  ROUTINE ice_fsd_add_newice  ***
      !!
      !! ** Purpose :   Add new ice growth in open water (not lateral growth
      !!                of existing ice) to the floe size distribution in the
      !!                appropriate floe size category.
      !!
      !! ** Method  :   New ice is added to the smallest floe size category.
      !!
      !! ** Input   :   pa_ifsd(nn_nfsd) : floe size distribution at one grid
      !!                                   point and for one thickness category
      !!                pa_newice        : area fraction of new ice formation
      !!                pa_i_before      : ice concentration *after* lateral growth
      !!                                   but *before* new ice growth at 1-D array
      !!                kcat             : floe size category index to add new ice to
      !!
      !! ** Note    :   This routine only updates the floe size distribution,
      !!                not ice concentration a_i, which is done in ice_thd_do
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(nn_nfsd), INTENT(inout) ::   pa_ifsd       ! FSD at one location, one thickness cat.
      REAL(wp)                    , INTENT(in)    ::   pa_newice     ! area fraction of new ice
      REAL(wp)                    , INTENT(in)    ::   pa_i_before   ! a_i after lat. growth of
      !                                                              ! existing ice but before addition
      !                                                              ! of pa_newice
      INTEGER                     , INTENT(in)    ::   kcat          ! FSD category index for new ice
      !
      INTEGER ::   jf   ! dummy loop index
      !
      !!-------------------------------------------------------------------

      IF( pa_newice > 0._wp ) THEN
         IF( SUM(pa_ifsd(:)) > epsi10 ) THEN
            !
            ! --- Add new ice to specified floe size category
            !
            ! The area fraction of ice in this floe size category, sk,
            ! and thickness category to which new ice is added, h, is:
            !
            !    [ L(sk,h)g(h)dsdh ]_before = pa_ifsd(sk) * pa_i_before
            !
            ! before addition of pa_newice. Then, after addition of new ice:
            !
            !    [ L(sk,h)g(h)dsdh ]_after = [ L(sk,h)g(h)dsdh ]_before + pa_newice
            !
            ! g(h) is already updated in ice_thd_do, but L(sk,h) needs updating
            ! too, achieved by rearranging the above. This is why it is necessary
            ! to pass pa_i_before to this routine rather than just using a_i_2d.
            !
            pa_ifsd(kcat) = (pa_ifsd(kcat)*pa_i_before + pa_newice) / (pa_i_before + pa_newice)

            ! --- Adjust other floe size categories
            !
            ! New ice area is only added to one floe size category, sk.
            ! So for the remaining floe size categories, s:
            !
            !    [ L(s,h)g(h)dsdh ]_after = [ L(s,h)g(h)dsdh ]_before
            !
            ! Since g(h)_before /= g(h)_after, L(s,h)_before /= L(s,h)_after.
            ! Rearranging gives L(s,h)_after and is thus updated:
            !
            DO jf = 1, nn_nfsd
               IF( jf /= kcat ) pa_ifsd(jf) = pa_ifsd(jf)*pa_i_before / (pa_i_before + pa_newice)
            ENDDO

         ELSE
            !
            ! --- Entirely new ice: put in specified floe size category:
            pa_ifsd(:)  = 0._wp
            pa_ifsd(kcat) = 1._wp

         ENDIF
      ENDIF

      CALL ice_fsd_cor( pa_ifsd )   ! small/negative value corrections, re-normalisation

   END SUBROUTINE ice_fsd_add_newice


   SUBROUTINE ice_fsd_weld( pa_ifsd, pa_i )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_fsd_weld  ***
      !!
      !! ** Purpose :   Evolve the floe size distribution subject to the
      !!                welding together of floes in freezing conditions
      !!
      !! ** Method  :   Floes are assumed to be placed randomly on the domain
      !!                (grid cell) and the probability of two floes overlapping
      !!                is described using a coagulation equation:
      !!
      !!                dN(x)/dt = 0.5*int[ K(x',x-x') dx'] - int[ K(x,x') dx' ]
      !!
      !!                   x    = floe area [m2]
      !!                   N(x) = number density of floes of area x [m-4]
      !!
      !!                K(x1,x2) = c_weld * x1 * x2 * N(x1) * N(x2)
      !!
      !!                   K      = coagulation kernel, the number of floe merging
      !!                            events per unit area of ocean, per unit x1, per
      !!                            unit x2, per unit time [m-6.s-1]
      !!                   c_weld = scale factor for welding. Can be interpreted as
      !!                            the total number of floes that weld with another
      !!                            per unit area of ocean per unit time, in the case
      !!                            of a fully ice-covered ocean [m-2.s-1]
      !!
      !!                See Roach et al. (2018a,b) for details of theory. The equation
      !!                is here solved in terms of area floe size distribution [rather
      !!                than N, using xN(x)dx = f(r)dr] and evolved using adaptive
      !!                time stepping.
      !!
      !! ** Input   :   pa_ifsd(nn_nfsd) : floe size distribution at one grid
      !!                                   point and for one thickness category
      !!                pa_i             : sea ice concentration in same thickness cat.
      !!
      !! ** Note    :   This routine does not check for local freezing conditions. It
      !!                does check input sea ice concentration is above a minimum
      !!                threshold set by namelist parameter rn_fsd_amin_weld. Welding is
      !!                considered unlikely below this threshold and in such cases this
      !!                routine does nothing.
      !!
      !!                The coefficient c_weld is set by namelist rn_fsd_c_weld, which
      !!                can be considered a tuning parameter.
      !!
      !!                This routine does not modify any other state variables.
      !!
      !! ** References
      !!    ----------
      !!    Roach, L. A., Smith, M. M., & Dean, S. M. (2018a).
      !!              Quantifying growth of pancake sea ice floes using images from drifting buoys
      !!              Journal of Geophysical Research: Oceans, 123(4), 2851-2866.
      !!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018b).
      !!              An emergent sea ice floe size distribution in a global coupled ocean-sea ice model
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(nn_nfsd), INTENT(inout) ::   pa_ifsd   ! FSD at one location, one thickness cat.
      REAL(wp)                    , INTENT(in)    ::   pa_i      ! ice conc. at one location, one thickness cat.
      !
      REAL(wp), DIMENSION(nn_nfsd) ::   zloss, zgain      ! exchange tendencies between FSD categories [1/s]
      REAL(wp)                     ::   zdfsd_weld        ! change in FSD due to a welding interaction [1/s]
      REAL(wp)                     ::   zt_elapsed        ! time elapsed during adaptive time stepping (units: s)
      INTEGER                      ::   isubt             ! number of iterations used in adaptive time stepping
      INTEGER                      ::   jf1, jf2, jf3     ! dummy loop indices
      !
      !!-------------------------------------------------------------------

      ! --- Additional conditions for floe welding (freezing conditions assumed):
      !        (1) ice concentration exceeds threshold (welding is unlikely
      !            with low sea ice concentrations)
      !        (2) must be some ice to weld in the first place (i.e., some
      !            ice in lower floe size categories)
      !
      IF( (pa_i > rn_fsd_amin_weld) .and. (SUM(pa_ifsd(1:nn_nfsd-1)) > epsi10) ) THEN

         ! --- Start adaptive time stepping
         zt_elapsed = 0._wp   ! time elapsed during adaptive time stepping
         isubt      = 0       ! number of sub time steps taken

         DO WHILE (zt_elapsed < rDt_ice)

            ! --- Calculate loss and gain rates of fractional area of floes
            !     in each floe size category due to welding
            zloss(:) = 0._wp
            zgain(:) = 0._wp   ! initialise

            DO jf1 = 1, nn_nfsd
               !
               ! --- This loop corresponds to calculation of loss in N(x)
               !     (here, FSD) for each category jf1, i.e., the second term
               !     of the dN/dt equation. Those losses are also counted as
               !     gains in other categories jf2 in next loop below. So, the
               !     gains in category jf1 are calculated indirectly by other
               !     iterations of this loop.
               !
               DO jf2 = 1, nn_nfsd
                  !
                  ! --- This loop corresponds to integral in coagulation equation,
                  !     i.e., considering interactions of floes in category jf1
                  !     with all other categories (jf2).
                  !
                  !     Calculate the loss from category jf1 and add it to
                  !     zloss(jf1), and add the same to the gain of whichever
                  !     category welded floes belong to (jf3).
                  !
                  !     Note corresponding loss from category jf2 is accounted for
                  !     when jf1 and jf2 are exchanged (i.e., outer loop).
                  !
                  !     If there can be no such welding, jf3 = 0 which is the
                  !     'missing value' in floe_iweld --> nothing happens.
                  !
                  !     Note lack of factor of 0.5 in equation because we just
                  !     calculate the losses/gains in area fraction directly, i.e.,
                  !     without explicitly calculating each of the two terms on the
                  !     right-hand side of the equation.
                  !
                  jf3 = floe_iweld(jf1,jf2)
                  !
                  IF( jf3 > jf1 ) THEN
                     zdfsd_weld = rn_fsd_c_weld * floe_ac(jf1) * pa_i * pa_ifsd(jf1) * pa_ifsd(jf2)
                     zloss(jf1) = zloss(jf1) + zdfsd_weld
                     zgain(jf3) = zgain(jf3) + zdfsd_weld
                  ENDIF
                  !
               ENDDO
            ENDDO

            ! Evolve pa_ifsd over maximum stable time step and increase zt_elapsed accordingly:
            CALL ice_fsd_tstep('ice_fsd_weld', pa_ifsd(:), zgain(:) - zloss(:), zt_elapsed, isubt)

            CALL ice_fsd_cor( pa_ifsd )   ! small/negative value corrections, re-normalisation

            ! --- Break adaptive time stepping loop if all ice is in
            !     the largest floe category (since all possible welding
            !     has occurred)
            IF( pa_ifsd(nn_nfsd) > (1._wp - epsi10)) EXIT

         ENDDO

         CALL ice_fsd_cor( pa_ifsd )   ! small/negative value corrections, re-normalisation

      ENDIF

   END SUBROUTINE ice_fsd_weld


   SUBROUTINE ice_fsd_thd( pa_ifsd, pG_s )
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_fsd_thd  ***
      !!
      !! ** Purpose :   Evolve the floe size distribution subject to lateral
      !!                growth/melt
      !!
      !! ** Method  :   dL(s,h)/dt = -G_s * div_s(L) + (2/s) * G_s * L(s,h)
      !!
      !!                where L(s,h) is the floe size (s) distribution at
      !!                             thickness h
      !!                      div_s  is divergence in s-space
      !!                      G_s    is the lateral growth/melt rate ds/dt, assumed
      !!                             to be independent of s and h, and G_s > 0
      !!                             implies growth
      !!
      !!                This equation is derived by Horvat and Tziperman (2015)
      !!                and adapted to the modified-areal floe size distribution
      !!                L(s,h) by Roach et al. (2018). This routine integrates
      !!                it forwards (for one thickness category) by one model
      !!                time step using adaptive time stepping (Horvat and
      !!                Tziperman, 2017; see ice_fsd_tstep).
      !!
      !! ** Input   :   pa_ifsd(nn_nfsd) : floe size distribution at one grid
      !!                                   point and for one thickness category
      !!                pG_s             : lateral growth/melt rate in m/s. Specifically ds/dt;
      !!                                   important as Horvat and Tziperman (2015) use 'radius'
      !!                                   whereas we have diameter for floe size (ds/dt = 2dr/dt).
      !!
      !! ** Note    :   This routine does not implement creation of new ice area or
      !!                loss of ice area due to complete melt of existing floes. Those
      !!                do affect the FSD but are handled in the separate routines
      !!                ice_fsd_add_newice and (indirectly) ice_thd_da.
      !!
      !! ** References
      !!    ----------
      !!    Horvat, C., & Tziperman, E. (2015).
      !!              A prognostic model of the sea-ice floe size and thickness distribution.
      !!              The Cryosphere, 9, 2119-2134.
      !!    Horvat, C., & Tziperman, E. (2017).
      !!              The evolution of scaling laws in the sea ice floe size distribution.
      !!              Journal of Geophysical Research: Oceans, 122(9), 7630-7650.
      !!    Roach, L. A., Horvat, C., Dean, S. M., & Bitz, C. M. (2018).
      !!              An emergent sea ice floe size distribution in a global coupled ocean-sea ice model
      !!              Journal of Geophysical Research: Oceans, 123(6), 4322-4337.
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(nn_nfsd), INTENT(inout) ::   pa_ifsd   ! FSD at one location, one thickness cat.
      REAL(wp),                     INTENT(in)    ::   pG_s      ! lateral growth/melt rate (ds/dt; m/s)
      !
      REAL(wp), DIMENSION(nn_nfsd) ::   za_ifsd_tend      ! FSD tendency (left side of eq. above)
      REAL(wp), DIMENSION(nn_nfsd) ::   zdiv_fsd          ! divergence term in equation
      REAL(wp)                     ::   zfsd_cor          ! correction factor to ensure area conservation
      REAL(wp)                     ::   zt_elapsed        ! time elapsed during adaptive time stepping (units: s)
      INTEGER                      ::   isubt             ! number of iterations used in adaptive time stepping
      INTEGER                      ::   jf                ! dummy loop index
      CHARACTER(len=1)             ::   cln               ! string for warning print to indicate ice_thd_da vs ice_thd_do
      !
      !!-------------------------------------------------------------------

      ! --- Start adaptive time stepping
      zt_elapsed = 0._wp   ! time elapsed during adaptive time stepping
      isubt      = 0       ! number of sub time steps taken

      DO WHILE (zt_elapsed < rDt_ice)

         za_ifsd_tend(:) = 0._wp   ! initialise (or reset with loop iteration)
         zdiv_fsd    (:) = 0._wp

         ! --- Calculate the divergence term [div(L), without the -G_s factor]
         !     of the FSD thermodynamic tendency equation using the divergence
         !     theorem.
         !
         ! The divergence in floe category jf equals the net 'flux' of floes
         ! out of that category. Since only growth or melt occurs at once, the
         ! array indices are different depending on the sign of pG_s. In
         ! growth, floes move from smaller to larger floe size categories only,
         ! so the 'flux' of floes from category jf is directed into category
         ! jf+1, while for melt it is into category jf-1.
         !
         IF( pG_s > 0._wp ) THEN   ! lateral growth

            DO jf = 2, nn_nfsd-1
               zdiv_fsd(jf) = (   (pa_ifsd(jf  ) / floe_ds(jf  ) )     &
                  &             - (pa_ifsd(jf-1) / floe_ds(jf-1) ) )
            ENDDO

            ! Smallest category: no 'floe flux' from smaller category:
            zdiv_fsd(1) = pa_ifsd(1) / floe_ds(1)

            ! Largest category: no 'floe flux' leaving this category:
            zdiv_fsd(nn_nfsd) = -pa_ifsd(nn_nfsd-1) / floe_ds(nn_nfsd-1)

            cln = 'o'   ! for warning print, to indicate ice_thd_do is calling

         ELSE   ! pG_r < 0; lateral melt
            !
            ! Note 'flux' of floes from category jf+1 goes into category jf
            ! which is a convergence in category jf, so need a minus sign
            ! for that term. But that minus sign is already provided by
            ! pG_s < 0, so need to add an extra minus sign to this and
            ! other 'flux' terms throughout so it cancels out later.
            !
            ! ToDo: may be clearer to write these fluxes with zG_s
            !       included in both growth and melt cases?
            !
            DO jf = 2, nn_nfsd-1
               zdiv_fsd(jf) = (   (pa_ifsd(jf+1) / floe_ds(jf+1) )     &
                  &             - (pa_ifsd(jf  ) / floe_ds(jf  ) ) )
            ENDDO

            ! Smallest category: there is a 'floe flux' leaving this category,
            ! but it represents complete melt of smallest floes and results in
            ! ice area loss. So that flux, which would be pa_ifsd(1) / floe_ds(1),
            ! is not here because this routine is just shifting ice between floe
            ! size categories, but the term appears directly in routine ice_thd_da.
            !
            ! Meanwhile, here we just have the 'floe flux' from category 2:
            !
            zdiv_fsd(1) = pa_ifsd(2) / floe_ds(2)

            ! Largest category: no 'floe flux' from larger category:
            zdiv_fsd(nn_nfsd) = -pa_ifsd(nn_nfsd) / floe_ds(nn_nfsd)

            cln = 'a'   ! for warning print, to indicate ice_thd_da is calling

         ENDIF

         ! --- Correction term
         !
         ! Sum over all floe size categories of the tendency equation must (in
         ! theory) be zero, because int(L dr) = 1 by definition, and so
         ! d/dt( int(L ds)) = 0. The divergence term also integrates to zero:
         ! indeed all elements of zdiv_fsd computed above cancel out when summed.
         ! Therefore, second term on RHS should sum to zero. In case of noise,
         ! which would manifest as spurious ice area, compute its integral,
         ! zfsd_cor, and subtract it from the actual tendency in each category
         ! weighted by that category's area fraction.
         !
         zfsd_cor = 2._wp * pG_s * SUM( pa_ifsd(:) / floe_sc(:) )

         ! --- Compute rate of change of FSD in each floe size category:
         DO jf = 1, nn_nfsd
            za_ifsd_tend(jf) = -pG_s * zdiv_fsd(jf)                                   &
               &               + 2._wp * pG_s * pa_ifsd(jf) * (1._wp / floe_sc(jf))   &
               &               - pa_ifsd(jf) * zfsd_cor
         ENDDO

         ! Evolve pa_ifsd over maximum stable time step and increase zt_elapsed accordingly:
         CALL ice_fsd_tstep( 'ice_thd_d'//cln//' -> ice_fsd_thd',           &
            &                pa_ifsd(:), za_ifsd_tend(:), zt_elapsed, isubt )

      ENDDO

      CALL ice_fsd_cor( pa_ifsd(:) )   ! small/negative value corrections, re-normalisation

   END SUBROUTINE ice_fsd_thd


   SUBROUTINE ice_fsd_wri( kt )
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE ice_fsd_wri  ***
      !!
      !! ** Purpose :   Writes output fields related to the FSD.
      !!
      !! ** Method  :   Calculates metrics if requested for output and
      !!                writes using iom routines.
      !!
      !!-------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt                    ! ocean time step index
      !
      REAL(wp), DIMENSION(A2D(0))     ::   zsavg       ! mean floe size, grid cell (m)
      REAL(wp), DIMENSION(A2D(0))     ::   zperi       ! perimeter density, grid cell (m.m-2)
      REAL(wp), DIMENSION(A2D(0))     ::   zseff       ! effective floe size, grid cell (m)
      REAL(wp), DIMENSION(A2D(0))     ::   zmsk00      ! 0% conc. mask, grid cell
      REAL(wp), DIMENSION(A2D(0))     ::   zmsk1_ati   ! mask = 1/at_i (ice) or 0 (no ice)
      REAL(wp), DIMENSION(A2D(0),jpl) ::   zperi_cat   ! perimeter density, each ITD category (m.m-2)
      REAL(wp), DIMENSION(A2D(0),jpl) ::   zseff_cat   ! effective floe size, each ITD category (m)
      REAL(wp), DIMENSION(A2D(0),jpl) ::   zmsk00c     ! 0% conc. mask, each ITD category
      !
      REAL(wp), DIMENSION(A2D(0),nn_nfsd,jpl) :: zmsk00fc         ! 0% conc. mask, each ITD and FSD category
      REAL(wp), DIMENSION(A2D(0),nn_nfsd)     :: zmsk00f          ! 0% conc. mask, each FSD category
      REAL(wp), DIMENSION(A2D(0),nn_nfsd)     :: zfsd             ! FSD integrated over ITD
      REAL(wp), DIMENSION(A2D(0),nn_nfsd)     :: zpdd             ! Perimeter density distribution
      INTEGER                                 :: ji, jj, jl, jf   ! dummy loop indices
      !
      !!-------------------------------------------------------------------

      ! --- Calculate sea ice threshold masks for outputs (as in subroutine ice_wri)
      zmsk00 (:,:)   = MERGE( 1._wp, 0._wp, at_i(A2D(0))  >= epsi06  )
      zmsk00c(:,:,:) = MERGE( 1._wp, 0._wp, a_i(A2D(0),:) >= epsi06  )

      ! --- Analogous masks including FSD dimension
      DO jf = 1, nn_nfsd
         zmsk00f (:,:,jf)   = MERGE( 1._wp, 0._wp, at_i(A2D(0))  >= epsi06 )
         zmsk00fc(:,:,jf,:) = MERGE( 1._wp, 0._wp, a_i(A2D(0),:) >= epsi06 )
      ENDDO

      ! --- Calculate new mask = 1/at_i or 0 if at_i too small:
      zmsk1_ati(A2D(0)) = 0._wp
      WHERE( at_i(A2D(0)) >= epsi06 ) zmsk1_ati(A2D(0)) = 1._wp / at_i(A2D(0))

      ! --- Calculate outputs
      !
      zseff_cat(A2D(0),:) = 0._wp   ! initialise
      zseff    (A2D(0))   = 0._wp
      zsavg    (A2D(0))   = 0._wp
      !
      DO_2D(0, 0, 0, 0)
         !
         ! Floe size distribution and perimeter density distributions:
         zfsd(ji,jj,:) = floe_size_dist( a_ifsd(ji,jj,:,:), a_i(ji,jj,:) )
         zpdd(ji,jj,:) = peri_dens_dist( a_ifsd(ji,jj,:,:), a_i(ji,jj,:) )
         !
         ! Area-weighted mean floe size:
         DO jf = 1, nn_nfsd
            zsavg(ji,jj) = zsavg(ji,jj) + floe_sc(jf) * zfsd(ji,jj,jf)
         ENDDO
         !
         ! Perimeter density and effective floe size, per ITD category:
         DO jl = 1, jpl
            zperi_cat(ji,jj,jl) = fsd_peri_dens( a_ifsd(ji,jj,:,jl) )
            zseff_cat(ji,jj,jl) = fsd_eff_size(  a_ifsd(ji,jj,:,jl) )
         ENDDO
         !
         ! Perimeter density and effective floe size, for all ice:
         zperi(ji,jj) = fsd_peri_dens( zfsd(ji,jj,:) )
         zseff(ji,jj) = fsd_eff_size(  zfsd(ji,jj,:) )
         !
      END_2D

      ! --- Write constant fields to output (if requested, case-by-case)
      IF(iom_use( 'icefsd_sl' )) CALL iom_put( 'icefsd_sl' , floe_sl(:) )
      IF(iom_use( 'icefsd_sc' )) CALL iom_put( 'icefsd_sc' , floe_sc(:) )
      IF(iom_use( 'icefsd_su' )) CALL iom_put( 'icefsd_su' , floe_su(:) )
      IF(iom_use( 'icefsd_al' )) CALL iom_put( 'icefsd_al' , floe_al(:) )
      IF(iom_use( 'icefsd_ac' )) CALL iom_put( 'icefsd_ac' , floe_ac(:) )
      IF(iom_use( 'icefsd_au' )) CALL iom_put( 'icefsd_au' , floe_au(:) )
      IF(iom_use( 'icefsd_ds' )) CALL iom_put( 'icefsd_ds' , floe_ds(:) )

      ! --- Write variable fields to output (if requested, case-by-case)
      IF(iom_use( 'icefsd_cat'     )) CALL iom_put( 'icefsd_cat'     , a_ifsd   (A2D(0),:,:) * zmsk00fc )
      IF(iom_use( 'icefsd'         )) CALL iom_put( 'icefsd'         , zfsd     (A2D(0),:)   * zmsk00f  )
      IF(iom_use( 'icepdd'         )) CALL iom_put( 'icepdd'         , zpdd     (A2D(0),:)   * zmsk00f  )
      IF(iom_use( 'icefsdperi_cat' )) CALL iom_put( 'icefsdperi_cat' , zperi_cat(A2D(0),:)   * zmsk00c  )
      IF(iom_use( 'icefsdseff_cat' )) CALL iom_put( 'icefsdseff_cat' , zseff_cat(A2D(0),:)   * zmsk00c  )
      IF(iom_use( 'icefsdperi'     )) CALL iom_put( 'icefsdperi'     , zperi    (A2D(0))     * zmsk00   )
      IF(iom_use( 'icefsdseff'     )) CALL iom_put( 'icefsdseff'     , zseff    (A2D(0))     * zmsk00   )
      IF(iom_use( 'icefsdsavg'     )) CALL iom_put( 'icefsdsavg'     , zsavg    (A2D(0))     * zmsk00   )

   END SUBROUTINE ice_fsd_wri


   SUBROUTINE ice_fsd_dia( cd_dia, pa_ifsdb, pa_ifsda, pa_ib, pa_ia )
      !!-------------------------------------------------------------------
      !!                 ***    ROUTINE ice_fsd_dia    ***
      !!
      !! ** Purpose :   Calculate and write FSD tendency diagnostics
      !!
      !! ** Method  :   The change in floe size-thickness distribution, FSTD = a_ifsd*a_i,
      !!                is calculated as (FSTD_a - FSTD_b) / rDt_ice, where '_a' and '_b' refer to
      !!                after and before the process of which the tendency is computed. This is done
      !!                similarly for other FSD-related diagnostics such as mean floe size.
      !!                Diagnostics are sent to IOM as required.
      !!
      !! ** Inputs  :   Length-3 character name of process for diagnostic suffix (e.g., 'lam' for lateral melt)
      !!                Prognostic FSTD and ice concentration (cat.) variables before and after process,
      !!                each on the inner domain only [i.e., send a_ifsd(A2D(0),:,:)].
      !!
      !!-------------------------------------------------------------------
      !
      CHARACTER(len=3)                          , INTENT(in) ::   cd_dia     ! process label (lam, lag, etc.)
      REAL(wp)   , DIMENSION(A2D(0),nn_nfsd,jpl), INTENT(in) ::   pa_ifsdb   ! FSTD before process (inner domain)
      REAL(wp)   , DIMENSION(A2D(0),nn_nfsd,jpl), INTENT(in) ::   pa_ifsda   ! FSTD after process (inner domain)
      REAL(wp)   , DIMENSION(A2D(0),jpl)        , INTENT(in) ::   pa_ib      ! a_i before process (inner domain)
      REAL(wp)   , DIMENSION(A2D(0),jpl)        , INTENT(in) ::   pa_ia      ! a_i after process (inner domain)
      !
      CHARACTER(len=25) ::   cl_ref   ! output field reference (whole name including suffix)
      CHARACTER(len=4)  ::   cl_sfx   ! output field reference (suffix)
      !
      REAL(wp), DIMENSION(A2D(0),nn_nfsd,jpl) ::   zmsk00fc     ! Ice present mask (2D + FSD and ITD dimensions)
      REAL(wp), DIMENSION(A2D(0),nn_nfsd)     ::   zmsk00f      ! Ice present mask (2D + FSD dimension)
      REAL(wp), DIMENSION(A2D(0))             ::   zmsk00       ! Ice present mask (2D)
      REAL(wp), DIMENSION(A2D(0),nn_nfsd,jpl) ::   zdfstd       ! Tendency of FSTD
      REAL(wp), DIMENSION(A2D(0),nn_nfsd)     ::   zdfsd        ! Tendency of FSD (FSTD integrated over ITD)
      REAL(wp), DIMENSION(A2D(0))             ::   zdsavg       ! Tendency of mean floe size (m/s)
      REAL(wp), DIMENSION(A2D(0))             ::   zat_ia       ! Total ice concentration (after process)
      INTEGER                                 ::   ji, jj, jf   ! dummy loop indices
      !
      !!-------------------------------------------------------------------

      zat_ia(:,:) = SUM( pa_ia(:,:,:), DIM=3 )   ! total ice conc. after

      ! --- Calculate sea ice threshold masks for outputs (as in subroutine ice_wri)
      zmsk00 (:,:)   = MERGE( 1._wp, 0._wp, zat_ia(:,:)  >= epsi06 )

      ! --- Analogous masks including FSD dimension
      DO jf = 1, nn_nfsd
         zmsk00f (:,:,jf)   = MERGE( 1._wp, 0._wp, zat_ia(:,:)  >= epsi06 )
         zmsk00fc(:,:,jf,:) = MERGE( 1._wp, 0._wp, pa_ia(:,:,:) >= epsi06 )
      ENDDO

      ! Calculate tendency diagnostics:
      !
      zdsavg(:,:) = 0._wp  ! initialise
      !
      DO_2D(0, 0, 0, 0)
         !
         ! Full FSTD tendency:
         zdfstd(ji,jj,:,:) = r1_Dt_ice * ( pa_ifsda(ji,jj,:,:) - pa_ifsdb(ji,jj,:,:) )
         !
         ! Floe size distribution tendency:
         zdfsd(ji,jj,:)    = r1_Dt_ice * (  floe_size_dist( pa_ifsda(ji,jj,:,:), pa_ia(ji,jj,:) )   &
            &                             - floe_size_dist( pa_ifsdb(ji,jj,:,:), pa_ib(ji,jj,:) )   )
         !
         ! Mean floe size tendency from integrating FSD, above, which already has 1/dt factor:
         DO jf = 1, nn_nfsd
            zdsavg(ji,jj) = zdsavg(ji,jj) + floe_sc(jf) * zdfsd(ji,jj,jf)
         ENDDO
         !
      END_2D

      ! Determine suffix for field references. If it is total (tendency across whole time step
      ! i.e. all processes) then we do not add a suffix, otherwise it is the 3-char. input:
      IF( TRIM(cd_dia) == 'tot' ) THEN
         cl_sfx = ''
      ELSE
         cl_sfx = '_'//TRIM(cd_dia)
      ENDIF

      ! Write diagnostics:
      cl_ref = 'icefsd_cat_tend'//TRIM(cl_sfx)
      IF( iom_use( cl_ref ) )   CALL iom_put( cl_ref, zdfstd * zmsk00fc )

      cl_ref = 'icefsd_tend'//TRIM(cl_sfx)
      IF( iom_use( cl_ref ) )   CALL iom_put( cl_ref, zdfsd  * zmsk00f  )

      cl_ref = 'icefsdsavg_tend'//TRIM(cl_sfx)
      IF( iom_use( cl_ref ) )   CALL iom_put( cl_ref, zdsavg * zmsk00   )

   END SUBROUTINE ice_fsd_dia


   SUBROUTINE fsd_initbounds
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE fsd_init_bounds  ***
      !!
      !! ** Purpose :   Calculate or read FSD category boundaries and related arrays
      !!
      !! ** Method  :   Select method to determine category limits from namelist parameter
      !!                nn_fsd_catini:
      !!
      !!                   0 = read nn_nfsd+1 directly from namelist parameter rn_fsd_catbnd
      !!                   1 = compute uniformly-spaced bounds
      !!                   2 = compute bounds with increasing spacing following Gaussian profile
      !!                   3 = compute bounds with exponentially-increasing spacing
      !!
      !!                For 1-3, bounds are placed between a minimum and maximum floe size (caliper diameter)
      !!                set via namelist parameters rn_fsd_smin and rn_fsd_smax. For 2-3, an additional
      !!                parameter rn_fsd_spc controls the degree of curvature/non-linearity in the
      !!                Gaussian or exponential curve.
      !!
      !!                nn_fsd_catini = 2 (Gaussian spacing); limits L(j) are computed as:
      !!
      !!                      L(j) = L(j-1) + k * [1 - EXP( -( (j-1)/(sigma*n) )^2 )]   for   j = 2..(n+1)
      !!
      !!                   where n = nn_nfsd, sigma = rn_fsd_spc, k is calculated to ensure that
      !!                   L(n+1) = smax, and L(1) is defined to be smin. The exponent includes a
      !!                   factor of n so that the overall shape is not affected by changing smin or
      !!                   smax and to make sigma a 'scaling' parameter rather than depending on choice of n.
      !!
      !!                nn_fsd_catini = 3 (exponentially-increasing spacing); limits L(j) are computed as:
      !!
      !!                      L(j) = L(j-1) + k * EXP( 10*sigma*(j-1)/n )   for   j = 2..(n+1)
      !!
      !!                   with parameters defined similarly to the Gaussian case. Here an ad-hoc factor
      !!                   of 10 is included to set an appropriate degree of non-linearity with default
      !!                   parameters. Particularly, increasing sigma much beyond 1 here can make spacing so
      !!                   small (at lower j) that it cannot be resolved. In all cases, a warning is thus
      !!                   written if any category width is below 1cm (arbitrarily).
      !!
      !!                Category limits are printed in ocean.output. The limits are then used to calculate
      !!                other related constant arrays, including the floe areas, welding array (floe_iweld),
      !!                and gradient in log space (for subroutine ice_fsd_brit).
      !!
      !!-------------------------------------------------------------------
      !
      REAL(wp), DIMENSION(nn_nfsd+1) ::   zlims   ! floe size category limits
      !
      REAL(wp) ::   znfsd           ! number of FSD categories as REAL
      REAL(wp) ::   zk              ! spacing scale factor for Gaussian/exponential limits case
      REAL(wp) ::   zfloe_aweld     ! area of two welded floes (for computing floe_iweld)
      INTEGER  ::   jf1, jf2, jf3   ! dummy loop indices
      INTEGER  ::   ierr            ! allocate status return value
      !
      !!-------------------------------------------------------------------

      ALLOCATE(floe_sl(nn_nfsd), floe_sc(nn_nfsd), floe_su(nn_nfsd), floe_ds(nn_nfsd),   &
         &     floe_al(nn_nfsd), floe_ac(nn_nfsd), floe_au(nn_nfsd),                     &
         &     floe_dlog_sc(nn_nfsd), floe_iweld(nn_nfsd, nn_nfsd), STAT=ierr)

      IF (ierr /= 0) CALL ctl_stop('fsd_init_bounds: could not allocate FSD size/area arrays')

      znfsd = REAL(nn_nfsd, KIND=wp)   ! for some computation of category limits below

      SELECT CASE( nn_fsd_catini )
            !
         CASE( 0 )   ! === Read from namelist === !
            !
            IF(lwp) WRITE(numout,*) 'nn_fsd_catini = 0  ==>>  FSD category limits written in namelist:'
            !
            zlims(:) = rn_fsd_catbnd(1:nn_nfsd+1)
            !
            ! These should NOT be used anywhere outside this routine, but just in case:
            rn_fsd_smin = zlims(1)
            rn_fsd_smax = zlims(nn_nfsd+1)
            !
         CASE( 1 )   ! === Uniformly-spaced bounds === !
            !
            IF(lwp) WRITE(numout,*) 'nn_fsd_catini = 1  ==>>  FSD category limits are uniformly spaced:'
            !
            zlims(1)         = rn_fsd_smin
            zlims(nn_nfsd+1) = rn_fsd_smax
            !
            DO jf1 = 2, nn_nfsd
               zlims(jf1) = rn_fsd_smin + (rn_fsd_smax - rn_fsd_smin) * REAL(jf1 - 1, KIND=wp) / znfsd
            ENDDO
            !
         CASE( 2 )   ! === Gaussian-spaced bounds === !
            !
            IF(lwp) WRITE(numout,*) 'nn_fsd_catini = 2  ==>>  FSD category limits are Gaussian spaced:'
            !
            ! Determine multiplier k:
            zk = 0._wp
            DO jf1 = 1, nn_nfsd
               zk = zk + 1._wp - EXP( -( REAL(nn_nfsd - jf1 + 1, KIND=wp) / (rn_fsd_spc * znfsd) )**2 )
            ENDDO
            zk = (rn_fsd_smax - rn_fsd_smin) / zk
            !
            zlims(1) = rn_fsd_smin
            DO jf1 = 2, nn_nfsd + 1
               zlims(jf1) = zlims(jf1-1) + zk * ( 1._wp - EXP( -(REAL(jf1 - 1, KIND=wp) / (rn_fsd_spc * znfsd) )**2) )
            ENDDO
            !
         CASE( 3 )   ! === Exponentially-spaced bounds === !
            !
            IF(lwp) WRITE(numout,*) 'nn_fsd_catini = 3  ==>>  FSD category spacing increases exponentially:'
            !
            ! Determine multiplier k:
            zk = 0._wp
            DO jf1 = 2, nn_nfsd + 1
               zk = zk + EXP( 10._wp * rn_fsd_spc * REAL(jf1 - 1, KIND=wp) / znfsd )
            ENDDO
            zk = (rn_fsd_smax - rn_fsd_smin) / zk
            !
            zlims(1) = rn_fsd_smin
            DO jf1 = 2, nn_nfsd + 1
               zlims(jf1) = zlims(jf1-1) + zk * EXP( 10._wp * rn_fsd_spc * REAL(jf1 - 1, KIND=wp) / znfsd )
            ENDDO
            !
         CASE DEFAULT
            !
            CALL ctl_stop('fsd_init_bounds: must choose nn_fsd_catini = 0, 1, 2, or 3')
            !
      ENDSELECT

      floe_sl = zlims(1:nn_nfsd)
      floe_su = zlims(2:nn_nfsd+1)
      floe_sc = .5_wp * (floe_su + floe_sl)

      floe_ds = floe_su - floe_sl

      ! Write FSD bounds in ocean.output (continuing from control print in ice_fsd_init)
      IF(lwp) THEN
         WRITE(numout,*)
         DO jf1 = 1, nn_nfsd
            WRITE(numout,'(A,F12.5,A,I2,A,F12.5,A)') '                         ',   &
               &    floe_sl(jf1), ' m <= category ', jf1, ' < ', floe_su(jf1), ' m'
         ENDDO
         WRITE(numout,*)
         !
         ! Write uniform or min./max. category width(s):
         IF( nn_fsd_catini == 1 ) THEN
            WRITE(numout,'(A,A,F12.5,A)') '                      ',   &
               &   '==>>> Uniform categories of width: ', floe_ds(1), ' m'
         ELSE
            WRITE(numout,'(A,A,F12.5,A)') '           ',   &
               &   '==>>> Non-uniform categories, smallest width: ', MINVAL(floe_ds(:)), ' m'
            WRITE(numout,'(A,A,F12.5,A)') '           ',   &
               &   '                               largest width: ', MAXVAL(floe_ds(:)), ' m'
         ENDIF
         WRITE(numout,*) ''
      ENDIF

      ! Sometimes automatic category spacing is too small, particularly in exponential case
      ! Check for small category widths and warn with suggested changes in each case:
      IF( ANY( ABS(floe_sl(:)) < 1.e-2 ) ) THEN
         CALL ctl_warn('fsd_init_bounds: some FSD categories are very small, < 1cm width; consider:'   ,   &
               &       '                 nn_fsd_catini = 0  : making your categories wider'            ,   &
               &       '                 nn_fsd_catini = 1-2: (in/de)creasing rn_fsd_smin/rn_fsd_smax)',   &
               &       '                 nn_fsd_catini = 2-3: decreasing rn_fsd_spc  (recommend <= 1)'     )
      ENDIF

      floe_al = rn_floeshape * floe_sl ** 2
      floe_ac = rn_floeshape * floe_sc ** 2
      floe_au = rn_floeshape * floe_su ** 2

      ! --- Calculate category index of default new ice floe size set in namelist
      nf_newice = nn_nfsd
      DO jf1 = nn_nfsd-1, 1, -1
         IF( (rn_fsd_s_newice >= floe_sl(jf1)) .AND. (rn_fsd_s_newice < floe_su(jf1)) ) THEN
            nf_newice = jf1
            EXIT
         ENDIF
      ENDDO

      ! --- Calculate floe welding array, floe_iweld
      ! floe_iweld(jf1,jf2) = index of FSD category that floes in category jf1,
      ! when welded with floes in category jf2, subsequently belong to
      !
      floe_iweld(:,:) = 0   ! 'missing' value (if no category for welding)
      !
      DO jf1 = 1, nn_nfsd
         DO jf2 = 1, nn_nfsd
            !
            ! --- If floes from centers of cat jf1 and jf2 weld, their new area is:
            zfloe_aweld = floe_ac(jf1) + floe_ac(jf2)
            !
            ! --- Find FSD category that fits into
            !     Check each floe size category; only one can be true:
            DO jf3 = 1, nn_nfsd-1
               IF( (zfloe_aweld >= floe_al(jf3)) .and. (zfloe_aweld < floe_au(jf3))) THEN
                  floe_iweld(jf1,jf2) = jf3
               ENDIF
            ENDDO
            ! --- Separate check for largest category:
            IF( zfloe_aweld >= floe_al(nn_nfsd)) floe_iweld(jf1,jf2) = nn_nfsd
         ENDDO
      ENDDO

      ! --- Calculate category spacing in log(s) space (for FSD restoring routine)
      !
      floe_dlog_sc(:) = 0._wp   ! initialise
      !
      DO jf1 = 2, nn_nfsd
         floe_dlog_sc(jf1) = LOG(floe_sc(jf1)) - LOG(floe_sc(jf1-1))
      ENDDO

   END SUBROUTINE fsd_initbounds


   SUBROUTINE ice_fsd_istate
      !!-------------------------------------------------------------------
      !!                 ***  ROUTINE ice_fsd_istate  ***
      !!
      !! ** Purpose :   Set initial values of floe size distribution
      !!
      !! ** Method  :   Set values based on namelist (namfsd) nn_fsd_ini:
      !!                   0 = no initialisation (i.e., all FSD values = 0)
      !!                   1 = all ice in largest floe size category
      !!                   2 = set all grid points, all ice thickness categories
      !!                       to have an imposed power law distribution. In this
      !!                       case the number density distribution exponent can
      !!                       be changed via namelist (namfsd) rn_fsd_ini_alpha
      !!                       (default = 2.1 as in Perovich and Jones, 2014).
      !!
      !! ** Note    :   Default nn_fsd_ini = 2. If general ice initialisation
      !!                flag, ln_iceini, is set to false, then nn_fsd_ini is
      !!                treated as in case 0 regardless of namelist value, i.e.,
      !!                no initialisation of the FSD is done. This allows the
      !!                FSD to 'emerge' from physical processes.
      !!
      !! ** References
      !!    ----------
      !!    Perovich, D. K. & Jones, K. F. (2014).
      !!              The seasonal evolution of sea ice floe size distribution.
      !!              Journal of Geophysical Research: Oceans, 119(12), 8767-8777.
      !!-------------------------------------------------------------------
      !
      LOGICAL  ::   llfsdini   ! condition whether to initialise FSD (T) or set to 0 (F)
      REAL(wp) ::   ztotfrac   ! for normalising
      INTEGER  ::   jf, jl     ! dummy variables for loop indices
      !
      !!-------------------------------------------------------------------

      ! === Determine whether to initialise or not === !
      !
      ! This routine is either called from ice_istate or from ice_rst_read.
      !
      ! If we are here and (ln_rstart = T OR nn_iceini_file == 2), this indicates restart read was
      ! attempted in the latter routine, but the restart file was found to have no FSD variables or
      ! wrong number of floe size categories and so was bypassed. But other variables *were* read
      ! from the restart file, are so are non-zero initialised. Therefore, FSD should be initialised.
      !
      ! If not (ln_rstart = F AND nn_iceini_file /= 2), we are here from ice_istate and so whether to
      ! initialise FSD or not is based on ln_iceini:
      !
      IF( ln_rstart .OR. (nn_iceini_file == 2) ) THEN
         llfsdini = .TRUE.      ! => here because we bypassed restart
      ELSE
         llfsdini = ln_iceini   ! => general initialisation case
      ENDIF

      ! === Warnings / Checks === !
      !
      ! We have no specific treatment for FSD if reading from a 'single category file' (nn_iceini_file == 1)
      ! If user wishes to start FSD from file, it must be a restart file, which is done in ice_rst_read
      ! for cases ln_restart = T .OR. (ln_iceini = T and nn_iceini_file == 2)
      !
      ! NOTE: value of nn_iceini_file only relevant when ln_iceini = T AND NOT ln_rstart
      ! Important to add conditions on the latter, otherwise irrelevant warning is raised
      !
      IF( nn_iceini_file == 1 .AND. ln_iceini .AND. (.NOT. ln_rstart) ) THEN
         CALL ctl_warn( 'ice_fsd_istate ===>>> : Single-category file read (nn_iceini_file == 1) not possible for FSD', &
            &           'we initialise FSD internally (i.e., NOT from file) according to nn_fsd_ini')
         llfsdini = .TRUE.  ! should be covered above, but does not hurt
      ENDIF

      ! === Initialise FSD values === !
      !
      IF( llfsdini ) THEN
         IF( nn_fsd_ini == 1 ) THEN
            IF(lwp) WRITE(numout,*) '   ice_fsd_istate   ==>>   floes initially all in largest category'
            !
            a_ifsd(:,:,nn_nfsd,:) = 1._wp
            !
         ELSE  ! >= 2
            IF(lwp) WRITE(numout,*) '   ice_fsd_istate   ==>>   imposed power law for initial FSD everywhere'
            !
            ztotfrac = 0._wp
            !
            ! Initial FSD is the same for each ice thickness category
            ! Calculate for first category:
            DO jf = 1, nn_nfsd
               ! Calculate power law FSD number distribution based on Perovich
               ! and Jones (2014) and convert to area fraction distribution:
               a_ifsd(:,:,jf,1) = floe_sc(jf) ** (-rn_fsd_ini_alpha - 1._wp) * floe_ac(jf) * floe_ds(jf)

               ztotfrac = ztotfrac + a_ifsd(1,1,jf,1)
            ENDDO
            !
            a_ifsd(:,:,:,1) = a_ifsd(:,:,:,1) / ztotfrac   ! normalise
            !
            ! Assign same initial FSD to remaining thickness categories:
            DO jl = 2, jpl
               a_ifsd(:,:,:,jl) = a_ifsd(:,:,:,1)
            ENDDO
            !
         ENDIF
      ELSE
         IF(lwp) WRITE(numout,*) '   ice_fsd_istate   ==>>   initial FSD = 0 (no initialisation)'
         !
         a_ifsd(:,:,:,:) = 0._wp
         !
      ENDIF

      IF(lwp) WRITE(numout,*) ''

   END SUBROUTINE ice_fsd_istate


   SUBROUTINE ice_fsd_init
      !!-------------------------------------------------------------------
      !!                  ***  ROUTINE ice_fsd_init   ***
      !! 
      !! ** Purpose :   Check whether FSD is to be activated, and if so carry
      !!                out initialisation of FSD, printing parameter values
      !!                to STDOUT, and call other FSD initialisation routines
      !! 
      !! ** Method  :   Read the namfsd namelist, call other initialisation
      !!                subroutines in the module if FSD is activated.
      !! 
      !! ** input   :   Namelist namfsd
      !!-------------------------------------------------------------------
      INTEGER ::   jf            ! Local loop index for FSD categories
      INTEGER ::   ios, ioptio   ! Local integer output status for namelist read
      INTEGER ::   ierr          ! Local integer allocate status
      !!
      NAMELIST/namfsd/ ln_fsd          , nn_fsd_catini   , nn_nfsd        , rn_fsd_smin     ,   &
         &             rn_fsd_smax     , rn_fsd_spc      , rn_fsd_catbnd  , rn_floeshape    ,   &
         &             nn_fsd_ini      , rn_fsd_ini_alpha, rn_fsd_s_newice                  ,   &
         &             ln_fsd_brit     , rn_fsd_brit_grad, rn_fsd_brit_tres                 ,   &
         &             rn_fsd_amin_weld, rn_fsd_c_weld
      !!-------------------------------------------------------------------
      !
      READ_NML_REF(numnam_ice, namfsd)
      READ_NML_CFG(numnam_ice, namfsd)
      IF(lwm) WRITE(numoni, namfsd)
      !
      IF(lwp) THEN   ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_fsd_init: ice parameters for floe size distribution'
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namfsd:'
         WRITE(numout,*) '      Floe size distribution activated or not                    ln_fsd = ', ln_fsd
         WRITE(numout,*) '         FSD category initialisation                      nn_fsd_catini = ', nn_fsd_catini
         WRITE(numout,*) '            Number of floe size categories                      nn_nfsd = ', nn_nfsd
         WRITE(numout,*) '            Minimum floe size     (nn_fsd_catini /= 0  )    rn_fsd_smin = ', rn_fsd_smin
         WRITE(numout,*) '            Maximum floe size     (nn_fsd_catini /= 0  )    rn_fsd_smax = ', rn_fsd_smax
         WRITE(numout,*) '            Spacing non-linearity (nn_fsd_catini  = 2,3)    rn_fsd_spc  = ', rn_fsd_spc
         WRITE(numout,*) '            Floe shape parameter, to determine floe areas  rn_floeshape = ', rn_floeshape
         WRITE(numout,*) '         FSD initialisation case (ln_iceini = T)             nn_fsd_ini = ', nn_fsd_ini
         WRITE(numout,*) '            Power law exponent  (nn_fsd_ini = 2)       rn_fsd_ini_alpha = ', rn_fsd_ini_alpha
         WRITE(numout,*) '         Floe size of new ice (in absence of waves)    rn_fsd_s_newice  = ', rn_fsd_s_newice
         WRITE(numout,*) '         Floe welding minimum sea ice concentration    rn_fsd_amin_weld = ', rn_fsd_amin_weld
         WRITE(numout,*) '         Floe welding coefficient                         rn_fsd_c_weld = ', rn_fsd_c_weld
         WRITE(numout,*) '         Activate brittle fracture scheme or not            ln_fsd_brit = ', ln_fsd_brit
         WRITE(numout,*) '            Max. gradient of number-density FSD        rn_fsd_brit_grad = ', rn_fsd_brit_grad
         WRITE(numout,*) '            Restoring time scale                       rn_fsd_brit_tres = ', rn_fsd_brit_tres
         WRITE(numout,*) ''
      ENDIF

      IF(ln_fsd) THEN

         ALLOCATE(a_ifsd   (jpi,jpj,nn_nfsd,jpl), a_ifsd_b(jpi,jpj,nn_nfsd,jpl),   &
            &     a_ifsd_b0(jpi,jpj,nn_nfsd,jpl), a_i_b0  (jpi,jpj,jpl),           &
            &     STAT=ierr                                                        )

         IF( ierr /= 0 )   CALL ctl_stop('ice_fsd_init: could not allocate arrays')

         CALL fsd_initbounds

      ELSE
         ! Set FSD-related logicals to false to avoid issues:
         ln_fsd_brit = .FALSE.
      ENDIF

   END SUBROUTINE ice_fsd_init

#else
   !!----------------------------------------------------------------------
   !!   Default option          Empty module          NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE icefsd
