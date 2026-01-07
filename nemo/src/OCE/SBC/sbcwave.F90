MODULE sbcwave
   !!======================================================================
   !!                       ***  MODULE  sbcwave  ***
   !! Wave module
   !!======================================================================
   !! History :  3.3  !  2011-09  (M. Adani)  Original code: Drag Coefficient
   !!         :  3.4  !  2012-10  (M. Adani)  Stokes Drift
   !!            3.6  !  2014-09  (E. Clementi,P. Oddo) New Stokes Drift Computation
   !!             -   !  2016-12  (G. Madec, E. Clementi) update Stoke drift computation
   !!                                                    + add sbc_wave_ini routine
   !!            4.2  !  2020-12  (G. Madec, E. Clementi) updates, new Stoke drift computation
   !!                                                    according to Couvelard et al.,2019
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_stokes    : calculate 3D Stokes-drift velocities
   !!   sbc_wave      : wave data from wave model: forced (netcdf files) or coupled mode
   !!   sbc_wave_init : initialisation fo surface waves
   !!----------------------------------------------------------------------
   USE phycst         ! physical constants
   USE oce            ! ocean variables
   USE dom_oce        ! ocean domain variables
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE bdy_oce        ! open boundary condition variables
   USE zdf_oce,  ONLY : ln_zdfswm ! Qiao wave enhanced mixing
   USE par_ice,  ONLY : ln_ice_wav, ln_ice_wav_attn, ln_ice_wav_spec ! sea-ice parameters
   !
   USE iom            ! I/O manager library
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE fldread        ! read input fields

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_stokes      ! routine called in sbccpl
   PUBLIC   sbc_wave        ! routine called in sbcmod
   PUBLIC   sbc_wave_init   ! routine called in sbcmod

   ! Variables checking if the wave parameters are coupled (if not, they are read from file)
   LOGICAL, PUBLIC ::   cpl_hsig          = .FALSE.
   LOGICAL, PUBLIC ::   cpl_phioc         = .FALSE.
   LOGICAL, PUBLIC ::   cpl_sdrft         = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wper          = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wnum          = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wstrf         = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wdrag         = .FALSE.
   LOGICAL, PUBLIC ::   cpl_charn         = .FALSE.
   LOGICAL, PUBLIC ::   cpl_taw           = .FALSE.
   LOGICAL, PUBLIC ::   cpl_bhd           = .FALSE.
   LOGICAL, PUBLIC ::   cpl_tsd           = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wpf           = .FALSE.
   LOGICAL, PUBLIC ::   cpl_wspec         = .FALSE.

   INTEGER ::   jpfld    ! number of files to read for stokes drift
   INTEGER ::   jp_usd   ! index of stokes drift  (i-component) (m/s)    at T-point
   INTEGER ::   jp_vsd   ! index of stokes drift  (j-component) (m/s)    at T-point
   INTEGER ::   jp_wmp   ! index of mean wave period            (s)      at T-point

   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_cdg     ! structure of input fields (file informations, fields read) Drag Coefficient
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_sd      ! structure of input fields (file informations, fields read) Stokes Drift
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_wnum    ! structure of input fields (file informations, fields read) wave number for Qiao
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_tauoc   ! structure of input fields (file informations, fields read) normalized wave stress into the ocean
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_hsw     ! structure of input fields (file informations, fields read) Significant Wave Height
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_wpf     ! structure of input fields (file informations, fields read) Wave Peak Frequency
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf_wspec   ! structure of input fields (file informations, fields read) Wave Energy Spectrum

   ! dynamics (Stokes)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   hsw             !: Significant Wave Height at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   wmp             !: Wave Mean Period at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   wpf             !: Wave Peak Frequency at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tsd2d           !: Surface Stokes Drift module at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   div_sd          !: barotropic stokes drift divergence
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   ut0sd, vt0sd    !: surface Stokes drift velocities at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tusd, tvsd      !: Stokes drift transport
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   usd, vsd, wsd   !: Stokes drift velocities at u-, v- & w-points, resp.u
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   bhd_wave        !: Bernoulli head. wave induce pression
   ! stress
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tauoc_wave      !: stress reduction factor  at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   cdn_wave        !: Neutral drag coefficient at t-point
!
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:)     ::   wfreq_l         !: Wave frequencies for discretised energy spectrum (lower limits of bins)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:)     ::   wfreq           !: Wave frequencies for discretised energy spectrum ('center' of bins)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:)     ::   wfreq_u         !: Wave frequencies for discretised energy spectrum (upper limits of bins)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:)     ::   wdfreq          !: Wave frequencies for discretised energy spectrum (bin widths)
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:,:) ::   wspec           !: Wave energy spectrum (function of frequency) at t-point
!
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   charn           !: charnock coefficient at t-point
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tawx            !: Net wave-supported stress, u
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tawy            !: Net wave-supported stress, v
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   twox            !: wave-ocean momentum flux, u
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   twoy            !: wave-ocean momentum flux, v
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tauoc_wavex     !: stress reduction factor  at, u component
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   tauoc_wavey     !: stress reduction factor  at, v component
   ! mixing
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   phioc           !: tke flux from wave model
   REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:)   ::   wnum            !: Wave Number at t-point
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_stokes( Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_stokes  ***
      !!
      !! ** Purpose :   compute the 3d Stokes Drift according to Breivik et al.,
      !!                2014 (DOI: 10.1175/JPO-D-14-0020.1)
      !!
      !! ** Method  : - Calculate the horizontal Stokes drift velocity (Breivik et al. 2014)
      !!              - Calculate its horizontal divergence
      !!              - Calculate the vertical Stokes drift velocity
      !!              - Calculate the barotropic Stokes drift divergence
      !!
      !! ** action  : - tsd2d         : module of the surface Stokes drift velocity
      !!              - usd, vsd, wsd : 3 components of the Stokes drift velocity
      !!              - div_sd        : barotropic Stokes drift divergence
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: Kmm ! ocean time level index
      INTEGER  ::   jj, ji, jk   ! dummy loop argument
      INTEGER  ::   ik           ! local integer
      REAL(wp) ::  ztransp, zfac, ztemp, zsp0, zsqrt, zbreiv16_w
      REAL(wp) ::  zdep_u, zdep_v, zkh_u, zkh_v, zda_u, zda_v, sdtrp, zInt_w0, zInt_w1
#if ! defined key_PSYCLONE_2p5p0
      REAL(wp), DIMENSION(:,:)  , ALLOCATABLE ::   zk_t, zk_u, zk_v, zu0_sd, zv0_sd ! 2D workspace
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   ze3divh                          ! 3D workspace
#else
      REAL(wp), DIMENSION(A2D(1))        ::   zk_t, zk_u, zk_v, zu0_sd, zv0_sd ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpkm1) ::   ze3divh                          ! 3D workspace
#endif
      !!---------------------------------------------------------------------
      !
#if ! defined key_PSYCLONE_2p5p0
      ALLOCATE( ze3divh(jpi,jpj,jpkm1) ) ! jpkm1 -> avoid lbc_lnk on jpk that is not defined
      ALLOCATE( zk_t(A2D(1)), zk_u(A2D(1)), zk_v(A2D(1)), zu0_sd(A2D(1)), zv0_sd(A2D(1)) )
#endif
      zk_t    (:,:) = 0._wp
      zk_u    (:,:) = 0._wp
      zk_v    (:,:) = 0._wp
      zu0_sd  (:,:) = 0._wp
      zv0_sd  (:,:) = 0._wp
      ze3divh (:,:,:) = 0._wp

      !
      ! select parameterization for the calculation of vertical Stokes drift
      ! exp. wave number at t-point
      IF( ln_breivikFV_2016 ) THEN
      ! Assumptions :  ut0sd and vt0sd are surface Stokes drift at T-points
      !                sdtrp is the norm of Stokes transport
      !
         zfac = 0.166666666667_wp
         DO_2D( 1, 1, 1, 1 ) ! In the deep-water limit we have ke = ||ust0||/( 6 * ||transport|| )
            zsp0          = SQRT( ut0sd(ji,jj)*ut0sd(ji,jj) + vt0sd(ji,jj)*vt0sd(ji,jj) ) !<-- norm of Surface Stokes drift
            tsd2d(ji,jj)  = zsp0
            IF( cpl_tsd ) THEN  !stokes transport is provided in coupled mode
               sdtrp      = SQRT( tusd(ji,jj)*tusd(ji,jj) + tvsd(ji,jj)*tvsd(ji,jj) )  !<-- norm of Surface Stokes drift transport
            ELSE
               ! Stokes drift transport estimated from Hs and Tmean
               sdtrp      = 2.0_wp * rpi / 16.0_wp *                             &
                   &        hsw(ji,jj)*hsw(ji,jj) / MAX( wmp(ji,jj), 0.0000001_wp )
            ENDIF
            zk_t (ji,jj)  = zfac * zsp0 / MAX ( sdtrp, 0.0000001_wp ) !<-- ke = ||ust0||/( 6 * ||transport|| )
         END_2D
         !
         DO jk = 1, jpkm1
            zfac = 0.166666666667_wp
            DO_2D( 1, 1, 1, 1 ) !++ Compute the FV Breivik 2016 function at T-points
               ! zInt at jk
               zfac       = - 2._wp * zk_t (ji,jj) * gdepw(ji,jj,jk,Kmm)  !<-- zfac should be negative definite
               ztemp      = EXP ( zfac )
               zsqrt      = SQRT( -zfac )
               zbreiv16_w = ztemp - SQRT(rpi)*zsqrt*ERFC(zsqrt) !Eq. 16 Breivik 2016
               zInt_w0    = ztemp - 4._wp * zk_t (ji,jj) * gdepw(ji,jj,jk,Kmm) * zbreiv16_w
               ! zInt at jk+1
               zfac       = - 2._wp * zk_t (ji,jj) * gdepw(ji,jj,jk+1,Kmm)  !<-- zfac should be negative definite
               ztemp      = EXP ( zfac )
               zsqrt      = SQRT( -zfac )
               zbreiv16_w = ztemp - SQRT(rpi)*zsqrt*ERFC(zsqrt) !Eq. 16 Breivik 2016
               zInt_w1    = ztemp - 4._wp * zk_t (ji,jj) * gdepw(ji,jj,jk+1,Kmm) * zbreiv16_w
               !
               !
               zsp0          = zfac / MAX(zk_t (ji,jj),0.0000001_wp)
               ztemp         = zInt_w0 - zInt_w1
               zu0_sd(ji,jj) = ut0sd(ji,jj) * zsp0 * ztemp * tmask(ji,jj,jk)
               zv0_sd(ji,jj) = vt0sd(ji,jj) * zsp0 * ztemp * tmask(ji,jj,jk)
            END_2D
            DO_2D( 1, 0, 1, 0 ) ! ++ Interpolate at U/V points
               zfac          =  1.0_wp / e3u(ji  ,jj,jk,Kmm)
               usd(ji,jj,jk) =  0.5_wp * zfac * ( zu0_sd(ji,jj)+zu0_sd(ji+1,jj) ) * umask(ji,jj,jk)
               zfac          =  1.0_wp / e3v(ji  ,jj,jk,Kmm)
               vsd(ji,jj,jk) =  0.5_wp * zfac * ( zv0_sd(ji,jj)+zv0_sd(ji,jj+1) ) * vmask(ji,jj,jk)
            END_2D
         ENDDO
         !
      ELSE
         zfac = 2.0_wp * rpi / 16.0_wp
         DO_2D( 1, 1, 1, 1 )
            ! Stokes drift velocity estimated from Hs and Tmean
            ztransp = zfac * hsw(ji,jj)*hsw(ji,jj) / MAX( wmp(ji,jj), 0.0000001_wp )
            ! Stokes surface speed
            tsd2d(ji,jj) = SQRT( ut0sd(ji,jj)*ut0sd(ji,jj) + vt0sd(ji,jj)*vt0sd(ji,jj))
            ! Wavenumber scale
            zk_t(ji,jj) = ABS( tsd2d(ji,jj) ) / MAX( ABS( 5.97_wp*ztransp ), 0.0000001_wp )
         END_2D
         DO_2D( 1, 0, 1, 0 )          ! exp. wave number & Stokes drift velocity at u- & v-points
            zk_u(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji+1,jj) )
            zk_v(ji,jj) = 0.5_wp * ( zk_t(ji,jj) + zk_t(ji,jj+1) )
            !
            zu0_sd(ji,jj) = 0.5_wp * ( ut0sd(ji,jj) + ut0sd(ji+1,jj) )
            zv0_sd(ji,jj) = 0.5_wp * ( vt0sd(ji,jj) + vt0sd(ji,jj+1) )
         END_2D

      !                       !==  horizontal Stokes Drift 3D velocity  ==!

         DO_3D( 1, 0, 1, 0, 1, jpkm1 )
            zdep_u = 0.5_wp * ( gdept(ji,jj,jk,Kmm) + gdept(ji+1,jj,jk,Kmm) )
            zdep_v = 0.5_wp * ( gdept(ji,jj,jk,Kmm) + gdept(ji,jj+1,jk,Kmm) )
            !
            zkh_u = zk_u(ji,jj) * zdep_u     ! k * depth
            zkh_v = zk_v(ji,jj) * zdep_v
            !                                ! Depth attenuation
            zda_u = EXP( -2.0_wp*zkh_u ) / ( 1.0_wp + 8.0_wp*zkh_u )
            zda_v = EXP( -2.0_wp*zkh_v ) / ( 1.0_wp + 8.0_wp*zkh_v )
            !
            usd(ji,jj,jk) = zda_u * zu0_sd(ji,jj) * umask(ji,jj,jk)
            vsd(ji,jj,jk) = zda_v * zv0_sd(ji,jj) * vmask(ji,jj,jk)
         END_3D
      ENDIF
      !
      !                       !==  vertical Stokes Drift 3D velocity  ==!
      !
      DO_3D( 0, 0, 0, 0, 1, jpkm1 )    ! Horizontal e3*divergence
         ze3divh(ji,jj,jk) = (  ( e2u(ji  ,jj) * e3u(ji  ,jj,jk,Kmm) * usd(ji  ,jj,jk)     &   ! add () for NP repro
            &                   - e2u(ji-1,jj) * e3u(ji-1,jj,jk,Kmm) * usd(ji-1,jj,jk) )   &
            &                 + ( e1v(ji,jj  ) * e3v(ji,jj  ,jk,Kmm) * vsd(ji,jj  ,jk)     &
            &                   - e1v(ji,jj-1) * e3v(ji,jj-1,jk,Kmm) * vsd(ji,jj-1,jk) ) ) * r1_e1e2t(ji,jj)
      END_3D
      !
      CALL lbc_lnk( 'sbcwave', ze3divh, 'T', 1.0_wp, usd, 'U', -1.0_wp, vsd, 'V', -1.0_wp )
      !
      IF( lk_linssh ) THEN   ;   ik = 1   ! none zero velocity through the sea surface
      ELSE                   ;   ik = 2   ! w=0 at the surface (set one for all in sbc_wave_init)
      ENDIF
      DO jk = jpkm1, ik, -1          ! integrate from the bottom the hor. divergence (NB: at k=jpk w is always zero)
         wsd(:,:,jk) = wsd(:,:,jk+1) - ze3divh(:,:,jk)
      END DO
      !
      IF( ln_bdy ) THEN
         DO jk = 1, jpkm1
            wsd(:,:,jk) = wsd(:,:,jk) * bdytmask(:,:)
         END DO
      ENDIF
      !                       !==  Horizontal divergence of barotropic Stokes transport  ==!
      div_sd(:,:) = 0._wp
      DO jk = 1, jpkm1                                 !
        div_sd(:,:) = div_sd(:,:) + ze3divh(:,:,jk)
      END DO
      !
      IF( iom_use( "ustokes" ) ) CALL iom_put( "ustokes", usd  )
      IF( iom_use( "vstokes" ) ) CALL iom_put( "vstokes", vsd  )
      IF( iom_use( "wstokes" ) ) CALL iom_put( "wstokes", wsd  )
      !
#if ! defined key_PSYCLONE_2p5p0
      DEALLOCATE( ze3divh )
      DEALLOCATE( zk_t, zk_u, zk_v, zu0_sd, zv0_sd )
#endif
      !
   END SUBROUTINE sbc_stokes
!
!
   SUBROUTINE sbc_wave( kt, Kmm )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_wave  ***
      !!
      !! ** Purpose :   read wave parameters from wave model in netcdf files
      !!                or from a coupled wave mdoel
      !!
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt   ! ocean time step
      INTEGER, INTENT(in   ) ::   Kmm  ! ocean time index
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_wave : update the read waves fields'
         WRITE(numout,*) '~~~~~~~~ '
      ENDIF
      !
      IF( ln_cdgw .AND. .NOT. cpl_wdrag ) THEN     !==  Neutral drag coefficient  ==!
         CALL fld_read( kt, nn_fsbc, sf_cdg )            ! read from external forcing
         cdn_wave(:,:) = sf_cdg(1)%fnow(:,:,1) * smask0(:,:)
      ENDIF
      !
      IF( ln_tauoc .AND. .NOT. cpl_wstrf ) THEN    !==  Wave induced stress  ==!
         CALL fld_read( kt, nn_fsbc, sf_tauoc )          ! read stress reduction factor due to wave from external forcing
         tauoc_wave(:,:) = sf_tauoc(1)%fnow(:,:,1) * smask0(:,:)
      ELSEIF ( ln_taw .AND. cpl_taw ) THEN
         IF (kt < 1) THEN ! The first fields gave by OASIS have very high erroneous values ....
            twox(:,:)=0._wp
            twoy(:,:)=0._wp
            tawx(:,:)=0._wp
            tawy(:,:)=0._wp
            tauoc_wavex(:,:) = 1._wp
            tauoc_wavey(:,:) = 1._wp
         ELSE
            tauoc_wavex(:,:) = abs(twox(:,:)/tawx(:,:))
            tauoc_wavey(:,:) = abs(twoy(:,:)/tawy(:,:))
         ENDIF
      ENDIF

      IF( .NOT. ln_wave_test ) THEN
         IF( ln_wave_spec .AND. .NOT. cpl_wspec ) THEN  !== Wave energy spectrum ==!
            CALL fld_read( kt, nn_fsbc, sf_wspec )
            wspec(:,:,:) = sf_wspec(1)%fnow(:,:,:) * SPREAD(tmask(:,:,1), 3, nn_nwfreq)  ! wave energy spectrum at T point (function of frequency)
         ENDIF

         IF( ln_ice_wav .AND. .NOT. cpl_wpf ) THEN      !== Peak Frequency ==!
            CALL fld_read( kt, nn_fsbc, sf_wpf )
            wpf(:,:) = sf_wpf(1)%fnow(:,:,1) * tmask(:,:,1)  ! wave peak frequency at T point
         ENDIF

         IF( (ln_sdw .OR. ln_ice_wav) .AND. .NOT. cpl_hsig ) THEN   !== Significant Wave Height ==!
            CALL fld_read( kt, nn_fsbc, sf_hsw )
            hsw(:,:) = sf_hsw(1)%fnow(:,:,1) * tmask(:,:,1)  ! significant wave height at T point
         ENDIF
      ENDIF
      !
      IF( ln_sdw .AND. .NOT. cpl_sdrft)  THEN       !==  Computation of the 3d Stokes Drift  ==!
         !
         IF( jpfld > 0 ) THEN                            ! Read from file only if the field is not coupled
            CALL fld_read( kt, nn_fsbc, sf_sd )          ! read wave parameters from external forcing
            !                                            ! NB: test case mode, not read as jpfld=0
            !                                            ! NB: hsw now read separately, above, as also needed by icewav
            IF( jp_wmp > 0 )   wmp  (:,:) = sf_sd(jp_wmp)%fnow(:,:,1) * tmask(:,:,1)  ! wave mean period
            IF( jp_usd > 0 )   ut0sd(:,:) = sf_sd(jp_usd)%fnow(:,:,1) * tmask(:,:,1)  ! 2D zonal Stokes Drift at T point
            IF( jp_vsd > 0 )   vt0sd(:,:) = sf_sd(jp_vsd)%fnow(:,:,1) * tmask(:,:,1)  ! 2D meridional Stokes Drift at T point
         ENDIF

         IF( ((jpfld == 3) .AND. .NOT. cpl_hsig) .OR. ln_wave_test )   &
            &      CALL sbc_stokes( Kmm )                ! Calculate only if all required fields are read
            !                                            ! or in wave test case
            !                                            ! In coupled case the call is done after (in sbc_cpl)
      ENDIF
      !
      IF( ln_zdfswm .AND. .NOT. cpl_wnum ) THEN     !== wavenumber ==!
         CALL fld_read( kt, nn_fsbc, sf_wnum )             ! read wave parameters from external forcing
         wnum(:,:) = sf_wnum(1)%fnow(:,:,1) * smask0(:,:)
      ENDIF
      !
   END SUBROUTINE sbc_wave


   SUBROUTINE sbc_wave_init
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_wave_init  ***
      !!
      !! ** Purpose :   Initialisation fo surface waves
      !!
      !! ** Method  : - Read namelist namsbc_wave
      !!              - create the structure used to read required wave fields
      !!                (its size depends on namelist options)
      !! ** action
      !!---------------------------------------------------------------------
      INTEGER ::   ierror, ios   ! local integer
      INTEGER ::   ifpr
      INTEGER ::   ji, jj, jf    ! dummy loop indices
      !!
      CHARACTER(len=100)     ::  cn_dir                            ! Root directory for location of drag coefficient files
      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::   slf_i            ! array of namelist informations on the fields to read
      TYPE(FLD_N)            ::  sn_cdg, sn_usd, sn_vsd,  &
                             &   sn_hsw, sn_wmp, sn_wnum, sn_tauoc, &
                             &   sn_wpf, sn_wspec                   ! informations about the fields to be read
      !
      NAMELIST/namsbc_wave/ cn_dir, sn_cdg, sn_usd, sn_vsd, sn_hsw, sn_wmp, sn_wnum, sn_tauoc,   &
         &                  sn_wpf, sn_wspec,                                                    &
         &                  ln_sdw, ln_wave_test, ln_breivikFV_2016, ln_vortex_force, ln_bern_srfc, &
         &                  ln_cdgw, ln_charn, ln_tauoc, ln_taw, ln_phioc, ln_stshear, &
         &                  ln_wave_spec, ln_wfreq_usr, ln_wfreq_usr_exp, rn_wfreq_usr,     &
         &                  nn_nwfreq, rn_wfreq_0, rn_wfreq_k
      !!---------------------------------------------------------------------
      IF( .NOT. ln_wave ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) '   No surface waves : all wave related logical set to false'
         ln_sdw   = .FALSE.
         ln_wave_test = .FALSE.
         ln_breivikFV_2016 = .FALSE.
         ln_vortex_force = .FALSE.
         ln_bern_srfc = .FALSE.
         ln_cdgw  = .FALSE.
         ln_charn = .FALSE.
         ln_tauoc = .FALSE.
         ln_taw   = .FALSE.
         ln_phioc = .FALSE.
         ln_stshear = .FALSE.
         RETURN
      ENDIF
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_wave_init : surface waves in the system'
         WRITE(numout,*) '~~~~~~~~~~~~~ '
      ENDIF
      !
      READ_NML_REF(numnam,namsbc_wave)
      READ_NML_CFG(numnam,namsbc_wave)
      IF(lwm) WRITE ( numond, namsbc_wave )
      !
      IF(lwp) THEN
         WRITE(numout,*) '   Namelist namsbc_wave'
         WRITE(numout,*) '      Stokes drift + Stokes Coriolis + tracer advection ln_sdw = ', ln_sdw
         WRITE(numout,*) '      Breivik 2016                       ln_breivikFV_2016 = ', ln_breivikFV_2016
         WRITE(numout,*) '      Vortex Force                         ln_vortex_force = ', ln_vortex_force
         WRITE(numout,*) '      Bernouilli Head Pressure                ln_bern_srfc = ', ln_bern_srfc
         WRITE(numout,*) '      Stress modificated by waves (forced case)   ln_tauoc = ', ln_tauoc
         WRITE(numout,*) '      Stress modificated by waves (coupled case)    ln_taw = ', ln_taw
         WRITE(numout,*) '      Neutral drag coefficient                     ln_cdgw = ', ln_cdgw
         WRITE(numout,*) '      Charnock coefficient                        ln_charn = ', ln_charn
         WRITE(numout,*) '      TKE flux from wave                          ln_phioc = ', ln_phioc
         WRITE(numout,*) '      Surface shear with Stokes drift           ln_stshear = ', ln_stshear
         WRITE(numout,*) '      Test with constant wave fields          ln_wave_test = ', ln_wave_test
         WRITE(numout,*) '      Read wave energy spectrum               ln_wave_spec = ', ln_wave_spec
         WRITE(numout,*) '         Number of frequency classes             nn_nwfreq = ', nn_nwfreq
         WRITE(numout,*) '         Frequencies defined by user          ln_wfreq_usr = ', ln_wfreq_usr
         WRITE(numout,*) '            Exponentially-spaced limits   ln_wfreq_usr_exp = ', ln_wfreq_usr_exp
         WRITE(numout,*) '         Lowest frequency (ln_wfreq_usr=F)      rn_wfreq_0 = ', rn_wfreq_0
         WRITE(numout,*) '         Interval scale factor (ln_wfreq_usr=F) rn_wfreq_k = ', rn_wfreq_k
      ENDIF

      !                                ! option check
      IF( .NOT.( ln_cdgw .OR. ln_sdw .OR. ln_tauoc .OR. ln_charn .OR. ln_ice_wav ) )   &
         &     CALL ctl_warn( 'Wave coupling/forcing activated but ln_cdgw=F, ln_sdw=F, ln_tauoc=F, ln_ice_wav=F')
      IF( ln_cdgw .AND. ln_blk )   &
         &     CALL ctl_warn( 'drag coefficient read from wave model available ONLY with ln_NCAR and ln_MFS aerobulk options')
      IF( ln_cdgw .AND. ln_charn )   &
         &     CALL ctl_stop( 'ln_cdgw=T and ln_charn=T, but only one option can be choosen')
      IF( ln_tauoc .AND. ln_taw )   &
         &     CALL ctl_stop( 'ln_tauoc=T and ln_taw=T, but only one option can be choosen')
      IF( ln_stshear .AND. .NOT.ln_sdw )   &
         &     CALL ctl_stop( 'Stokes shear term calculated only if activated Stokes Drift ln_sdw=T')
      IF( ln_bern_srfc .AND. .NOT.cpl_bhd )   &
         &     CALL ctl_stop( 'ln_bern_srfc option only works in coupled mode')

      !            !== Check options for wave/sea-ice interactions ==!
      !
      ! Unavoidable to do this here due to order of initialisation routines and hence namelist reads.
      ! Sea-ice namelists -- specifically, namwav for wave-ice interaction module, icewav -- is read
      ! before namsbc_wave but must check namwav options w.r.t. choice of ln_wave_spec in namsbc_wave.
      !
      IF( ln_ice_wav ) THEN
         !
         IF( ln_ice_wav_spec .AND. .NOT.ln_wave_spec )   &
            &   CALL ctl_stop( 'ln_ice_wav_spec=T but spectrum NOT read because ln_wave_spec=F')
         !
         ! Warn if spectrum is being read but not used by icewav module (possible to do so, but strange to
         ! read full wave spectrum and *not* use it for wave breakup, so possibly a user oversight):
         IF( .NOT.ln_ice_wav_spec .AND. ln_wave_spec )   &
            &   CALL ctl_warn('ln_ice_wav_spec=F but spectrum IS being read in (ln_wave_spec=T); intentional?')
         !
      ENDIF

      !            !== Check options for setting frequencies for wave spectrum ==!
      !
      IF( rn_wfreq_0 <= 0._wp ) CALL ctl_stop( 'Smallest wave frequency rn_wfreq_0 cannot be negative'  )
      IF( rn_wfreq_k <= 1._wp ) CALL ctl_stop( 'Wave frequency increment factor rn_wfreq_k must be > 1' )

      !                !==  Set frequency arrays for wave energy spectrum  ==!
      !
      ALLOCATE( wfreq_l(nn_nwfreq), wfreq(nn_nwfreq), wfreq_u(nn_nwfreq), wdfreq(nn_nwfreq) )
      !
      IF( ln_wfreq_usr ) THEN
         !
         ! !== Set bounds from namelist ==!
         !
         wfreq_l(:) = rn_wfreq_usr(1:nn_nwfreq  )   ! lower limit of frequency bins
         wfreq_u(:) = rn_wfreq_usr(2:nn_nwfreq+1)   ! upper limit of frequency bins
         !
         IF( ln_wfreq_usr_exp ) THEN
            ! Put frequencies at fixed-ratio spacing assuming bin limits are exponentially spaced:
            DO jf = 1, nn_nwfreq
               wfreq(jf) = SQRT( wfreq_u(jf) / wfreq_l(jf) ) * wfreq_l(jf)
            ENDDO
         ELSE
            ! Put frequencies at mid-points of bin limits:
            wfreq(:) = .5_wp * (wfreq_l(:) + wfreq_u(:))
         ENDIF
         !
      ELSE
         ! !== Calculate exponentially-spaced frequency bin intervals such that f(n) = k*f(n-1) ==!
         !
         ! In general, frequency fn = f0 * k^n, for n = 0, 1, ..., corresponding to the bin boundaries,
         ! and n = 1/2, 3/2, ... corresponding to the 'central' values used in calculations. This is all
         ! determined by f0, k, and nn_wfreq (rn_wfreq_0, rn_wfreq_k, and nn_nwfreq, respectively)
         !
         wfreq_l(1) =                    rn_wfreq_0   ! lower limit of lowest frequency bin
         wfreq  (1) = SQRT(rn_wfreq_k) * rn_wfreq_0   ! 'central' frequency of lowest frequency bin
         wfreq_u(1) =       rn_wfreq_k * rn_wfreq_0   ! upper limit of lowest frequency bin
         !
         DO jf = 2, nn_nwfreq
            wfreq_l(jf) = rn_wfreq_k * wfreq_l(jf-1)
            wfreq  (jf) = rn_wfreq_k * wfreq  (jf-1)
            wfreq_u(jf) = rn_wfreq_k * wfreq_u(jf-1)
         ENDDO
         !
      ENDIF

      wdfreq(:) = wfreq_u(:) - wfreq_l(:)   ! width of frequency bins

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) '      Frequency (Hz) bins for discretised wave spectrum (n / lower bound / frequency / upper bound):'
         DO jf = 1, nn_nwfreq
            WRITE(numout,*) '      ', jf, wfreq_l(jf), wfreq(jf), wfreq_u(jf)
         ENDDO
         WRITE(numout,*)
      ENDIF

      !                             !==  Allocate wave arrays  ==!
      IF( ln_sdw ) THEN
         ALLOCATE( ut0sd (jpi,jpj), vt0sd (jpi,jpj) )
         ALLOCATE( wmp   (jpi,jpj) )
         ALLOCATE( tsd2d (jpi,jpj), div_sd(jpi,jpj) )
         ALLOCATE( tusd  (jpi,jpj), tvsd  (jpi,jpj) )
         ALLOCATE( usd   (jpi,jpj,jpk), vsd(jpi,jpj,jpk), wsd(jpi,jpj,jpk) )
         usd   (:,:,:) = 0._wp
         vsd   (:,:,:) = 0._wp
         wsd   (:,:,:) = 0._wp
         wmp     (:,:) = 0._wp
         ut0sd   (:,:) = 0._wp
         vt0sd   (:,:) = 0._wp
         tusd    (:,:) = 0._wp
         tvsd    (:,:) = 0._wp
      ENDIF
      IF( ln_sdw .OR. ln_ice_wav ) THEN
         ALLOCATE( hsw(jpi,jpj) )
         hsw(:,:) = 0._wp
      ENDIF
      IF( ln_ice_wav ) THEN
         ALLOCATE( wpf(jpi,jpj), wspec(jpi,jpj,nn_nwfreq) )
         wpf  (:,:)   = 0._wp
         wspec(:,:,:) = 0._wp
      ENDIF
      IF( ln_zdfswm    .AND. .NOT. cpl_wnum ) ALLOCATE( wnum    (A2D(0) ) )
      IF( ln_bern_srfc .AND.       cpl_bhd  ) ALLOCATE( bhd_wave(jpi,jpj) )
      !
      IF( ln_cdgw ) THEN                       ! wave drag from wave model
         IF( .NOT. cpl_wdrag ) THEN
            ALLOCATE( sf_cdg(1), STAT=ierror )               !* allocate and fill sf_wave with sn_cdg
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_cdg structure' )
                                   ALLOCATE( sf_cdg(1)%fnow(A2D(0),1)   )
            IF( sn_cdg%ln_tint )   ALLOCATE( sf_cdg(1)%fdta(A2D(0),1,2) )
            CALL fld_fill( sf_cdg, (/ sn_cdg /), cn_dir, 'sbc_wave_init', 'Wave module', 'namsbc_wave' )
         ENDIF
         ALLOCATE( cdn_wave(A2D(0)) )
         cdn_wave(:,:) = 0._wp
      ELSE IF ( ln_charn ) THEN                     ! wave drag from Charnok coefficient
         IF( .NOT. cpl_charn )   CALL ctl_stop( 'STOP', 'Charnock based wind stress can be used in coupled mode only' )
         ALLOCATE( charn(A2D(0)) )
         charn(:,:) = 0._wp
      ENDIF

      IF( ln_taw ) THEN                     ! wind/wave stress x/y components from wave model (coupled only)
         IF( .NOT. cpl_taw )   CALL ctl_stop( 'STOP', 'wind stress from wave model can be used in coupled mode only, use ln_cdgw instead' )
         ALLOCATE( tawx(A2D(0)), tawy(A2D(0)) )
         ALLOCATE( twox(A2D(0)), twoy(A2D(0)) )
         ALLOCATE( tauoc_wavex(A2D(0)), tauoc_wavey(A2D(0)) )
         tawx(:,:) = 0._wp
         tawy(:,:) = 0._wp
         twox(:,:) = 0._wp
         twoy(:,:) = 0._wp
         tauoc_wavex(:,:) = 1._wp
         tauoc_wavey(:,:) = 1._wp
      ENDIF

      IF( ln_phioc ) THEN                                  ! normalized wave energy flux
         IF( .NOT. cpl_phioc )   CALL ctl_stop( 'STOP', 'phioc can be used in coupled mode only' )
         ALLOCATE( phioc(A2D(0)) )
         phioc(:,:) = 0._wp
      ENDIF

      IF( ln_tauoc ) THEN                    ! normalized wave stress into the ocean
         IF( .NOT. cpl_wstrf ) THEN
            ALLOCATE( sf_tauoc(1), STAT=ierror )           !* allocate and fill sf_wave with sn_tauoc
            IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_tauoc structure' )
                                    ALLOCATE( sf_tauoc(1)%fnow(A2D(0),1)   )
            IF( sn_tauoc%ln_tint )  ALLOCATE( sf_tauoc(1)%fdta(A2D(0),1,2) )
            CALL fld_fill( sf_tauoc, (/ sn_tauoc /), cn_dir, 'sbc_wave_init', 'Wave module', 'namsbc_wave' )
         ENDIF
         ALLOCATE( tauoc_wave(A2D(0)) )
         tauoc_wave(:,:) = 0._wp
      ENDIF

      IF( ln_wave_test ) THEN       !==  Wave TEST case  ==!   set uniform waves fields
         !
         jpfld = 0   ! No field read
         !
         IF( ln_sdw ) THEN
            ut0sd(:,:) = 0.13_wp * tmask(:,:,1)   ! m/s
            vt0sd(:,:) = 0.00_wp                  ! m/s
            wmp  (:,:) = 8.00_wp                  ! seconds
         ENDIF
         IF( ln_sdw .OR. ln_ice_wav ) THEN
            hsw  (:,:) = 2.80_wp                  ! meters
         ENDIF
         IF( ln_ice_wav ) THEN
            wpf  (:,:) = 0.08_wp                  ! Hz
            !
            ! For wave spectrum test case calculate Bretschneider formula, which depends on wpf and hsw,
            ! using test case values defined above. See module icewav, function wav_spec_bret, for
            ! explanation of expression below (cannot use same function here due to circular dependence;
            ! may make more sense to move that function to sbcwave, i.e., this, module?)
            !
            DO_2D(0,0,0,0)
               wspec(ji,jj,:) = .3125_wp * hsw(ji,jj)**2 * (wpf(ji,jj)**4 / wfreq(:)**5)   &
                  &           * EXP( -1.25_wp * ( wpf(ji,jj) / wfreq(:) )**4 )
            END_2D
            !
         ENDIF
         !
      ELSE                          !==  create the structure associated with fields to be read  ==!
         
         IF( ln_ice_wav_spec ) THEN             ! wave energy spectrum (currently, needed for wave-ice interactions only)
            IF( .NOT. cpl_wspec ) THEN
               ALLOCATE( sf_wspec(1), STAT=ierror )         !* allocate and fill sf_wspec with sn_wpsec
               IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wspec structure' )
                                      ALLOCATE( sf_wspec(1)%fnow(jpi,jpj,nn_nwfreq) )
               IF( sn_wspec%ln_tint ) ALLOCATE( sf_wspec(1)%fdta(jpi,jpj,nn_nwfreq,2) )
               CALL fld_fill( sf_wspec, (/ sn_wspec /), cn_dir, 'sbc_wave', 'Wave module', 'namsbc_wave' )
            ENDIF
         ENDIF

         IF( ln_ice_wav ) THEN                  ! wave peak frequency (currently, needed for wave-ice interactions only)
            IF( .NOT. cpl_wpf ) THEN
               ALLOCATE( sf_wpf(1), STAT=ierror )           !* allocate and fill sf_wpf with sn_wpf
               IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wpf structure' )
                                    ALLOCATE( sf_wpf(1)%fnow(jpi,jpj,1) )
               IF( sn_wpf%ln_tint ) ALLOCATE( sf_wpf(1)%fdta(jpi,jpj,1,2) )
               CALL fld_fill( sf_wpf, (/ sn_wpf /), cn_dir, 'sbc_wave', 'Wave module', 'namsbc_wave' )
            ENDIF
         ENDIF

         IF( ln_sdw .OR. ln_ice_wav ) THEN      ! significant wave height
            IF( .NOT. cpl_hsig ) THEN
               ALLOCATE( sf_hsw(1), STAT=ierror )           !* allocate and fill sf_hsw with sn_hsw
               IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_hsw structure' )
                                    ALLOCATE( sf_hsw(1)%fnow(jpi,jpj,1) )
               IF( sn_hsw%ln_tint ) ALLOCATE( sf_hsw(1)%fdta(jpi,jpj,1,2) )
               CALL fld_fill( sf_hsw, (/ sn_hsw /), cn_dir, 'sbc_wave', 'Wave module', 'namsbc_wave' )
            ENDIF
         ENDIF

         IF( ln_sdw ) THEN                      ! Stokes drift
            ! 1. Find out how many fields have to be read from file if not coupled
            !    (NB. hsw now done separately, above, as it is also needed by icewav, not just sdw)
            jpfld=0
            jp_usd=0   ;   jp_vsd=0   ;   jp_wmp=0
            IF( .NOT. cpl_sdrft ) THEN
               jpfld  = jpfld + 1
               jp_usd = jpfld
               jpfld  = jpfld + 1
               jp_vsd = jpfld
            ENDIF
            IF( .NOT. cpl_wper ) THEN
               jpfld  = jpfld + 1
               jp_wmp = jpfld
            ENDIF
            ! 2. Read from file only the non-coupled fields
            IF( jpfld > 0 ) THEN
               ALLOCATE( slf_i(jpfld) )
               IF( jp_usd > 0 )   slf_i(jp_usd) = sn_usd
               IF( jp_vsd > 0 )   slf_i(jp_vsd) = sn_vsd
               IF( jp_wmp > 0 )   slf_i(jp_wmp) = sn_wmp
               ALLOCATE( sf_sd(jpfld), STAT=ierror )   !* allocate and fill sf_sd with stokes drift
               IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wave structure' )
               !
               DO ifpr= 1, jpfld
                                              ALLOCATE( sf_sd(ifpr)%fnow(jpi,jpj,1)   )
                  IF( slf_i(ifpr)%ln_tint )   ALLOCATE( sf_sd(ifpr)%fdta(jpi,jpj,1,2) )
               END DO
               !
               CALL fld_fill( sf_sd, slf_i, cn_dir, 'sbc_wave_init', 'Wave module', 'namsbc_wave' )
               sf_sd(jp_usd)%zsgn = -1._wp   ;  sf_sd(jp_vsd)%zsgn = -1._wp   ! vector field at T point: overwrite default definition of zsgn
               DEALLOCATE( slf_i )
            ENDIF
            !
            ! 3. Wave number (only needed for Qiao parametrisation, ln_zdfswm=T)
            IF( ln_zdfswm .AND. .NOT. cpl_wnum ) THEN
               ALLOCATE( sf_wnum(1), STAT=ierror )           !* allocate and fill sf_wave with sn_wnum
               IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_wave_init: unable to allocate sf_wnum structure' )
                                      ALLOCATE( sf_wnum(1)%fnow(A2D(0),1)   )
               IF( sn_wnum%ln_tint )  ALLOCATE( sf_wnum(1)%fdta(A2D(0),1,2) )
               CALL fld_fill( sf_wnum, (/ sn_wnum /), cn_dir, 'sbc_wave', 'Wave module', 'namsbc_wave' )
            ENDIF
            !
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE sbc_wave_init

   !!======================================================================
END MODULE sbcwave
