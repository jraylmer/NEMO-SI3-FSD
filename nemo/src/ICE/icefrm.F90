MODULE icefrm
   !!======================================================================
   !!                       ***  MODULE  icefrm  ***
   !! Sea-ice form drag param from CICE
   !!======================================================================
   !! History :  4.2        ! 2023  (D. Schroeder)       Form drag from CICE5
   !!            5.0        ! 2025  (J.R. Aylmer)        Add FSD dependence
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   USE par_ice
   USE ice   , ONLY : a_i, at_i, at_ip, vt_i, vt_s, drag_io, drag_ia           ! sea-ice variables
   USE icefsd, ONLY : a_ifsd, floe_size_dist, fsd_eff_size                     ! floe size distribution
   USE in_out_manager       ! I/O manager
   USE iom                  ! for iom_put
   USE timing               ! timing
   USE phycst                ! physical constants
   USE sbcblk , ONLY : nn_frm, rn_Cs_io, rn_Cs_ia, rn_Cr_ia, rn_Cr_io, rn_Cf_ia, rn_Cf_io   ! Form Drag params

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_frm        ! called by icestp.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_frm( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ice_frm  ***
      !!                   
      !! ** Purpose :   Neutral drag coefficients for turbulent momentum
      !!                exchange between ocean - sea-ice and atmosphere sea-ice 
      !!                calculated from ridge height, distance, floe size and pond
      !!                fraction based upon Tsamados et al. (2014), 
      !!                JPO, DOI: 10.1175/JPO-D-13-0215.1.
      !!
      !! ** Action  :   Calculates drag_io and drag_ia 
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER ::   ji, jj           ! dummy loop indices
      !
      ! Sail/keel parameters    
      REAL(wp), PARAMETER :: zazlpha   = 0._wp       ! weight functions for area of ridged ice 
      REAL(wp), PARAMETER :: zbeta2    = 0.75_wp       
      REAL(wp), PARAMETER :: zbeta     = 0.5_wp      ! power exponent appearing in astar and 
      ! L=Lmin(A*/(A*-A))**beta [0,1]
      REAL(wp), PARAMETER :: ztanar    = 0.4_wp      ! sail slope 
      REAL(wp), PARAMETER :: ztanak    = 0.4_wp      ! keel slope
      REAL(wp), PARAMETER :: zhkoverhr = 4._wp       ! hkeel/hridge ratio
      REAL(wp), PARAMETER :: zdkoverdr = 1._wp       ! dkeel/distrdg ratio
      REAL(wp), PARAMETER :: zphir     = 0.8_wp      ! porosity of ridges
      REAL(wp), PARAMETER :: zphik     = 0.8_wp      ! porosity of keels
      ! Floe parameters   
      REAL(wp), PARAMETER :: zLmin     = 8._wp       ! min length of floe (m) [5,100]
      REAL(wp), PARAMETER :: zLmax     = 300._wp     ! max length of floe (m) [30,3000]
      REAL(wp), PARAMETER :: zLmoy     = 300._wp     ! average length of floe (m) [30,1000]
      ! Melt pond parameters   
      REAL(wp), PARAMETER :: zlpmin    = 2.26_wp     ! min pond length (m) see Eq. 17 [1,10]
      REAL(wp), PARAMETER :: zlpmax    = 24.63_wp    ! max pond length (m) see Eq. 17 [10,100]
      REAL(wp), PARAMETER :: zCpond    = 0.2_wp      ! ratio of local form drag over geometrical parameter [0,1]
      ! Limits for stability
      REAL(wp), PARAMETER :: zcamax    = 0.0112_wp   ! Maximum for atmospheric drag
      REAL(wp), PARAMETER :: zcomax    = 0.04_wp     ! Maximum for ocean drag
      REAL(wp), PARAMETER :: zcamin    = 0.0007_wp   ! Minimum for atmospheric drag
      REAL(wp), PARAMETER :: zcomin    = 0.0025_wp   ! Minimum for ocean drag
      ! Other parameters   
      REAL(wp), PARAMETER :: zsHGB     = 0.18_wp     ! attenuation parameter for the sheltering function
      REAL(wp), PARAMETER :: zmrdgi    = 20._wp      ! screening effect see Lu2011 [5,50]
      REAL(wp), PARAMETER :: zmrdgo    = 10._wp      ! screening effect see Lu2011 [5,50]
      REAL(wp), PARAMETER :: zref      = 10._wp      ! reference height
      REAL(wp), PARAMETER :: zsl       = 22._wp      ! Sheltering parameter Lupkes2012 [10,30]
      REAL(wp), PARAMETER :: zorough   = 0.000327_wp ! ocean surface roughness (m)
      REAL(wp), PARAMETER :: zirough   = 0.0005_wp   ! ice surface roughness (m)
      ! Local variables
      REAL(wp) :: zhfreebd     ! freeboard (m)
      REAL(wp) :: zhdraft      ! draft of ice + snow column (Stoessel1993)
      REAL(wp) :: zhridge      ! ridge height
      REAL(wp) :: zdistrdg     ! distance between ridges
      REAL(wp) :: zhkeel       ! keel depth
      REAL(wp) :: zdkeel       ! distance between keels
      REAL(wp) :: zlfloe       ! floe length
      !
      REAL(wp) :: zai, zaip    ! ice area 
      REAL(wp) :: z1_ai        ! ice area inverse
      REAL(wp) :: zdha, zdho   ! temporary value for ridges and keels conditions
      REAL(wp) :: zsc          ! ocean and wind attenuation function
      REAL(wp) :: zastar       ! new constant for form drag
      REAL(wp) :: zlpond       ! pond length (m)
      !
      REAL(wp) :: zlog, z1_tanak, z1_tanar, z1_dkoverdr, z1_ilog, z1_olog
      !
      REAL(wp), PARAMETER ::   zardgi  = -0.004216_wp ! parameter used to estimate ardg value (intercept)
      REAL(wp), PARAMETER ::   zardgi1 =  0.157815_wp ! parameter used to estimate ardg value (fist degree)
      REAL(wp), PARAMETER ::   zardgi2 = -0.005002_wp ! parameter used to estimate ardg value (second degree)
      REAL(wp), PARAMETER ::   zvrdgi  = -0.059959_wp ! parameter used to estimate vrdg value (intercept)
      REAL(wp), PARAMETER ::   zvrdgi1 =  0.501150_wp ! parameter used to estimate vrdg value (first degree)
      REAL(wp), PARAMETER ::   zvrdgi2 =  0.077504_wp ! parameter used to estimate vrdg value (second degree)
      !
      REAL(wp), DIMENSION(A2D(nn_hls)) ::   zardg_drag    ! ridged ice concentration
      REAL(wp), DIMENSION(A2D(nn_hls)) ::   zvrdg_drag    ! ridged ice thickness
      REAL(wp), DIMENSION(A2D(nn_hls)) ::   zdrag_ia_skin ! neutral skin drag coefficient
      REAL(wp), DIMENSION(A2D(nn_hls)) ::   zdrag_ia_floe ! neutral floe edge drag coefficient
      REAL(wp), DIMENSION(A2D(nn_hls)) ::   zdrag_ia_pond ! neutral pond edge drag coefficient
      REAL(wp), DIMENSION(A2D(nn_hls)) ::   zdrag_ia_rdg  ! neutral ridge drag coefficient
      REAL(wp), DIMENSION(A2D(nn_hls)) ::   zdrag_io_skin ! skin drag coefficient
      REAL(wp), DIMENSION(A2D(nn_hls)) ::   zdrag_io_floe ! floe edge drag coefficient
      REAL(wp), DIMENSION(A2D(nn_hls)) ::   zdrag_io_keel ! keel drag coefficient
      REAL(wp), DIMENSION(A2D(0)) ::   zmsk00
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('ice_frm')
      !
      IF( kt == nit000 ) THEN   ! at first time-step
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ice_frm: sea-ice form drag param'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         IF(lwp) WRITE(numout,*)
      ENDIF
      !
      zmsk00(:,:) = MERGE( 1._wp, 0._wp, at_i(A2D(0)) >= epsi06  )
      !
      zastar = 1._wp / ( 1._wp - ( zLmin / zLmax ) ** ( 1._wp / zbeta ) )
      !
      z1_tanak = 1._wp / ztanak
      z1_tanar = 1._wp / ztanar
      z1_dkoverdr = zdkoverdr

      z1_ilog = 1._wp / LOG( zref / zirough )
      z1_olog = 1._wp / LOG( zref / zorough )
 
      !----------------------------------------------------------!
      !   Evaluation of ridged ice concentration and thickness   !
      !----------------------------------------------------------!
      ! Note: Currently there is no tracer for level or ridged-ice in SI3.
      ! Therefore, we preliminary use this approximation from Félicien Tournay
      ! and Antoine Barthélemy in 2017, at UCLouvain. They used model outputs
      ! from CICE to fit polynomials giving the volume and concentration of
      ! deformed ice as functions of sea ice properties. 
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         zardg_drag(ji,jj) = zardgi2 * ( vt_i(ji,jj) * vt_i(ji,jj) ) + zardgi1 * vt_i(ji,jj) + zardgi
         zvrdg_drag(ji,jj) = zvrdgi2 * ( vt_i(ji,jj) * vt_i(ji,jj) ) + zvrdgi1 * vt_i(ji,jj) + zvrdgi
         zardg_drag(ji,jj) = MIN( MAX( 0._wp , zardg_drag(ji,jj) ) , at_i(ji,jj) )
         zvrdg_drag(ji,jj) = MIN( MAX( 0._wp , zvrdg_drag(ji,jj) ) , vt_i(ji,jj) )

         ! For small values of ice volume, vrdg_drag can be zero and ardg_drag very small.
         ! To avoid problems with skin drags, ardg_drag is set to zero when vrdg_drag is.
         IF( zvrdg_drag(ji,jj) < epsi10 )   zardg_drag(ji,jj) = 0._wp
      END_2D

      !---------------!
      !   Main loop   !
      !---------------!
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )

         ! default values      
         zdrag_ia_skin(ji,jj) = rn_Cs_ia
         zdrag_io_skin(ji,jj) = rn_Cs_io
         zdrag_ia_rdg (ji,jj) = 0._wp
         zdrag_io_keel(ji,jj) = 0._wp
         zdrag_ia_floe(ji,jj) = 0._wp
         zdrag_io_floe(ji,jj) = 0._wp
         zdrag_ia_pond(ji,jj) = 0._wp

         ! start calculations
         IF( at_i(ji,jj) > 0.01_wp ) THEN
            
            zaip  = at_ip(ji,jj)
            zai   = at_i (ji,jj)
            z1_ai = 1._wp / zai
            
            IF( ln_fsd ) THEN
               !
               ! Calculate 'effective floe size' from FSD (FSTD integrated over all ice thicknesses)
               ! This replaces zlfloe without needing to change other parts of calculation, such that
               ! (zai / zlfloe) correctly gives the total cross-wind (or cross-current) floe edge per
               ! unit domain area in either case; see Tsamados et al. (2014) and docs for details.
               !
               ! Note: changing FSD floe shape parameter (rn_floeshape) will not consistently
               ! modify form drag due to floe edges since default Lupkes et al. (2012)
               ! formulation includes shape factor absorbed into coefficients (namelist
               ! parameters rn_Cf_i{a,o}), so cannot easily be accounted for at present.
               !
               zlfloe = fsd_eff_size( floe_size_dist(a_ifsd(ji,jj,:,:), a_i(ji,jj,:)) )
               !
            ELSE
               ! Original non-FSD implementation
               ! floe size parameterization see Eq. 13
               zlfloe = zLmin * ( zastar / ( zastar - zai ) ) ** zbeta
            ENDIF

            ! pond length
            zlpond = zlpmin * ( 1._wp - zaip ) + zlpmax * zaip

            ! draft and freeboard
            zhdraft  = ( rhoi * vt_i(ji,jj) + rhos * vt_s(ji,jj) ) * z1_ai * r1_rho0 ! without ponds
            zhfreebd = (        vt_i(ji,jj) +        vt_s(ji,jj) ) * z1_ai - zhdraft
            ! Do not allow draft larger than ice thickness (see Eq. 28)
            !    Clem: this condition should rarely be fulfilled 
            IF( zhdraft >= vt_i(ji,jj) * z1_ai ) THEN          ! replace excess snow with ice so hi ~= zhdraft
               zhfreebd = zhdraft * ( 1._wp - rhoi * r1_rho0 ) &
                  &     + ( vt_s(ji,jj) - ( vt_i(ji,jj) - zhdraft * zai ) * rhoi*r1_rhos ) * ( 1._wp - rhos*r1_rho0 ) * z1_ai ! Stoessel1993  
            ENDIF
            !
            IF( zardg_drag(ji,jj) > 0.01_wp ) THEN

               ! ridges and keels
               !    => hridge, hkeel, distrdg and dkeel estimates for simple triangular geometry as in the CICE code
               zhridge = 2._wp * zvrdg_drag(ji,jj) / zardg_drag(ji,jj)                                       &
                  &            * ( zazlpha + zbeta2 * zhkoverhr * z1_dkoverdr * ztanar * z1_tanak )          &
                  &            / ( zphir + zphik * ztanar * z1_tanak * zhkoverhr * zhkoverhr * z1_dkoverdr )
               zdistrdg = 2._wp * zhridge * zai / zardg_drag(ji,jj)                                          &
                  &             * ( zazlpha * z1_tanar + zbeta2 * z1_tanak * zhkoverhr * z1_dkoverdr )
               zhkeel = zhkoverhr * zhridge
               zdkeel = zdkoverdr * zdistrdg
               !
               !
               zdha = MAX( 0._wp , zhridge - zhfreebd )    ! Use the height of ridges relative to the mean freeboard of the pack
               zdho = MAX( 0._wp , zhkeel  - zhdraft  )    ! Use the height of keels  relative to the mean draft     of the pack
               !---------------!
               !   Skin drag   !
               !---------------!
               zdrag_ia_skin(ji,jj) = MAX( 0._wp, &
                  &                   MIN( rn_Cs_ia * ( 1._wp - zmrdgi * ( zdha / zdistrdg ) ) , zcamax ) )
               zdrag_io_skin(ji,jj) = MAX( 0._wp, &
                  &                   MIN( rn_Cs_io * ( 1._wp - zmrdgo * ( zdho / zdkeel   ) ) , zcomax ) )

               !----------------!
               !   Ridge drag   !
               !----------------!
               IF( zdha > epsi10 ) THEN
                  zsc = 1._wp - EXP( -zsHGB * zdistrdg / zdha )
                  zlog =  LOG( zdha / zirough ) 
                  zdrag_ia_rdg(ji,jj) = MAX( 0._wp, &
                     &                  MIN( rn_Cr_ia * 0.5_wp * zdha / zdistrdg * zsc * ( zlog * z1_ilog ) &
                     &                                                                 * ( zlog * z1_ilog ) , zcamax ) )
               ENDIF
               
               !---------------!
               !   Keel drag   !
               !---------------!
               IF( zdho > epsi10 ) THEN
                  zsc = 1._wp - EXP( -zsHGB * zdkeel / zdho )
                  zlog =  LOG( zdho / zirough )
                  zdrag_io_keel(ji,jj) = MAX( 0._wp, &
                     &                   MIN( rn_Cr_io * 0.5_wp * zdho / zdkeel * zsc * ( zlog * z1_ilog ) &
                     &                                                                * ( zlog * z1_ilog ) , zcomax ) )
               ENDIF

            ENDIF   ! ardg_drag > 0.01_wp

            !-----------------!
            ! Floe edge drag  !
            !-----------------!
            zsc =  1._wp - EXP( -zsl * zbeta * ( 1._wp - zai ) )
            IF( zhfreebd > epsi10 ) THEN
               zlog = LOG( zhfreebd / zorough )
               zdrag_ia_floe(ji,jj) = MAX( 0._wp, &
                  &                   MIN( rn_Cf_ia * 0.5_wp * zhfreebd / zlfloe * zsc * ( zlog * z1_olog ) &
                  &                                                                    * ( zlog * z1_olog ), zcamax ) )
            ENDIF
            IF( zhdraft > epsi10 ) THEN
               zlog = LOG( zhdraft / zorough )
               zdrag_io_floe(ji,jj) = MAX( 0._wp, &
                  &                   MIN( rn_Cf_io * 0.5_wp * zhdraft / zlfloe * zsc * ( zlog * z1_olog ) &
                  &                                                                   * ( zlog * z1_olog ), zcomax ) )
            ENDIF

            !----------------!
            ! Pond edge drag !
            !----------------!
            IF( zhfreebd > epsi10 ) THEN
               zsc = zaip ** ( 1._wp + 1._wp / ( zref * zbeta ) )
               zlog = LOG ( zhfreebd / zorough )
               zdrag_ia_pond(ji,jj) = MAX( 0._wp, &
                  &                   MIN( zCpond * 0.5_wp * zhfreebd / zlpond * zsc * ( zlog * z1_olog ) &
                  &                                                                  * ( zlog * z1_olog ), zcamax ) ) 
            ENDIF

         ENDIF
         
         !-------------------------!
         ! Total drag coefficients !
         !-------------------------!
         drag_ia(ji,jj) = MAX( zcamin, &
            &             MIN( zdrag_ia_skin(ji,jj) + zdrag_ia_floe(ji,jj) + zdrag_ia_rdg(ji,jj) + zdrag_ia_pond(ji,jj), zcamax ) ) &
            &             * tmask(ji,jj,1)
         drag_io(ji,jj) = MAX( zcomin, &
            &             MIN( zdrag_io_skin(ji,jj) + zdrag_io_floe(ji,jj) + zdrag_io_keel(ji,jj), zcomax ) ) &
            &             * tmask(ji,jj,1)
      END_2D

      ! outputs
      CALL iom_put( 'drag_io'     ,  drag_io     (A2D(0)) * zmsk00 )
      CALL iom_put( 'drag_io_skin', zdrag_io_skin(A2D(0)) * zmsk00 )
      CALL iom_put( 'drag_io_keel', zdrag_io_keel(A2D(0)) * zmsk00 )
      CALL iom_put( 'drag_io_floe', zdrag_io_floe(A2D(0)) * zmsk00 )
      !
      CALL iom_put( 'drag_ia'     ,  drag_ia     (A2D(0)) * zmsk00 )
      CALL iom_put( 'drag_ia_skin', zdrag_ia_skin(A2D(0)) * zmsk00 )
      CALL iom_put( 'drag_ia_rdg' , zdrag_ia_rdg (A2D(0)) * zmsk00 )
      CALL iom_put( 'drag_ia_floe', zdrag_ia_floe(A2D(0)) * zmsk00 )
      CALL iom_put( 'drag_ia_pond', zdrag_ia_pond(A2D(0)) * zmsk00 )
      !
      CALL iom_put( 'ardg_drag', zardg_drag(A2D(0)) * zmsk00 )
      CALL iom_put( 'vrdg_drag', zvrdg_drag(A2D(0)) * zmsk00 )
      !
      !
      IF( ln_timing )   CALL timing_stop('ice_frm')
      !
   END SUBROUTINE ice_frm

#else
    !!----------------------------------------------------------------------
    !!   Default option         Empty Module                No sea-ice model
    !!----------------------------------------------------------------------
CONTAINS
    SUBROUTINE ice_frm( kt )         ! Empty routine
       WRITE(*,*) 'ice_frm: You should not have seen this print! error?', kt
    END SUBROUTINE ice_frm
#endif


END MODULE icefrm
