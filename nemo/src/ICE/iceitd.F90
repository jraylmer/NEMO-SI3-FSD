MODULE iceitd
   !!======================================================================
   !!                       ***  MODULE iceitd ***
   !!   sea-ice : ice thickness distribution
   !!======================================================================
   !! History :  3.0  !  2005-12  (M. Vancoppenolle) original code (based on CICE)
   !!            4.0  !  2018     (many people)      SI3 [aka Sea Ice cube]
   !!----------------------------------------------------------------------
#if defined key_si3
   !!----------------------------------------------------------------------
   !!   'key_si3'                                       SI3 sea-ice model
   !!----------------------------------------------------------------------
   !!   ice_itd_rem   : redistribute ice thicknesses after thermo growth and melt
   !!   itd_glinear   : build g(h) satisfying area and volume constraints
   !!   itd_shiftice  : shift ice across category boundaries, conserving everything
   !!   ice_itd_reb   : rebin ice thicknesses into bounded categories
   !!   ice_itd_init  : read ice thicknesses mean and min from namelist
   !!----------------------------------------------------------------------
   USE par_ice        ! SI3 parameters
   USE phycst         ! physical constants
   USE ice            ! sea-ice: variables
   USE icevar  , ONLY : ice_var_roundoff
   USE icectl         ! sea-ice: conservation tests
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ice_itd_init  ! called in icestp
   PUBLIC   ice_itd_rem   ! called in icethd
   PUBLIC   ice_itd_reb   ! called in icecor

   INTEGER            ::   nice_catbnd     ! choice of the type of ice category function
   !                                       ! associated indices:
   INTEGER, PARAMETER ::   np_cathfn = 1   ! categories defined by a function
   INTEGER, PARAMETER ::   np_catusr = 2   ! categories defined by the user
   !
   !                                           !! ** namelist (namitd) **
   LOGICAL                    ::   ln_cat_hfn   ! ice categories are defined by function like rn_himean**(-0.05)
   REAL(wp)                   ::   rn_himean    ! mean thickness of the domain
   LOGICAL                    ::   ln_cat_usr   ! ice categories are defined by rn_catbnd
   REAL(wp), DIMENSION(0:100) ::   rn_catbnd    ! ice categories bounds
   REAL(wp)                   ::   rn_himax     ! maximum ice thickness allowed
   !
   !! * Substitutions
#  include "do_loop_substitute.h90"
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/ICE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ice_itd_rem( kt )
      !!------------------------------------------------------------------
      !!                ***  ROUTINE ice_itd_rem ***
      !!
      !! ** Purpose :   computes the redistribution of ice thickness
      !!                after thermodynamic growth of ice thickness
      !!
      !! ** Method  :   Linear remapping
      !!
      !! References :   W.H. Lipscomb, JGR 2001
      !!------------------------------------------------------------------
      INTEGER , INTENT (in) ::   kt      ! Ocean time step
      !
      INTEGER  ::   ji, jj, jl, jcat !, npti ! dummy loop index
      REAL(wp) ::   zx1, zwk1, zdh0, zetamin, zdamax   ! local scalars
      REAL(wp) ::   zx2, zwk2, zda0, zetamax           !   -      -
      REAL(wp) ::   zx3
      REAL(wp) ::   zslope          ! used to compute local thermodynamic "speeds"
      !
      LOGICAL , DIMENSION(A2D(0))       ::   iptidx          !
      INTEGER , DIMENSION(A2D(0),jpl-1) ::   jdonor          ! donor category index
      REAL(wp), DIMENSION(A2D(0),jpl)   ::   zdhice          ! ice thickness increment
      REAL(wp), DIMENSION(A2D(0),jpl)   ::   g0, g1          ! coefficients for fitting the line of the ITD
      REAL(wp), DIMENSION(A2D(0),jpl)   ::   hL, hR          ! left and right boundary for the ITD for each thickness
      REAL(wp), DIMENSION(A2D(0),jpl-1) ::   zdaice, zdvice  ! local increment of ice area and volume
      REAL(wp), DIMENSION(A2D(0))       ::   zhb0, zhb1      ! category boundaries for thinnes categories
      REAL(wp), DIMENSION(A2D(0),0:jpl) ::   zhbnew          ! new boundaries of ice categories
      !!------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('iceitd_rem')

      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_itd_rem: remapping ice thickness distribution'

      IF( ln_icediachk )   CALL ice_cons_hsm(0, 'iceitd_rem', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      IF( ln_icediachk )   CALL ice_cons2D  (0, 'iceitd_rem',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)

      !-----------------------------------------------------------------------------------------------
      !  1) Identify grid cells with ice
      !-----------------------------------------------------------------------------------------------
      at_i(A2D(0)) = SUM( a_i(A2D(0),:), dim=3 )
      !
      iptidx(:,:) = .false.
      DO_2D( 0, 0, 0, 0 )
         IF ( at_i(ji,jj) > epsi10 ) THEN
             iptidx(ji,jj) = .true.
         ENDIF
      END_2D
      !-----------------------------------------------------------------------------------------------
      !  2) Compute new category boundaries
      !-----------------------------------------------------------------------------------------------
      IF( ANY (iptidx(:,:) ) ) THEN
         !
         zdhice(:,:,:) = 0._wp
         zhbnew(:,:,:) = 0._wp
         !
         DO jl = 1, jpl
            ! Compute thickness change in each ice category
            DO_2D(0, 0, 0, 0)
               IF( a_i(ji,jj,jl) > epsi10 )   zdhice(ji,jj,jl) = h_i(ji,jj,jl) - h_i_b(ji,jj,jl)
            END_2D
         END DO
         !
         ! --- New boundaries for category 1:jpl-1 --- !
         DO jl = 1, jpl - 1
            !
            DO_2D(0, 0, 0, 0)
               !
               IF( iptidx(ji,jj) ) THEN
                  !
                  ! --- New boundary: Hn* = Hn + Fn*dt --- !
                  !     Fn*dt = ( fn + (fn+1 - fn)/(hn+1 - hn) * (Hn - hn) ) * dt = zdhice + zslope * (Hmax - h_i_b)
                  !
                  IF    ( a_i_b(ji,jj,jl) >  epsi10 .AND. a_i_b(ji,jj,jl+1) >  epsi10 ) THEN   ! a(jl+1) & a(jl) /= 0
                     zslope        = ( zdhice(ji,jj,jl+1) - zdhice(ji,jj,jl) ) / ( h_i_b(ji,jj,jl+1) - h_i_b(ji,jj,jl) )
                     zhbnew(ji,jj,jl) = hi_max(jl) + zdhice(ji,jj,jl) + zslope * ( hi_max(jl) - h_i_b(ji,jj,jl) )
                  ELSEIF( a_i_b(ji,jj,jl) >  epsi10 .AND. a_i_b(ji,jj,jl+1) <= epsi10 ) THEN   ! a(jl+1)=0 => Hn* = Hn + fn*dt
                     zhbnew(ji,jj,jl) = hi_max(jl) + zdhice(ji,jj,jl)
                  ELSEIF( a_i_b(ji,jj,jl) <= epsi10 .AND. a_i_b(ji,jj,jl+1) >  epsi10 ) THEN   ! a(jl)=0 => Hn* = Hn + fn+1*dt
                     zhbnew(ji,jj,jl) = hi_max(jl) + zdhice(ji,jj,jl+1)
                  ELSE                                                                       ! a(jl+1) & a(jl) = 0
                     zhbnew(ji,jj,jl) = hi_max(jl)
                  ENDIF
                  !
                  ! --- 2 conditions for remapping --- !
                  ! 1) hn(t+1)+espi < Hn* < hn+1(t+1)-epsi
                  !    Note: hn(t+1) must not be too close to either HR or HL otherwise a division by nearly 0 is possible
                  !          in itd_glinear in the case (HR-HL) = 3(Hice - HL) or = 3(HR - Hice)
# if defined key_single
                  IF( a_i(ji,jj,jl  ) > epsi10 .AND. h_i(ji,jj,jl  ) > ( zhbnew(ji,jj,jl) - epsi06 ) )   iptidx(ji,jj) = .false.
                  IF( a_i(ji,jj,jl+1) > epsi10 .AND. h_i(ji,jj,jl+1) < ( zhbnew(ji,jj,jl) + epsi06 ) )   iptidx(ji,jj) = .false.
# else
                  IF( a_i(ji,jj,jl  ) > epsi10 .AND. h_i(ji,jj,jl  ) > ( zhbnew(ji,jj,jl) - epsi10 ) )   iptidx(ji,jj) = .false.
                  IF( a_i(ji,jj,jl+1) > epsi10 .AND. h_i(ji,jj,jl+1) < ( zhbnew(ji,jj,jl) + epsi10 ) )   iptidx(ji,jj) = .false.
# endif
                  !
                  ! 2) Hn-1 < Hn* < Hn+1
                  IF( zhbnew(ji,jj,jl) < hi_max(jl-1) )   iptidx(ji,jj) = .false.
                  IF( zhbnew(ji,jj,jl) > hi_max(jl+1) )   iptidx(ji,jj) = .false.
                  !
               END IF ! iptidx == 1
            END_2D
         END DO
         !
         ! --- New boundaries for category jpl --- !
         DO_2D(0, 0, 0, 0)
            IF (iptidx(ji,jj) ) THEN
               IF( a_i(ji,jj,jpl) > epsi10 ) THEN
                  zhbnew(ji,jj,jpl) = MAX( hi_max(jpl-1), 3._wp * h_i(ji,jj,jpl) - 2._wp * zhbnew(ji,jj,jpl-1) )
               ELSE
                  zhbnew(ji,jj,jpl) = hi_max(jpl)
               ENDIF
               !
               ! --- 1 additional condition for remapping (1st category) --- !
               ! H0+epsi < h1(t) < H1-epsi
               !    h1(t) must not be too close to either HR or HL otherwise a division by nearly 0 is possible
               !    in itd_glinear in the case (HR-HL) = 3(Hice - HL) or = 3(HR - Hice)
# if defined key_single
               IF( h_i_b(ji,jj,1) < ( hi_max(0) + epsi06 ) )   iptidx(ji,jj) = .false.
               IF( h_i_b(ji,jj,1) > ( hi_max(1) - epsi06 ) )   iptidx(ji,jj) = .false.
# else
               IF( h_i_b(ji,jj,1) < ( hi_max(0) + epsi10 ) )   iptidx(ji,jj) = .false.
               IF( h_i_b(ji,jj,1) > ( hi_max(1) - epsi10 ) )   iptidx(ji,jj) = .false.
# endif
            ENDIF
         END_2D
         !
      ENDIF

      !-----------------------------------------------------------------------------------------------
      !  3) Compute g(h)
      !-----------------------------------------------------------------------------------------------
      IF( ANY (iptidx(:,:) ) ) THEN
         !
         zhb0(:,:) = hi_max(0)   ;   zhb1(:,:) = hi_max(1)
         g0(:,:,:) = 0._wp       ;   g1(:,:,:) = 0._wp
         hL(:,:,:) = 0._wp       ;   hR(:,:,:) = 0._wp
         !
         DO jl = 1, jpl
            !
            IF( jl == 1 ) THEN
               !
               ! --- g(h) for category 1 --- !
               CALL itd_glinear( iptidx(:,:), zhb0(:,:)  , zhb1(:,:)  , h_i_b(:,:,1)  , a_i(A2D(0),1)  ,  &   ! in
                  &              g0  (:,:,1), g1  (:,:,1), hL     (:,:,1), hR    (:,:,1)   )   ! out
               !
               ! Area lost due to melting of thin ice
               DO_2D(0, 0, 0, 0)
                  !
                  IF( ( iptidx(ji,jj) ) .AND. a_i(ji,jj,jl) > epsi10 ) THEN
                     !
                     zdh0 =  h_i(ji,jj,1) - h_i_b(ji,jj,1)
                     IF( zdh0 < 0.0 ) THEN      ! remove area from category 1
                        zdh0 = MIN( -zdh0, hi_max(1) )
                        !Integrate g(1) from 0 to dh0 to estimate area melted
                        zetamax = MIN( zdh0, hR(ji,jj,1) ) - hL(ji,jj,1)
                        !
                        IF( zetamax > 0.0 ) THEN
                           zx1    = zetamax
                           zx2    = 0.5 * zetamax * zetamax
                           zda0   = g1(ji,jj,1) * zx2 + g0(ji,jj,1) * zx1                ! ice area removed
                           zdamax = a_i(ji,jj,1) * (1.0 - h_i(ji,jj,1) / h_i_b(ji,jj,1) ) ! Constrain new thickness <= h_i
                           zda0   = MIN( zda0, zdamax )                            ! ice area lost due to melting of thin ice (zdamax > 0)
                           ! Remove area, conserving volume
                           h_i(ji,jj,1) = h_i(ji,jj,1) * a_i(ji,jj,1) / ( a_i(ji,jj,1) - zda0 )
                           a_i(ji,jj,1) = a_i(ji,jj,1) - zda0
                           v_i(ji,jj,1) = a_i(ji,jj,1) * h_i(ji,jj,1) ! useless ?
                        ENDIF
                        !
                     ELSE ! if ice accretion zdh0 > 0
                        ! zhbnew was 0, and is shifted to the right to account for thin ice growth in openwater (F0 = f1)
                        zhbnew(ji,jj,0) = MIN( zdh0, hi_max(1) )
                     ENDIF
                     !
                  ENDIF
                  !
               END_2D
               !
            ENDIF ! jl=1
            !
            ! --- g(h) for each thickness category --- !
            CALL itd_glinear( iptidx(:,:), zhbnew(:,:,jl-1), zhbnew(:,:,jl), h_i(A2D(0),jl)   , a_i(A2D(0),jl)   ,  &   ! in
               &              g0    (:,:,jl  ), g1    (:,:,jl), hL    (:,:,jl), hR    (:,:,jl)   )   ! out
            !
         END DO

         !-----------------------------------------------------------------------------------------------
         !  4) Compute area and volume to be shifted across each boundary (Eq. 18)
         !-----------------------------------------------------------------------------------------------
         DO jl = 1, jpl - 1
            !
            DO_2D(0, 0, 0, 0)
               !
               IF (iptidx(ji,jj) ) THEN
                  ! left and right integration limits in eta space
                  IF (zhbnew(ji,jj,jl) > hi_max(jl)) THEN ! Hn* > Hn => transfer from jl to jl+1
                     zetamin = MAX( hi_max(jl)   , hL(ji,jj,jl) ) - hL(ji,jj,jl)   ! hi_max(jl) - hL
                     zetamax = MIN( zhbnew(ji,jj,jl), hR(ji,jj,jl) ) - hL(ji,jj,jl)   ! hR - hL
                     jdonor(ji,jj,jl) = jl
                  ELSE                                 ! Hn* <= Hn => transfer from jl+1 to jl
                     zetamin = 0.0
                     zetamax = MIN( hi_max(jl), hR(ji,jj,jl+1) ) - hL(ji,jj,jl+1)  ! hi_max(jl) - hL
                     jdonor(ji,jj,jl) = jl + 1
                  ENDIF
                  zetamax = MAX( zetamax, zetamin ) ! no transfer if etamax < etamin
                  !
                  zx1  = zetamax - zetamin
                  zwk1 = zetamin * zetamin
                  zwk2 = zetamax * zetamax
                  zx2  = 0.5 * ( zwk2 - zwk1 )
                  zwk1 = zwk1 * zetamin
                  zwk2 = zwk2 * zetamax
                  zx3  = 1.0 / 3.0 * ( zwk2 - zwk1 )
                  jcat = jdonor(ji,jj,jl)
                  zdaice(ji,jj,jl) = g1(ji,jj,jcat)*zx2 + g0(ji,jj,jcat)*zx1
                  zdvice(ji,jj,jl) = g1(ji,jj,jcat)*zx3 + g0(ji,jj,jcat)*zx2 + zdaice(ji,jj,jl)*hL(ji,jj,jcat)
               !
               ENDIF
            END_2D
         END DO

         !----------------------------------------------------------------------------------------------
         ! 5) Shift ice between categories
         !----------------------------------------------------------------------------------------------
         CALL itd_shiftice ( iptidx, jdonor, zdaice, zdvice )
         !
      ENDIF ! ANY iptidx
      !
      ! the following fields need to be updated in the halos (done afterwards):
      ! a_i, v_i, v_s, sv_i, oa_i, h_i, a_ip, v_ip, v_il, t_su, e_i, e_s
      !
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'iceitd_rem', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      IF( ln_icediachk )   CALL ice_cons2D  (1, 'iceitd_rem',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)
      IF( ln_timing    )   CALL timing_stop ('iceitd_rem')
      !
   END SUBROUTINE ice_itd_rem


   SUBROUTINE itd_glinear( iptidx, HbL, Hbr, phice, paice, pg0, pg1, phL, phR )
      !!------------------------------------------------------------------
      !!                ***  ROUTINE itd_glinear ***
      !!
      !! ** Purpose :   build g(h) satisfying area and volume constraints (Eq. 6 and 9)
      !!
      !! ** Method  :   g(h) is linear and written as: g(eta) = g1(eta) + g0
      !!                with eta = h - HL
      !!------------------------------------------------------------------
      LOGICAL, DIMENSION(A2D(0)), INTENT(in   ) ::   iptidx
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   HbL, HbR      ! left and right category boundaries
      REAL(wp), DIMENSION(A2D(0)), INTENT(in   ) ::   phice, paice  ! ice thickness and concentration
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) ::   pg0, pg1      ! coefficients in linear equation for g(eta)
      REAL(wp), DIMENSION(A2D(0)), INTENT(inout) ::   phL, phR      ! min and max value of range over which g(h) > 0
      !
      INTEGER  ::   ji,jj        ! horizontal indices
      REAL(wp) ::   z1_3 , z2_3  ! 1/3 , 2/3
      REAL(wp) ::   zh13         ! HbL + 1/3 * (HbR - HbL)
      REAL(wp) ::   zh23         ! HbL + 2/3 * (HbR - HbL)
      REAL(wp) ::   zdhr         ! 1 / (hR - hL)
      REAL(wp) ::   zwk1, zwk2   ! temporary variables
      !!------------------------------------------------------------------
      !
      z1_3 = 1._wp / 3._wp
      z2_3 = 2._wp / 3._wp
      !
      DO_2D(0, 0, 0, 0)
         !
         IF (iptidx(ji,jj) ) THEN

            IF( paice(ji,jj) > epsi10  .AND. phice(ji,jj) > epsi10 )  THEN
               !
               ! Initialize hL and hR
               phL(ji,jj) = HbL(ji,jj)
               phR(ji,jj) = HbR(ji,jj)
               !
               ! Change hL or hR if hice falls outside central third of range,
               ! so that hice is in the central third of the range [HL HR]
               zh13 = z1_3 * ( 2._wp * phL(ji,jj) +         phR(ji,jj) )
               zh23 = z1_3 * (         phL(ji,jj) + 2._wp * phR(ji,jj) )
               !
               IF    ( phice(ji,jj) < zh13 ) THEN   ;   phR(ji,jj) = 3._wp * phice(ji,jj) - 2._wp * phL(ji,jj) ! move HR to the left
               ELSEIF( phice(ji,jj) > zh23 ) THEN   ;   phL(ji,jj) = 3._wp * phice(ji,jj) - 2._wp * phR(ji,jj) ! move HL to the right
               ENDIF
               !
               ! Compute coefficients of g(eta) = g0 + g1*eta
               IF( phR(ji,jj) > phL(ji,jj) ) THEN   ;   zdhr = 1._wp / (phR(ji,jj) - phL(ji,jj))
               ELSE                           ;   zdhr = 0._wp ! if hR=hL=hice => no remapping
               ENDIF
              !!zdhr = 1._wp / (phR(ii) - phL(ii))
               zwk1 = 6._wp * paice(ji,jj) * zdhr
               zwk2 = ( phice(ji,jj) - phL(ji,jj) ) * zdhr
               pg0(ji,jj) = zwk1 * ( z2_3 - zwk2 )                    ! Eq. 14
               pg1(ji,jj) = 2._wp * zdhr * zwk1 * ( zwk2 - 0.5_wp )   ! Eq. 14
               !
            ELSE  ! remap_flag = .false. or a_i < epsi10
               phL(ji,jj) = 0._wp
               phR(ji,jj) = 0._wp
               pg0(ji,jj) = 0._wp
               pg1(ji,jj) = 0._wp
            ENDIF
            !
         ENDIF ! iptidx
      END_2D
      !
   END SUBROUTINE itd_glinear


   SUBROUTINE itd_shiftice( iptidx, kdonor, pdaice, pdvice )
      !!------------------------------------------------------------------
      !!                ***  ROUTINE itd_shiftice ***
      !!
      !! ** Purpose :   shift ice across category boundaries, conserving everything
      !!              ( area, volume, energy, age*vol, and mass of salt )
      !!------------------------------------------------------------------
      LOGICAL , DIMENSION(A2D(0)), INTENT(in) ::   iptidx   !
      INTEGER , DIMENSION(A2D(0),jpl-1), INTENT(in) ::   kdonor   ! donor category index
      REAL(wp), DIMENSION(A2D(0),jpl-1), INTENT(in) ::   pdaice   ! ice area transferred across boundary
      REAL(wp), DIMENSION(A2D(0),jpl-1), INTENT(in) ::   pdvice   ! ice volume transferred across boundary
      !
      INTEGER  ::   ji, jj, jl, jk         ! dummy loop indices
      INTEGER  ::   jl2, jl1           ! local integers
      REAL(wp) ::   zworka, zworkv, ztrans ! ice/snow transferred
      REAL(wp), DIMENSION(A2D(0))            ::   ztmp             ! workspace
      REAL(wp), DIMENSION(A2D(0),jpl)        ::   zaTsfn           !  -    -
      !!------------------------------------------------------------------

      !----------------------------------------------------------------------------------------------
      ! 1) Define a variable equal to a_i*T_su
      !----------------------------------------------------------------------------------------------
      DO jl = 1, jpl
         DO_2D( 0, 0, 0, 0 )
            IF( iptidx(ji,jj) ) zaTsfn(ji,jj,jl) = a_i(ji,jj,jl) * t_su(ji,jj,jl)
         END_2D
      END DO

      !-------------------------------------------------------------------------------
      ! 2) Transfer volume and energy between categories
      !-------------------------------------------------------------------------------
      DO jl = 1, jpl - 1
         DO_2D(0, 0, 0, 0)
            !
            IF( iptidx(ji,jj) ) THEN
               jl1 = kdonor(ji,jj,jl)
               !
               IF( jl1 > 0 ) THEN
                  !
                  IF ( jl1 == jl  ) THEN   ;   jl2 = jl1+1
                  ELSE                     ;   jl2 = jl
                  ENDIF
                  !
                  IF( v_i(ji,jj,jl1) >= epsi10 ) THEN   ;   zworkv = pdvice(ji,jj,jl) / v_i(ji,jj,jl1)
                  ELSE                                  ;   zworkv = 0._wp
                  ENDIF
                  IF( a_i(ji,jj,jl1) >= epsi10 ) THEN   ;   zworka = pdaice(ji,jj,jl) / a_i(ji,jj,jl1)
                  ELSE                                  ;   zworka = 0._wp
                  ENDIF
                  !
                  a_i(ji,jj,jl1) = a_i(ji,jj,jl1) - pdaice(ji,jj,jl)       ! Ice areas
                  a_i(ji,jj,jl2) = a_i(ji,jj,jl2) + pdaice(ji,jj,jl)
                  !
                  v_i(ji,jj,jl1) = v_i(ji,jj,jl1) - pdvice(ji,jj,jl)       ! Ice volumes
                  v_i(ji,jj,jl2) = v_i(ji,jj,jl2) + pdvice(ji,jj,jl)
                  !
                  ztrans         = v_s(ji,jj,jl1) * zworkv              ! Snow volumes
                  v_s(ji,jj,jl1) = v_s(ji,jj,jl1) - ztrans
                  v_s(ji,jj,jl2) = v_s(ji,jj,jl2) + ztrans
                  !
                  ztrans          = oa_i(ji,jj,jl1) * zworka            ! Ice age
                  oa_i(ji,jj,jl1) = oa_i(ji,jj,jl1) - ztrans
                  oa_i(ji,jj,jl2) = oa_i(ji,jj,jl2) + ztrans
                  !
                  ztrans          = zaTsfn(ji,jj,jl1) * zworka             ! Surface temperature
                  zaTsfn(ji,jj,jl1)  = zaTsfn(ji,jj,jl1) - ztrans
                  zaTsfn(ji,jj,jl2)  = zaTsfn(ji,jj,jl2) + ztrans
                  !
                  IF ( ln_pnd_LEV .OR. ln_pnd_TOPO ) THEN
                     ztrans          = a_ip(ji,jj,jl1) * zworka         ! Pond fraction
                     a_ip(ji,jj,jl1) = a_ip(ji,jj,jl1) - ztrans
                     a_ip(ji,jj,jl2) = a_ip(ji,jj,jl2) + ztrans
                     !
                     ztrans          = v_ip(ji,jj,jl1) * zworkv         ! Pond volume
                     v_ip(ji,jj,jl1) = v_ip(ji,jj,jl1) - ztrans
                     v_ip(ji,jj,jl2) = v_ip(ji,jj,jl2) + ztrans
                     !
                     IF ( ln_pnd_lids ) THEN                            ! Pond lid volume
                        ztrans          = v_il(ji,jj,jl1) * zworkv
                        v_il(ji,jj,jl1) = v_il(ji,jj,jl1) - ztrans
                        v_il(ji,jj,jl2) = v_il(ji,jj,jl2) + ztrans
                     ENDIF
                  ENDIF
                  !
                  DO jk = 1, nlay_s                                     ! Snow heat content
                     ztrans            = e_s(ji,jj,jk,jl1) * zworkv
                     e_s(ji,jj,jk,jl1) = e_s(ji,jj,jk,jl1) - ztrans
                     e_s(ji,jj,jk,jl2) = e_s(ji,jj,jk,jl2) + ztrans
                  ENDDO
                  DO jk = 1, nlay_i                                     ! Ice heat content
                     ztrans            = e_i(ji,jj,jk,jl1) * zworkv
                     e_i(ji,jj,jk,jl1) = e_i(ji,jj,jk,jl1) - ztrans
                     e_i(ji,jj,jk,jl2) = e_i(ji,jj,jk,jl2) + ztrans
                  ENDDO
                  !                                                     ! Ice salinity
                  IF( nn_icesal == 4 ) THEN
                     DO jk = 1, nlay_i
                        ztrans              = szv_i(ji,jj,jk,jl1) * zworkv
                        szv_i(ji,jj,jk,jl1) = szv_i(ji,jj,jk,jl1) - ztrans
                        szv_i(ji,jj,jk,jl2) = szv_i(ji,jj,jk,jl2) + ztrans
                     ENDDO
                  ELSE
                     ztrans          = sv_i(ji,jj,jl1) * zworkv
                     sv_i(ji,jj,jl1) = sv_i(ji,jj,jl1) - ztrans
                     sv_i(ji,jj,jl2) = sv_i(ji,jj,jl2) + ztrans
                  ENDIF
                  !
               ENDIF   ! jl1 >0
            ENDIF ! iptidx true
         END_2D
            !
      END DO                   ! boundaries, 1 to jpl-1

      !-------------------
      ! 3) roundoff errors
      !-------------------
      ! clem: The transfer between one category to another can lead to very small negative values (-1.e-20)
      !       because of truncation error ( i.e. 1. - 1. /= 0 )

      CALL ice_var_roundoff( a_i, v_i, v_s, sv_i, oa_i, a_ip, v_ip, v_il, e_s, e_i, szv_i, iptidx )
      
      ! at_i must be <= rn_amax
      DO_2D( 0, 0, 0, 0 )
         IF( iptidx(ji,jj) ) THEN
            ztmp(ji,jj) = SUM( a_i(ji,jj,:) )
         ENDIF
      END_2D

      DO jl  = 1, jpl
         DO_2D(0, 0, 0, 0)
            IF( iptidx(ji,jj) .AND. ztmp(ji,jj) > rn_amax_2d(ji,jj) )   &
               &   a_i(ji,jj,jl) = a_i(ji,jj,jl) * rn_amax_2d(ji,jj) / ztmp(ji,jj)
         END_2D
      END DO

      !-------------------------------------------------------------------------------
      ! 4) Update ice thickness and temperature
      !-------------------------------------------------------------------------------
      DO jl = 1, jpl
         DO_2D(0, 0, 0, 0)
            IF( iptidx(ji,jj) ) THEN
# if defined key_single
               IF ( a_i(ji,jj,jl) >= epsi06 ) THEN
# else
               IF ( a_i(ji,jj,jl) >= epsi20 ) THEN
# endif
                  h_i (ji,jj,jl)  =  v_i(ji,jj,jl) / a_i(ji,jj,jl)
                  t_su(ji,jj,jl)  =  zaTsfn(ji,jj,jl) / a_i(ji,jj,jl)
               ELSE
                  h_i (ji,jj,jl)  = 0._wp
                  t_su(ji,jj,jl)  = rt0
               ENDIF
            ENDIF
         END_2D
      END DO
      !
      !
   END SUBROUTINE itd_shiftice


   SUBROUTINE ice_itd_reb( kt )
      !!------------------------------------------------------------------
      !!                ***  ROUTINE ice_itd_reb ***
      !!
      !! ** Purpose : rebin - rebins thicknesses into defined categories
      !!
      !! ** Method  : If a category thickness is out of bounds, shift part (for down to top)
      !!              or entire (for top to down) area, volume, and energy
      !!              to the neighboring category
      !!------------------------------------------------------------------
      LOGICAL , DIMENSION(A2D(0))       ::   iptidx          !
      INTEGER , INTENT (in) ::   kt      ! Ocean time step
      INTEGER ::   ji, jj, jl !, npti   ! dummy loop indices
      !
      INTEGER , DIMENSION(A2D(0),jpl-1) ::   jdonor           ! donor category index
      REAL(wp), DIMENSION(A2D(0),jpl-1) ::   zdaice, zdvice   ! ice area and volume transferred
      !!------------------------------------------------------------------
      IF( ln_timing )   CALL timing_start('iceitd_reb')
      !
      IF( kt == nit000 .AND. lwp )   WRITE(numout,*) '-- ice_itd_reb: rebining ice thickness distribution'
      !
      IF( ln_icediachk )   CALL ice_cons_hsm(0, 'iceitd_reb', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      IF( ln_icediachk )   CALL ice_cons2D  (0, 'iceitd_reb',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)
      !
      jdonor(:,:,:) = 0
      zdaice(:,:,:) = 0._wp
      zdvice(:,:,:) = 0._wp
      !
      !                       !---------------------------------------
      DO jl = 1, jpl-1        ! identify thicknesses that are too big
         !                    !---------------------------------------
         DO_2D( 0, 0, 0, 0 )
            IF( a_i(ji,jj,jl) > 0._wp .AND. v_i(ji,jj,jl) > (a_i(ji,jj,jl) * hi_max(jl)) ) THEN
               iptidx(ji,jj) = .true.
            ELSE
               iptidx(ji,jj) = .false.
            ENDIF
         END_2D
         !
         IF( ANY (iptidx(:,:) ) ) THEN
            !
            DO_2D(0, 0, 0, 0)
             IF (iptidx(ji,jj)) THEN
               jdonor(ji,jj,jl)  = jl
               ! how much of a_i you send in cat sup is somewhat arbitrary
               ! these are from CICE => transfer everything
               ! these are from LLN => transfer only half of the category
               zdaice(ji,jj,jl)  =                       0.5_wp  * a_i(ji,jj,jl)
               zdvice(ji,jj,jl)  = v_i(ji,jj,jl) - (1._wp - 0.5_wp) * a_i(ji,jj,jl) * hi_mean(jl)
             ENDIF
            END_2D
            !
            CALL itd_shiftice( iptidx, jdonor, zdaice, zdvice )  ! Shift jl=>jl+1
            ! Reset shift parameters
            jdonor(:,:,jl) = 0
            zdaice(:,:,jl) = 0._wp
            zdvice(:,:,jl) = 0._wp
         ENDIF
         !
      END DO
      !                       !-----------------------------------------
      DO jl = jpl-1, 1, -1    ! Identify thicknesses that are too small
         !                    !-----------------------------------------
         DO_2D( 0, 0, 0, 0 )
            IF( a_i(ji,jj,jl+1) > 0._wp .AND. v_i(ji,jj,jl+1) <= (a_i(ji,jj,jl+1) * hi_max(jl)) ) THEN
               iptidx(ji,jj) = .true.
            ELSE
               iptidx(ji,jj) = .false.
            ENDIF
         END_2D
         !
         IF( ANY (iptidx(:,:) ) ) THEN
            DO_2D(0, 0, 0, 0)
             IF (iptidx(ji,jj)) THEN
               jdonor(ji,jj,jl) = jl + 1
               zdaice(ji,jj,jl) = a_i(ji,jj,jl+1)
               zdvice(ji,jj,jl) = v_i(ji,jj,jl+1)
             ENDIF
            END_2D
            !
            CALL itd_shiftice( iptidx, jdonor, zdaice, zdvice )  ! Shift jl+1=>jl
            ! Reset shift parameters
            jdonor(:,:,jl) = 0
            zdaice(:,:,jl) = 0._wp
            zdvice(:,:,jl) = 0._wp
         ENDIF
         !
      END DO
      !
      ! clem: those fields must be updated on the halos: a_i, v_i, v_s, sv_i, oa_i, h_i, t_su, a_ip, v_ip, v_il, e_i, e_s
      !       note: ice_itd_reb is called in icedyn
      !             and in icethd (but once the arrays are already updated on the boundaries)
      !
      IF( ln_icediachk )   CALL ice_cons_hsm(1, 'iceitd_reb', rdiag_v, rdiag_s, rdiag_t, rdiag_fv, rdiag_fs, rdiag_ft)
      IF( ln_icediachk )   CALL ice_cons2D  (1, 'iceitd_reb',  diag_v,  diag_s,  diag_t,  diag_fv,  diag_fs,  diag_ft)
      IF( ln_timing    )   CALL timing_stop ('iceitd_reb')
      !
   END SUBROUTINE ice_itd_reb


   SUBROUTINE ice_itd_init
      !!------------------------------------------------------------------
      !!                ***  ROUTINE ice_itd_init ***
      !!
      !! ** Purpose :   Initializes the ice thickness distribution
      !! ** Method  :   ...
      !! ** input   :   Namelist namitd
      !!-------------------------------------------------------------------
      INTEGER  ::   jl            ! dummy loop index
      INTEGER  ::   ios, ioptio   ! Local integer output status for namelist read
      REAL(wp) ::   zhmax, znum, zden, zalpha   !   -      -
      !
      NAMELIST/namitd/ ln_cat_hfn, rn_himean, ln_cat_usr, rn_catbnd, rn_himin, rn_himax
      !!------------------------------------------------------------------
      !
      rn_catbnd(:) =  0._wp ! Circumvent possible initialization by compiler
                            ! to prevent from errors when writing output
      READ_NML_REF(numnam_ice,namitd)
      READ_NML_CFG(numnam_ice,namitd)
      IF(lwm) WRITE( numoni, namitd )
      !
      IF(lwp) THEN                  ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ice_itd_init: Initialization of ice cat distribution '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namitd: '
         WRITE(numout,*) '      Ice categories are defined by a function of rn_himean**(-0.05)    ln_cat_hfn = ', ln_cat_hfn
         WRITE(numout,*) '         mean ice thickness in the domain                               rn_himean  = ', rn_himean
         WRITE(numout,*) '      Ice categories are defined by rn_catbnd                           ln_cat_usr = ', ln_cat_usr
         WRITE(numout,*) '      minimum ice thickness allowed                                     rn_himin   = ', rn_himin
         WRITE(numout,*) '      maximum ice thickness allowed                                     rn_himax   = ', rn_himax
      ENDIF
      !
      !-----------------------------------!
      !  Thickness categories boundaries  !
      !-----------------------------------!
      !                             !== set the choice of ice categories ==!
      ioptio = 0
      IF( ln_cat_hfn ) THEN   ;   ioptio = ioptio + 1   ;   nice_catbnd = np_cathfn    ;   ENDIF
      IF( ln_cat_usr ) THEN   ;   ioptio = ioptio + 1   ;   nice_catbnd = np_catusr    ;   ENDIF
      IF( ioptio /= 1 )   CALL ctl_stop( 'ice_itd_init: choose one and only one ice categories boundaries' )
      !
      SELECT CASE( nice_catbnd )
      !                                !------------------------!
      CASE( np_cathfn )                ! h^(-alpha) function
         !                             !------------------------!
         zalpha = 0.05_wp
         zhmax  = 3._wp * rn_himean
         hi_max(0) = 0._wp
         DO jl = 1, jpl
            znum = jpl * ( zhmax+1 )**zalpha
            zden = REAL( jpl-jl , wp ) * ( zhmax + 1._wp )**zalpha + REAL( jl , wp )
            hi_max(jl) = ( znum / zden )**(1./zalpha) - 1
         END DO
         !                             !------------------------!
      CASE( np_catusr )                ! user defined
         !                             !------------------------!
         DO jl = 0, jpl
            hi_max(jl) = rn_catbnd(jl)
         END DO
         !
      END SELECT
      !
      DO jl = 1, jpl                ! mean thickness by category
         hi_mean(jl) = ( hi_max(jl) + hi_max(jl-1) ) * 0.5_wp
      END DO
      !
      hi_max(jpl) = rn_himax        ! set to a big value to ensure that all ice is thinner than hi_max(jpl)
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) '   ===>>>   resulting thickness category boundaries :'
      IF(lwp) WRITE(numout,*) '            hi_max(:)= ', hi_max(0:jpl)
      !
      IF( hi_max(1) < rn_himin )   CALL ctl_stop('ice_itd_init: the upper bound of the 1st category must be bigger than rn_himin')
      !
   END SUBROUTINE ice_itd_init

#else
   !!----------------------------------------------------------------------
   !!   Default option :         Empty module         NO SI3 sea-ice model
   !!----------------------------------------------------------------------
#endif

   !!======================================================================
END MODULE iceitd
